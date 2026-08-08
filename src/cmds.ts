import { Delaunay } from "d3-delaunay";
import type { Bounds, Point } from "./types.js";

/** A rectangle as opposing corners — the form the clamp and bound tests want. */
interface Corners {
    x1: number;
    x2: number;
    y1: number;
    y2: number;
}

/** Clamps a point back inside the boundary `Gamma`. */
function gamma([px, py]: Point, { x1, x2, y1, y2 }: Corners): Point {
    let [nx, ny] = [px, py];
    if (px < x1) nx = x1;
    else if (px > x2) nx = x2;
    if (py < y1) ny = y1;
    else if (py > y2) ny = y2;
    return [nx, ny];
}

function distance([ax, ay]: Point, [bx, by]: Point): number {
    return Math.hypot(bx - ax, by - ay);
}

/** Overlap between two glyphs of side `s` centred at `p` and `q`. */
function compute_overlap([px, py]: Point, [qx, qy]: Point, s: number): number {
    const dx = Math.abs(px - qx);
    const dy = Math.abs(py - qy);
    if (dx < s || dy < s) {
        const ox = Math.abs(Math.min(dx - s, 0));
        const oy = Math.abs(Math.min(dy - s, 0));
        const dd = distance([px, py], [qx, qy]);
        if (ox < oy) {
            return dx === 0 ? (dd * oy) / dy : (dd * ox) / dx;
        } else {
            return dy === 0 ? (dd * ox) / dx : (dd * oy) / dy;
        }
    }
    return 0;
}

/**
 * Overlap below this fraction of a node counts as none. Used both to pick
 * separation edges and to decide when the layout is done; the two must agree.
 */
const OVERLAP_EPS = 1e-3;

/** An edge of the proximity graph, carrying its ideal separation `delta`. */
export interface ProximityEdge {
    source: number;
    target: number;
    /** Distance between the two points in the original layout, `l_ij`. */
    distance: number;
    /** Overlap between the two glyphs, `delta_ij`. */
    overlap: number;
    /** Ideal distance, `d_ij = l_ij + delta_ij`. */
    delta: number;
}

/**
 * Builds the Delaunay proximity graph over the current layout `A`, with ideal
 * distances measured against the original layout `B`.
 *
 * @param augment - Also keeps overlapping pairs that are not Delaunay neighbours.
 */
function create_proximity_graph(
    A: readonly Point[],
    B: readonly Point[],
    s: number,
    augment = false,
): ProximityEdge[] {
    const n = A.length;
    const delaunay = Delaunay.from(A as Point[]);
    const edges: ProximityEdge[] = [];
    const eps = s * OVERLAP_EPS;
    const is_neighbour = new Uint8Array(n);

    const add = (i: number, j: number, overlap: number) => {
        const odistance = distance(B[i], B[j]);
        edges.push({
            source: i,
            target: j,
            distance: odistance,
            overlap,
            delta: odistance + overlap,
        });
    };

    for (let i = 0; i < n; ++i) {
        const pa = A[i];
        let touched: number[] = [];
        for (const j of delaunay.neighbors(i)) {
            is_neighbour[j] = 1;
            touched.push(j);
            add(i, j, compute_overlap(pa, A[j], s));
        }
        if (augment) {
            // Naive O(n^2) collision detection, as the paper prescribes, but
            // allocating only the edges it keeps. The threshold must match the
            // end condition's, or a pair can sit permanently between the two.
            for (let j = 0; j < n; ++j) {
                if (j === i || is_neighbour[j]) continue;
                const overlap = compute_overlap(pa, A[j], s);
                if (overlap > eps) add(i, j, overlap);
            }
        }
        for (const j of touched) is_neighbour[j] = 0;
        touched = [];
    }
    return edges;
}

/**
 * Algorithm 1 terminates when "all data items are inside Γ with no overlap".
 * Both halves are per-item, not aggregate.
 *
 * `Gamma` is already inset by half a node (§4.2), so a centre inside it means
 * the whole node is displayed. The epsilon absorbs the drift the weighted
 * average in (4) leaves behind; an exact comparison rejects layouts that
 * actually succeeded.
 */
function end_condition(
    proximity_graph: readonly ProximityEdge[],
    G: readonly Point[],
    { x1, x2, y1, y2 }: Corners,
    size: number,
): boolean {
    // A thousandth of a node: visually nothing, but comfortably above the
    // few-1e-4-of-a-node drift the weighted average in (4) leaves behind.
    const eps = size * OVERLAP_EPS;
    for (const e of proximity_graph) {
        if (e.overlap > eps) return false;
    }
    for (const [px, py] of G) {
        if (px < x1 - eps || px > x2 + eps || py < y1 - eps || py > y2 + eps) return false;
    }
    return true;
}

/**
 * The `m × n` grid of §4.2: small enough that every cell is used, large enough
 * to hold every item. The paper's conditions are `N ≤ m·n`, `m·(n−1) ≤ N` and
 * `(m−1)·n ≤ N`. A cell is also the node size.
 */
export function cmds_grid(width: number, height: number, count: number): [number, number] {
    const aspect = height > 0 && Number.isFinite(width / height) ? width / height : 1;
    let m = Math.max(1, Math.ceil(Math.sqrt(count * aspect)));
    let n = Math.max(1, Math.ceil(count / m));
    // Trim any row or column that is not needed; dropping one keeps m·n ≥ N.
    while (m > 1 && (m - 1) * n >= count) --m;
    while (n > 1 && m * (n - 1) >= count) --n;
    return [m, n];
}

/** Parameters for {@link gridify_cmds}. */
export interface CMDSParameters {
    /** Boundary the layout is confined to, as returned by `get_bounds`. */
    Gamma: Bounds;
    /**
     * Starting node size, also the grid pitch the result is snapped to.
     *
     * Defaults to a cell of the `m × n` grid of §4.2. If the layout does not
     * converge this shrinks, so it is a starting point rather than a fixed
     * value — see {@link gridify_cmds}.
     */
    size?: number;
    /** Weight of the boundary constraint. Defaults to `0.1`, per §4.2. */
    alpha?: number;
    /** Inner iterations per round, `N` in Algorithm 1. Defaults to `20`, per §4.2. */
    inner_iterations?: number;
    /**
     * Rounds of Algorithm 1 before the grid is grown. Defaults to `150`.
     *
     * The quality/time dial: overlaps fall slowly but steadily at a fixed node
     * size, so more rounds keep the nodes larger — roughly 2% of the plot area
     * at 10 rounds, 33% at 200, 71% at 800. Cost is linear in it.
     */
    max_iter?: number;
    /** How many times the grid may grow before giving up. Defaults to `8`. */
    max_rescales?: number;
}

/**
 * One step of the optimisation: the current layout, its proximity graph, and
 * the node size in force. The size changes when the grid is grown.
 */
export type CMDSStep = [Point[], ProximityEdge[], number];

/**
 * Removes overlaps with constrained MDS (Liu, Hu, North & Shen, *Computer
 * Graphics Forum* 37(1), 2018), then snaps the result to a grid.
 *
 * Alternates stress-majorization updates with rebuilding the Delaunay
 * proximity graph, augmented with edges for overlapping pairs, until every
 * item is inside `Gamma` with no overlap. If a run does not converge the grid
 * grows by a row or column and the nodes shrink (§4.2).
 *
 * Yields each intermediate layout, so callers can animate the optimisation;
 * the last yield is the grid-snapped result. Note that a run which stops with
 * overlaps left can snap two points into one cell, and that coincident input
 * points produce non-finite output.
 *
 * @param D - Input points.
 */
export function* gridify_cmds(
    D: readonly Point[],
    parameters: CMDSParameters,
): Generator<CMDSStep, void, undefined> {
    const N = D.length;
    const Gamma = parameters.Gamma;
    const alpha = parameters.alpha ?? 0.1;
    const inner = parameters.inner_iterations ?? 20;
    const max_iter = parameters.max_iter ?? 150;
    const max_rescales = parameters.max_rescales ?? 8;

    const W = Gamma.width;
    const H = Gamma.height;

    // §4.2: an m × n grid that fits Gamma and holds every item; a cell is the
    // node size. An explicit `size` just picks the starting grid instead.
    let [m, n] =
        parameters.size !== undefined && parameters.size > 0
            ? [
                  Math.max(1, Math.round(W / parameters.size)),
                  Math.max(1, Math.round(H / parameters.size)),
              ]
            : cmds_grid(W, H, N);

    // §4.2: the input is uniformly scaled to fit inside the inset Gamma, or
    // its extreme points start outside it and the boundary term can never be
    // satisfied. Fitting the tightest inset — the first, largest node size —
    // suffices, since later rescales only shrink nodes and widen the inset.
    const size0 = Math.min(W / m, H / n);
    const source = (() => {
        let min_x = Infinity, max_x = -Infinity, min_y = Infinity, max_y = -Infinity;
        for (const [x, y] of D) {
            if (x < min_x) min_x = x;
            if (x > max_x) max_x = x;
            if (y < min_y) min_y = y;
            if (y > max_y) max_y = y;
        }
        const dw = max_x - min_x;
        const dh = max_y - min_y;
        const iw = Math.max(W - size0, 0);
        const ih = Math.max(H - size0, 0);
        const k = Math.min(dw > 0 ? iw / dw : Infinity, dh > 0 ? ih / dh : Infinity);
        const scale = Number.isFinite(k) && k > 0 ? k : 1;
        // Centre the scaled layout inside Gamma.
        const cx = Gamma.x + W / 2 - ((min_x + max_x) / 2) * scale;
        const cy = Gamma.y + H / 2 - ((min_y + max_y) / 2) * scale;
        return D.map(([x, y]) => [x * scale + cx, y * scale + cy] as Point);
    })();

    let G: Point[] = source.map(([x, y]) => [x, y] as Point);
    let size = size0;
    let proximity_graph = create_proximity_graph(G, source, size, false);

    for (let rescale = 0; ; ++rescale) {
        size = Math.min(W / m, H / n);
        // Inset so a node whose centre is inside is drawn completely.
        const inset: Corners = {
            x1: Gamma.x + size / 2,
            x2: Gamma.x + W - size / 2,
            y1: Gamma.y + size / 2,
            y2: Gamma.y + H - size / 2,
        };

        proximity_graph = create_proximity_graph(G, source, size, false);
        let converged = false;

        for (let round = 0; round < max_iter; ++round) {
            for (let it = 0; it < inner; ++it) {
                // update step according to (4)
                const nominator: Point[] = Array.from({ length: N }, () => [0, 0]);
                const denominator: number[] = Array.from({ length: N }, () => 0);

                proximity_graph.forEach(({ source: i, target: j, delta }) => {
                    const pi = G[i];
                    const pj = G[j];
                    const wij = 1 / Math.pow(delta, 2);
                    const dist = distance(pi, pj);
                    nominator[i][0] += wij * (pj[0] + delta * ((pi[0] - pj[0]) / dist));
                    nominator[i][1] += wij * (pj[1] + delta * ((pi[1] - pj[1]) / dist));
                    denominator[i] += wij;
                });

                G = G.map((p, i) => {
                    const [gamma_px, gamma_py] = gamma(p, inset);
                    return [
                        (nominator[i][0] + alpha * gamma_px) / (denominator[i] + alpha),
                        (nominator[i][1] + alpha * gamma_py) / (denominator[i] + alpha),
                    ] as Point;
                });
            }
            // regenerate the graph, augmented with edges for overlapping nodes
            proximity_graph = create_proximity_graph(G, source, size, true);

            yield [G, proximity_graph, size];

            if (end_condition(proximity_graph, G, inset, size)) {
                converged = true;
                break;
            }
        }

        if (converged || rescale >= max_rescales) break;

        // §4.2: still overlapping, so grow the grid by one row or column and
        // rescale the nodes. Nodes are square with side min(W/m, H/n), so grow
        // whichever dimension currently sets that minimum — growing the other
        // one leaves the node size unchanged and the rescale achieves nothing.
        if (W / m <= H / n) ++m;
        else ++n;
    }

    // Grid alignment, (6).
    G = G.map((d) => d.map((v) => Math.round((v - size / 2) / size) * size + size / 2) as Point);

    proximity_graph = create_proximity_graph(G, source, size, true);
    yield [G, proximity_graph, size];
}
