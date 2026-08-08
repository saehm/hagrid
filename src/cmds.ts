import { Delaunay } from "d3-delaunay";
import type { Point, Region } from "./types.js";

/** Clamps a point back inside the boundary `Gamma`. */
function gamma([px, py]: Point, { x1, x2, y1, y2 }: Region): Point {
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
    const delaunay = Delaunay.from(A as Point[]);
    return A.map((_, i) => {
        const neighbors = [...delaunay.neighbors(i)];
        return A.map((_, j) => {
            const pa = A[i];
            const qa = A[j];
            const po = B[i];
            const qo = B[j];
            const odistance = distance(po, qo);
            const overlap = compute_overlap(pa, qa, s);
            return {
                source: i,
                target: j,
                distance: odistance,
                overlap,
                delta: odistance + overlap,
            };
        }).filter(({ target, overlap }) => {
            const keep_overlap = overlap > 1e-2;
            const keep_neighbor = neighbors.indexOf(target) >= 0;
            return augment ? keep_overlap || keep_neighbor : keep_neighbor;
        });
    }).flat();
}

/** True once no overlaps remain and every point sits inside `Gamma`. */
function end_condition(
    proximity_graph: readonly ProximityEdge[],
    G: readonly Point[],
    { x1, x2, y1, y2 }: Region,
): boolean {
    const tol = 1;
    const sum_overlaps = proximity_graph.reduce((a, b) => a + b.overlap, 0);
    const inside_Gamma = G.map(
        ([px, py]) => px >= x1 && px <= x2 && py >= y1 && py <= y2,
    ).reduce((a, b) => a && b, true);
    return sum_overlaps < tol && inside_Gamma;
}

/** Parameters for {@link gridify_cmds}. */
export interface CMDSParameters {
    /** Boundary the layout is confined to. */
    Gamma: Region;
    /** Glyph size; also the grid pitch the result is snapped to. */
    size: number;
    /** Weight of the boundary constraint. Defaults to `0.1`. */
    alpha?: number;
    /** Maximum outer iterations. Defaults to `50`. */
    max_iter?: number;
}

/** One step of the optimisation: the current layout and its proximity graph. */
export type CMDSStep = [Point[], ProximityEdge[]];

/**
 * Removes overlaps with constrained MDS: stress majorization toward ideal
 * distances that include the glyph overlap, with the layout pulled back into
 * `Gamma`, then snapped to a grid of pitch `size`.
 *
 * Yields each intermediate layout, so callers can animate the optimisation.
 * The final yield is the grid-snapped result.
 *
 * @param D - Input points.
 */
export async function* gridify_cmds(
    D: readonly Point[],
    parameters: CMDSParameters,
): AsyncGenerator<CMDSStep, void, undefined> {
    const N = D.length;
    const Gamma = parameters.Gamma;
    const s = parameters.size;
    const alpha = parameters.alpha ?? 0.1;
    let max_iter = parameters.max_iter ?? 50;

    // copy D
    let G: Point[] = D.map(([x, y]) => [x, y] as Point);
    // construct proximity graph and calculate the ideal distances according to (3)
    let proximity_graph = create_proximity_graph(G, D, s, false);

    do {
        for (let it = 0; it < 20; ++it) {
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
                const [gamma_px, gamma_py] = gamma(p, Gamma);
                return [
                    (nominator[i][0] + alpha * gamma_px) / (denominator[i] + alpha),
                    (nominator[i][1] + alpha * gamma_py) / (denominator[i] + alpha),
                ] as Point;
            });
        }
        // reconstruct the proximity graph with the new layout
        proximity_graph = create_proximity_graph(G, D, s, true);

        yield [G, proximity_graph];
    } while (--max_iter > 0 && !end_condition(proximity_graph, G, Gamma));

    // round to grid
    G = G.map((d) => d.map((v) => Math.round((v - s / 2) / s) * s + s / 2) as Point);

    proximity_graph = create_proximity_graph(G, D, s, true);
    yield [G, proximity_graph];
}
