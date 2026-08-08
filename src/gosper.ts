import { distance, extent, get_scales } from "./utils.js";
import type { Point } from "./types.js";

/**
 * One refinement level of the Gosper L-system: the segment scale `s`, the
 * production symbols `t`, and the turn directions `d` in sixths of a turn.
 */
interface GosperLevel {
    s: number;
    t: ("a" | "b")[];
    d: number[];
}

/**
 * Expands the Gosper curve L-system up to `level` refinements.
 *
 * The two productions `a -> abbaaab` and `b -> abbbaab` (with their associated
 * turn sequences) are the standard Gosper rewriting rules.
 */
function create_gosper_fractal(level: number): GosperLevel[] {
    const t1: ("a" | "b")[] = ["a", "b", "b", "a", "a", "a", "b"];
    const d1 = [0, 5, 3, 4, 0, 0, 1];

    const t2: ("a" | "b")[] = ["a", "b", "b", "b", "a", "a", "b"];
    const d2 = [1, 0, 0, 4, 3, 5, 0];

    const SQRT7 = Math.sqrt(7);

    const result: GosperLevel[] = [{ s: SQRT7, t: ["a"], d: [0] }];

    for (let l = 1; l < level; ++l) {
        const previous_result = result[l - 1];
        result.push({
            s: previous_result.s * (1 / SQRT7),
            t: [],
            d: [],
        });

        const t = previous_result.t;
        const d = previous_result.d;

        for (let i = 0; i < t.length; ++i) {
            if (t[i] === "a") {
                result[l].t.push(...t1);
                result[l].d.push(...add_mod_6(d[i], d1));
            } else {
                result[l].t.push(...t2);
                result[l].d.push(...add_mod_6(d[i], d2));
            }
        }
    }

    return result;
}

function add_mod_6(m: number, d: readonly number[]): number[] {
    return d.map((e) => (m + e) % 6);
}

/** Walks a level's turn sequence to produce the curve's vertex positions. */
function generate_level(level: GosperLevel): Point[] {
    const k1 = 0.5;
    const k2 = Math.sqrt(3) / 2;
    const d_cos = [1, k1, -k1, -1, -k1, k1];
    const d_sin = [0, k2, k2, 0, -k2, -k2];
    const scale = level.s;
    const n = level.d.length + 1;
    const V: Point[] = [[0, 0]];
    for (let i = 1; i < n - 1; ++i) {
        const d = level.d[i];
        const [px, py] = V[i - 1];
        V.push([px + scale * d_cos[d], py + scale * d_sin[d]]);
    }
    return V;
}

/**
 * Maps a point to the index of the nearest cell on the Gosper curve.
 *
 * Materializes the whole curve on every call, so prefer {@link gosper_curve}
 * plus a reused lookup when encoding more than one point.
 */
export function gosper_encode(p: Point, level: number): number {
    const gosper_fractal = create_gosper_fractal(level + 1);
    const V = generate_level(gosper_fractal[level]);
    return find_nearest(p, V);
}

/** Maps a curve index back to its cell centre — the inverse of {@link gosper_encode}. */
export function gosper_decode(i: number, level: number): Point {
    const gosper_fractal = create_gosper_fractal(level + 1);
    const V = generate_level(gosper_fractal[level]);
    return V[i];
}

/** Returns every vertex of the Gosper curve at `level`, in curve order. */
export function gosper_curve(level: number): Point[] {
    const gosper_fractal = create_gosper_fractal(level + 1);
    return generate_level(gosper_fractal[level]);
}

/**
 * Places point `i` when its nearest cell `p` is taken, by walking outward along
 * the curve in both directions and claiming the nearer free cell.
 */
function gosper_collision(
    P: Map<number, number>,
    p: number,
    d: Point,
    i: number,
    V: readonly Point[],
): void {
    const e = V.length - 1;
    const valid = (q: number) => q >= 0 && q <= e;
    let pl = p;
    let pr = p;
    for (;;) {
        ++pl;
        --pr;
        const el = !P.has(pl);
        const er = !P.has(pr);
        const vl = valid(pl);
        const vr = valid(pr);
        if (vl && el && !er) {
            P.set(pl, i);
            return;
        } else if (!el && vr && er) {
            P.set(pr, i);
            return;
        } else if (el && er) {
            if (vl && vr) {
                const dl = distance(d, V[pl]);
                const dr = distance(d, V[pr]);
                P.set(dl < dr ? pl : pr, i);
                return;
            } else if (vl) {
                P.set(pl, i);
                return;
            } else if (vr) {
                P.set(pr, i);
                return;
            }
        }
    }
}

/** Parameters for {@link gridify_gosper}. */
export interface GosperParameters {
    /** Curve level; the curve has `7^level` cells. */
    level: number;
    /** Shrinks the target extent so the data sits inside the Gosper island. */
    scale_factor?: number;
}

/**
 * Removes overlaps by snapping each point to a cell of a Gosper curve,
 * resolving contested cells along the curve. Produces a hexagonal grid.
 *
 * Runs in `O(N · 7^level)`: every point is matched against the materialized
 * curve by linear scan. {@link gridify_gosper_uher} computes the index
 * directly instead, but does not yet produce correct results.
 *
 * @param D - Input points.
 * @returns One cell centre per input point, in input order.
 */
export function gridify_gosper(
    D: readonly Point[],
    { level, scale_factor = 0.8 }: GosperParameters,
): Point[] {
    const N = D.length;
    const P = new Map<number, number>();
    const Y: Point[] = new Array(N).fill(0);
    const V = gosper_curve(level);

    const grid_extent: [Point, Point] = [
        extent(V, (d) => d[0] * scale_factor),
        extent(V, (d) => d[1] * scale_factor),
    ];
    const scales = get_scales(D, grid_extent, { round: false });

    D.forEach((d, i) => {
        const [x, y] = scales.map((s) => s(d)) as Point;
        const p = find_nearest([x, y], V);
        if (P.has(p)) {
            gosper_collision(P, p, [x, y], i, V);
        } else {
            P.set(p, i);
        }
    });

    for (const [p, i] of P.entries()) {
        Y[i] = V[p];
    }

    return Y;
}

/** Index of the entry in `list` closest to `p`, or `-1` if `list` is empty. */
function find_nearest(p: Point, list: readonly Point[]): number {
    let nearest_index = -1;
    let nearest_distance = Infinity;

    list.forEach((q, i) => {
        const dist = distance(p, q);
        if (dist < nearest_distance) {
            nearest_index = i;
            nearest_distance = dist;
        }
    });

    return nearest_index;
}
