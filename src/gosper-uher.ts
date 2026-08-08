/**
 * Direct hierarchical indexing of the Gosper curve, after
 *
 *   Uher, Gajdoš, Snášel, Lai & Radecký,
 *   "Hierarchical Hexagonal Clustering and Indexing",
 *   Symmetry 11(6):731, 2019. https://doi.org/10.3390/sym11060731
 *
 * Unlike {@link gridify_gosper}, which scans a materialized curve, this
 * computes a point's curve index in `O(level)` — the approach the paper
 * promises. It is **not yet correct**: the cell coordinates it produces do not
 * match the reference implementation. Prefer {@link gridify_gosper} until the
 * discrepancy is resolved.
 */
import type { Point } from "./types.js";

/** A cube coordinate on the hexagonal lattice, where `x + y + z === 0`. */
type Cube = [number, number, number];

const REC_SQRT7 = 0.3779644730092272; // 1 / Math.sqrt(7)
const SQRT3_DIV_3 = 0.5773502691896257; // Math.sqrt(3) / 3
const ALPHA = 0.3334731722518321; // Math.asin(Math.sqrt(3) / (2 * Math.sqrt(7)))
const SQRT7 = 2.6457513110645907; // Math.sqrt(7)

/** Maps a raw hexagon index onto its position in the Gosper node pattern. */
const transform_index = [4, 0, 1, 2, 3, 6, 5];

/** Radii of the circles inscribed in the Gosper island, indexed by level. */
const inscribed_circles = [
    0.755928946, 0.755928946, 0.750121467308, 0.746782631146, 0.746782631146,
    0.746577727521, 0.746348363909, 0.746348363909, 0.746344578768, 0.746327538283,
    0.746327538283, 0.746327538283, 0.746326555879, 0.746326555879, 0.746326555879,
    0.746326510616, 0.746326510616, 0.746326510616, 0.746326508597, 0.746326508597,
    0.746326508597,
];

function xy2cube([x, y]: Point, size: number): Cube {
    const cx = (x * SQRT3_DIV_3 - (y * 1) / 3) / size;
    const cz = ((y * 2) / 3) / size;
    const cy = -cx - cz;
    return [cx, cy, cz];
}

function cube2xy([cx, , cz]: Cube, size: number): Point {
    const x = (3 / 2) * cx * size;
    const y = Math.sqrt(3) * (cx / 2 + cz) * size;
    return [x, y];
}

/**
 * Classifies which of the six neighbouring hexagons a fractional cube offset
 * points at, or `0` when the offset is negligible.
 *
 * Note: the original returned `undefined` on an exact tie between two dominant
 * axes, which propagated `NaN` through the index. Ties now fall through to `0`.
 */
function B([x, y, z]: Cube, thd = 10e-3): number {
    const [ax, ay, az] = [x, y, z].map((v) => Math.abs(v));
    if (ax < thd && ay < thd && az < thd) {
        return 0;
    } else if (ax > ay && ax > az) {
        return x > 0 ? 2 : 5;
    } else if (ay > ax && ay > az) {
        return y > 0 ? 4 : 1;
    } else if (az > ax && az > ay) {
        return z > 0 ? 6 : 3;
    }
    return 0;
}

/** Rounds a fractional cube coordinate to the nearest lattice hexagon. */
function cube_round([x, y, z]: Cube): Cube {
    let [l, m, n] = [x, y, z].map((d) => Math.round(d)) as Cube;
    const [dx, dy, dz] = [l - x, m - y, n - z].map((d) => Math.abs(d));
    if (dx > dy && dx > dz) {
        l = -m - n;
    } else if (dy > dz) {
        m = -l - n;
    } else {
        n = -l - m;
    }
    return [l, m, n];
}

/** Expands one level of the hierarchy, the inverse of the contraction in encode. */
function inverse_transform([x, , z]: Cube): Cube {
    const nx = 3 * x + z;
    const nz = -x + 2 * z;
    return [nx, -nx - nz, nz];
}

/** Side length of a level-0 hexagon such that `D` fits the Gosper island. */
export function gosper_uher_size(D: readonly Point[], level: number): number {
    let min_x = Infinity;
    let max_x = -Infinity;
    let min_y = Infinity;
    let max_y = -Infinity;
    D.forEach(([x, y]) => {
        min_x = x < min_x ? x : min_x;
        max_x = x > max_x ? x : max_x;
        min_y = y < min_y ? y : min_y;
        max_y = y > max_y ? y : max_y;
    });

    const half_diagonal =
        Math.sqrt(Math.pow(max_x - min_x, 2) + Math.pow(max_y - min_y, 2)) / 2;
    const s = half_diagonal / inscribed_circles[level];
    return s * Math.pow(REC_SQRT7, level);
}

// Rotations applied when descending into a child node, in radians.
const f_yArr = [-2.0943951023931953, 0, 0, -2.0943951023931953, 0, 2.0943951023931953, 0];
const gammaArr = [
    0, -1.0471975511965976, -2.0943951023931953, -3.141592653589793, 0,
    1.0471975511965976, 2.0943951023931953,
];

function rotate([x, y]: Point, alpha: number): Point {
    const sin_alpha = Math.sin(alpha);
    const cos_alpha = Math.cos(alpha);
    return [cos_alpha * x - sin_alpha * y, sin_alpha * x + cos_alpha * y];
}

function scale([x, y]: Point, s: number): Point {
    return [x * s, y * s];
}

/**
 * Maps a point to its position along the Gosper curve.
 *
 * @param point - The point to index.
 * @param level - Curve level; the curve has `7^level` cells.
 * @param size - Level-0 hexagon size, from {@link gosper_uher_size}.
 * @returns The curve index, in `[0, 7^level)`.
 */
export function gosper_uher_encode([xp, yp]: Point, level: number, size: number): number {
    const b: number[] = new Array(level).fill(0);
    let cube = xy2cube([xp, yp], size);
    let rounded = cube_round(cube);
    let p = cube2xy(rounded, size);
    for (let i = 1; i <= level; ++i) {
        p = scale(p, Math.pow(REC_SQRT7, i));
        p = rotate(p, -ALPHA * i);
        cube = xy2cube(p, size);
        rounded = cube_round(cube);
        b[level - i] = B([
            cube[0] - rounded[0],
            cube[1] - rounded[1],
            cube[2] - rounded[2],
        ]);
        if (i === level) continue;
        p = cube2xy(rounded, size);
        p = scale(p, Math.pow(SQRT7, i));
        p = rotate(p, ALPHA * i);
    }

    // Rewrite the raw hexagon sequence into Gosper node order, tracking the
    // rotation and traversal direction each node imposes on its children.
    let rotDir = 0;
    let order = true;
    const k: number[] = new Array(level).fill(0);
    for (let i = 0; i < level; ++i) {
        let bi = b[i];
        if (bi !== 0 && rotDir !== 0) {
            bi += 2 * rotDir;
            if (bi < 1) {
                bi += 6;
            } else if (bi > 6) {
                bi -= 6;
            }
        }
        const ki = transform_index[bi];
        if (ki === 0 || ki === 3) {
            --rotDir;
            if (rotDir < -1) rotDir = 1;
        } else if (ki === 5) {
            ++rotDir;
            if (rotDir > 1) rotDir = -1;
        }
        k[i] = order ? ki : 6 - ki;
        if (ki === 0 || ki === 4 || ki === 5) {
            order = !order;
        }
    }
    return k.reduce((a, digit) => a * 7 + digit, 0);
}

/**
 * Maps a position along the Gosper curve back to its cell centre — the inverse
 * of {@link gosper_uher_encode}.
 *
 * @param t - Curve index, in `[0, 7^level)`.
 * @param level - Curve level.
 * @param size - Level-0 hexagon size, from {@link gosper_uher_size}.
 */
export function gosper_uher_decode(t: number, level: number, size: number): Point {
    const k: number[] = new Array(level).fill(0);
    for (let i = level - 1; i >= 0; --i) {
        const mod_t = t % 7;
        k[i] = mod_t;
        t = (t - mod_t) / 7;
    }

    let d: Cube = [-size, 0, 0];
    for (let l = 0; l < level - 1; ++l) {
        d = inverse_transform(d);
    }
    let step = cube2xy(d, 1);

    let c: Point = [0, 0];
    let ord = true;
    for (let i = 0; i < level; ++i) {
        const ki = i > 0 && !ord ? 6 - k[i] : k[i];
        const dk = rotate(step, gammaArr[ki]);
        if (ki !== 4) {
            c = [c[0] + dk[0], c[1] + dk[1]];
        }
        step = rotate(step, ALPHA + f_yArr[ki]);
        step = scale(step, REC_SQRT7);
        if (ki === 0 || ki === 4 || ki === 5) {
            ord = !ord;
        }
    }
    return c;
}

function distance([ax, ay]: Point, [bx, by]: Point): number {
    return Math.hypot(bx - ax, by - ay);
}

/** Resolves a contested cell by walking outward along the curve. */
function collision(
    P: Map<number, number>,
    p: number,
    d: Point,
    i: number,
    level: number,
    size: number,
): void {
    const e = 7 ** level;
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
                const dl = distance(d, gosper_uher_decode(pl, level, size));
                const dr = distance(d, gosper_uher_decode(pr, level, size));
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

/** Parameters for {@link gridify_gosper_uher}. */
export interface GosperUherParameters {
    /** Curve level; the curve has `7^level` cells. */
    level: number;
}

/**
 * Removes overlaps using the direct Gosper indexing of Uher et al.
 *
 * **Known broken** — retained for the Phase 2 work of reconciling it with
 * {@link gridify_gosper}, which is correct but scans the whole curve.
 *
 * @param D - Input points.
 * @returns One cell centre per input point, in input order.
 */
export function gridify_gosper_uher(
    D: readonly Point[],
    { level }: GosperUherParameters,
): Point[] {
    const size = gosper_uher_size(D, level);
    const N = D.length;
    const P = new Map<number, number>();
    const Y: Point[] = new Array(N).fill(0);

    D.forEach((d, i) => {
        const [x, y] = d.map(Math.round) as Point;
        const p = gosper_uher_encode([x, y], level, size);
        if (P.has(p)) {
            collision(P, p, d, i, level, size);
        } else {
            P.set(p, i);
        }
    });

    for (const [p, i] of P.entries()) {
        Y[i] = gosper_uher_decode(p, level, size);
    }
    return Y;
}
