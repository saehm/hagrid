import { distance, get_scales } from "./utils.js";
import type { Point } from "./types.js";

/**
 * Reflects `[px, py]` across the diagonal of a `size`-wide quadrant, the
 * rotation step that keeps successive Hilbert quadrants connected.
 */
function rotate(size: number, px: number, py: number, rx: number, ry: number): Point {
    if (ry === 0) {
        if (rx === 1) {
            px = size - 1 - px;
            py = size - 1 - py;
        }
        return [py, px];
    }
    return [px, py];
}

/**
 * Maps a grid cell to its position along the Hilbert curve.
 *
 * @param point - Integer grid cell, both components in `[0, size)`.
 * @param size - Side length of the grid; must be a power of two.
 * @returns The curve index, in `[0, size²)`.
 */
export function hilbert_encode([px, py]: Point, size: number): number {
    let n = 0;
    for (let s = size / 2; s > 0; s /= 2) {
        const rx = (px & s) > 0 ? 1 : 0;
        const ry = (py & s) > 0 ? 1 : 0;
        n += s ** 2 * ((3 * rx) ^ ry);
        [px, py] = rotate(size, px, py, rx, ry);
    }
    return n;
}

/**
 * Maps a position along the Hilbert curve back to its grid cell — the inverse
 * of {@link hilbert_encode}.
 *
 * @param n - Curve index, in `[0, size²)`.
 * @param size - Side length of the grid; must be a power of two.
 */
export function hilbert_decode(n: number, size: number): Point {
    let t = n;
    let [px, py] = [0, 0];
    for (let s = 1; s < size; s *= 2) {
        const rx = 1 & (t / 2);
        const ry = 1 & (t ^ rx);
        [px, py] = rotate(s, px, py, rx, ry);
        px += s * rx;
        py += s * ry;
        t = t >> 2;
    }
    return [px, py];
}

/**
 * Places point `i` when its preferred cell `p` is taken, by walking outward
 * along the curve in both directions and claiming the nearer free cell.
 */
function hilbert_collision(
    P: Map<number, number>,
    p: number,
    d: Point,
    i: number,
    size: number,
): void {
    const e = size ** 2;
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
                const dl = distance(d, hilbert_decode(pl, size));
                const dr = distance(d, hilbert_decode(pr, size));
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

/** Parameters for {@link gridify_hilbert}. */
export interface HilbertParameters {
    /** Curve level; the grid is `2^level` cells on a side. */
    level: number;
    /** Pads the shorter axis so the data's aspect ratio survives the mapping. */
    keep_aspect_ratio?: boolean;
}

/**
 * Removes overlaps by snapping each point to a cell of a Hilbert curve,
 * resolving contested cells along the curve. Produces a square grid.
 *
 * @param D - Input points.
 * @returns One grid cell per input point, in input order.
 */
export function gridify_hilbert(
    D: readonly Point[],
    { level, keep_aspect_ratio = false }: HilbertParameters,
): Point[] {
    const size = 2 ** level;
    const curveLength = 4 ** level;
    const N = D.length;
    const P = new Map<number, number>();
    const Y: Point[] = new Array(N).fill(0);
    const scales = get_scales(D, [
        [0, size - 1],
        [0, size - 1],
    ], { keep_aspect_ratio, round: true });
    D.forEach((d, i) => {
        const [x, y] = scales.map((s) => s(d)) as Point;
        const p = hilbert_encode([x, y], size);
        if (P.has(p)) {
            // collision detected
            hilbert_collision(P, p, [x, y], i, curveLength);
        } else {
            P.set(p, i);
        }
    });

    for (const [p, i] of P.entries()) {
        Y[i] = hilbert_decode(p, size);
    }

    return Y;
}
