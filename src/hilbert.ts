import { distance, get_scales } from "./utils.js";
import type { Point } from "./types.js";

/**
 * Largest curve level this implementation can index exactly.
 *
 * The index is a double, exact up to 2^53, and a level-`L` curve has `4^L`
 * cells. Beyond level 26 indices would start colliding silently, so they are
 * rejected instead.
 */
export const HILBERT_MAX_LEVEL = 26;

/** Hoisted so the guard costs a comparison, not a `Math.pow`, per call. */
const HILBERT_MAX_SIZE = 2 ** HILBERT_MAX_LEVEL;

function check_size(size: number): void {
    if (!(size > 0) || (size & (size - 1)) !== 0) {
        throw new Error(`hilbert: size must be a power of two, got ${size}`);
    }
    if (size > HILBERT_MAX_SIZE) {
        throw new Error(
            `hilbert: level above ${HILBERT_MAX_LEVEL} cannot be indexed exactly ` +
                `(size ${size} needs more than 2^53 distinct indices)`,
        );
    }
}

/**
 * Maps a grid cell to its position along the Hilbert curve.
 *
 * The quadrant reflection is written out inline rather than delegated to a
 * helper returning a tuple: allocating one array per level dominated
 * everything else and made this ~70x slower than it needed to be.
 *
 * @param point - Integer grid cell, both components in `[0, size)`.
 * @param size - Side length of the grid; must be a power of two.
 * @returns The curve index, in `[0, size²)`.
 */
export function hilbert_encode([px, py]: Point, size: number): number {
    check_size(size);
    let n = 0;
    for (let s = size >>> 1; s > 0; s >>>= 1) {
        const rx = (px & s) > 0 ? 1 : 0;
        const ry = (py & s) > 0 ? 1 : 0;
        n += s * s * ((3 * rx) ^ ry);
        // rotate: reflect across the quadrant diagonal, then swap
        if (ry === 0) {
            if (rx === 1) {
                px = size - 1 - px;
                py = size - 1 - py;
            }
            const t = px;
            px = py;
            py = t;
        }
    }
    return n;
}

/**
 * Maps a position along the Hilbert curve back to its grid cell — the inverse
 * of {@link hilbert_encode}.
 *
 * The index is kept in floating point rather than passed through `>>` and `&`,
 * which coerce to 32 bits: above level 15 that silently returned wrong
 * coordinates. Only the low bit is ever combined, and 0/1 operands are safe.
 *
 * @param n - Curve index, in `[0, size²)`.
 * @param size - Side length of the grid; must be a power of two.
 */
export function hilbert_decode(n: number, size: number): Point {
    check_size(size);
    let t = n;
    let px = 0;
    let py = 0;
    for (let s = 1; s < size; s *= 2) {
        const rx = Math.floor(t / 2) % 2;
        const ry = (t % 2) ^ rx;
        if (ry === 0) {
            if (rx === 1) {
                px = s - 1 - px;
                py = s - 1 - py;
            }
            const tmp = px;
            px = py;
            py = tmp;
        }
        px += s * rx;
        py += s * ry;
        t = Math.floor(t / 4);
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
        // Both probes have run off the ends, so every cell is taken. Without
        // this the walk spins forever, which is what a level too small for
        // the data used to produce.
        if (!vl && !vr) {
            throw new Error(
                "hilbert: more points than cells — raise the level (the grid has 4^level cells)",
            );
        }
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
