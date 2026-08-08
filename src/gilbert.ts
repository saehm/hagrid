/*
 * Generalized Hilbert ("gilbert") curve for arbitrary rectangles.
 *
 * This file is part of hagrid and is distributed under the LGPL-3.0-or-later,
 * like the rest of the library.
 *
 * The algorithm is Jakub Červený's, published at
 * https://github.com/jakubcerveny/gilbert under the BSD 2-Clause licence, and
 * this implementation follows it closely enough to be treated as derived from
 * it. Its notice is therefore reproduced, as that licence requires:
 *
 *   Copyright (c) 2018, Jakub Červený
 *
 *   Redistribution and use in source and binary forms, with or without
 *   modification, are permitted provided that the following conditions are met:
 *
 *   1. Redistributions of source code must retain the above copyright notice,
 *      this list of conditions and the following disclaimer.
 *
 *   2. Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *   POSSIBILITY OF SUCH DAMAGE.
 *
 * Indexing is `O(log(max(W, H)))` — the same shape of cost as
 * {@link hilbert_encode}, with a larger constant.
 */
import { get_bounds } from "./utils.js";
import type { Point } from "./types.js";

const sign = (v: number): number => (v > 0 ? 1 : v < 0 ? -1 : 0);

/** Does the rectangle spanned from `(x, y)` by vectors `a` and `b` hold `(px, py)`? */
function inside(
    px: number,
    py: number,
    x: number,
    y: number,
    ax: number,
    ay: number,
    bx: number,
    by: number,
): boolean {
    const x1 = x + (ax - sign(ax)) + (bx - sign(bx));
    const y1 = y + (ay - sign(ay)) + (by - sign(by));
    return (
        px >= Math.min(x, x1) &&
        px <= Math.max(x, x1) &&
        py >= Math.min(y, y1) &&
        py <= Math.max(y, y1)
    );
}

/**
 * Maps a cell of a `W × H` grid to its position along the curve.
 *
 * Unlike {@link hilbert_encode} the grid need not be square nor a power of
 * two, so elongated data can be gridded at its own aspect ratio.
 *
 * @param px - Column, in `[0, W)`.
 * @param py - Row, in `[0, H)`.
 * @returns The curve index, in `[0, W·H)`.
 */
export function gilbert_xy2d(px: number, py: number, W: number, H: number): number {
    let dst = 0;
    let x = 0;
    let y = 0;
    let ax: number;
    let ay: number;
    let bx: number;
    let by: number;
    if (W >= H) {
        ax = W; ay = 0; bx = 0; by = H;
    } else {
        ax = 0; ay = H; bx = W; by = 0;
    }

    for (;;) {
        const w = Math.abs(ax + ay);
        const h = Math.abs(bx + by);
        const dax = sign(ax), day = sign(ay);
        const dbx = sign(bx), dby = sign(by);

        // A single row or column is walked directly.
        if (h === 1) return dst + (dax * (px - x) + day * (py - y));
        if (w === 1) return dst + (dbx * (px - x) + dby * (py - y));

        let ax2 = Math.trunc(ax / 2), ay2 = Math.trunc(ay / 2);
        let bx2 = Math.trunc(bx / 2), by2 = Math.trunc(by / 2);
        const w2 = Math.abs(ax2 + ay2);
        const h2 = Math.abs(bx2 + by2);

        if (2 * w > 3 * h) {
            // Elongated: split the long axis in two, keeping each half even so
            // the halves can be traversed without a seam.
            if (w2 % 2 && w > 2) { ax2 += dax; ay2 += day; }
            if (inside(px, py, x, y, ax2, ay2, bx, by)) {
                ax = ax2; ay = ay2;
                continue;
            }
            dst += Math.abs((ax2 + ay2) * (bx + by));
            x += ax2; y += ay2; ax -= ax2; ay -= ay2;
            continue;
        }

        // Squarish: a step across, a long run, then a step back.
        if (h2 % 2 && h > 2) { bx2 += dbx; by2 += dby; }

        if (inside(px, py, x, y, bx2, by2, ax2, ay2)) {
            const nax = bx2, nay = by2;
            bx = ax2; by = ay2; ax = nax; ay = nay;
            continue;
        }
        dst += Math.abs((bx2 + by2) * (ax2 + ay2));

        if (inside(px, py, x + bx2, y + by2, ax, ay, bx - bx2, by - by2)) {
            x += bx2; y += by2; bx -= bx2; by -= by2;
            continue;
        }
        dst += Math.abs((ax + ay) * ((bx - bx2) + (by - by2)));

        const nx = x + (ax - dax) + (bx2 - dbx);
        const ny = y + (ay - day) + (by2 - dby);
        const nbx = -(ax - ax2), nby = -(ay - ay2);
        x = nx; y = ny;
        ax = -bx2; ay = -by2; bx = nbx; by = nby;
    }
}

/**
 * Maps a position along the generalized Hilbert curve back to its cell — the
 * inverse of {@link gilbert_xy2d}.
 *
 * @param d - Curve index, in `[0, W·H)`.
 */
export function gilbert_d2xy(d: number, W: number, H: number): Point {
    let x = 0;
    let y = 0;
    let ax: number;
    let ay: number;
    let bx: number;
    let by: number;
    if (W >= H) {
        ax = W; ay = 0; bx = 0; by = H;
    } else {
        ax = 0; ay = H; bx = W; by = 0;
    }

    for (;;) {
        const w = Math.abs(ax + ay);
        const h = Math.abs(bx + by);
        const dax = sign(ax), day = sign(ay);
        const dbx = sign(bx), dby = sign(by);

        if (h === 1) return [x + dax * d, y + day * d];
        if (w === 1) return [x + dbx * d, y + dby * d];

        let ax2 = Math.trunc(ax / 2), ay2 = Math.trunc(ay / 2);
        let bx2 = Math.trunc(bx / 2), by2 = Math.trunc(by / 2);
        const w2 = Math.abs(ax2 + ay2);
        const h2 = Math.abs(bx2 + by2);

        if (2 * w > 3 * h) {
            if (w2 % 2 && w > 2) { ax2 += dax; ay2 += day; }
            const first = Math.abs((ax2 + ay2) * (bx + by));
            if (d < first) {
                ax = ax2; ay = ay2;
                continue;
            }
            d -= first;
            x += ax2; y += ay2; ax -= ax2; ay -= ay2;
            continue;
        }

        if (h2 % 2 && h > 2) { bx2 += dbx; by2 += dby; }

        const first = Math.abs((bx2 + by2) * (ax2 + ay2));
        if (d < first) {
            const nax = bx2, nay = by2;
            bx = ax2; by = ay2; ax = nax; ay = nay;
            continue;
        }
        d -= first;

        const second = Math.abs((ax + ay) * ((bx - bx2) + (by - by2)));
        if (d < second) {
            x += bx2; y += by2; bx -= bx2; by -= by2;
            continue;
        }
        d -= second;

        const nx = x + (ax - dax) + (bx2 - dbx);
        const ny = y + (ay - day) + (by2 - dby);
        const nbx = -(ax - ax2), nby = -(ay - ay2);
        x = nx; y = ny;
        ax = -bx2; ay = -by2; bx = nbx; by = nby;
    }
}

/** Parameters for {@link gridify_gilbert}. */
export interface GilbertParameters {
    /** Grid width in cells. Derived from the data's aspect ratio if omitted. */
    width?: number;
    /** Grid height in cells. Derived from the data's aspect ratio if omitted. */
    height?: number;
    /** Spare cells to allow beyond the point count. Defaults to `1` (none). */
    whitespace?: number;
}

/**
 * Picks a `width × height` grid of at least `n` cells whose shape is as close
 * as possible to `aspect`.
 *
 * Both dimensions are integers, so only certain shapes exist at a given cell
 * count and the best one is not always nearest `sqrt(n / aspect)`. A few
 * heights either side of the ideal are tried, scored on the log ratio so that
 * twice too wide and twice too tall count equally, with a light penalty for
 * overshooting `n`.
 */
export function gilbert_dimensions(n: number, aspect: number): [number, number] {
    const a = Number.isFinite(aspect) && aspect > 0 ? aspect : 1;
    const count = Math.max(1, Math.ceil(n));
    const ideal = Math.sqrt(count / a);

    let best: [number, number] = [count, 1];
    let best_score = Infinity;
    const lo = Math.max(1, Math.floor(ideal) - 3);
    const hi = Math.max(1, Math.ceil(ideal) + 3);
    for (let height = lo; height <= hi; ++height) {
        const [width, h] = continuous(Math.max(1, Math.ceil(count / height)), height);
        const shape = Math.abs(Math.log(width / h / a));
        const slack = Math.abs(Math.log((width * h) / count));
        const score = shape + 0.25 * slack;
        if (score < best_score) {
            best_score = score;
            best = [width, h];
        }
    }
    return best;
}

/**
 * Nudges a grid off the one shape family whose curve is discontinuous:
 * an even short side with an odd long side. Growing the long side by one
 * escapes it and never removes cells.
 *
 * A discontinuity does not break correctness — the mapping stays a bijection —
 * but consecutive indices stop being neighbours in the plane, which is the
 * property collision resolution relies on.
 */
function continuous(width: number, height: number): [number, number] {
    const short = Math.min(width, height);
    const long = Math.max(width, height);
    if (short % 2 === 0 && long % 2 === 1) {
        return width >= height ? [width + 1, height] : [width, height + 1];
    }
    return [width, height];
}

/**
 * Removes overlaps by snapping points onto a generalized Hilbert curve whose
 * grid matches the data's aspect ratio, rather than padding to a square as
 * `hilbert` must.
 *
 * @param D - Input points.
 * @returns One `[column, row]` cell per input point, in input order.
 * @throws If the grid has fewer cells than there are points.
 */
export function gridify_gilbert(
    D: readonly Point[],
    parameters: GilbertParameters = {},
): Point[] {
    const N = D.length;
    const { x, y, width: dw, height: dh } = get_bounds(D);

    let W: number;
    let H: number;
    if (parameters.width !== undefined && parameters.height !== undefined) {
        W = parameters.width;
        H = parameters.height;
    } else {
        const target = Math.max(1, Math.ceil(N * (parameters.whitespace ?? 1)));
        [W, H] = gilbert_dimensions(target, dh > 0 ? dw / dh : 1);
    }
    if (W * H < N) {
        throw new Error(
            `gridify_gilbert: ${W}x${H} grid has ${W * H} cells for ${N} points`,
        );
    }

    // Scale the data onto the grid. A zero-width or zero-height extent collapses
    // to a single column or row rather than dividing by zero.
    const sx = dw > 0 ? (W - 1) / dw : 0;
    const sy = dh > 0 ? (H - 1) / dh : 0;

    const P = new Map<number, number>();
    const Y: Point[] = new Array(N).fill(0);
    const last = W * H - 1;

    D.forEach((d, i) => {
        const cx = Math.round((d[0] - x) * sx);
        const cy = Math.round((d[1] - y) * sy);
        const p = gilbert_xy2d(cx, cy, W, H);
        if (!P.has(p)) {
            P.set(p, i);
            return;
        }
        // Contested cell: walk outward along the curve for the nearest free one.
        let pl = p;
        let pr = p;
        for (;;) {
            ++pl;
            --pr;
            const vl = pl <= last;
            const vr = pr >= 0;
            if (!vl && !vr) {
                throw new Error("gridify_gilbert: more points than cells");
            }
            if (vl && !P.has(pl)) {
                P.set(pl, i);
                return;
            }
            if (vr && !P.has(pr)) {
                P.set(pr, i);
                return;
            }
        }
    });

    for (const [p, i] of P.entries()) {
        Y[i] = gilbert_d2xy(p, W, H);
    }
    return Y;
}
