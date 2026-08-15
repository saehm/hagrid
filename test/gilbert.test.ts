/**
 * The generalized Hilbert curve must be a genuine space-filling curve on any
 * rectangle: a bijection onto `[0, W·H)`, covering every cell, moving one step
 * at a time, and invertible.
 */
import { describe, expect, it } from "vitest";

import {
    gilbert_d2xy,
    gilbert_dimensions,
    gilbert_xy2d,
    gridify_gilbert,
} from "../src/gilbert.js";
import type { Point } from "../src/types.js";

const SHAPES: [number, number][] = [
    [1, 1], [1, 9], [9, 1], [2, 2], [5, 5], [8, 8], [16, 16], [32, 32],
    [10, 7], [7, 10], [13, 5], [3, 100], [100, 3], [64, 12], [37, 91], [128, 9],
];

function rng(seed: number) {
    return () => {
        seed |= 0;
        seed = (seed + 0x6d2b79f5) | 0;
        let t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
        t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
}

describe("gilbert_xy2d is a space-filling curve on any rectangle", () => {
    for (const [W, H] of SHAPES) {
        it(`${W}x${H}: bijection onto [0, ${W * H})`, () => {
            const seen = new Set<number>();
            for (let x = 0; x < W; ++x) {
                for (let y = 0; y < H; ++y) {
                    const d = gilbert_xy2d(x, y, W, H);
                    expect(d).toBeGreaterThanOrEqual(0);
                    expect(d).toBeLessThan(W * H);
                    expect(seen.has(d)).toBe(false);
                    seen.add(d);
                }
            }
            expect(seen.size).toBe(W * H);
        });

        it(`${W}x${H}: consecutive indices are adjacent cells`, () => {
            const cells: Point[] = new Array(W * H);
            for (let x = 0; x < W; ++x) {
                for (let y = 0; y < H; ++y) cells[gilbert_xy2d(x, y, W, H)] = [x, y];
            }
            for (let d = 1; d < W * H; ++d) {
                const [x0, y0] = cells[d - 1];
                const [x1, y1] = cells[d];
                expect(Math.abs(x1 - x0) + Math.abs(y1 - y0)).toBe(1);
            }
        });

        it(`${W}x${H}: d2xy inverts xy2d`, () => {
            for (let x = 0; x < W; ++x) {
                for (let y = 0; y < H; ++y) {
                    expect(gilbert_d2xy(gilbert_xy2d(x, y, W, H), W, H)).toEqual([x, y]);
                }
            }
        });
    }
});

describe("the curve is continuous wherever the library can put it", () => {
    /** True when every consecutive pair of indices is orthogonally adjacent. */
    function continuous(W: number, H: number): boolean {
        const cells: Point[] = new Array(W * H);
        for (let x = 0; x < W; ++x) {
            for (let y = 0; y < H; ++y) cells[gilbert_xy2d(x, y, W, H)] = [x, y];
        }
        if (cells.some((c) => !c)) return false;
        for (let d = 1; d < cells.length; ++d) {
            const [x0, y0] = cells[d - 1];
            const [x1, y1] = cells[d];
            if (Math.abs(x1 - x0) + Math.abs(y1 - y0) !== 1) return false;
        }
        return true;
    }

    // The original 16 hand-picked shapes all happened to be continuous, which
    // hid that 506 of the 2304 shapes up to 48x48 are not: every one with an
    // even short side and an odd long side, and nothing outside that class.
    it("every shape gilbert_dimensions can return is continuous", () => {
        const broken: string[] = [];
        for (let n = 1; n <= 400; n += 3) {
            for (const aspect of [0.05, 0.2, 0.5, 1, 2, 5, 20, 100]) {
                const [w, h] = gilbert_dimensions(n, aspect);
                if (w * h < n) broken.push(`${n}@${aspect} undersized ${w}x${h}`);
                else if (!continuous(w, h)) broken.push(`${n}@${aspect} -> ${w}x${h}`);
            }
        }
        expect(broken).toEqual([]);
    });

    it("dodges the even-short/odd-long family rather than returning it", () => {
        // 4x5 is the smallest member and is discontinuous if asked for directly.
        expect(continuous(4, 5)).toBe(false);
        expect(continuous(4, 6)).toBe(true);
        for (let n = 1; n <= 300; ++n) {
            for (const aspect of [0.3, 1, 3, 30]) {
                const [w, h] = gilbert_dimensions(n, aspect);
                const short = Math.min(w, h);
                const long = Math.max(w, h);
                expect(short % 2 === 0 && long % 2 === 1, `${n}@${aspect} -> ${w}x${h}`).toBe(false);
            }
        }
    });
});

describe("gilbert_dimensions matches the requested aspect", () => {
    it("always yields room for every point", () => {
        for (const n of [1, 7, 100, 1000, 12345]) {
            for (const aspect of [0.01, 0.5, 1, 3, 100]) {
                const [w, h] = gilbert_dimensions(n, aspect);
                expect(w * h).toBeGreaterThanOrEqual(n);
                expect(w).toBeGreaterThanOrEqual(1);
                expect(h).toBeGreaterThanOrEqual(1);
            }
        }
    });

    it("tracks the aspect ratio rather than squaring up", () => {
        const [w, h] = gilbert_dimensions(10000, 100);
        expect(w / h).toBeGreaterThan(50);
        // A square grid for the same point count would be 100x100.
        expect(h).toBeLessThan(100);
    });

    it("falls back to square on a degenerate aspect", () => {
        for (const bad of [0, -1, NaN, Infinity]) {
            const [w, h] = gilbert_dimensions(100, bad);
            expect(w * h).toBeGreaterThanOrEqual(100);
        }
    });
});

describe("gridify_gilbert", () => {
    const r = rng(1234);

    it("gives every point a distinct cell", () => {
        for (const n of [1, 2, 50, 500, 2000]) {
            const data: Point[] = Array.from({ length: n }, () => [r() * 100, r() * 100]);
            const positions = gridify_gilbert(data);
            expect(positions).toHaveLength(n);
            expect(new Set(positions.map((p) => p.join(","))).size).toBe(n);
        }
    });

    it("keeps an elongated dataset elongated", () => {
        // The point of the method: hilbert would pad this to a square.
        const data: Point[] = Array.from({ length: 400 }, (_, i) => [
            (i * 997) % 1000,
            (i * 31) % 10,
        ]);
        const positions = gridify_gilbert(data);
        const xs = positions.map((p) => p[0]);
        const ys = positions.map((p) => p[1]);
        const aspect =
            (Math.max(...xs) - Math.min(...xs)) / (Math.max(...ys) - Math.min(...ys));
        expect(aspect).toBeGreaterThan(10);
    });

    it("handles degenerate extents", () => {
        const cases: Record<string, Point[]> = {
            collinear_x: Array.from({ length: 40 }, (_, i) => [i, 0]),
            collinear_y: Array.from({ length: 40 }, (_, i) => [0, i]),
            identical: Array.from({ length: 20 }, () => [5, 5]),
            single: [[1, 1]],
        };
        for (const [name, data] of Object.entries(cases)) {
            const positions = gridify_gilbert(data);
            expect(positions, name).toHaveLength(data.length);
            for (const p of positions) {
                expect(Number.isFinite(p[0]) && Number.isFinite(p[1]), name).toBe(true);
            }
            expect(new Set(positions.map((p) => p.join(","))).size, name).toBe(data.length);
        }
    });

    it("honours an explicit grid and rejects one that is too small", () => {
        const data: Point[] = Array.from({ length: 30 }, () => [r() * 10, r() * 10]);
        const positions = gridify_gilbert(data, { cols: 10, rows: 5 });
        for (const [x, y] of positions) {
            expect(x).toBeGreaterThanOrEqual(0);
            expect(x).toBeLessThan(10);
            expect(y).toBeGreaterThanOrEqual(0);
            expect(y).toBeLessThan(5);
        }
        expect(() => gridify_gilbert(data, { cols: 4, rows: 4 })).toThrow(/cells for/);
    });

    it("rejects a half-given grid rather than silently deriving one", () => {
        const data: Point[] = Array.from({ length: 30 }, () => [r() * 10, r() * 10]);
        expect(() => gridify_gilbert(data, { cols: 10 })).toThrow(/together/);
        expect(() => gridify_gilbert(data, { rows: 5 })).toThrow(/together/);
    });

    it("spreads points out when given whitespace", () => {
        const data: Point[] = Array.from({ length: 100 }, () => [r() * 100, r() * 100]);
        const tight = gridify_gilbert(data);
        const loose = gridify_gilbert(data, { whitespace: 4 });
        const span = (ps: Point[]) => {
            const xs = ps.map((p) => p[0]);
            return Math.max(...xs) - Math.min(...xs);
        };
        expect(span(loose)).toBeGreaterThan(span(tight));
    });
});
