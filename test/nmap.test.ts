/**
 * N-Map partitions the plotting area, so its rectangles must tile the bounding
 * box exactly: areas summing to the whole, and no two overlapping.
 *
 * `sizes()` used to return the first half's share twice. The left half is
 * pinned to the region's left edge and the right half to its right edge, so
 * shares that do not add up to the full size leave a gap down the middle —
 * up to 38% of the plotting area, while the outer extent still looked correct.
 */
import { describe, expect, it } from "vitest";

import { gridify_nmap, nmap_regions } from "../src/nmap.js";
import { get_bounds } from "../src/utils.js";
import type { NMapCell } from "../src/nmap.js";
import type { Point } from "../src/types.js";

function rng(seed: number) {
    return () => {
        seed |= 0;
        seed = (seed + 0x6d2b79f5) | 0;
        let t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
        t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
}

function uniform(n: number, seed: number): Point[] {
    const r = rng(seed);
    return Array.from({ length: n }, () => [r() * 100, r() * 100]);
}

function overlap(a: NMapCell, b: NMapCell): number {
    const w = Math.min(a.x + a.width, b.x + b.width) - Math.max(a.x, b.x);
    const h = Math.min(a.y + a.height, b.y + b.height) - Math.max(a.y, b.y);
    return w > 0 && h > 0 ? w * h : 0;
}

describe("nmap rectangles tile the bounding box", () => {
    // Odd counts are the ones the bug hit: an even split has pA === pB, so the
    // wrong share happened to be the right one. 8 stayed correct throughout.
    for (const n of [8, 9, 50, 101, 257, 500]) {
        it(`N=${n}: areas sum to the full bounding box`, () => {
            const data = uniform(n, n * 13);
            const cells = nmap_regions(data);
            const bb = get_bounds(data);
            const total = cells.reduce((a, c) => a + Math.abs(c.width * c.height), 0);
            expect(cells).toHaveLength(n);
            expect(total / (bb.width * bb.height)).toBeCloseTo(1, 6);
        });

        it(`N=${n}: no two rectangles overlap`, () => {
            const cells = nmap_regions(uniform(n, n * 13));
            let worst = 0;
            for (let i = 0; i < cells.length; ++i) {
                for (let j = i + 1; j < cells.length; ++j) {
                    worst = Math.max(worst, overlap(cells[i], cells[j]));
                }
            }
            expect(worst).toBeLessThan(1e-9);
        });

        it(`N=${n}: every rectangle stays inside the bounding box`, () => {
            const data = uniform(n, n * 13);
            const bb = get_bounds(data);
            for (const c of nmap_regions(data)) {
                expect(c.x).toBeGreaterThanOrEqual(bb.x - 1e-9);
                expect(c.y).toBeGreaterThanOrEqual(bb.y - 1e-9);
                expect(c.x + c.width).toBeLessThanOrEqual(bb.x + bb.width + 1e-9);
                expect(c.y + c.height).toBeLessThanOrEqual(bb.y + bb.height + 1e-9);
            }
        });
    }
});

describe("nmap centres", () => {
    it("sit inside their own rectangle", () => {
        const data = uniform(200, 77);
        const cells = nmap_regions(data);
        const centres = gridify_nmap(data);
        expect(centres).toHaveLength(cells.length);
        centres.forEach(([x, y], i) => {
            const c = cells[i];
            expect(x).toBeGreaterThanOrEqual(Math.min(c.x, c.x + c.width) - 1e-9);
            expect(x).toBeLessThanOrEqual(Math.max(c.x, c.x + c.width) + 1e-9);
            expect(y).toBeGreaterThanOrEqual(Math.min(c.y, c.y + c.height) - 1e-9);
            expect(y).toBeLessThanOrEqual(Math.max(c.y, c.y + c.height) + 1e-9);
        });
    });

    it("are all distinct and finite", () => {
        for (const n of [9, 101, 500]) {
            const centres = gridify_nmap(uniform(n, n));
            for (const [x, y] of centres) {
                expect(Number.isFinite(x) && Number.isFinite(y)).toBe(true);
            }
            expect(new Set(centres.map((p) => p.join(","))).size).toBe(n);
        }
    });

    it("respects an explicit bounding box", () => {
        const data = uniform(64, 5);
        const BB = { x: -10, y: 20, width: 4, height: 8 };
        for (const c of nmap_regions(data, { BB })) {
            expect(c.x).toBeGreaterThanOrEqual(BB.x - 1e-9);
            expect(c.y).toBeGreaterThanOrEqual(BB.y - 1e-9);
            expect(c.x + c.width).toBeLessThanOrEqual(BB.x + BB.width + 1e-9);
            expect(c.y + c.height).toBeLessThanOrEqual(BB.y + BB.height + 1e-9);
        }
    });
});

describe("nmap rejects degenerate subdivisions instead of returning NaN", () => {
    // A region whose points all share the split coordinate has no extent to
    // divide, so the affine scale becomes share/0. This used to propagate
    // silently and the layout came back full of nulls.
    const degenerate: Record<string, Point[]> = {
        "vertical line": [[5, 0], [5, 1], [5, 2], [5, 3]],
        "horizontal line": [[0, 5], [1, 5], [2, 5], [3, 5]],
        "all identical": [[7, 7], [7, 7], [7, 7]],
        "duplicate pair": [[0, 0], [0, 0], [9, 9]],
        "10x10 integer lattice": Array.from({ length: 100 }, (_, i) => [
            i % 10,
            Math.floor(i / 10),
        ]),
    };

    for (const [name, data] of Object.entries(degenerate)) {
        it(`throws on ${name}`, () => {
            expect(() => gridify_nmap(data)).toThrow(/cannot subdivide/);
            expect(() => nmap_regions(data)).toThrow(/cannot subdivide/);
        });
    }

    it("names the axis and suggests a way out", () => {
        try {
            gridify_nmap([[5, 0], [5, 1], [5, 2], [5, 3]]);
            expect.unreachable("should have thrown");
        } catch (e) {
            const msg = (e as Error).message;
            expect(msg).toMatch(/along x/);
            expect(msg).toMatch(/gilbert|hilbert/);
        }
    });

    // Distinctness is not the property that matters — these are all distinct.
    it("still accepts data with spread on both axes", () => {
        const diagonal: Point[] = Array.from({ length: 16 }, (_, i) => [i, i]);
        expect(() => gridify_nmap(diagonal)).not.toThrow();
        expect(() => gridify_nmap(uniform(300, 42))).not.toThrow();
    });

    it("never returns a non-finite position when it does not throw", () => {
        for (let seed = 0; seed < 40; ++seed) {
            const data = uniform(50 + seed, seed * 31 + 1);
            let out: Point[] | null = null;
            try {
                out = gridify_nmap(data);
            } catch {
                continue; // rejected up front, which is the point
            }
            for (const [x, y] of out) {
                expect(Number.isFinite(x) && Number.isFinite(y)).toBe(true);
            }
        }
    });
});
