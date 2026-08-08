/**
 * Gridfit, after Keim & Herrmann, IEEE Visualization '98.
 *
 * The invariant the whole algorithm rests on is `P(B) <= A(B)` (§3.2): every
 * subregion must have at least as many grid cells as it has points. If that
 * holds all the way down, every point gets its own cell — so the properties
 * worth pinning are that everything is placed, placed distinctly, and placed
 * close to where it started.
 */
import { describe, expect, it } from "vitest";

import { gridify_gridfit } from "../src/gridfit.js";
import { gridify } from "../src/gridify.js";
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

/** Integer grid coordinates, as the paper's pixel model assumes. */
function scatter(n: number, seed: number, spread = 2): Point[] {
    const r = rng(seed);
    const side = Math.ceil(Math.sqrt(n * spread));
    return Array.from({ length: n }, () => [
        Math.floor(r() * side),
        Math.floor(r() * side),
    ]);
}

const key = (p: Point) => p.join(",");

describe("gridfit places every point exactly once", () => {
    for (const n of [1, 2, 16, 100, 400, 1000, 5000]) {
        it(`N=${n}: one distinct cell per point`, () => {
            const data = scatter(n, n * 7);
            const cells = gridify_gridfit(data);
            expect(cells).toHaveLength(n);
            for (const [x, y] of cells) {
                expect(Number.isFinite(x) && Number.isFinite(y)).toBe(true);
            }
            expect(new Set(cells.map(key)).size).toBe(n);
        });
    }

    it("copes when the data has barely enough cells", () => {
        // Exactly N cells: the region cannot grow, so every cell must be used.
        const side = 10;
        const data: Point[] = Array.from({ length: side * side }, (_, i) => [
            i % side,
            Math.floor(i / side),
        ]);
        const cells = gridify_gridfit(data);
        expect(new Set(cells.map(key)).size).toBe(data.length);
    });

    it("handles degenerate inputs that defeat nmap", () => {
        const cases: Record<string, Point[]> = {
            identical: Array.from({ length: 50 }, () => [7, 7]),
            single_column: Array.from({ length: 50 }, (_, i) => [3, i]),
            single_row: Array.from({ length: 50 }, (_, i) => [i, 3]),
            heavy_duplicates: Array.from({ length: 200 }, (_, i) => [i % 5, Math.floor(i / 5) % 5]),
            negative: Array.from({ length: 80 }, (_, i) => [-50 + (i % 9), -30 + Math.floor(i / 9)]),
        };
        for (const [name, data] of Object.entries(cases)) {
            const cells = gridify_gridfit(data);
            expect(cells, name).toHaveLength(data.length);
            for (const [x, y] of cells) {
                expect(Number.isFinite(x) && Number.isFinite(y), name).toBe(true);
            }
            expect(new Set(cells.map(key)).size, name).toBe(data.length);
        }
    });

    it("accepts non-integer coordinates", () => {
        const r = rng(31);
        const data: Point[] = Array.from({ length: 200 }, () => [r() * 20, r() * 20]);
        const cells = gridify_gridfit(data);
        expect(new Set(cells.map(key)).size).toBe(data.length);
    });
});

describe("gridfit preserves position", () => {
    // Optimization goal (1) of §2.1 is to minimise the distance between a
    // point and its new position.
    for (const n of [100, 1000, 5000]) {
        it(`N=${n}: mean displacement stays below one cell`, () => {
            const data = scatter(n, n * 3);
            const cells = gridify_gridfit(data);
            const total = cells.reduce(
                (a, c, i) => a + Math.hypot(c[0] - data[i][0], c[1] - data[i][1]),
                0,
            );
            expect(total / n).toBeLessThan(1);
        });
    }

    it("keeps the layout inside a region with room for every point", () => {
        const n = 500;
        const data = scatter(n, 99);
        const cells = gridify_gridfit(data);
        const xs = cells.map((c) => c[0]);
        const ys = cells.map((c) => c[1]);
        const w = Math.max(...xs) - Math.min(...xs) + 1;
        const h = Math.max(...ys) - Math.min(...ys) + 1;
        expect(w * h).toBeGreaterThanOrEqual(n);
    });
});

describe("gridfit through gridify", () => {
    it("is reachable as a method and reports a runtime", () => {
        const data = scatter(300, 5);
        const { positions, runtime } = gridify(data, "gridfit");
        expect(positions).toHaveLength(data.length);
        expect(new Set(positions.map(key)).size).toBe(data.length);
        expect(runtime).toBeGreaterThanOrEqual(0);
    });
});
