/**
 * Uher et al.'s direct Gosper indexing, verified against its own decode.
 *
 * decode was always correct; encode was not. The two faults were that
 * `xy2cube` used a pointy-topped hexagon convention while `cube2xy` is
 * flat-topped, and that B() inverted the sign test on all three axes relative
 * to paper eq. (16). Either fix alone leaves the round trip at ~14% or 0%.
 */
import { describe, expect, it } from "vitest";

import {
    gosper_decode,
    gosper_encode,
    gosper_size,
    gridify_gosper,
} from "../src/gosper.js";
import type { Point } from "../src/types.js";

const SIZE_REF: Point[] = [
    [0, 0],
    [100, 100],
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

describe("gosper_uher encode inverts decode", () => {
    for (const level of [1, 2, 3, 4, 5]) {
        it(`level ${level}: every one of ${7 ** level} indices round-trips`, () => {
            const size = gosper_size(SIZE_REF, level);
            for (let t = 0; t < 7 ** level; ++t) {
                const p = gosper_decode(t, level, size);
                expect(gosper_encode(p, level, size)).toBe(t);
            }
        });
    }

    it("is stable under a hexagon-sized jitter around each cell centre", () => {
        // A point anywhere well inside a cell must index to that cell.
        const level = 3;
        const size = gosper_size(SIZE_REF, level);
        const r = rng(2024);
        // The step between consecutive cells is constant along the curve, so
        // one measurement gives the pitch for every cell.
        const c0 = gosper_decode(0, level, size);
        const c1 = gosper_decode(1, level, size);
        const pitch = Math.hypot(c1[0] - c0[0], c1[1] - c0[1]);

        for (let t = 0; t < 7 ** level; ++t) {
            const [x, y] = gosper_decode(t, level, size);
            const angle = r() * 2 * Math.PI;
            const radius = pitch * 0.25;
            const jittered: Point = [
                x + radius * Math.cos(angle),
                y + radius * Math.sin(angle),
            ];
            expect(gosper_encode(jittered, level, size)).toBe(t);
        }
    });
});

describe("gosper_uher decode geometry", () => {
    for (const level of [1, 2, 3, 4]) {
        it(`level ${level}: ${7 ** level} distinct cells, constant step`, () => {
            const size = gosper_size(SIZE_REF, level);
            const cells: Point[] = [];
            for (let t = 0; t < 7 ** level; ++t) {
                cells.push(gosper_decode(t, level, size));
            }
            expect(new Set(cells.map((p) => p.map((v) => v.toFixed(6)).join(","))).size)
                .toBe(7 ** level);

            const first = Math.hypot(cells[1][0] - cells[0][0], cells[1][1] - cells[0][1]);
            for (let i = 1; i < cells.length; ++i) {
                const d = Math.hypot(
                    cells[i][0] - cells[i - 1][0],
                    cells[i][1] - cells[i - 1][1],
                );
                expect(d).toBeCloseTo(first, 8);
            }
        });
    }
});

describe("points outside the island are rejected, not aliased", () => {
    for (const level of [2, 3, 4]) {
        it(`level ${level}: returns -1 well outside, a real index inside`, () => {
            const size = gosper_size(SIZE_REF, level);
            const cells: Point[] = [];
            for (let t = 0; t < 7 ** level; ++t) {
                cells.push(gosper_decode(t, level, size));
            }
            const xs = cells.map((c) => c[0]);
            const ys = cells.map((c) => c[1]);
            const w = Math.max(...xs) - Math.min(...xs);
            const hgt = Math.max(...ys) - Math.min(...ys);

            // Every genuine cell centre is inside.
            for (const c of cells) {
                expect(gosper_encode(c, level, size)).toBeGreaterThanOrEqual(0);
            }

            // A ring well beyond the island is outside. Before this check the
            // index aliased onto an arbitrary interior cell, often at the far
            // end of the curve.
            for (let i = 0; i < 180; ++i) {
                const a = (i / 180) * Math.PI * 2;
                for (const mult of [1.5, 3, 10]) {
                    const p: Point = [
                        Math.cos(a) * (w / 2) * mult,
                        Math.sin(a) * (hgt / 2) * mult,
                    ];
                    expect(gosper_encode(p, level, size)).toBe(-1);
                }
            }
        });
    }
});

describe("gridify_gosper", () => {
    it("centres the data on the island so nothing aliases", () => {
        // The island is centred on the origin. Data sitting anywhere else — the
        // usual case — used to leave over half its points off the tile, each
        // silently folded onto an unrelated cell.
        const r = rng(31337);
        const level = 4;
        for (const offset of [0, 50, -400, 1e4]) {
            const data: Point[] = Array.from({ length: 300 }, () => [
                offset + r() * 100,
                offset + r() * 100,
            ]);
            const xs = data.map((d) => d[0]);
            const ys = data.map((d) => d[1]);
            const mid: Point = [
                (Math.min(...xs) + Math.max(...xs)) / 2,
                (Math.min(...ys) + Math.max(...ys)) / 2,
            ];
            const size = gosper_size(data, level);
            for (const d of data) {
                const centred: Point = [d[0] - mid[0], d[1] - mid[1]];
                expect(gosper_encode(centred, level, size)).toBeGreaterThanOrEqual(0);
            }
            expect(() => gridify_gosper(data, { level })).not.toThrow();
        }
    });

    it("keeps its output near the input's own coordinates", () => {
        const r = rng(808);
        const data: Point[] = Array.from({ length: 200 }, () => [
            1000 + r() * 100,
            -500 + r() * 100,
        ]);
        const positions = gridify_gosper(data, { level: 4 });
        const xs = positions.map((p) => p[0]);
        const ys = positions.map((p) => p[1]);
        expect(Math.min(...xs)).toBeGreaterThan(900);
        expect(Math.max(...xs)).toBeLessThan(1160);
        expect(Math.min(...ys)).toBeGreaterThan(-600);
        expect(Math.max(...ys)).toBeLessThan(-340);
    });

    it("reports a too-small level instead of hanging", () => {
        // 400 points cannot fit in 7^3 = 343 cells. The collision walk used to
        // spin forever once every cell was taken.
        const r = rng(5);
        const data: Point[] = Array.from({ length: 400 }, () => [r() * 100, r() * 100]);
        expect(() => gridify_gosper(data, { level: 3 })).toThrow(
            /more points than cells/,
        );
    });

    it("gives every point a distinct cell of the curve", () => {
        const r = rng(555);
        const data: Point[] = Array.from({ length: 400 }, () => [r() * 100, r() * 100]);
        const level = 4;
        const size = gosper_size(data, level);
        // gridify reports cells in the input's own frame, so shift the island
        // by the same amount before comparing.
        const xs = data.map((d) => d[0]);
        const ys = data.map((d) => d[1]);
        const mid_x = (Math.min(...xs) + Math.max(...xs)) / 2;
        const mid_y = (Math.min(...ys) + Math.max(...ys)) / 2;
        const cells = new Set<string>();
        for (let t = 0; t < 7 ** level; ++t) {
            const [cx, cy] = gosper_decode(t, level, size);
            cells.add([cx + mid_x, cy + mid_y].map((v) => v.toFixed(6)).join(","));
        }

        const positions = gridify_gosper(data, { level });
        expect(positions).toHaveLength(data.length);
        expect(new Set(positions.map((p) => p.join(","))).size).toBe(data.length);
        for (const p of positions) {
            expect(cells.has(p.map((v) => v.toFixed(6)).join(","))).toBe(true);
            expect(Number.isFinite(p[0]) && Number.isFinite(p[1])).toBe(true);
        }
    });

    it("never materializes the curve, so level is not a memory cost", () => {
        // A level far beyond any practical dataset would be 7^12 ≈ 1.4e10
        // cells; this must still return promptly.
        const size = gosper_size(SIZE_REF, 12);
        const t = gosper_encode([37, 21], 12, size);
        expect(Number.isInteger(t)).toBe(true);
        expect(t).toBeGreaterThanOrEqual(0);
        expect(t).toBeLessThan(7 ** 12);
        expect(gosper_encode(gosper_decode(t, 12, size), 12, size)).toBe(t);
    });
});
