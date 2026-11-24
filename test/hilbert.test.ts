/**
 * The Hilbert index must be a genuine space-filling curve on a `2^level`
 * square: a bijection onto `[0, size²)`, covering every cell, moving one step
 * at a time, and invertible.
 */
import { describe, expect, it } from "vitest";

import {
    HILBERT_MAX_LEVEL,
    gridify_hilbert,
    hilbert_decode,
    hilbert_encode,
} from "../src/hilbert.js";
import { gridify } from "../src/gridify.js";

const LEVELS = [1, 2, 3, 4, 5];

describe("hilbert_encode and hilbert_decode are inverses", () => {
    for (const level of LEVELS) {
        const size = 2 ** level;

        it(`level ${level} round-trips every index`, () => {
            for (let n = 0; n < size ** 2; ++n) {
                const p = hilbert_decode(n, size);
                expect(p[0]).toBeGreaterThanOrEqual(0);
                expect(p[0]).toBeLessThan(size);
                expect(p[1]).toBeGreaterThanOrEqual(0);
                expect(p[1]).toBeLessThan(size);
                expect(hilbert_encode(p, size)).toBe(n);
            }
        });

        it(`level ${level} visits all ${size ** 2} cells exactly once`, () => {
            const seen = new Set<string>();
            for (let n = 0; n < size ** 2; ++n) {
                seen.add(hilbert_decode(n, size).join(","));
            }
            expect(seen.size).toBe(size ** 2);
        });

        it(`level ${level} moves one cell at a time`, () => {
            for (let n = 1; n < size ** 2; ++n) {
                const [x0, y0] = hilbert_decode(n - 1, size);
                const [x1, y1] = hilbert_decode(n, size);
                expect(Math.abs(x1 - x0) + Math.abs(y1 - y0)).toBe(1);
            }
        });
    }
});

/**
 * Every point gets its own cell whenever the grid has room, however contested
 * the preferred cells are. This is the property the collision walk exists for,
 * and it was broken for years by a units mismatch: the walk was handed the
 * curve length where it expected the grid's side, so it accepted indices past
 * the end of the curve, and those decoded back onto cells already taken.
 */
describe("gridify_hilbert seats every point when the grid has room", () => {
    /** Deterministic, so a failure is reproducible. */
    function rng(seed: number) {
        let a = seed;
        return () => {
            a |= 0;
            a = (a + 0x6d2b79f5) | 0;
            let t = Math.imul(a ^ (a >>> 15), 1 | a);
            t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
            return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
        };
    }

    const cases: [string, number, number][] = [
        // [name, points, distinct positions among them]
        ["all distinct", 800, 800],
        ["a few duplicates", 200, 195],
        ["heavy duplicates", 400, 300],
        ["mostly one point", 600, 100],
    ];

    for (const [name, n, distinct] of cases) {
        it(`${name}: ${n} points get ${n} distinct cells`, () => {
            const next = rng(1);
            const base = Array.from(
                { length: distinct },
                () => [next() * 100, next() * 100] as [number, number],
            );
            const points = Array.from({ length: n }, (_, i) =>
                i < distinct ? base[i] : base[Math.floor(next() * distinct)],
            );

            const level = Math.ceil(Math.log2(n) / Math.log2(4));
            const out = gridify_hilbert(points, { level });

            expect(out).toHaveLength(n);
            expect(new Set(out.map((p) => p.join(","))).size).toBe(n);
        });
    }

    it("indexes stay on the curve", () => {
        // A grid barely larger than the point count forces the walk to the ends,
        // which is where an off-by-one lands a point outside the lattice.
        const level = 4;
        const size = 2 ** level;
        const points = Array.from({ length: size ** 2 }, () => [0, 0] as [number, number]);
        const out = gridify_hilbert(points, { level });

        expect(new Set(out.map((p) => p.join(","))).size).toBe(size ** 2);
        for (const [x, y] of out) {
            expect(x).toBeGreaterThanOrEqual(0);
            expect(x).toBeLessThan(size);
            expect(y).toBeGreaterThanOrEqual(0);
            expect(y).toBeLessThan(size);
        }
    });

    it("throws rather than overflowing when there is no room", () => {
        const level = 2;
        const points = Array.from({ length: 4 ** level + 1 }, () => [0, 0] as [number, number]);
        expect(() => gridify_hilbert(points, { level })).toThrow(/more points than cells/);
    });
});

describe("hilbert rejects sizes it cannot index", () => {
    for (const size of [0, -4, 3, 6, 100]) {
        it(`size ${size} is not a power of two`, () => {
            expect(() => hilbert_encode([0, 0], size)).toThrow(/power of two/);
            expect(() => hilbert_decode(0, size)).toThrow(/power of two/);
        });
    }

    it(`level above ${HILBERT_MAX_LEVEL} cannot be indexed exactly`, () => {
        const size = 2 ** (HILBERT_MAX_LEVEL + 1);
        expect(() => hilbert_encode([0, 0], size)).toThrow(/indexed exactly/);
        expect(() => hilbert_decode(0, size)).toThrow(/indexed exactly/);
    });
});

describe("keep_aspect_ratio keeps a wide dataset wide", () => {
    // A 4:1 strip of points with plenty of free cells, so the collision
    // resolution has no reason to spread them over the whole square.
    const D: [number, number][] = [];
    for (let i = 0; i < 8; i++) for (let j = 0; j < 2; j++) D.push([i * 4, j]);
    const span = (P: readonly [number, number][]) =>
        Math.max(...P.map(([, y]) => y)) - Math.min(...P.map(([, y]) => y));

    it("gridify_hilbert pads the short axis instead of stretching it", () => {
        const level = 5;
        const stretched = gridify_hilbert(D, { level });
        const kept = gridify_hilbert(D, { level, keep_aspect_ratio: true });
        for (const [x, y] of kept) {
            expect(x).toBeGreaterThanOrEqual(0);
            expect(x).toBeLessThan(2 ** level);
            expect(y).toBeGreaterThanOrEqual(0);
            expect(y).toBeLessThan(2 ** level);
        }
        expect(span(kept)).toBeLessThan(span(stretched) / 2);
        // padding is split evenly, so the strip sits in the middle rows
        const mid = (2 ** level - 1) / 2;
        const mean_y = kept.reduce((s, [, y]) => s + y, 0) / kept.length;
        expect(Math.abs(mean_y - mid)).toBeLessThanOrEqual(1);
    });

    it("gridify passes the aspect ratio through to hilbert", () => {
        const off = gridify(D, "hilbert", { pluslevel: 2 }).positions;
        const on = gridify(D, "hilbert", { pluslevel: 2, keep_aspect_ratio: true }).positions;
        expect(span(on)).toBeLessThan(span(off) / 2);
    });
});

describe("hilbert_encode rejects points off the grid", () => {
    const size = 8;
    for (const p of [
        [-1, 0],
        [0, -1],
        [size, 0],
        [0, size],
        [1.5, 2],
        [NaN, 0],
    ] as [number, number][]) {
        it(`throws for [${p}]`, () => {
            expect(() => hilbert_encode(p, size)).toThrow(RangeError);
        });
    }
});

describe("gridify_hilbert handles data with no extent on an axis", () => {
    const level = 4;
    const mid = Math.round((2 ** level - 1) / 2);

    it("a horizontal line lands on the middle row", () => {
        const D: [number, number][] = [0, 1, 2, 3].map((x) => [x, 5]);
        const Y = gridify_hilbert(D, { level });
        expect(Y.map(([x]) => x)).toEqual([0, 5, 10, 15]);
        for (const [, y] of Y) expect(y).toBe(mid);
    });

    it("a vertical line lands on the middle column", () => {
        const D: [number, number][] = [0, 1, 2, 3].map((y) => [5, y]);
        const Y = gridify_hilbert(D, { level });
        expect(Y.map(([, y]) => y)).toEqual([0, 5, 10, 15]);
        for (const [x] of Y) expect(x).toBe(mid);
    });

    it("a single point lands in the centre", () => {
        for (const keep_aspect_ratio of [false, true]) {
            expect(gridify_hilbert([[3, 3]], { level, keep_aspect_ratio })).toEqual([[mid, mid]]);
        }
    });
});
