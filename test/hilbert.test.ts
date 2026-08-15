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
