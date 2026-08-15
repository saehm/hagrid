/**
 * The Hilbert index must be a genuine space-filling curve on a `2^level`
 * square: a bijection onto `[0, size²)`, covering every cell, moving one step
 * at a time, and invertible.
 */
import { describe, expect, it } from "vitest";

import { HILBERT_MAX_LEVEL, hilbert_decode, hilbert_encode } from "../src/hilbert.js";

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
