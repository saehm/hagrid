/**
 * Hexagon orientation is a 30° rotation of the lattice, nothing more.
 * `"pointy"` puts neighbours at multiples of 60°, `"flat"` at 30° off those,
 * and switching must leave the curve's ordering and step length untouched.
 */
import { describe, expect, it } from "vitest";

import {
    gosper_curve,
    gosper_decode,
    gosper_encode,
    gosper_size,
    gridify_gosper,
} from "../src/gosper.js";
import type { HexOrientation, Point } from "../src/types.js";

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

/** Angle to the nearest neighbouring cell, folded into one 60° sector. */
function neighbour_angle(cells: readonly Point[]): number {
    let best = Infinity;
    let angle = 0;
    for (let i = 1; i < Math.min(cells.length, 300); ++i) {
        const dx = cells[i][0] - cells[0][0];
        const dy = cells[i][1] - cells[0][1];
        const d = Math.hypot(dx, dy);
        if (d > 1e-12 && d < best - 1e-12) {
            best = d;
            angle = Math.atan2(dy, dx);
        }
    }
    return (((angle * 180) / Math.PI + 360) % 360) % 60;
}

describe("orientation names mean what they say", () => {
    it("pointy puts neighbours at multiples of 60°", () => {
        expect(neighbour_angle(gosper_curve(4, "pointy"))).toBeCloseTo(0, 6);
    });

    it("flat offsets them by 30°", () => {
        expect(neighbour_angle(gosper_curve(4, "flat"))).toBeCloseTo(30, 6);
    });

    it("defaults to the paper's own flat-topped lattice", () => {
        expect(neighbour_angle(gosper_curve(4))).toBeCloseTo(30, 6);
    });
});

describe("orientation is a rigid rotation, not a relabelling", () => {
    for (const level of [2, 3, 4]) {
        it(`level ${level}: cells are the pointy set turned -30°`, () => {
            const pointy = gosper_curve(level, "pointy");
            const flat = gosper_curve(level, "flat");
            expect(flat).toHaveLength(pointy.length);
            const c = Math.cos(-Math.PI / 6);
            const s = Math.sin(-Math.PI / 6);
            for (let i = 0; i < flat.length; ++i) {
                expect(c * flat[i][0] - s * flat[i][1]).toBeCloseTo(pointy[i][0], 8);
                expect(s * flat[i][0] + c * flat[i][1]).toBeCloseTo(pointy[i][1], 8);
            }
        });

        it(`level ${level}: step length is unchanged`, () => {
            const pointy = gosper_curve(level, "pointy");
            const flat = gosper_curve(level, "flat");
            for (let i = 1; i < pointy.length; ++i) {
                const a = Math.hypot(
                    pointy[i][0] - pointy[i - 1][0],
                    pointy[i][1] - pointy[i - 1][1],
                );
                const b = Math.hypot(
                    flat[i][0] - flat[i - 1][0],
                    flat[i][1] - flat[i - 1][1],
                );
                expect(b).toBeCloseTo(a, 9);
            }
        });
    }
});

describe("indexing round-trips in either orientation", () => {
    for (const orientation of ["pointy", "flat"] as const) {
        for (const level of [1, 2, 3, 4, 5]) {
            it(`${orientation} level ${level}: all ${7 ** level} indices`, () => {
                const size = gosper_size(SIZE_REF, level);
                for (let t = 0; t < 7 ** level; ++t) {
                    const p = gosper_decode(t, level, size, orientation);
                    expect(gosper_encode(p, level, size, orientation)).toBe(t);
                }
            });
        }
    }
});

describe("gridify works in either orientation", () => {
    const r = rng(9090);
    const data: Point[] = Array.from({ length: 300 }, () => [r() * 100, r() * 100]);
    const level = 4;

    /** The cells gridify can produce: the island at the fitted size, recentred. */
    function placed_cells(orientation: HexOrientation): Set<string> {
        const size = gosper_size(data, level);
        const xs = data.map((d) => d[0]);
        const ys = data.map((d) => d[1]);
        const mid_x = (Math.min(...xs) + Math.max(...xs)) / 2;
        const mid_y = (Math.min(...ys) + Math.max(...ys)) / 2;
        const out = new Set<string>();
        for (let t = 0; t < 7 ** level; ++t) {
            const [x, y] = gosper_decode(t, level, size, orientation);
            out.add([x + mid_x, y + mid_y].map((v) => v.toFixed(6)).join(","));
        }
        return out;
    }

    for (const orientation of ["pointy", "flat"] as const) {
        it(`places every point on a distinct cell — ${orientation}`, () => {
            const cells = placed_cells(orientation);
            const positions = gridify_gosper(data, { level, orientation });
            expect(positions).toHaveLength(data.length);
            expect(new Set(positions.map((p) => p.join(","))).size).toBe(data.length);
            for (const p of positions) {
                expect(cells.has(p.map((v) => v.toFixed(6)).join(","))).toBe(true);
            }
        });
    }

    it("turns the lattice under the data, so cells are reassigned", () => {
        // Rotating the lattice moves the data relative to it, which is the
        // point of the option — it does not carry each point to a matching cell.
        const a = gridify_gosper(data, { level, orientation: "pointy" });
        const b = gridify_gosper(data, { level, orientation: "flat" });
        const differing = a.filter((_, i) => a[i].join(",") !== b[i].join(",")).length;
        expect(differing).toBeGreaterThan(0);
    });
});
