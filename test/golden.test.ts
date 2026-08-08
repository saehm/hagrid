/**
 * Verifies the TypeScript port against outputs captured from the pre-migration
 * JavaScript bundle (v0.1.x). Any difference here is a regression in the port,
 * not a change in taste — deliberate behaviour changes must re-record the
 * fixtures and say why in the commit.
 *
 * `hilbert` and `dgrid` are still pinned to the original v0.1.1 output. The
 * `nmap` entries were re-recorded when `sizes()` was fixed: it had been giving
 * the first half's share to both halves, leaving up to 38% of the plotting
 * area as gaps. Those entries carry a `rerecorded` note, and the tiling
 * property that replaced them is checked directly in `nmap.test.ts`.
 */
import { describe, expect, it } from "vitest";

import { gridify, type GridifyMethod } from "../src/gridify.js";
import { hilbert_decode, hilbert_encode } from "../src/hilbert.js";
import * as utils from "../src/utils.js";
import type { Point } from "../src/types.js";

import datasets from "./__goldens__/datasets.json" with { type: "json" };
import gridify_golden from "./__goldens__/gridify.json" with { type: "json" };
import hilbert_codec from "./__goldens__/hilbert-codec.json" with { type: "json" };
import utils_golden from "./__goldens__/utils.json" with { type: "json" };

const DATA = datasets as Record<string, Point[]>;

describe("gridify matches the v0.1.x bundle", () => {
    const entries = Object.entries(
        gridify_golden as Record<string, { ok: boolean; value: Point[] }>,
    );

    it("covers every fixture", () => {
        expect(entries.length).toBe(24);
    });

    for (const [key, golden] of entries) {
        const [dataset, method, params] = key.split("::");
        it(`${method} ${params} on ${dataset}`, () => {
            const { positions } = gridify(
                DATA[dataset],
                method as Exclude<GridifyMethod, "cmds">,
                JSON.parse(params),
            );
            expect(positions).toEqual(golden.value);
        });
    }
});

describe("hilbert index round-trips", () => {
    for (const { level, size, points } of hilbert_codec as {
        level: number;
        size: number;
        points: { n: number; p: Point; roundtrip: number }[];
    }[]) {
        it(`level ${level} decodes and re-encodes every cell`, () => {
            for (const { n, p, roundtrip } of points) {
                expect(hilbert_decode(n, size)).toEqual(p);
                expect(hilbert_encode(p as Point, size)).toBe(roundtrip);
            }
        });

        it(`level ${level} visits all ${size ** 2} cells exactly once`, () => {
            const seen = new Set(
                points.map(({ n }) => hilbert_decode(n, size).join(",")),
            );
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

describe("utils match the v0.1.x bundle", () => {
    const g = utils_golden as {
        get_bounds: Record<string, utils.Bounds>;
        remap: Point[];
        extent: Record<string, { x: Point; y: Point }>;
        distance: number[];
    };

    it("get_bounds", () => {
        for (const [name, expected] of Object.entries(g.get_bounds)) {
            expect(utils.get_bounds(DATA[name])).toEqual(expected);
        }
    });

    it("remap", () => {
        expect(utils.remap([5, 5], [0, 0, 10, 10], [0, 0, 100, 100])).toEqual(g.remap[0]);
        expect(utils.remap([0, 10], [0, 0, 10, 10], [5, 5, 20, 20])).toEqual(g.remap[1]);
        expect(utils.remap([-3, 7], [-10, -10, 20, 20], [0, 0, 1, 1])).toEqual(g.remap[2]);
    });

    it("extent", () => {
        for (const [name, expected] of Object.entries(g.extent)) {
            expect(utils.extent(DATA[name], (p) => p[0])).toEqual(expected.x);
            expect(utils.extent(DATA[name], (p) => p[1])).toEqual(expected.y);
        }
    });

    it("distance", () => {
        expect(utils.distance([0, 0], [3, 4])).toBe(g.distance[0]);
        expect(utils.distance([-1, -1], [2, 3])).toBe(g.distance[1]);
    });
});
