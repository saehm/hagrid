/**
 * Behavioural tests for the parts that had no working baseline to capture:
 * cmds, gosper and gridfit were all broken in v0.1.x, so these assert
 * properties rather than recorded outputs.
 */
import { describe, expect, it } from "vitest";

import { gridify } from "../src/gridify.js";
import { gosper_curve, gridify_gosper } from "../src/gosper.js";
import { gridify_gridfit } from "../src/gridfit.js";
import { version } from "../src/main.js";
import type { CMDSStep } from "../src/cmds.js";
import type { Point } from "../src/types.js";

function ring(n: number, seed = 1): Point[] {
    return Array.from({ length: n }, (_, i) => {
        const t = (i * 2.399963 + seed) % (2 * Math.PI);
        const r = 20 + (i % 7) * 3;
        return [50 + r * Math.cos(t), 50 + r * Math.sin(t)] as Point;
    });
}

const unique = (ps: readonly Point[]) => new Set(ps.map((p) => p.join(","))).size;

describe("gridify contract", () => {
    it("reports a version", () => {
        expect(version).toMatch(/^\d+\.\d+\.\d+/);
    });

    it("returns positions and a runtime", () => {
        const result = gridify(ring(40), "hilbert");
        expect(result.positions).toHaveLength(40);
        expect(result.runtime).toBeGreaterThanOrEqual(0);
    });

    it("rejects an unknown method", () => {
        // @ts-expect-error — exercising the runtime guard for untyped callers.
        expect(() => gridify(ring(10), "nope")).toThrow(/unknown method/);
    });

    it("does not mutate the caller's data", () => {
        const data = ring(30);
        const snapshot = JSON.stringify(data);
        for (const method of ["hilbert", "dgrid", "nmap"] as const) {
            gridify(data, method);
        }
        expect(JSON.stringify(data)).toBe(snapshot);
    });
});

describe("hilbert and dgrid remove every overlap", () => {
    for (const n of [16, 64, 257]) {
        it(`assigns ${n} points to distinct cells`, () => {
            for (const method of ["hilbert", "dgrid"] as const) {
                const { positions } = gridify(ring(n), method);
                expect(positions).toHaveLength(n);
                expect(unique(positions)).toBe(n);
            }
        });
    }
});

describe("gosper", () => {
    it("builds a curve of 7^level cells", () => {
        for (const level of [1, 2, 3]) {
            expect(gosper_curve(level)).toHaveLength(7 ** level);
        }
    });

    it("steps a constant distance between consecutive cells", () => {
        const V = gosper_curve(3);
        const first = Math.hypot(V[1][0] - V[0][0], V[1][1] - V[0][1]);
        for (let i = 1; i < V.length; ++i) {
            const d = Math.hypot(V[i][0] - V[i - 1][0], V[i][1] - V[i - 1][1]);
            expect(d).toBeCloseTo(first, 10);
        }
    });

    it("assigns distinct cells and defaults pluslevel", () => {
        // Previously `gridify(data, "gosper")` read an absent pluslevel and
        // produced a NaN level, so this call could not work at all.
        const { positions } = gridify(ring(50), "gosper");
        expect(positions).toHaveLength(50);
        expect(unique(positions)).toBe(50);
    });

    it("places points on actual curve cells", () => {
        const cells = new Set(gosper_curve(4).map((p) => p.join(",")));
        const positions = gridify_gosper(ring(60), { level: 4 });
        for (const p of positions) {
            expect(cells.has(p.join(","))).toBe(true);
        }
    });
});

describe("cmds", () => {
    it("runs to a grid-snapped final layout without throwing", async () => {
        // v0.1.x always threw here: the final proximity graph was built with
        // the wrong arity, so the snapped layout never reached a consumer.
        const data = ring(40);
        const { steps } = gridify(data, "cmds", { max_iter: 5 });
        const yielded: CMDSStep[] = [];
        for await (const step of steps) {
            yielded.push(step);
        }
        expect(yielded.length).toBeGreaterThan(1);

        const [final] = yielded[yielded.length - 1];
        expect(final).toHaveLength(data.length);
        for (const [x, y] of final) {
            expect(Number.isFinite(x)).toBe(true);
            expect(Number.isFinite(y)).toBe(true);
        }
    });

    it("snaps the final layout to a regular lattice", async () => {
        const size = 8;
        const { steps } = gridify(ring(30), "cmds", { size, max_iter: 3 });
        let final: Point[] = [];
        for await (const [G] of steps) final = G;

        // Every coordinate must sit at a half-offset multiple of the pitch.
        for (const [x, y] of final) {
            expect(((x - size / 2) / size) % 1).toBeCloseTo(0, 9);
            expect(((y - size / 2) / size) % 1).toBeCloseTo(0, 9);
        }
    });
});

describe("gridfit", () => {
    it("runs instead of throwing on an undefined d3", () => {
        // v0.1.x referenced an unimported `d3`, so this threw a ReferenceError
        // before doing any work.
        const data: Point[] = ring(40).map(([x, y]) => [Math.round(x), Math.round(y)]);
        const result = gridify_gridfit(data);
        expect(result.nodes.length).toBeGreaterThan(0);
        expect(result.final.length + result.leftovers.length).toBeGreaterThan(0);
    });

    it("never places two points in the same cell", () => {
        const data: Point[] = ring(64, 2).map(([x, y]) => [Math.round(x), Math.round(y)]);
        const { final } = gridify_gridfit(data);
        const cells = final.map(({ cell }) => cell.join(","));
        expect(new Set(cells).size).toBe(cells.length);
    });
});
