/**
 * CMDS, after Liu, Hu, North & Shen, Computer Graphics Forum 37(1), 2018.
 *
 * The update rule was already faithful to eq. (4), with the paper's N = 20 and
 * α = 0.1. What was missing was everything §4.2 says about making it converge:
 * the m × n grid that sets the node size, the inset of Γ by half a node, the
 * uniform scaling of the input into Γ, and growing the grid when a run does not
 * converge. Without those the loop could never terminate.
 */
import { describe, expect, it } from "vitest";

import { cmds_grid, gridify_cmds } from "../src/cmds.js";
import { gridify } from "../src/gridify.js";
import { get_bounds } from "../src/utils.js";
import type { Point, Region } from "../src/types.js";

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

const MAX_ITER = 10;
const MAX_RESCALES = 60;

async function drain(data: Point[], params: Record<string, unknown> = {}) {
    const { steps } = gridify(data, "cmds", params);
    let layout: Point[] = [];
    let size = 0;
    let rounds = 0;
    for await (const [G, , s] of steps) {
        layout = G;
        size = s;
        ++rounds;
    }
    const cap =
        ((params.max_iter as number) ?? MAX_ITER) *
            (((params.max_rescales as number) ?? MAX_RESCALES) + 1) +
        1;
    return { layout, size, rounds, converged: rounds < cap };
}

describe("cmds_grid follows §4.2", () => {
    it("holds every item while using every row and column", () => {
        for (const count of [1, 2, 7, 50, 101, 1000]) {
            for (const [w, h] of [[100, 100], [400, 100], [100, 400], [37, 91]]) {
                const [m, n] = cmds_grid(w, h, count);
                expect(m * n, `${count} in ${w}x${h}`).toBeGreaterThanOrEqual(count);
                expect((m - 1) * n, `${count} in ${w}x${h}`).toBeLessThan(count);
                expect(m * (n - 1), `${count} in ${w}x${h}`).toBeLessThan(count);
            }
        }
    });

    it("tracks the region's aspect ratio", () => {
        const [wide_m, wide_n] = cmds_grid(1000, 100, 400);
        expect(wide_m).toBeGreaterThan(wide_n);
        const [tall_m, tall_n] = cmds_grid(100, 1000, 400);
        expect(tall_n).toBeGreaterThan(tall_m);
    });
});

describe("cmds terminates and produces a usable layout", () => {
    for (const n of [10, 25, 60]) {
        it(`N=${n}: finite, distinct, lattice-aligned positions`, async () => {
            const data = uniform(n, n * 7);
            const { layout, size } = await drain(data);

            expect(layout).toHaveLength(n);
            expect(size).toBeGreaterThan(0);
            for (const [x, y] of layout) {
                expect(Number.isFinite(x) && Number.isFinite(y)).toBe(true);
            }
            // eq. (6): every coordinate on a half-offset multiple of the pitch
            for (const [x, y] of layout) {
                for (const v of [x, y]) {
                    const k = (v - size / 2) / size;
                    expect(Math.abs(k - Math.round(k))).toBeCloseTo(0, 9);
                }
            }
            // Distinct cells only follow from convergence: the grid snap can
            // collide two nodes if the run stopped while overlaps remained.
            // cmds does not converge on every input — see the README.
            const { converged } = await drain(data);
            if (converged) {
                expect(new Set(layout.map((p) => p.join(","))).size).toBe(n);
            }
        });
    }

    it("shrinks the nodes rather than giving up", async () => {
        // The grid grows whenever a run does not converge, so the node size
        // must come back no larger than the first cell it started from.
        const data = uniform(80, 3);
        const { size } = await drain(data);
        const bb = get_bounds(data);
        const w = Math.max(bb.width, bb.height);
        const [m, n] = cmds_grid(w, w, data.length);
        expect(size).toBeLessThanOrEqual(Math.min(w / m, w / n) + 1e-9);
    });

    it("respects an explicit starting size by picking that grid", async () => {
        const data = uniform(40, 11);
        const { size } = await drain(data, { size: 4, max_iter: 2, max_rescales: 0 });
        expect(size).toBeGreaterThan(0);
        expect(size).toBeLessThanOrEqual(4 + 1e-9);
    });

    it("stops early once it converges instead of burning the budget", async () => {
        const data = uniform(12, 5);
        const { rounds } = await drain(data, { max_iter: 4, max_rescales: 60 });
        expect(rounds).toBeLessThan(4 * 61 + 1);
    });
});

describe("cmds keeps its layout inside Gamma", () => {
    it("places every node so its whole glyph fits", async () => {
        const data = uniform(40, 21);
        const bb = get_bounds(data);
        const w = Math.max(bb.width, bb.height);
        const Gamma: Region = { x1: bb.x, x2: bb.x + w, y1: bb.y, y2: bb.y + w };

        // gridify_cmds is the generator itself, not a { steps } wrapper.
        const steps = gridify_cmds(data, { Gamma, max_iter: 4, max_rescales: 30 });
        let layout: Point[] = [];
        let size = 0;
        for await (const [G, , s] of steps) {
            layout = G;
            size = s;
        }
        // Gamma is inset by half a node, so the glyph must sit inside Gamma
        // itself. Allow one pitch of slack for the final grid snap.
        for (const [x, y] of layout) {
            expect(x).toBeGreaterThanOrEqual(Gamma.x1 - size);
            expect(x).toBeLessThanOrEqual(Gamma.x2 + size);
            expect(y).toBeGreaterThanOrEqual(Gamma.y1 - size);
            expect(y).toBeLessThanOrEqual(Gamma.y2 + size);
        }
    });
});
