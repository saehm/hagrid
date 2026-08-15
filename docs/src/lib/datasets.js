import { gaussian, mulberry32 } from "./rng.js";

/**
 * Synthetic datasets for the demos.
 *
 * These are not stand-ins for real projections — they are chosen so that each
 * one provokes a documented behaviour of some method: the elongated set shows
 * why `gilbert` exists, the lattice set is what makes `nmap` throw, and the
 * duplicates set is what `cmds` cannot lay out at all. A demo you cannot make
 * fail teaches nothing about when to reach for the alternative.
 *
 * Every generator takes the seeded RNG, so the same seed is the same figure on
 * every reload and in every browser.
 */

const make = (n, fn) => Array.from({ length: n }, fn);

/** @type {Record<string, {label: string, note: string, make: (n: number, r: () => number) => [number, number][]}>} */
export const datasets = {
    uniform: {
        label: "Uniform random",
        note: "Uniform random in a square. No structure to preserve.",
        make: (n, r) => make(n, () => [r() * 100, r() * 100]),
    },
    clusters: {
        label: "Gaussian clusters",
        note: "Six gaussian blobs. The ordinary case.",
        make: (n, r) => {
            const k = 6;
            const c = make(k, () => [12 + r() * 76, 12 + r() * 76]);
            return make(n, (_, i) => [
                c[i % k][0] + gaussian(r) * 6,
                c[i % k][1] + gaussian(r) * 6,
            ]);
        },
    },
    ring: {
        label: "Ring",
        note: "A hollow annulus. The empty middle is what the grid has to decide about.",
        make: (n, r) =>
            make(n, () => {
                const a = r() * Math.PI * 2;
                const rad = 32 + gaussian(r) * 3;
                return [50 + rad * Math.cos(a), 50 + rad * Math.sin(a)];
            }),
    },
    moons: {
        label: "Two moons",
        note: "Two interleaved arcs. Watch whether the two stay separable after gridding.",
        make: (n, r) =>
            make(n, (_, i) => {
                const t = r() * Math.PI;
                const x = Math.cos(t) * 30;
                const y = Math.sin(t) * 30;
                return i % 2 === 0
                    ? [50 + x - 8 + gaussian(r) * 2.5, 42 + y + gaussian(r) * 2.5]
                    : [50 - x + 8 + gaussian(r) * 2.5, 62 - y + gaussian(r) * 2.5];
            }),
    },
    spiral: {
        label: "Spiral",
        note: "Three turns. Colour makes it obvious when the winding order survives.",
        make: (n, r) =>
            make(n, (_, i) => {
                const t = (i / n) * Math.PI * 6;
                const rad = 4 + (i / n) * 40;
                return [
                    50 + rad * Math.cos(t) + gaussian(r) * 1.2,
                    50 + rad * Math.sin(t) + gaussian(r) * 1.2,
                ];
            }),
    },
    lattice: {
        label: "Integer lattice",
        note: "All points distinct, all on a regular grid. Quantised data — nmap refuses it.",
        make: (n) => {
            const side = Math.ceil(Math.sqrt(n));
            return make(n, (_, i) => [
                (i % side) * (100 / side),
                Math.floor(i / side) * (100 / side),
            ]);
        },
    },
    elongated: {
        label: "Elongated (100:1)",
        note: "A hundred times wider than tall. Punishes methods that insist on a square.",
        make: (n, r) => make(n, () => [r() * 100, 49.5 + r() * 1]),
    },
    duplicates: {
        label: "Heavy duplicates",
        note: "Snapped to a coarse grid, so most positions repeat many times over.",
        make: (n, r) => {
            const k = Math.max(4, Math.floor(Math.sqrt(n) / 2));
            return make(n, () => [
                Math.floor(r() * k) * (100 / k),
                Math.floor(r() * k) * (100 / k),
            ]);
        },
    },
    line: {
        label: "Single column",
        note: "Every point shares an x. Degenerate on one axis outright.",
        make: (n, r) => make(n, () => [50, r() * 100]),
    },
};

export const dataset_keys = Object.keys(datasets);

/** Builds a dataset by key, from a seed. */
export function build(key, n, seed) {
    // The seed is spread through a large odd multiplier so that adjacent slider
    // positions give unrelated point sets rather than near-identical ones.
    return datasets[key].make(n, mulberry32(seed * 7919));
}

// The point encodings live in ./colour.js.
