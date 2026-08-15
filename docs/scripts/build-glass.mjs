/**
 * Builds the glass showcase.
 *
 * The oxides are percentages of a whole — they sum to 100 — so euclidean
 * distance on them counts the closure as signal. The Aitchison geometry is the
 * right one, but only four of the eight oxides are always present: Ba is zero
 * in 176 of 214 samples, Fe in 144. Those are absences, not small numbers, and
 * replacing them to satisfy a logarithm would invent the very structure the
 * projection then shows.
 *
 * So this compares the honest options and lets the labels decide:
 *
 *   standardised  every column z-scored, euclidean — ignores the closure
 *   aitchison     CLR over the four always-present oxides, the traces kept as
 *                 z-scored amounts, refractive index z-scored alongside
 *   clr-all       CLR over all eight with the zeros replaced — included to
 *                 show what the replacement does, not because it is sound
 *
 *   node scripts/build-glass.mjs [--sweep]
 */
import { writeFileSync, mkdirSync } from "node:fs";
import { dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";

import * as druid from "@saehrimnir/druidjs";
import { GLASS } from "@saehrimnir/mistle";

import { zscore, clr, take, concat, label_agreement } from "./representations.mjs";
import { candidates, evaluate, report } from "./projection-quality.mjs";

const here = dirname(fileURLToPath(import.meta.url));
const OUT = resolve(here, "../src/data/glass-2d.json");
const SWEEP = process.argv.includes("--sweep");

const RI = 0;
const OXIDES = [1, 2, 3, 4, 5, 6, 7, 8]; // Na Mg Al Si K Ca Ba Fe
const ALWAYS = [1, 3, 4, 6]; //            Na Al Si Ca — no zeros
const TRACE = [2, 5, 7, 8]; //             Mg K  Ba Fe — mostly zero

const values = GLASS.values;
const labels = GLASS.labels;

const REPRESENTATIONS = {
    standardised: () => zscore(values),
    aitchison: () =>
        concat(clr(take(values, ALWAYS)), zscore(take(values, TRACE)), zscore(take(values, [RI]))),
    "clr-all": () => {
        // Multiplicative replacement: zeros become a fraction of the smallest
        // positive value in their column, the rest rescaled to keep the sum.
        const parts = take(values, OXIDES);
        const floors = OXIDES.map((_, j) => {
            const positive = parts.map((r) => r[j]).filter((v) => v > 0);
            return Math.min(...positive) * 0.65;
        });
        const filled = parts.map((r) => {
            const patched = r.map((v, j) => (v > 0 ? v : floors[j]));
            const total = patched.reduce((a, b) => a + b, 0);
            return patched.map((v) => (v / total) * 100);
        });
        return concat(clr(filled), zscore(take(values, [RI])));
    },
};

console.log("representation      dims   label agreement (k=10)");
const spaces = {};
for (const [name, build] of Object.entries(REPRESENTATIONS)) {
    const rows = build();
    spaces[name] = rows;
    console.log(
        `${name.padEnd(18)}  ${String(rows[0].length).padStart(4)}   ${(label_agreement(rows, labels) * 100).toFixed(1)}%`,
    );
}

/*
 * Measured: standardised 58.6%, aitchison 59.3%, clr-all 60.5%. Two points
 * separates them, which is the actual finding — the closure is not what makes
 * this data hard, and no metric here is a breakthrough.
 *
 * `clr-all` edges it, but its margin is bought with 176 invented values for
 * absent barium, so a fraction of a percent is not enough to prefer it. The
 * choice is `aitchison`: it is the right geometry for the parts that really are
 * parts, it invents nothing, and it is within noise of the best score.
 */
const chosen_space = ["aitchison", spaces.aitchison];
console.log(`\nusing: ${chosen_space[0]} (principled; the three are within two points)`);
const high = chosen_space[1];
const Z = druid.Matrix.from(high);

function project(method, params) {
    const Y = new druid[method](Z, params).transform();
    const pts = values.map((_, i) => [Y.entry(i, 0), Y.entry(i, 1)]);
    return pts.every(([x, y]) => Number.isFinite(x) && Number.isFinite(y)) ? pts : null;
}

// The winner of `--sweep`, pinned so a plain run reproduces the shipped file
// rather than a near-miss. Re-run the sweep if the data or the representation
// changes.
let chosen = { method: "TSNE", params: { seed: 1, perplexity: 40, epsilon: 5 } };
let points = project(chosen.method, chosen.params);

if (SWEEP) {
    const scored = [];
    for (const c of candidates()) {
        const pts = project(c.method, c.params);
        if (!pts) continue;
        scored.push({
            ...c,
            points: pts,
            scores: evaluate(high, pts, { digit: labels, speaker: labels }),
        });
    }
    const best = report(scored);
    chosen = { method: best.method, params: best.params };
    points = best.points;
}

const ranges = GLASS.columns.map((_, j) => [
    Math.min(...values.map((r) => r[j])),
    Math.max(...values.map((r) => r[j])),
]);

mkdirSync(dirname(OUT), { recursive: true });
writeFileSync(
    OUT,
    JSON.stringify(
        {
            name: "Glass",
            source: "GLASS (mistle)",
            representation: chosen_space[0],
            method: chosen.method,
            params: chosen.params,
            note: "Generated by scripts/build-glass.mjs — do not edit by hand.",
            columns: GLASS.columns,
            oxides: OXIDES,
            ranges,
            labels,
            values,
            points,
        },
        null,
        1,
    ),
);
console.log(
    `\n${values.length} samples, ${chosen_space[0]} space, ${chosen.method} ${JSON.stringify(chosen.params)} → ${OUT}`,
);
