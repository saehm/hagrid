/**
 * Builds the heartbeat showcase from ECG5000.
 *
 * The OpenML copies of this set carry a binarised target; the original five
 * classes only survive in the UCR archive, so the data comes from there.
 *
 * Each beat is 140 samples, already z-normalised and aligned on the R peak —
 * which matters for the metric. The usual reason to reach for dynamic time
 * warping is phase drift, and there is none left to absorb here, so the choice
 * between plain euclidean and DTW is measured rather than assumed.
 *
 *   node scripts/build-ecg.mjs [--sweep]
 */
import { writeFileSync, mkdirSync, existsSync, readFileSync } from "node:fs";
import { dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";
import { inflateRawSync } from "node:zlib";

import * as druid from "@saehrimnir/druidjs";
import { label_agreement } from "./representations.mjs";
import { candidates, evaluate, report } from "./projection-quality.mjs";

const here = dirname(fileURLToPath(import.meta.url));
const OUT = resolve(here, "../src/data/ecg-2d.json");
const CACHE = resolve(here, "../.cache/ecg5000.json");
const SWEEP = process.argv.includes("--sweep");

const ZIP = "https://www.timeseriesclassification.com/aeon-toolkit/ECG5000.zip";
const PER_CLASS = 60; // capped by whatever a class actually has

/** 1 is a normal beat; the rest are the ways it goes wrong. */
const CLASS_NAMES = {
    1: "Normal",
    2: "R-on-T PVC",
    3: "PVC",
    4: "Supraventricular",
    5: "Unclassifiable",
};

/* ---------------------------------------------------------------- zip ---- */

/**
 * Just enough ZIP to pull the two ARFFs out. Node ships the decompressor but no
 * container reader, and this beats a dependency for two files.
 *
 * Read via the central directory rather than by walking local headers: this
 * archive sets the data-descriptor flag, which leaves the sizes in the local
 * headers as zero and the real ones in a trailer. The central directory always
 * has them.
 */
function* zip_entries(buf) {
    const eocd = buf.lastIndexOf(Buffer.from([0x50, 0x4b, 0x05, 0x06]));
    if (eocd < 0) throw new Error("not a zip: no end-of-central-directory");
    const count = buf.readUInt16LE(eocd + 10);

    let at = buf.readUInt32LE(eocd + 16);
    for (let i = 0; i < count; ++i) {
        if (buf.readUInt32LE(at) !== 0x02014b50) throw new Error("bad central directory");
        const method = buf.readUInt16LE(at + 10);
        const compressed = buf.readUInt32LE(at + 20);
        const name_len = buf.readUInt16LE(at + 28);
        const extra_len = buf.readUInt16LE(at + 30);
        const comment_len = buf.readUInt16LE(at + 32);
        const local = buf.readUInt32LE(at + 42);
        const name = buf.toString("utf8", at + 46, at + 46 + name_len);

        yield {
            name,
            data() {
                // The local header repeats the name and extra field, and its
                // extra field is not always the same length as this one.
                const body =
                    local + 30 + buf.readUInt16LE(local + 26) + buf.readUInt16LE(local + 28);
                const raw = buf.subarray(body, body + compressed);
                return method === 0 ? raw : inflateRawSync(raw);
            },
        };
        at += 46 + name_len + extra_len + comment_len;
    }
}

/** Rows of `[...140 samples, class]` from an ARFF. */
function parse_arff(text) {
    const lines = text.split(/\r?\n/);
    const start = lines.findIndex((l) => /^@data/i.test(l)) + 1;
    return lines
        .slice(start)
        .filter((l) => l.trim())
        .map((l) => l.split(",").map(Number));
}

async function load() {
    if (existsSync(CACHE)) {
        console.log("beats from cache");
        return JSON.parse(readFileSync(CACHE, "utf8"));
    }
    console.log(`fetching ${ZIP} …`);
    const res = await fetch(ZIP);
    if (!res.ok) throw new Error(`${res.status} fetching the archive`);
    const buf = Buffer.from(await res.arrayBuffer());

    const rows = [];
    for (const entry of zip_entries(buf)) {
        if (!/ECG5000_(TRAIN|TEST)\.arff$/i.test(entry.name)) continue;
        rows.push(...parse_arff(entry.data().toString("utf8")));
    }
    if (rows.length < 4000) throw new Error(`only ${rows.length} beats found`);

    const beats = rows.map((r) => ({
        series: r.slice(0, r.length - 1),
        label: CLASS_NAMES[r[r.length - 1]] ?? String(r[r.length - 1]),
    }));
    mkdirSync(dirname(CACHE), { recursive: true });
    writeFileSync(CACHE, JSON.stringify(beats));
    return beats;
}

/* -------------------------------------------------------------- sample ---- */

const all = await load();
console.log(`${all.length} beats, ${all[0].series.length} samples each`);

const by_class = new Map();
for (const b of all) {
    if (!by_class.has(b.label)) by_class.set(b.label, []);
    by_class.get(b.label).push(b);
}
// Stratified, in file order so the choice is not another random seed to tune.
const beats = [];
for (const [name, list] of by_class) {
    const take = list.slice(0, Math.min(PER_CLASS, list.length));
    beats.push(...take);
    console.log(`  ${name.padEnd(17)} ${String(take.length).padStart(3)} of ${list.length}`);
}
/**
 * Roll every beat by half its length.
 *
 * The segmentation cuts each beat so the QRS complex lands late, which puts the
 * interesting part against the right edge of a small glyph and the flat stretch
 * where the eye lands first. Rotating by half brings the spike to the front.
 *
 * Rolling all of them by the same amount is a relabelling of the sample axis,
 * so every pairwise distance is unchanged and the projection is exactly the one
 * the unrolled data would have produced.
 */
const roll = (s) => {
    const half = Math.floor(s.length / 2);
    return [...s.slice(half), ...s.slice(0, half)];
};

const series = beats.map((b) => roll(b.series));
const labels = beats.map((b) => b.label);

/* -------------------------------------------------------------- metric ---- */

/**
 * Dynamic time warping with a Sakoe-Chiba band: the classic time-series
 * distance, which lets two beats line up even if one runs slightly early.
 */
function dtw(a, b, band = 10) {
    const n = a.length;
    const INF = Infinity;
    let prev = new Float64Array(n + 1).fill(INF);
    let curr = new Float64Array(n + 1).fill(INF);
    prev[0] = 0;
    for (let i = 1; i <= n; ++i) {
        curr.fill(INF);
        const lo = Math.max(1, i - band);
        const hi = Math.min(n, i + band);
        for (let j = lo; j <= hi; ++j) {
            const cost = (a[i - 1] - b[j - 1]) ** 2;
            curr[j] = cost + Math.min(prev[j], curr[j - 1], prev[j - 1]);
        }
        [prev, curr] = [curr, prev];
    }
    return Math.sqrt(prev[n]);
}

function dtw_embedding(rows, dims = 40) {
    const n = rows.length;
    const D = Array.from({ length: n }, () => new Float64Array(n));
    for (let i = 0; i < n; ++i) {
        for (let j = i + 1; j < n; ++j) {
            const d = dtw(rows[i], rows[j]);
            D[i][j] = d;
            D[j][i] = d;
        }
    }
    // Classical scaling turns the distances back into vectors, which is what
    // both the projections and the quality metrics need.
    const Y = new druid.MDS(druid.Matrix.from(D.map((r) => Array.from(r))), {
        d: dims,
        metric: "precomputed",
    }).transform();
    return Array.from({ length: n }, (_, i) =>
        Array.from({ length: dims }, (_, j) => Y.entry(i, j)),
    );
}

console.log("\nrepresentation      label agreement (k=10)");
const SPACES = { euclidean: series, dtw: dtw_embedding(series) };
for (const [name, rows] of Object.entries(SPACES)) {
    console.log(`${name.padEnd(18)}  ${(label_agreement(rows, labels) * 100).toFixed(1)}%`);
}

/*
 * The beats are already z-normalised and R-peak aligned, so DTW has little
 * misalignment left to absorb. Whichever wins, it wins on the numbers above.
 */
const chosen_space = Object.entries(SPACES).sort(
    (a, b) => label_agreement(b[1], labels) - label_agreement(a[1], labels),
)[0];
console.log(`\nusing: ${chosen_space[0]}`);
const high = chosen_space[1];
const Z = druid.Matrix.from(high.map((r) => Array.from(r)));

/* ------------------------------------------------------------- project ---- */

function project(method, params) {
    const Y = new druid[method](Z, params).transform();
    const pts = beats.map((_, i) => [Y.entry(i, 0), Y.entry(i, 1)]);
    return pts.every(([x, y]) => Number.isFinite(x) && Number.isFinite(y)) ? pts : null;
}

// The winner of `--sweep`, pinned so a plain run reproduces the shipped file
// rather than a near-miss. Re-run the sweep if the sample or the representation
// changes.
let chosen = { method: "TSNE", params: { seed: 1, perplexity: 20, epsilon: 5 } };
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

mkdirSync(dirname(OUT), { recursive: true });
writeFileSync(
    OUT,
    JSON.stringify({
        name: "Heartbeats",
        source: "ECG5000, UCR Time Series Archive",
        representation: chosen_space[0],
        method: chosen.method,
        params: chosen.params,
        note: "Generated by scripts/build-ecg.mjs — do not edit by hand.",
        classes: [...by_class.keys()],
        labels,
        series: series.map((s) => s.map((v) => Math.round(v * 1000) / 1000)),
        points,
    }),
);
console.log(
    `\n${beats.length} beats, ${chosen_space[0]} space, ${chosen.method} ${JSON.stringify(chosen.params)} → ${OUT}`,
);
