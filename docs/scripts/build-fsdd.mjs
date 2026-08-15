/**
 * Builds the spoken-digit showcase from the Free Spoken Digit Dataset.
 *
 * Downloads a sample of recordings, computes a log-mel spectrogram for each,
 * projects those with PaCMAP, and writes the result. The audio itself is *not*
 * copied here — the page streams each clip from the dataset's own repository
 * when you click it, so nothing is bundled and nothing loads up front.
 *
 *   node scripts/build-fsdd.mjs                 # defaults below
 *   node scripts/build-fsdd.mjs 20 PaCMAP 4     # per speaker+digit, method, seed
 *
 * FSDD is CC BY-SA 4.0 — Jakobovski/free-spoken-digit-dataset.
 */
import { writeFileSync, mkdirSync, existsSync, readFileSync } from "node:fs";
import { dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";
import { deflateSync } from "node:zlib";

import * as druid from "@saehrimnir/druidjs";

const here = dirname(fileURLToPath(import.meta.url));
const OUT = resolve(here, "../src/data/fsdd-2d.json");

const RAW = "https://raw.githubusercontent.com/Jakobovski/free-spoken-digit-dataset/master/recordings";
const SPEAKERS = ["george", "jackson", "lucas", "nicolas"];
const DIGITS = [0, 2, 4, 6, 8];

const args = process.argv.slice(2).filter((a) => a !== "--sweep");
const SWEEP = process.argv.includes("--sweep");

const PER = Number(args[0] ?? 15); // recordings per speaker and digit
// Defaults are the winner of `--sweep`, so a plain run reproduces the shipped
// file rather than a near-miss. Re-run the sweep if the sample changes.
const METHOD = args[1] ?? "TSNE";
const SEED = Number(args[2] ?? 2);
const PERPLEXITY = Number(args[3] ?? 40);

// Spectrogram shape. Ships in the page's JSON as a PNG per clip, so it wants
// to stay modest — but 16x16 was mush at the sizes the grid affords.
const BANDS = 36;
const FRAMES = 36;
const FFT = 256; // 32 ms at 8 kHz
const HOP = 48; // enough columns to fill 36 frames without inventing any

/* ---------------------------------------------------------------- wav ---- */

/** FSDD is 8 kHz 16-bit mono PCM; this reads just enough RIFF for that. */
function decode_wav(buffer) {
    const view = new DataView(buffer);
    if (view.getUint32(0, false) !== 0x52494646) throw new Error("not RIFF");

    let offset = 12;
    let rate = 8000;
    let samples = null;
    while (offset + 8 <= view.byteLength) {
        const id = String.fromCharCode(
            view.getUint8(offset), view.getUint8(offset + 1),
            view.getUint8(offset + 2), view.getUint8(offset + 3),
        );
        const size = view.getUint32(offset + 4, true);
        const body = offset + 8;
        if (id === "fmt ") {
            rate = view.getUint32(body + 4, true);
        } else if (id === "data") {
            const count = Math.floor(size / 2);
            samples = new Float32Array(count);
            for (let i = 0; i < count; ++i) samples[i] = view.getInt16(body + 2 * i, true) / 32768;
        }
        offset = body + size + (size % 2);
    }
    if (!samples) throw new Error("no data chunk");
    return { samples, rate };
}

/* ---------------------------------------------------------------- fft ---- */

/** In-place iterative radix-2 FFT. `re`/`im` are power-of-two length. */
function fft(re, im) {
    const n = re.length;
    for (let i = 1, j = 0; i < n; ++i) {
        let bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) {
            [re[i], re[j]] = [re[j], re[i]];
            [im[i], im[j]] = [im[j], im[i]];
        }
    }
    for (let len = 2; len <= n; len <<= 1) {
        const ang = (-2 * Math.PI) / len;
        const wr = Math.cos(ang), wi = Math.sin(ang);
        for (let i = 0; i < n; i += len) {
            let cr = 1, ci = 0;
            for (let k = 0; k < len / 2; ++k) {
                const ar = re[i + k], ai = im[i + k];
                const br = re[i + k + len / 2] * cr - im[i + k + len / 2] * ci;
                const bi = re[i + k + len / 2] * ci + im[i + k + len / 2] * cr;
                re[i + k] = ar + br; im[i + k] = ai + bi;
                re[i + k + len / 2] = ar - br; im[i + k + len / 2] = ai - bi;
                const ncr = cr * wr - ci * wi;
                ci = cr * wi + ci * wr;
                cr = ncr;
            }
        }
    }
}

const hann = Array.from({ length: FFT }, (_, i) => 0.5 - 0.5 * Math.cos((2 * Math.PI * i) / (FFT - 1)));

/** Mel-ish band edges over the FFT bins, so low frequencies get more detail. */
const mel = (f) => 2595 * Math.log10(1 + f / 700);
const unmel = (m) => 700 * (10 ** (m / 2595) - 1);
function band_edges(rate) {
    const lo = mel(80), hi = mel(rate / 2);
    return Array.from({ length: BANDS + 2 }, (_, i) =>
        Math.round((unmel(lo + ((hi - lo) * i) / (BANDS + 1)) / (rate / 2)) * (FFT / 2)),
    );
}

/**
 * A log-mel spectrogram, resampled to a fixed `BANDS × FRAMES` grid so every
 * clip is comparable regardless of how long it was.
 *
 * Returned row-major, band 0 (lowest) first, each value 0..1.
 */
function spectrogram({ samples, rate }) {
    // Trim leading and trailing silence — the recordings have plenty, and it
    // would otherwise dominate the time axis.
    let peak = 0;
    for (const s of samples) peak = Math.max(peak, Math.abs(s));
    const gate = peak * 0.04;
    let a = 0, b = samples.length - 1;
    while (a < b && Math.abs(samples[a]) < gate) ++a;
    while (b > a && Math.abs(samples[b]) < gate) --b;
    const voiced = samples.subarray(Math.max(0, a - HOP), Math.min(samples.length, b + HOP));

    const edges = band_edges(rate);
    const columns = [];
    for (let start = 0; start + FFT <= voiced.length; start += HOP) {
        const re = new Float64Array(FFT);
        const im = new Float64Array(FFT);
        for (let i = 0; i < FFT; ++i) re[i] = voiced[start + i] * hann[i];
        fft(re, im);

        const col = new Float64Array(BANDS);
        for (let bnd = 0; bnd < BANDS; ++bnd) {
            let sum = 0;
            for (let k = edges[bnd]; k < Math.max(edges[bnd] + 1, edges[bnd + 2]); ++k) {
                sum += re[k] * re[k] + im[k] * im[k];
            }
            col[bnd] = Math.log10(1 + sum * 1e3);
        }
        columns.push(col);
    }
    if (columns.length === 0) return new Array(BANDS * FRAMES).fill(0);

    // Resample the time axis to FRAMES, then normalise the whole clip.
    const grid = new Float64Array(BANDS * FRAMES);
    let hi = 0;
    for (let t = 0; t < FRAMES; ++t) {
        const at = (t / Math.max(1, FRAMES - 1)) * (columns.length - 1);
        const j = Math.min(columns.length - 2, Math.floor(at));
        const frac = columns.length > 1 ? at - j : 0;
        for (let bnd = 0; bnd < BANDS; ++bnd) {
            const v = columns.length > 1
                ? columns[j][bnd] * (1 - frac) + columns[j + 1][bnd] * frac
                : columns[0][bnd];
            grid[bnd * FRAMES + t] = v;
            hi = Math.max(hi, v);
        }
    }
    return Array.from(grid, (v) => (hi > 0 ? Math.round((v / hi) * 255) / 255 : 0));
}

/* ---------------------------------------------------------------- png ---- */

/**
 * A spectrogram is 16x16 = 256 cells, and 300 of them would be 76,800 SVG
 * rectangles. Encoding each as a tiny greyscale PNG instead makes every glyph
 * one `<image>`, which the browser rasterises without any DOM cost — and it
 * still renders server-side, which a canvas would not.
 *
 * Colour is left to the page: the spectrogram stays grey and the class is shown
 * as a frame around it, so the same image serves whichever labelling is chosen.
 */
const CRC_TABLE = Array.from({ length: 256 }, (_, n) => {
    let c = n;
    for (let k = 0; k < 8; ++k) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1;
    return c >>> 0;
});
function crc32(buf) {
    let c = 0xffffffff;
    for (const byte of buf) c = CRC_TABLE[(c ^ byte) & 0xff] ^ (c >>> 8);
    return (c ^ 0xffffffff) >>> 0;
}
function chunk(type, data) {
    const head = Buffer.alloc(8);
    head.writeUInt32BE(data.length, 0);
    head.write(type, 4, "ascii");
    const crc = Buffer.alloc(4);
    crc.writeUInt32BE(crc32(Buffer.concat([head.subarray(4), data])), 0);
    return Buffer.concat([head, data, crc]);
}
function png_grey(width, height, pixels) {
    const raw = Buffer.alloc((width + 1) * height);
    for (let y = 0; y < height; ++y) {
        raw[y * (width + 1)] = 0; // filter: none
        for (let x = 0; x < width; ++x) raw[y * (width + 1) + 1 + x] = pixels[y * width + x];
    }
    const ihdr = Buffer.alloc(13);
    ihdr.writeUInt32BE(width, 0);
    ihdr.writeUInt32BE(height, 4);
    ihdr[8] = 8; // bit depth
    ihdr[9] = 0; // greyscale
    return Buffer.concat([
        Buffer.from([0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a]),
        chunk("IHDR", ihdr),
        chunk("IDAT", deflateSync(raw, { level: 9 })),
        chunk("IEND", Buffer.alloc(0)),
    ]);
}

/** Row-major `[band][frame]` to an image with low frequencies at the bottom. */
function spectrogram_png(grid) {
    const pixels = new Uint8Array(BANDS * FRAMES);
    for (let b = 0; b < BANDS; ++b) {
        for (let t = 0; t < FRAMES; ++t) {
            pixels[(BANDS - 1 - b) * FRAMES + t] = Math.round(grid[b * FRAMES + t] * 255);
        }
    }
    return `data:image/png;base64,${png_grey(FRAMES, BANDS, pixels).toString("base64")}`;
}

/* --------------------------------------------------------------- build ---- */

async function grab(name) {
    const res = await fetch(`${RAW}/${name}`);
    if (!res.ok) throw new Error(`${res.status} ${name}`);
    return decode_wav(await res.arrayBuffer());
}

const wanted = [];
for (const speaker of SPEAKERS) {
    for (const digit of DIGITS) {
        for (let i = 0; i < PER; ++i) wanted.push({ file: `${digit}_${speaker}_${i}.wav`, digit, speaker });
    }
}

// Decoding is deterministic and the download is the slow part, so cache it —
// sweeping projections should not re-fetch 300 files each time.
const CACHE = resolve(here, `../.cache/fsdd-${PER}-${BANDS}x${FRAMES}.json`);
let items = [];
if (existsSync(CACHE)) {
    items = JSON.parse(readFileSync(CACHE, "utf8"));
    console.log(`${items.length} recordings from cache`);
} else {
    console.log(`fetching ${wanted.length} recordings…`);
    const CONCURRENCY = 12;
    for (let i = 0; i < wanted.length; i += CONCURRENCY) {
        const batch = wanted.slice(i, i + CONCURRENCY);
        const done = await Promise.all(
            batch.map(async (w) => {
                try {
                    return { ...w, spectrogram: spectrogram(await grab(w.file)) };
                } catch (e) {
                    console.warn(`  skipped ${w.file}: ${e.message}`);
                    return null;
                }
            }),
        );
        items.push(...done.filter(Boolean));
        process.stdout.write(`  ${items.length}/${wanted.length}\r`);
    }
    console.log(`\n${items.length} recordings decoded`);
    mkdirSync(dirname(CACHE), { recursive: true });
    writeFileSync(CACHE, JSON.stringify(items));
}
if (items.length < 10) throw new Error("too few recordings; is the network up?");

// Project the spectrograms themselves — standardised per feature, as elsewhere.
const X = items.map((it) => it.spectrogram);
const d = X[0].length;
const mean = Array.from({ length: d }, (_, j) => X.reduce((s, r) => s + r[j], 0) / X.length);
const std = Array.from({ length: d }, (_, j) =>
    Math.sqrt(X.reduce((s, r) => s + (r[j] - mean[j]) ** 2, 0) / X.length) || 1,
);
const standardised = X.map((r) => r.map((v, j) => (v - mean[j]) / std[j]));

// A 36x36 spectrogram is 1296 numbers, most of them redundant. PCA down to 50
// first: it is what these methods expect to be given, it costs almost nothing
// in fidelity here, and it turns a sweep from unusable into a minute's work.
// The metrics still score against the full-resolution features.
const PCA_DIMS = 50;
const P = new druid.PCA(druid.Matrix.from(standardised), { d: PCA_DIMS }).transform();
const Z = druid.Matrix.from(
    Array.from({ length: items.length }, (_, i) =>
        Array.from({ length: PCA_DIMS }, (_, j) => P.entry(i, j)),
    ),
);

function project(method, params) {
    if (!(method in druid)) throw new Error(`unknown projection ${method}`);
    const Y = new druid[method](Z, params).transform();
    const pts = items.map((_, i) => [Y.entry(i, 0), Y.entry(i, 1)]);
    if (!pts.every(([x, y]) => Number.isFinite(x) && Number.isFinite(y))) return null;
    return pts;
}

// t-SNE takes perplexity and epsilon; the neighbour-graph methods do not, so
// they are handed the seed alone rather than parameters they would ignore.
const DEFAULTS =
    METHOD === "TSNE"
        ? { seed: SEED, perplexity: PERPLEXITY, epsilon: 5 }
        : { seed: SEED };

let chosen = { method: METHOD, params: DEFAULTS };
let points = project(METHOD, DEFAULTS);

if (SWEEP) {
    const { candidates, evaluate, report } = await import("./projection-quality.mjs");
    const scored = [];
    for (const c of candidates()) {
        const pts = project(c.method, c.params);
        if (!pts) {
            console.warn(`  ${c.method} ${JSON.stringify(c.params)} → non-finite, skipped`);
            continue;
        }
        scored.push({ ...c, points: pts, scores: evaluate(standardised, pts, labels_for(items)) });
    }
    const best = report(scored);
    chosen = { method: best.method, params: best.params };
    points = best.points;
}

function labels_for(list) {
    return {
        digit: list.map((it) => String(it.digit)),
        speaker: list.map((it) => it.speaker),
    };
}
if (!points) throw new Error(`${METHOD} produced non-finite coordinates`);

mkdirSync(dirname(OUT), { recursive: true });
writeFileSync(
    OUT,
    JSON.stringify({
        name: "Spoken digits",
        source: "https://github.com/Jakobovski/free-spoken-digit-dataset",
        licence: "CC BY-SA 4.0",
        audio_base: RAW,
        method: chosen.method,
        params: chosen.params,
        note: "Generated by scripts/build-fsdd.mjs — do not edit by hand.",
        bands: BANDS,
        frames: FRAMES,
        digits: DIGITS,
        speakers: SPEAKERS,
        files: items.map((it) => it.file),
        digit: items.map((it) => it.digit),
        speaker: items.map((it) => it.speaker),
        images: items.map((it) => spectrogram_png(it.spectrogram)),
        points,
    }),
);
console.log(
    `${items.length} clips, ${BANDS}x${FRAMES} spectrograms, ` +
        `${chosen.method} ${JSON.stringify(chosen.params)} → ${OUT}`,
);
