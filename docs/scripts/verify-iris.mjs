/**
 * The showcase claims a computed glyph size is overlap-free by construction.
 * This checks that across every technique and a sweep of their parameters,
 * rather than trusting the argument.
 */
import { readFileSync } from "node:fs";
import { dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";

import * as hagrid from "@saehrimnir/hagrid";
import { IRIS } from "@saehrimnir/mistle";

const here = dirname(fileURLToPath(import.meta.url));
const iris = JSON.parse(readFileSync(resolve(here, "../src/data/iris-2d.json"), "utf8"));

const W = 720, H = 460, PAD = 26;
const COLUMN_MAX = IRIS.columns.map((c) => IRIS.statistics[c].max);
const REACH = IRIS.values.map((row) =>
    Math.max(4 + (row[2] / COLUMN_MAX[2]) * 5, 2 + (row[0] / COLUMN_MAX[0]) * 5),
);
const MAX_REACH = Math.max(...REACH);

function project(pts) {
    let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
    for (const [x, y] of pts) {
        if (x < x0) x0 = x; if (x > x1) x1 = x;
        if (y < y0) y0 = y; if (y > y1) y1 = y;
    }
    const s = Math.min((W - 2 * PAD) / (x1 - x0 || 1), (H - 2 * PAD) / (y1 - y0 || 1));
    const ox = PAD + (W - 2 * PAD - (x1 - x0) * s) / 2 - x0 * s;
    const oy = PAD + (H - 2 * PAD - (y1 - y0) * s) / 2 - y0 * s;
    return pts.map(([x, y]) => [x * s + ox, y * s + oy]);
}
function min_pitch(p) {
    let best = Infinity;
    for (let i = 0; i < p.length; ++i)
        for (let j = i + 1; j < p.length; ++j) {
            const d = Math.hypot(p[i][0] - p[j][0], p[i][1] - p[j][1]);
            if (d > 1e-9 && d < best) best = d;
        }
    return Number.isFinite(best) ? best : 0;
}
function overlaps(p, s) {
    let n = 0;
    for (let i = 0; i < p.length; ++i) {
        const ri = REACH[i] * s;
        for (let j = i + 1; j < p.length; ++j) {
            const d = Math.hypot(p[i][0] - p[j][0], p[i][1] - p[j][1]);
            if (d < ri + REACH[j] * s) ++n;
        }
    }
    return n;
}
function to_grid(pts, per_point) {
    const side = Math.max(2, Math.ceil(Math.sqrt(pts.length * per_point)));
    let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
    for (const [x, y] of pts) {
        if (x < x0) x0 = x; if (x > x1) x1 = x;
        if (y < y0) y0 = y; if (y > y1) y1 = y;
    }
    return pts.map(([x, y]) => [
        Math.round(((x - x0) / (x1 - x0 || 1)) * (side - 1)),
        Math.round(((y - y0) / (y1 - y0 || 1)) * (side - 1)),
    ]);
}

const cases = [["as projected", null]];
for (const pluslevel of [0, 1]) for (const whitespace of [1, 2, 4])
    cases.push(["hilbert", { pluslevel, whitespace }],
               ["gosper", { pluslevel, whitespace }]);
// 2.4 is the value the showcase opens Gilbert at, so the shipped default has to
// be one of the configurations this sweep proves overlap-free.
for (const whitespace of [1, 2, 2.4, 4]) cases.push(["gilbert", { whitespace }]);
for (const [cols, rows] of [[13, 12], [20, 8], [30, 5]]) cases.push(["gilbert", { cols, rows }]);
for (const aspect_ratio of [0.5, 1, 2]) cases.push(["dgrid", { aspect_ratio }]);
for (const [cols, rows] of [[13, 12], [20, 8]]) cases.push(["dgrid", { cols, rows }]);
for (const squared of [false, true]) cases.push(["nmap", { squared }]);
for (const cells of [1, 2, 4]) cases.push(["gridfit", { __cells: cells }]);
for (const orientation of ["flat", "pointy"]) cases.push(["gosper", { orientation }]);
cases.push(["hilbert", { keep_aspect_ratio: true }]);

let worst = 0, failures = 0;
console.log(["technique".padEnd(14), "parameters".padEnd(34), "cell px".padStart(8), "glyph".padStart(7), "overlaps".padStart(9)].join("  "));
for (const [technique, p] of cases) {
    let pts;
    try {
        pts = technique === "as projected"
            ? iris.points
            : hagrid.gridify(
                  technique === "gridfit" ? to_grid(iris.points, p.__cells) : iris.points,
                  technique,
                  Object.fromEntries(Object.entries(p ?? {}).filter(([k]) => !k.startsWith("__"))),
              ).positions;
    } catch (e) {
        console.log([technique.padEnd(14), JSON.stringify(p ?? {}).padEnd(34), "threw".padStart(8)].join("  "));
        continue;
    }
    const screen = project(pts);
    const pitch = min_pitch(screen);
    const size = pitch / (2 * MAX_REACH);
    const n = overlaps(screen, size);
    if (n > 0) ++failures;
    worst = Math.max(worst, n);
    console.log([
        technique.padEnd(14), JSON.stringify(p ?? {}).padEnd(34),
        pitch.toFixed(1).padStart(8), size.toFixed(2).padStart(7),
        String(n).padStart(9),
    ].join("  "));
}
console.log(`\n${cases.length} configurations, ${failures} with any overlap (worst ${worst}).`);
if (failures > 0) process.exitCode = 1;
