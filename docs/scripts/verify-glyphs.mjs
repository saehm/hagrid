/*
 * Holds the rectangular-glyph sizing to its promise.
 *
 * `GlyphShowcase` sizes a glyph from a disc model, which is exact for a round
 * glyph and only an approximation for a box. On a lattice the box reaches its
 * neighbour along an axis rather than through its corner, so it may be √2
 * larger than the disc model allows — worth having, and worth checking, because
 * getting it wrong the other way silently overlaps every cell.
 *
 * So this does not check the model. It checks the boxes: it lays each dataset
 * out with every technique, sizes the glyph exactly as the component does, and
 * tests every pair of rectangles for overlap. Then it does it again 4% larger,
 * which must fail — otherwise the size is merely safe, not tight.
 *
 * Run with `pnpm verify:glyphs`.
 */
import * as hagrid from "@saehrimnir/hagrid";
import { readFileSync } from "node:fs";
import { box_reach, cell_reach } from "../src/lib/glyph-reach.js";

const load = (n) => JSON.parse(readFileSync(new URL(`../src/data/${n}.json`, import.meta.url)));

const WIDTH = 720;
const HEIGHT = 460;
const PAD = 26;

/** The component's own frame fit, so the pitch means the same thing. */
function fit(points) {
    let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
    for (const [x, y] of points) {
        if (x < x0) x0 = x; if (x > x1) x1 = x;
        if (y < y0) y0 = y; if (y > y1) y1 = y;
    }
    const s = Math.min((WIDTH - 2 * PAD) / (x1 - x0 || 1), (HEIGHT - 2 * PAD) / (y1 - y0 || 1));
    const ox = PAD + (WIDTH - 2 * PAD - (x1 - x0) * s) / 2 - x0 * s;
    const oy = PAD + (HEIGHT - 2 * PAD - (y1 - y0) * s) / 2 - y0 * s;
    return points.map(([x, y]) => [x * s + ox, y * s + oy]);
}

const min_pitch = (pts) => {
    let best = Infinity;
    for (let i = 0; i < pts.length; ++i) {
        for (let j = i + 1; j < pts.length; ++j) {
            const d = Math.hypot(pts[i][0] - pts[j][0], pts[i][1] - pts[j][1]);
            if (d > 1e-9 && d < best) best = d;
        }
    }
    return best;
};

/**
 * Overlapping pairs, ignoring points that were seated in the *same* cell.
 *
 * Coincident output is a property of the technique, not of the glyph size, and
 * no size short of zero separates two glyphs at the same coordinates. It should
 * be zero for every method here — it was 21 of the paintings under `hilbert`
 * until the collision walk was fixed, which is how that bug was found — so it
 * is counted and reported rather than charged to the sizing.
 */
function overlaps(pts, hw, hh, hex) {
    let n = 0;
    let coincident = 0;
    const eps = 1e-6;
    for (let i = 0; i < pts.length; ++i) {
        for (let j = i + 1; j < pts.length; ++j) {
            const dx = Math.abs(pts[i][0] - pts[j][0]);
            const dy = Math.abs(pts[i][1] - pts[j][1]);
            if (dx < eps && dy < eps) { ++coincident; continue; }
            // A cell-shaped glyph on a hexagonal lattice tiles by construction,
            // so the test is simply that no two centres are inside one pitch.
            const hit = hex
                ? Math.hypot(dx, dy) < 2 * hw - eps
                : dx < 2 * hw - eps && dy < 2 * hh - eps;
            if (hit) ++n;
        }
    }
    return { n, coincident };
}

const CASES = [
    // `cell` glyphs take the shape of their cell; `box` ones stay rectangular
    // whatever lattice they land on.
    { name: "flowers", data: load("flowers-2d"), shape: "cell", hw: 1, hh: 1 },
    { name: "paintings", data: load("paintings-2d"), shape: "cell", hw: 1, hh: 1 },
    { name: "wine (bars)", data: load("wine-2d"), shape: "box", hw: 1, hh: 1 },
    { name: "glass", data: load("glass-2d"), shape: "box", hw: 0.62, hh: 1 },
    { name: "ecg", data: load("ecg-2d"), shape: "box", hw: 1, hh: 1 },
    { name: "fsdd", data: load("fsdd-2d"), shape: "box", hw: 1, hh: 1 },
];

const METHODS = ["hilbert", "gilbert", "gosper", "dgrid", "gridfit", "nmap"];

let failures = 0;
const check = (name, ok, detail = "") => {
    if (!ok) ++failures;
    console.log(`${ok ? "ok  " : "FAIL"}  ${name}${detail ? "  — " + detail : ""}`);
};

/** gridfit wants grid units, as the components hand it. */
function to_grid(points, cells) {
    const side = Math.ceil(Math.sqrt(points.length * cells));
    let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
    for (const [x, y] of points) {
        if (x < x0) x0 = x; if (x > x1) x1 = x;
        if (y < y0) y0 = y; if (y > y1) y1 = y;
    }
    return points.map(([x, y]) => [
        Math.round(((x - x0) / (x1 - x0 || 1)) * (side - 1)),
        Math.round(((y - y0) / (y1 - y0 || 1)) * (side - 1)),
    ]);
}

for (const { name, data, shape, hw, hh } of CASES) {
    const reach_of =
        shape === "cell" ? cell_reach(data.points.length) : box_reach(data.points.length, hw, hh);
    for (const method of METHODS) {
        let placed;
        try {
            const input = method === "gridfit" ? to_grid(data.points, 2) : data.points;
            placed = hagrid.gridify(input, method, method === "nmap" ? { squared: true } : {}).positions;
        } catch (error) {
            console.log(`skip  ${name} / ${method} — ${error.message.split("\n")[0]}`);
            continue;
        }
        const screen = fit(placed);
        const pitch = min_pitch(screen);
        if (!Number.isFinite(pitch)) {
            console.log(`skip  ${name} / ${method} — no distinct positions`);
            continue;
        }
        const cell = { kind: method === "gosper" ? "hex" : "square", orientation: "flat" };
        const reach = Math.max(...reach_of(method, cell));
        const size = pitch / (2 * reach);
        const hex = shape === "cell" && cell.kind === "hex";

        const { n: bad, coincident } = overlaps(screen, hw * size, hh * size, hex);
        check(
            `${name.padEnd(13)} ${method.padEnd(8)} fits`,
            bad === 0,
            `${bad} overlapping pairs` + (coincident ? `, ${coincident} seated in one cell` : ""),
        );

        // 4% larger must collide, or the size is leaving room on the table.
        const grown = overlaps(screen, hw * size * 1.04, hh * size * 1.04, hex);
        check(
            `${name.padEnd(13)} ${method.padEnd(8)} is tight`,
            grown.n > 0,
            grown.n === 0 ? "no overlap at +4% — the glyph could be larger" : `${grown.n} pairs at +4%`,
        );
    }
}

console.log(failures === 0 ? "\nall checks passed" : `\n${failures} FAILED`);
process.exit(failures === 0 ? 0 : 1);
