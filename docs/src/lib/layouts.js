import {
    gridify,
    gridify_cmds,
    cmds_grid,
    gilbert_dimensions,
    gosper_curve,
    gosper_size,
    nmap_regions,
    utils,
} from "@saehrimnir/hagrid";

const { get_bounds, get_scales } = utils;

/**
 * Runs one method and returns everything a drawing needs, in a single
 * coordinate space.
 *
 * The reason this module exists at all is that the methods do not agree on what
 * space they answer in. `gosper`, `nmap` and `cmds` return the input's own
 * coordinates; `hilbert`, `gilbert` and `dgrid` return integer cell indices;
 * `gridfit` returns grid units and must be *given* grid units. So the original
 * points and the layout are only sometimes comparable, and normalising each by
 * its own bounding box — the obvious thing — invents displacement that is not
 * there, because the two bounding boxes differ.
 *
 * Every entry below therefore returns `source`: the input points projected into
 * that method's own output space, using the method's own scaling. A point that
 * did not move then draws exactly on its cell, and what is left on screen is
 * real displacement.
 *
 * Results are shaped as a list of `frames` even though only `cmds` has more
 * than one, so that the scrubber and the readouts have a single shape to read.
 */

/** Cell centres sit on integers, so cell *boundaries* are the half-integers. */
const integer_grid = (cols, rows) => ({
    type: "rect",
    pitch: 1,
    x0: -0.5,
    y0: -0.5,
    cols,
    rows,
});

/** Projects `data` onto `[0, cols-1] x [0, rows-1]`, one axis at a time. */
function to_lattice(data, cols, rows) {
    const { x, y, width, height } = get_bounds(data);
    const sx = width > 0 ? (cols - 1) / width : 0;
    const sy = height > 0 ? (rows - 1) / height : 0;
    return data.map(([px, py]) => [(px - x) * sx, (py - y) * sy]);
}

/**
 * The centre-to-centre spacing of a grid, so displacement can be reported in
 * cell pitches for every method rather than in whatever units it answers in.
 * That is the unit the library's own figures use ("mean displacement 0.73
 * cells"), and unlike a fraction of the extent it does not change meaning when
 * the grid is resized.
 */
export function pitch_of(grid) {
    if (!grid) return 1;
    if (grid.type === "rect") return grid.pitch;
    // Hexagon centres are sqrt(3) * circumradius apart.
    if (grid.type === "hex") return grid.radius * Math.sqrt(3);
    if (grid.type === "regions" && grid.cells.length) {
        const area =
            grid.cells.reduce((s, c) => s + c.width * c.height, 0) / grid.cells.length;
        return Math.sqrt(area) || 1;
    }
    return 1;
}

const single = (source, positions, grid, runtime, facts) => ({
    source,
    frames: [{ positions, grid }],
    runtime,
    facts,
});

function hilbert(data, params) {
    const { pluslevel = 0, whitespace = 1, keep_aspect_ratio = false } = params;
    // Mirrors gridify's own level choice; if this drifts, the grid drawn is not
    // the grid the points were placed on.
    const level = Math.ceil(Math.log2(data.length * whitespace) / Math.log2(4)) + pluslevel;
    const size = 2 ** level;
    const { positions, runtime } = gridify(data, "hilbert", params);
    // gridify_hilbert scales through get_scales, so use it here too rather than
    // a hand-rolled affine — keep_aspect_ratio changes the mapping.
    const scales = get_scales(
        data,
        [
            [0, size - 1],
            [0, size - 1],
        ],
        { keep_aspect_ratio, round: false },
    );
    return single(
        data.map((d) => [scales[0](d), scales[1](d)]),
        positions,
        integer_grid(size, size),
        runtime,
        [
            ["level", level],
            ["grid", `${size} x ${size}`],
            ["cells", size * size],
        ],
    );
}

/**
 * The grid `gilbert` picks for itself at a given whitespace.
 *
 * Exposed so the demo can show the answer in its columns and rows fields rather
 * than leaving them blank: a reader should not have to guess what the library
 * chose before they can adjust it.
 */
export function auto_grid(data, whitespace = 1) {
    const { width, height } = get_bounds(data);
    const aspect = height > 0 ? width / height : 1;
    const target = Math.max(1, Math.ceil(data.length * whitespace));
    const [cols, rows] = gilbert_dimensions(target, aspect);
    return { cols, rows };
}

function gilbert(data, params) {
    // Columns and rows win when set; the demo keeps them filled with whatever
    // the whitespace dial last proposed, so editing one pins the grid.
    const auto = auto_grid(data, params.whitespace ?? 1);
    const W = params.cols || auto.cols;
    const H = params.rows || auto.rows;
    const pinned = W !== auto.cols || H !== auto.rows;
    const { positions, runtime } = gridify(data, "gilbert", { cols: W, rows: H });
    return single(to_lattice(data, W, H), positions, integer_grid(W, H), runtime, [
        ["grid", `${W} x ${H}`],
        ["cells", W * H],
        ["sized by", pinned ? "columns and rows" : "whitespace"],
    ]);
}

function dgrid(data, params) {
    const n = data.length;
    const aspect_ratio = params.aspect_ratio ?? 1;
    let rows = params.rows;
    let cols = params.cols;
    if (rows === undefined || cols === undefined) {
        rows = Math.max(1, Math.floor(Math.sqrt(n * aspect_ratio)));
        cols = Math.ceil(n / rows);
    }
    const { positions, runtime } = gridify(data, "dgrid", params);
    return {
        ...single(
            to_lattice(data, cols, rows),
            positions,
            integer_grid(cols, rows),
            runtime,
            [
                ["grid", `${cols} x ${rows}`],
                ["cells", cols * rows],
            ],
        ),
        // dgrid assigns by rank, not by an affine map, so the projection above
        // is a fair reference frame rather than an exact preimage.
        approximate_source: true,
    };
}

function gosper(data, params) {
    const { pluslevel = 0, whitespace = 1, orientation = "flat" } = params;
    const level = Math.ceil(Math.log2(data.length * whitespace) / Math.log2(7)) + pluslevel;
    const { positions, runtime } = gridify(data, "gosper", params);

    // Rebuild the island the layout was placed on: same size, same centring.
    const size = gosper_size(data, level);
    const { x, y, width, height } = get_bounds(data);
    const mid = [x + width / 2, y + height / 2];
    const centres = gosper_curve(level, orientation, size).map(([cx, cy]) => [
        cx + mid[0],
        cy + mid[1],
    ]);

    return single(data, positions, hex_grid(centres), runtime, [
        ["level", level],
        ["cells", 7 ** level],
    ]);
}

/**
 * Hexagons have no axis-aligned lattice to step, so they are drawn per cell.
 * Radius and rotation are measured off the decoded centres rather than derived
 * from the orientation constant — that way the drawing cannot disagree with the
 * lattice, whichever way up it is.
 */
function hex_grid(centres) {
    if (centres.length < 2) return { type: "hex", centres, radius: 1, angle: 0 };
    const [ax, ay] = centres[0];
    let best = Infinity;
    let angle = 0;
    for (let i = 1; i < centres.length; ++i) {
        const dx = centres[i][0] - ax;
        const dy = centres[i][1] - ay;
        const d = Math.hypot(dx, dy);
        if (d > 1e-9 && d < best) {
            best = d;
            angle = Math.atan2(dy, dx);
        }
    }
    // Neighbouring centres of a hex lattice are sqrt(3) * circumradius apart,
    // and a vertex sits 30 degrees off the direction to a neighbour.
    return {
        type: "hex",
        centres,
        radius: best / Math.sqrt(3),
        angle: angle + Math.PI / 6,
    };
}

function nmap(data, params) {
    // Always the squared variant: padding the point set to a square number gives
    // the treemap a full grid to subdivide, which is the form worth comparing
    // against the curve methods.
    const p = { ...params, squared: true };
    const { positions, runtime } = gridify(data, "nmap", p);
    const cells = nmap_regions(data, p);
    return single(data, positions, { type: "regions", cells }, runtime, [
        ["rectangles", cells.length],
    ]);
}

function gridfit(data, params) {
    // gridfit is the one method whose *input* is grid units, so the data is
    // scaled before it is handed over, not after.
    const n = data.length;
    const whitespace = params.whitespace ?? 2.5;
    const { width, height } = get_bounds(data);
    const aspect = height > 0 ? width / height : 1;
    const target = Math.max(1, Math.ceil(n * whitespace));
    const rows = Math.max(1, Math.round(Math.sqrt(target / Math.max(aspect, 1e-9))));
    const cols = Math.max(1, Math.ceil(target / rows));
    const scaled = to_lattice(data, cols, rows);
    const { positions, runtime } = gridify(scaled, "gridfit", params);
    return single(scaled, positions, integer_grid(cols, rows), runtime, [
        ["input grid", `${cols} x ${rows}`],
        ["cells", cols * rows],
    ]);
}

/** @type {Record<string, (data: any, params: any) => any>} */
const runners = { hilbert, gilbert, gosper, dgrid, nmap, gridfit };

/**
 * Runs `method` and returns `{source, frames, runtime, facts}`, or `{error}` if
 * the method rejected the input. Several methods throw by design — `nmap` on
 * quantised data, `hilbert` when the grid is too small — and the demos are meant
 * to show that happening rather than hide it.
 */
export function layout(method, data, params = {}) {
    const run = runners[method];
    if (!run) throw new Error(`layout: no runner for ${method}`);
    try {
        const result = run(data, params);
        const last = result.frames[result.frames.length - 1].positions;
        const bad = last.filter(
            (p) => !p || !Number.isFinite(p[0]) || !Number.isFinite(p[1]),
        ).length;
        if (bad) {
            return {
                error: `${bad} of ${last.length} positions came back non-finite — this method cannot lay out this dataset.`,
            };
        }
        return { ...result, error: null };
    } catch (error) {
        return { error: error instanceof Error ? error.message : String(error) };
    }
}

/**
 * Drives `cmds` to completion, collecting every intermediate layout so the
 * result can be scrubbed and replayed.
 *
 * The generator is synchronous and there is no way to make it otherwise, so this
 * blocks for as long as the optimisation takes — seconds, on a few hundred
 * points. That is why it runs in a worker: see `layout-worker.js`. `on_progress`
 * reports the round count as it goes, since a run long enough to need a worker
 * is long enough to need a progress readout.
 */
export function run_cmds(data, params = {}, on_progress = () => {}) {
    const { x, y, width, height } = get_bounds(data);
    const w = Math.max(width, height);
    const Gamma = { x, y, width: w, height: w };

    // Replays the scaling gridify_cmds applies on the way in (§4.2: the input is
    // uniformly scaled into the inset Gamma). Without this the "before" points
    // are drawn somewhere the optimisation never saw them.
    const [m, k] = cmds_grid(Gamma.width, Gamma.height, data.length);
    const size0 = Math.min(Gamma.width / m, Gamma.height / k);
    const iw = Math.max(Gamma.width - size0, 0);
    const ih = Math.max(Gamma.height - size0, 0);
    const kx = width > 0 ? iw / width : Infinity;
    const ky = height > 0 ? ih / height : Infinity;
    const s = Math.min(kx, ky);
    const scale = Number.isFinite(s) && s > 0 ? s : 1;
    const cx = Gamma.x + Gamma.width / 2 - (x + width / 2) * scale;
    const cy = Gamma.y + Gamma.height / 2 - (y + height / 2) * scale;
    const source = data.map(([px, py]) => [px * scale + cx, py * scale + cy]);

    const started = performance.now();
    const frames = [];
    try {
        const steps = gridify_cmds(data, { Gamma, ...params });
        for (const [positions, , size] of steps) {
            frames.push({
                positions: positions.map(([px, py]) => [px, py]),
                // Cell centres land on (k + 1/2) * size, so boundaries are the
                // multiples of size — in absolute coordinates, not Gamma's.
                grid: {
                    type: "rect",
                    pitch: size,
                    x0: Math.floor(Gamma.x / size) * size,
                    y0: Math.floor(Gamma.y / size) * size,
                    cols: Math.ceil(Gamma.width / size) + 1,
                    rows: Math.ceil(Gamma.height / size) + 1,
                },
            });
            if (frames.length % 25 === 0) on_progress(frames.length);
        }
    } catch (error) {
        return { error: error instanceof Error ? error.message : String(error) };
    }

    if (!frames.length) return { error: "cmds produced no layout." };
    const last = frames[frames.length - 1].positions;
    const bad = last.filter((p) => !p || !Number.isFinite(p[0]) || !Number.isFinite(p[1]))
        .length;
    if (bad) {
        // Coincident points have an ideal distance of zero, so the weight 1/d^2
        // is infinite and the whole layout goes non-finite.
        return {
            error: `${bad} of ${last.length} positions came back non-finite — cmds cannot lay out data with exact duplicates.`,
        };
    }

    return {
        source,
        frames,
        runtime: performance.now() - started,
        facts: [
            ["rounds", frames.length],
            ["node size", frames[frames.length - 1].grid.pitch.toFixed(2)],
        ],
        error: null,
    };
}

/** Points sharing a cell with another point. Should be zero. */
export function count_distinct(positions) {
    const seen = new Set();
    for (const [x, y] of positions) seen.add(`${x},${y}`);
    return seen.size;
}

/** Per-point distance from where it started, in cell pitches. */
export function displacements(source, positions, pitch = 1) {
    const n = Math.min(source.length, positions.length);
    const out = new Array(n);
    for (let i = 0; i < n; ++i) {
        out[i] =
            Math.hypot(positions[i][0] - source[i][0], positions[i][1] - source[i][1]) /
            (pitch || 1);
    }
    return out;
}

/** Five-number summary plus the mean and Tukey fences, for the box plot. */
export function summarise(d) {
    const empty = {
        mean: 0, median: 0, p95: 0, max: 0, min: 0,
        q1: 0, q3: 0, lo_fence: 0, hi_fence: 0, outliers: [],
    };
    if (!d.length) return empty;
    const sorted = [...d].sort((a, b) => a - b);
    const at = (p) => sorted[Math.min(sorted.length - 1, Math.floor(sorted.length * p))];
    const q1 = at(0.25);
    const q3 = at(0.75);
    const iqr = q3 - q1;
    // Whiskers reach the furthest point still inside 1.5 IQR, not the fence
    // itself, so a tight distribution does not draw whiskers into empty space.
    const lo_fence = sorted.find((v) => v >= q1 - 1.5 * iqr) ?? sorted[0];
    let hi_fence = sorted[0];
    for (const v of sorted) if (v <= q3 + 1.5 * iqr) hi_fence = v;
    return {
        mean: d.reduce((s, v) => s + v, 0) / d.length,
        median: at(0.5),
        p95: at(0.95),
        max: sorted[sorted.length - 1],
        min: sorted[0],
        q1,
        q3,
        lo_fence,
        hi_fence,
        outliers: sorted.filter((v) => v < lo_fence || v > hi_fence),
    };
}
