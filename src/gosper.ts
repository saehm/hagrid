/**
 * Gosper (Node-Gosper) space-filling curve indexing, after
 *
 *   Uher, Gajdoš, Snášel, Lai & Radecký,
 *   "Hierarchical Hexagonal Clustering and Indexing",
 *   Symmetry 11(6):731, 2019. https://doi.org/10.3390/sym11060731
 *
 * A point's curve index is computed arithmetically in `O(level)`; the curve is
 * never materialized, so memory is constant in the level and a level of 12
 * (7^12 ≈ 1.4e10 cells) costs no more than a level of 3.
 *
 * Encoding runs bottom-up (Algorithm 1) because Gosper islands are non-convex
 * fractal tiles whose jagged edges overstep their parent hexagon; a top-down
 * descent misassigns points near those edges. Hexagon centres always lie
 * inside their parent, so climbing the hierarchy is exact. The resulting
 * indices are then rewritten top-down into Node-Gosper order (Algorithm 2),
 * which is where pattern rotation and traversal order are resolved.
 *
 * This replaced an earlier implementation that built the whole curve as an
 * L-system and matched points against it. That one was ~110x slower end to
 * end, needed `O(7^level)` memory, and threw on collinear input.
 */
import type { HexOrientation, Point } from "./types.js";


const REC_SQRT7 = 0.3779644730092272; // 1 / Math.sqrt(7)
const SQRT3_DIV_3 = 0.5773502691896257; // Math.sqrt(3) / 3
const ALPHA = 0.3334731722518321; // Math.asin(Math.sqrt(3) / (2 * Math.sqrt(7)))

/**
 * Transfer function E(P2, ·) from the paper's Table 1, mapping a default
 * pattern index onto its position in the Node-Gosper pattern.
 */
const transform_index = [4, 0, 1, 2, 3, 6, 5];

/**
 * Paper eq. (13): `T = C⁻¹ · G · C = [[2/7, -1/7], [1/7, 3/7]]`, the level
 * `i` to level `i-1` transform on axial coordinates. Stored row-major.
 */
const T = [2 / 7, -1 / 7, 1 / 7, 3 / 7];

/** Radii of the circles inscribed in the Gosper island, indexed by level. */
const inscribed_circles = [
    0.755928946, 0.755928946, 0.750121467308, 0.746782631146, 0.746782631146,
    0.746577727521, 0.746348363909, 0.746348363909, 0.746344578768, 0.746327538283,
    0.746327538283, 0.746327538283, 0.746326555879, 0.746326555879, 0.746326555879,
    0.746326510616, 0.746326510616, 0.746326510616, 0.746326508597, 0.746326508597,
    0.746326508597,
];

/**
 * Lattice convention, paper eqs. (11) and (12).
 *
 * Hexagons are laid out **flat-topped**: axial `(q, l)` sits at
 * `x = 3/2·q·size`, `y = √3·(q/2 + l)·size`, which puts axial `(1, 0)` at 30°.
 * The inverse is `q = (2x/3)/size`, `l = (-x/3 + y·√3/3)/size`.
 *
 * These two must agree. The original paired the flat-topped forward map with a
 * *pointy-topped* inverse, so the two disagreed about the lattice itself and
 * every index above level 0 came out wrong. Both directions are now inlined
 * into encode and decode, which is why there are no helpers here to drift
 * apart again.
 *
 * Cube coordinates add `y = -x - z` so the three axes sum to zero; rounding
 * restores that invariant by rebuilding whichever axis rounded least
 * confidently.
 */

/**
 * Classifies which of the six neighbouring hexagons a fractional cube offset
 * points at, or `0` when the offset is negligible. Paper eq. (16).
 *
 * The original inverted the sign test on all three axes and used a threshold
 * of `10e-3` where the paper specifies `1e-5`. It also fell off the end
 * returning `undefined` on a tie between two dominant axes, which propagated
 * `NaN` through the index; ties now resolve to the first dominant axis.
 */
function B(xd: number, yd: number, zd: number, thd = 1e-5): number {
    const ax = Math.abs(xd);
    const ay = Math.abs(yd);
    const az = Math.abs(zd);
    if (ax < thd && ay < thd && az < thd) {
        return 0;
    } else if (ax >= ay && ax >= az) {
        return xd < 0 ? 2 : 5;
    } else if (ay >= ax && ay >= az) {
        return yd < 0 ? 4 : 1;
    }
    return zd < 0 ? 6 : 3;
}

/** Side length of a level-0 hexagon such that `D` fits the Gosper island. */
export function gosper_size(D: readonly Point[], level: number): number {
    let min_x = Infinity;
    let max_x = -Infinity;
    let min_y = Infinity;
    let max_y = -Infinity;
    D.forEach(([x, y]) => {
        min_x = x < min_x ? x : min_x;
        max_x = x > max_x ? x : max_x;
        min_y = y < min_y ? y : min_y;
        max_y = y > max_y ? y : max_y;
    });

    const half_diagonal =
        Math.sqrt(Math.pow(max_x - min_x, 2) + Math.pow(max_y - min_y, 2)) / 2;
    const s = half_diagonal / inscribed_circles[level];
    return s * Math.pow(REC_SQRT7, level);
}

// Rotations applied when descending into a child node, in radians.
const f_yArr = [-2.0943951023931953, 0, 0, -2.0943951023931953, 0, 2.0943951023931953, 0];
const gammaArr = [
    0, -1.0471975511965976, -2.0943951023931953, -3.141592653589793, 0,
    1.0471975511965976, 2.0943951023931953,
];

// The per-level rotations only ever come from these seven fixed angles, so
// their sines and cosines are tabulated instead of recomputed per cell.
const COS_GAMMA = gammaArr.map(Math.cos);
const SIN_GAMMA = gammaArr.map(Math.sin);
const COS_STEP = f_yArr.map((f) => Math.cos(ALPHA + f));
const SIN_STEP = f_yArr.map((f) => Math.sin(ALPHA + f));

const COS_30 = Math.cos(Math.PI / 6);
const SIN_30 = Math.sin(Math.PI / 6);

/**
 * Scratch buffer for the per-level digits.
 *
 * Both encode and decode need every digit before they can use any of them, but
 * allocating an array per call dominated the arithmetic. Neither function
 * yields, recurses or calls out, so a single shared buffer is safe.
 */
let digit_scratch = new Int8Array(32);
function digits(level: number): Int8Array {
    if (digit_scratch.length < level) digit_scratch = new Int8Array(level * 2);
    return digit_scratch;
}

/**
 * This construction is natively flat-topped, because {@link cube2xy} sends
 * axial `(1, 0)` to 30°. The pointy-topped lattice is the same lattice turned
 * by -30°, so a point is rotated into the native frame on the way in and the
 * cell is rotated back out again. Rotation leaves the index ordering alone.
 */
function to_native(p: Point, orientation: HexOrientation): Point {
    if (orientation === "flat") return p;
    return [COS_30 * p[0] - SIN_30 * p[1], SIN_30 * p[0] + COS_30 * p[1]];
}

function from_native(p: Point, orientation: HexOrientation): Point {
    if (orientation === "flat") return p;
    return [COS_30 * p[0] + SIN_30 * p[1], -SIN_30 * p[0] + COS_30 * p[1]];
}

/**
 * Maps a point to its position along the Gosper curve.
 *
 * Returns `-1` when the point lies outside the level-`level` Gosper island.
 * The climb only ever records a hexagon's position *relative* to its parent,
 * so without this check an outside point aliases onto an arbitrary interior
 * cell — frequently one at the far end of the curve — rather than being
 * rejected. A contained point always lands on the root hexagon `(0, 0)`.
 *
 * @param point - The point to index.
 * @param level - Curve level; the curve has `7^level` cells.
 * @param size - Level-0 hexagon size, from {@link gosper_size}.
 * @param orientation - Which way up the hexagons sit. Defaults to `"flat"`.
 * @returns The curve index in `[0, 7^level)`, or `-1` if outside the island.
 */
export function gosper_encode(
    point: Point,
    level: number,
    size: number,
    orientation: HexOrientation = "flat",
): number {
    let xp = point[0];
    let yp = point[1];
    if (orientation !== "flat") {
        const rx = COS_30 * xp - SIN_30 * yp;
        yp = SIN_30 * xp + COS_30 * yp;
        xp = rx;
    }
    // Algorithm 1: bottom-up hierarchical localization. Mapping must go
    // bottom-up because Gosper islands are non-convex fractal tiles whose
    // jagged edges overstep their parent hexagon — a top-down descent
    // misassigns points near those edges. Hexagon *centres*, by contrast,
    // always lie inside the parent hexagon, so climbing is exact.
    //
    // The original walked xy space applying cumulative `REC_SQRT7^i` scales
    // and `-ALPHA*i` rotations. The paper instead applies one constant matrix
    // to axial coordinates, repeatedly.
    const b = digits(level);

    // xy2cube, eq. (12), inlined.
    const cx = ((2 * xp) / 3) / size;
    const cz = (-xp / 3 + yp * SQRT3_DIV_3) / size;
    const cy = -cx - cz;

    // cube_round, inlined: round each axis, then rebuild the least confident
    // one from the other two so the x + y + z = 0 invariant survives.
    let rq = Math.round(cx);
    let rm = Math.round(cy);
    let rl = Math.round(cz);
    const ex = Math.abs(rq - cx);
    const ey = Math.abs(rm - cy);
    const ez = Math.abs(rl - cz);
    if (ex > ey && ex > ez) rq = -rm - rl;
    else if (ey > ez) rm = -rq - rl;
    else rl = -rq - rm;

    let q = rq;
    let l = rl;

    for (let i = 1; i <= level; ++i) {
        // Paper eq. (14): a_{i-1} = T · a_i, with T from eq. (13).
        const nq = T[0] * q + T[1] * l;
        const nl = T[2] * q + T[3] * l;
        const vy = -nq - nl;
        q = Math.round(nq);
        l = Math.round(nl);
        // Paper eq. (15): the fractional part locates the child within the parent.
        b[level - i] = B(nq - q, vy - Math.round(vy), nl - l);
    }

    // The climb ends on the root hexagon for any point the island contains.
    // Anything else started in a neighbouring island, and the digits collected
    // above describe a cell of *that* tile — which, folded back into this one,
    // lands somewhere unrelated. Report no containing cell instead.
    if (q !== 0 || l !== 0) return -1;

    // Algorithm 2. Rewrite the raw hexagon sequence into Gosper node order,
    // tracking the rotation and traversal direction each node imposes on its
    // children.
    let rotDir = 0;
    let order = true;
    const k: number[] = new Array(level).fill(0);
    for (let i = 0; i < level; ++i) {
        let bi = b[i];
        if (bi !== 0 && rotDir !== 0) {
            bi += 2 * rotDir;
            if (bi < 1) {
                bi += 6;
            } else if (bi > 6) {
                bi -= 6;
            }
        }
        const ki = transform_index[bi];
        if (ki === 0 || ki === 3) {
            --rotDir;
            if (rotDir < -1) rotDir = 1;
        } else if (ki === 5) {
            ++rotDir;
            if (rotDir > 1) rotDir = -1;
        }
        k[i] = order ? ki : 6 - ki;
        if (ki === 0 || ki === 4 || ki === 5) {
            order = !order;
        }
    }
    return k.reduce((a, digit) => a * 7 + digit, 0);
}

/**
 * Maps a position along the Gosper curve back to its cell centre — the inverse
 * of {@link gosper_encode}. Paper Section 3.3.
 *
 * @param t - Curve index, in `[0, 7^level)`.
 * @param level - Curve level.
 * @param size - Level-0 hexagon size, from {@link gosper_size}.
 * @param orientation - Which way up the hexagons sit. Defaults to `"flat"`.
 */
export function gosper_decode(
    t: number,
    level: number,
    size: number,
    orientation: HexOrientation = "flat",
): Point {
    const k = digits(level);
    for (let i = level - 1; i >= 0; --i) {
        const mod_t = t % 7;
        k[i] = mod_t;
        t = (t - mod_t) / 7;
    }

    // inverse_transform, inlined: expand the unit step out to level 0.
    let dx = -size;
    let dz = 0;
    for (let lv = 0; lv < level - 1; ++lv) {
        const nx = 3 * dx + dz;
        dz = -dx + 2 * dz;
        dx = nx;
    }
    // cube2xy at unit size.
    let sx = (3 / 2) * dx;
    let sy = Math.sqrt(3) * (dx / 2 + dz);

    let cxp = 0;
    let cyp = 0;
    let ord = true;
    for (let i = 0; i < level; ++i) {
        const ki = i > 0 && !ord ? 6 - k[i] : k[i];
        if (ki !== 4) {
            cxp += COS_GAMMA[ki] * sx - SIN_GAMMA[ki] * sy;
            cyp += SIN_GAMMA[ki] * sx + COS_GAMMA[ki] * sy;
        }
        const nsx = COS_STEP[ki] * sx - SIN_STEP[ki] * sy;
        const nsy = SIN_STEP[ki] * sx + COS_STEP[ki] * sy;
        sx = nsx * REC_SQRT7;
        sy = nsy * REC_SQRT7;
        if (ki === 0 || ki === 4 || ki === 5) {
            ord = !ord;
        }
    }
    if (orientation !== "flat") {
        const rx = COS_30 * cxp + SIN_30 * cyp;
        cyp = -SIN_30 * cxp + COS_30 * cyp;
        cxp = rx;
    }
    return [cxp, cyp];
}

function distance([ax, ay]: Point, [bx, by]: Point): number {
    return Math.hypot(bx - ax, by - ay);
}

/** Resolves a contested cell by walking outward along the curve. */
function collision(
    P: Map<number, number>,
    p: number,
    d: Point,
    i: number,
    level: number,
    size: number,
): void {
    const e = 7 ** level;
    const valid = (q: number) => q >= 0 && q <= e;
    let pl = p;
    let pr = p;
    for (;;) {
        ++pl;
        --pr;
        const el = !P.has(pl);
        const er = !P.has(pr);
        const vl = valid(pl);
        const vr = valid(pr);
        // Both probes have run off the ends, so every cell is taken. Without
        // this the walk spins forever, which is what a level too small for
        // the data used to produce.
        if (!vl && !vr) {
            throw new Error(
                "gosper_uher: more points than cells — raise the level (the curve has 7^level cells)",
            );
        }
        if (vl && el && !er) {
            P.set(pl, i);
            return;
        } else if (!el && vr && er) {
            P.set(pr, i);
            return;
        } else if (el && er) {
            if (vl && vr) {
                const dl = distance(d, gosper_decode(pl, level, size));
                const dr = distance(d, gosper_decode(pr, level, size));
                P.set(dl < dr ? pl : pr, i);
                return;
            } else if (vl) {
                P.set(pl, i);
                return;
            } else if (vr) {
                P.set(pr, i);
                return;
            }
        }
    }
}

/**
 * Every cell centre of the Gosper curve at `level`, in curve order.
 *
 * Useful for drawing the grid. Nothing else needs it — indexing goes straight
 * through {@link gosper_encode} — so this is the only place the curve is ever
 * materialized, and the caller decides when to pay for it.
 *
 * @param level - Curve level; the curve has `7^level` cells.
 * @param orientation - Which way up the hexagons sit. Defaults to `"flat"`.
 * @param size - Level-0 hexagon size. Defaults to `1`.
 */
export function gosper_curve(
    level: number,
    orientation: HexOrientation = "flat",
    size = 1,
): Point[] {
    const n = 7 ** level;
    const V: Point[] = new Array(n);
    for (let t = 0; t < n; ++t) V[t] = gosper_decode(t, level, size, orientation);
    return V;
}

/** Parameters for {@link gridify_gosper}. */
export interface GosperParameters {
    /** Curve level; the curve has `7^level` cells. */
    level: number;
    /** Which way up the hexagons sit. Defaults to `"flat"`. */
    orientation?: HexOrientation;
}

/**
 * Removes overlaps using the direct Gosper indexing of Uher et al.
 *
 * Produces a hexagonal grid like {@link gridify_gosper}, but without building
 * the curve, so it stays fast as the level grows. Cell coordinates are in the
 * paper's own frame and differ in scale and orientation from
 * {@link gridify_gosper}'s.
 *
 * @param D - Input points.
 * @returns One cell centre per input point, in input order.
 */
export function gridify_gosper(
    D: readonly Point[],
    { level, orientation = "flat" }: GosperParameters,
): Point[] {
    const size = gosper_size(D, level);
    const N = D.length;
    const P = new Map<number, number>();
    const Y: Point[] = new Array(N).fill(0);

    // The island is centred on the origin, so the data has to be moved there
    // before indexing. Without this the points sit off the tile entirely and
    // over half of them alias onto unrelated cells. `gosper_size` scales
    // the island so the data's bounding circle lands inside its inscribed
    // circle, which only holds about a shared centre.
    let min_x = Infinity;
    let max_x = -Infinity;
    let min_y = Infinity;
    let max_y = -Infinity;
    D.forEach(([x, y]) => {
        if (x < min_x) min_x = x;
        if (x > max_x) max_x = x;
        if (y < min_y) min_y = y;
        if (y > max_y) max_y = y;
    });
    const mid_x = (min_x + max_x) / 2;
    const mid_y = (min_y + max_y) / 2;

    D.forEach((d, i) => {
        // Note: no rounding. The original quantized inputs to integers here,
        // which collapsed distinct points on any small-scale dataset.
        // Everything below stays in the native (flat) frame; the requested
        // orientation is applied once on the way in and once on the way out,
        // so collision resolution never has to know about it.
        const centred = to_native([d[0] - mid_x, d[1] - mid_y], orientation);
        const p = gosper_encode(centred, level, size);
        if (p < 0) {
            throw new Error(
                "gridify_gosper: point fell outside the Gosper island; " +
                    "the level is too small for this data",
            );
        }
        if (P.has(p)) {
            collision(P, p, centred, i, level, size);
        } else {
            P.set(p, i);
        }
    });

    for (const [p, i] of P.entries()) {
        const [x, y] = from_native(gosper_decode(p, level, size), orientation);
        Y[i] = [x + mid_x, y + mid_y];
    }
    return Y;
}
