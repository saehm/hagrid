import { get_bounds } from "./utils.js";
import type { Bounds, Point } from "./types.js";

/**
 * A point that may be a synthetic filler added by {@link squared}, in which
 * case its index slot is `undefined`. `dist` caches the nearest-neighbour
 * distance used to rank candidates.
 */
type FillerPoint = [number, number, number | undefined] & { dist?: number };

/** A rectangle in the subdivision, once assigned, carries its point. */
type Tile = Bounds & { d?: FillerPoint };

/** A 2×3 affine matrix, row-major, as `[[a, b, tx], [c, d, ty]]`. */
type Affine = [[number, number, number], [number, number, number]];

/**
 * Rejects a subdivision that cannot be scaled.
 *
 * The split boundary is the midpoint between the two halves along the split
 * axis. When every point in the region shares that coordinate the boundary
 * lands on the region's own edge, one side gets zero extent, and the affine
 * scale becomes `share / 0`.
 *
 * Distinct points are not enough to avoid this: four points on a vertical line
 * all differ yet share an `x`, and an integer lattice trips it too.
 */
function assert_splittable(axis: "x" | "y", a: number, b: number, count: number): void {
    // Only an exact zero is fatal. Extents are routinely negative here — the
    // points are sorted descending and the affine transforms can flip a
    // rectangle — and a negative extent still divides perfectly well.
    if (a !== 0 && b !== 0 && Number.isFinite(a) && Number.isFinite(b)) return;
    throw new Error(
        `gridify_nmap: cannot subdivide ${count} points along ${axis} — every point in ` +
            `this region shares the same ${axis} coordinate, so there is no extent to ` +
            `divide. nmap needs each subdivision to have spread on the axis it is cut ` +
            `along; integer or rounded coordinates commonly break this. Jitter the input ` +
            `slightly, or use another method such as "gilbert" or "hilbert".`,
    );
}

/**
 * Recursively splits the region alongside the point set, alternating between
 * vertical and horizontal cuts, until each region holds a single point.
 */
function nmap(G: Tile[], P: FillerPoint[], R: Tile, direction = true): void {
    const N = P.length;
    const N_half = Math.floor(N / 2);
    if (N === 1) {
        R.d = P[0];
        G.push(R);
    } else {
        const dim = direction ? 0 : 1;
        P = P.sort((a, b) => b[dim] - a[dim]);
        const [Da, Db] = split(P, N_half);
        const Ra: Tile = Object.assign({}, R);
        const Rb: Tile = Object.assign({}, R);
        if (direction) {
            const [wRa, wRb] = sizes(Da, Db, R.width);
            const bh = (Da[Da.length - 1][0] + Db[0][0]) / 2;
            Ra.width = bh - R.x;
            Rb.x += Ra.width;
            Rb.width = R.width - Ra.width;
            assert_splittable("x", Ra.width, Rb.width, N);
            const HRa: Affine = [
                [wRa / Ra.width, 0, R.x * (1 - wRa / Ra.width)],
                [0, 1, 0],
            ];
            const HRb: Affine = [
                [wRb / Rb.width, 0, (R.x + R.width) * (1 - wRb / Rb.width)],
                [0, 1, 0],
            ];
            nmap(G, affine_transform_points(Da, HRa), affine_transform_rect(Ra, HRa), !direction);
            nmap(G, affine_transform_points(Db, HRb), affine_transform_rect(Rb, HRb), !direction);
        } else {
            const [hRa, hRb] = sizes(Da, Db, R.height);
            const bv = (Da[Da.length - 1][1] + Db[0][1]) / 2;
            Ra.height = bv - R.y;
            Rb.y = R.y + Ra.height;
            Rb.height = R.height - Ra.height;
            assert_splittable("y", Ra.height, Rb.height, N);
            const VRa: Affine = [
                [1, 0, 0],
                [0, hRa / Ra.height, R.y * (1 - hRa / Ra.height)],
            ];
            const VRb: Affine = [
                [1, 0, 0],
                [0, hRb / Rb.height, (R.y + R.height) * (1 - hRb / Rb.height)],
            ];
            nmap(G, affine_transform_points(Da, VRa), affine_transform_rect(Ra, VRa), !direction);
            nmap(G, affine_transform_points(Db, VRb), affine_transform_rect(Rb, VRb), !direction);
        }
    }
}

function affine_transform_points(
    D: FillerPoint[],
    [[scale_x, shear_x, translate_x], [scale_y, shear_y, translate_y]]: Affine,
): FillerPoint[] {
    return D.map(
        ([x, y, i]) =>
            [
                x * scale_x + y * shear_x + translate_x,
                x * scale_y + y * shear_y + translate_y,
                i,
            ] as FillerPoint,
    );
}

function affine_transform_rect(
    { x, y, width, height }: Tile,
    [[scale_x, shear_x, translate_x], [scale_y, shear_y, translate_y]]: Affine,
): Tile {
    const nx = x * scale_x + y * shear_x + translate_x;
    const ny = x * scale_y + y * shear_y + translate_y;
    const nwidth = (x + width) * scale_x + (y + height) * shear_x + translate_x - nx;
    const nheight = (x + width) * scale_y + (y + height) * shear_y + translate_y - ny;
    return { x: nx, y: ny, width: nwidth, height: nheight };
}

function split(P: FillerPoint[], pos: number): [FillerPoint[], FillerPoint[]] {
    const N = P.length;
    return [P.slice(0, pos), P.slice(pos, N)];
}

/**
 * Divides `size` between the two halves in proportion to their point counts.
 *
 * The two shares must add back up to `size`: the left half is pinned to the
 * region's left edge and the right half to its right edge, so shares that fall
 * short leave a gap down the middle that nothing fills.
 */
function sizes(Da: FillerPoint[], Db: FillerPoint[], size: number): Point {
    const pA = Da.length;
    const pB = Db.length;
    const total = pA + pB;
    return [(pA / total) * size, (pB / total) * size];
}

/** The region a run should fill: the data's bounding box unless one is given. */
function region_for(D: readonly Point[], parameters: NMapParameters): Bounds {
    return parameters.BB ?? get_bounds(D);
}

/**
 * mulberry32, seeded per call.
 *
 * The reference jitters its candidate positions with `Math.random`, which makes
 * a run unreproducible — the same input gave three different layouts across
 * five calls. The placement rule is unchanged; only the source of randomness
 * is. This is a deliberate deviation, so that a measurement can be repeated.
 */
function jitter(seed: number): () => number {
    let s = seed >>> 0;
    return () => {
        s = (s + 0x6d2b79f5) | 0;
        let t = Math.imul(s ^ (s >>> 15), 1 | s);
        t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
}

/**
 * Pads the point set with synthetic border points until it fills a square
 * grid, after https://github.com/sebastian-meier/nmap-squared.js
 *
 * This is the reference's "border" method: normalise into the region keeping
 * the data's aspect, lay candidates along the perimeter, and repeatedly take
 * whichever sits furthest from anything already placed. The fillers carry no
 * index and are dropped from the output, so they show as empty cells.
 *
 * Note that this pads the point *count* to a square number. It does not make
 * the cells square — their aspect follows the region, which keeps the data's
 * shape. See the README.
 */
function squared(D: FillerPoint[], { x: x0, y: y0, width, height }: Bounds): FillerPoint[] {
    const random = jitter(D.length * 2654435761);
    const BB2 = get_bounds(D as unknown as Point[]);
    let sx = width / BB2.width;
    let sy = height / BB2.height;
    if (sx < sy) {
        sy = sx;
    } else {
        sx = sy;
    }
    D = D.map(([x, y, i]) => [(x - x0) * sx + 1, (y - y0) * sy + 1, i] as FillerPoint);
    const grid_size = Math.ceil(Math.ceil(Math.sqrt(D.length)) / 4) * 4;
    const sq_amount = Math.pow(grid_size, 2);
    let sq_missing = sq_amount - D.length;
    const ep_num = 1 / sq_missing;

    // border method
    let extra: FillerPoint[] = [];
    // top and bottom border
    for (let x = 0; x <= grid_size - ep_num; x += ep_num) {
        extra.push([
            ((width * sx) / grid_size) * x + ep_num * random(),
            ep_num * random(),
            undefined,
        ] as FillerPoint);
        extra.push([
            ((width * sx) / grid_size) * x + ep_num * random(),
            height * sy + 1 + ep_num * random(),
            undefined,
        ] as FillerPoint);
    }

    // left and right border
    for (let y = 0; y <= grid_size - ep_num; y += ep_num) {
        extra.push([
            ep_num * random(),
            ((height * sy) / grid_size) * y + ep_num * random(),
            undefined,
        ] as FillerPoint);
        extra.push([
            width * sx + 1 + ep_num * random(),
            ((height * sy) / grid_size) * y + ep_num * random(),
            undefined,
        ] as FillerPoint);
    }

    while (sq_missing > 0 && extra.length > 0) {
        extra = calcDist(extra, D);
        // Take the candidate in the emptiest space, and take it *out* of the
        // pool: peeking would append the same point every iteration.
        D.push(extra.pop() as FillerPoint);
        --sq_missing;
    }

    return D;
}

/**
 * Ranks candidates in `A` by distance to the nearest point of `B`, farthest
 * last — so the tail of the list is the emptiest space.
 *
 * `dist` carries across calls on purpose: `B` only grows, so the nearest
 * distance only shrinks, and the running minimum stays correct.
 */
function calcDist(A: FillerPoint[], B: FillerPoint[]): FillerPoint[] {
    for (let i = 0; i < A.length; ++i) {
        const a = A[i];
        if (a.dist === undefined) a.dist = Infinity;
        for (let j = 0; j < B.length; ++j) {
            const b = B[j];
            let t_dist = Math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2);
            if (b[2] === undefined) {
                t_dist *= 2;
            }
            if (t_dist < a.dist) {
                a.dist = t_dist;
            }
        }
    }
    return A.sort((a, b) => (a.dist as number) - (b.dist as number));
}

/** Parameters for {@link gridify_nmap}. */
export interface NMapParameters {
    /** Region to fill. Defaults to the bounding box of the input. */
    BB?: Bounds | null;
    /** Pads the input with border points so it fills a square grid. */
    squared?: boolean;
}

/** A rectangle owned by one input point. */
export interface NMapCell extends Bounds {
    /** Index into the original input. */
    index: number;
}

/**
 * The rectangles N-Map assigns, one per input point, in input order. They tile
 * the region exactly. {@link gridify_nmap} reduces these to their centres.
 *
 * @throws If any subdivision has no spread along the axis it is cut on.
 *   Distinct points are not enough — collinear runs and integer lattices
 *   both trip it.
 */
export function nmap_regions(
    D: readonly Point[],
    parameters: NMapParameters = {},
): NMapCell[] {
    const G: Tile[] = [];
    const BB = region_for(D, parameters);
    let added_index: FillerPoint[] = D.map(([x, y], i) => [x, y, i] as FillerPoint);
    if (parameters.squared === true) {
        added_index = squared(added_index, BB);
    }
    nmap(G, added_index, { ...BB }, false);
    return G.filter((d) => d.d![2] !== undefined)
        .sort((a, b) => (a.d![2] as number) - (b.d![2] as number))
        .map(({ x, y, width, height }, i) => ({ x, y, width, height, index: i }));
}

/**
 * Removes overlaps with the space-filling N-Map treemap approach, recursively
 * subdividing the plotting area so each point owns a rectangle.
 *
 * @throws If any subdivision has no spread along the axis it is cut on.
 *   Distinct points are not enough — collinear runs and integer lattices
 *   both trip it.
 * * @param D - Input points.
 * @returns The centre of each point's rectangle, in input order.
 */
export function gridify_nmap(
    D: readonly Point[],
    parameters: NMapParameters = {},
): Point[] {
    const G: Tile[] = [];
    const BB = region_for(D, parameters);
    let added_index: FillerPoint[] = D.map(([x, y], i) => [x, y, i] as FillerPoint);
    if (parameters.squared === true) {
        added_index = squared(added_index, BB);
    }
    nmap(G, added_index, { ...BB }, false);
    return G.filter((d) => d.d![2] !== undefined)
        .sort((a, b) => (a.d![2] as number) - (b.d![2] as number))
        .map(({ x, y, width, height }) => [x + width / 2, y + height / 2] as Point);
}
