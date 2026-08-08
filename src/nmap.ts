import { get_bounds } from "./utils.js";
import type { Bounds, IndexedPoint, Point } from "./types.js";

/**
 * A point that may be a synthetic filler added by {@link squared}, in which
 * case its index slot is `undefined`. `dist` caches the nearest-neighbour
 * distance used to rank candidates.
 */
type FillerPoint = [number, number, number | undefined] & { dist?: number };

/** A rectangle in the subdivision, once assigned, carries its point. */
type Region = Bounds & { d?: FillerPoint };

/** A 2×3 affine matrix, row-major, as `[[a, b, tx], [c, d, ty]]`. */
type Affine = [[number, number, number], [number, number, number]];

/**
 * Rejects a subdivision that cannot be scaled.
 *
 * The split boundary is the midpoint between the two halves along the split
 * axis. When every point in the region shares that coordinate the boundary
 * lands on the region's own edge, one side gets zero extent, and the affine
 * scale becomes `share / 0` — `NaN` if the share is zero too. That used to
 * propagate silently through every child region, so the layout came back full
 * of `null`s with no indication of why.
 *
 * Note this is not prevented by the points being distinct: four points on a
 * vertical line are all different yet share an `x`, and an integer lattice
 * trips it too.
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
function nmap(G: Region[], P: FillerPoint[], R: Region, direction = true): void {
    const N = P.length;
    const N_half = Math.floor(N / 2);
    if (N === 1) {
        R.d = P[0];
        G.push(R);
    } else {
        const dim = direction ? 0 : 1;
        P = P.sort((a, b) => b[dim] - a[dim]);
        const [Da, Db] = split(P, N_half);
        const Ra: Region = Object.assign({}, R);
        const Rb: Region = Object.assign({}, R);
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
    { x, y, width, height }: Region,
    [[scale_x, shear_x, translate_x], [scale_y, shear_y, translate_y]]: Affine,
): Region {
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
 * The two shares must add back up to `size`, because the left half is pinned to
 * the region's left edge and the right half to its right edge. The original
 * returned `pA`'s share twice, so whenever the halves were uneven — every split
 * of an odd number of points — they summed to less than `size` and left a gap
 * down the middle that nothing ever filled. The outer extent looked fine, which
 * is why this went unnoticed; the interior did not tile.
 */
function sizes(Da: FillerPoint[], Db: FillerPoint[], size: number): Point {
    const pA = Da.length;
    const pB = Db.length;
    const total = pA + pB;
    return [(pA / total) * size, (pB / total) * size];
}

/**
 * The region a run should fill: the data's bounding box, or the smallest
 * square containing it when `squared` is asked for.
 *
 * Padding the point set alone does not square anything — the layout's shape
 * follows the region it is given, so the region has to be square too.
 */
function region_for(D: readonly Point[], parameters: NMapParameters): Bounds {
    if (parameters.BB) return parameters.BB;
    const bb = get_bounds(D);
    if (parameters.squared !== true) return bb;
    const side = Math.max(bb.width, bb.height);
    return {
        x: bb.x + bb.width / 2 - side / 2,
        y: bb.y + bb.height / 2 - side / 2,
        width: side,
        height: side,
    };
}

/**
 * Pads the point set with synthetic border points until it fills a square
 * grid, after https://github.com/sebastian-meier/nmap-squared.js
 *
 * Note: uses `Math.random`, so results are not reproducible run to run.
 */
function squared(D: FillerPoint[], { x: x0, y: y0, width, height }: Bounds): FillerPoint[] {
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
            ((width * sx) / grid_size) * x + ep_num * Math.random(),
            ep_num * Math.random(),
            undefined,
        ] as FillerPoint);
        extra.push([
            ((width * sx) / grid_size) * x + ep_num * Math.random(),
            height * sy + 1 + ep_num * Math.random(),
            undefined,
        ] as FillerPoint);
    }

    // left and right border
    for (let y = 0; y <= grid_size - ep_num; y += ep_num) {
        extra.push([
            ep_num * Math.random(),
            ((height * sy) / grid_size) * y + ep_num * Math.random(),
            undefined,
        ] as FillerPoint);
        extra.push([
            width * sx + 1 + ep_num * Math.random(),
            ((height * sy) / grid_size) * y + ep_num * Math.random(),
            undefined,
        ] as FillerPoint);
    }

    while (sq_missing > 0 && extra.length > 0) {
        extra = calcDist(extra, D);
        // Take the candidate sitting in the emptiest space, and take it *out*
        // of the pool. Peeking at it instead — which is what this used to do —
        // appended the same point `sq_missing` times, and a pile of identical
        // points has no extent to subdivide, so the layout then failed.
        D.push(extra.pop() as FillerPoint);
        --sq_missing;
    }

    return D;
}

/**
 * Ranks candidates in `A` by distance to the nearest point of `B`, farthest
 * last — so the tail of the list is the emptiest space.
 *
 * `dist` starts at `Infinity`. Leaving it `undefined`, as this used to, made
 * every `t_dist < a.dist` comparison false, so no distance was ever recorded
 * and the comparator returned `NaN` — the sort did nothing at all.
 *
 * The value carries across calls on purpose: `B` only ever grows, so the
 * nearest distance can only shrink, and keeping the running minimum is both
 * correct and cheaper than starting over.
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
 * The rectangles N-Map assigns, one per input point, in input order.
 *
 * {@link gridify_nmap} reduces these to their centres; this exposes the cells
 * themselves, which is what you want to draw a treemap — and what makes the
 * tiling verifiable.
 *
 * Throws if any subdivision has no extent along the axis it is being cut on —
 * see {@link assert_splittable}. Distinct points are not enough to guarantee
 * this: collinear runs and integer lattices both trip it.
 */
export function nmap_regions(
    D: readonly Point[],
    parameters: NMapParameters = {},
): NMapCell[] {
    const G: Region[] = [];
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
 * Throws if any subdivision has no extent along the axis it is being cut on —
 * see {@link assert_splittable}. Distinct points are not enough to guarantee
 * this: collinear runs and integer lattices both trip it.
 * * @param D - Input points.
 * @returns The centre of each point's rectangle, in input order.
 */
export function gridify_nmap(
    D: readonly Point[],
    parameters: NMapParameters = {},
): Point[] {
    const G: Region[] = [];
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
