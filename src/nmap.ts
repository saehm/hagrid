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
 * BUG (pre-existing, preserved for now): the second element uses `pA` where it
 * should use `pB`, so the split is not proportional unless the halves are
 * equal. Fixing it changes every nmap output, so it is deferred rather than
 * folded into the TypeScript port. See the Phase 2 notes.
 */
function sizes(Da: FillerPoint[], Db: FillerPoint[], size: number): Point {
    const pA = Da.length;
    const pB = Db.length;
    return [(pA / (pA + pB)) * size, (pA / (pA + pB)) * size];
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

    while (sq_missing > 0) {
        extra = calcDist(extra, D);
        D.push(extra[extra.length - 1]);
        --sq_missing;
    }

    return D;
}

/** Ranks candidates in `A` by distance to the nearest point of `B`. */
function calcDist(A: FillerPoint[], B: FillerPoint[]): FillerPoint[] {
    for (let i = 0; i < A.length; ++i) {
        const a = A[i];
        for (let j = 0; j < B.length; ++j) {
            const b = B[j];
            let t_dist = Math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2);
            if (b[2] === undefined) {
                t_dist *= 2;
            }
            if (t_dist < (a.dist as number)) {
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

/**
 * Removes overlaps with the space-filling N-Map treemap approach, recursively
 * subdividing the plotting area so each point owns a rectangle.
 *
 * @param D - Input points.
 * @returns The centre of each point's rectangle, in input order.
 */
export function gridify_nmap(
    D: readonly Point[],
    parameters: NMapParameters = {},
): Point[] {
    const G: Region[] = [];
    const BB = parameters.BB ?? get_bounds(D);
    let added_index: FillerPoint[] = D.map(([x, y], i) => [x, y, i] as FillerPoint);
    if (parameters.squared === true) {
        added_index = squared(added_index, BB);
    }
    nmap(G, added_index, { ...BB }, false);
    return G.filter((d) => d.d![2] !== undefined)
        .sort((a, b) => (a.d![2] as number) - (b.d![2] as number))
        .map(({ x, y, width, height }) => [x + width / 2, y + height / 2] as Point);
}
