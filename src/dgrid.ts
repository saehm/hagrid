import type { IndexedPoint, Point } from "./types.js";

/** A point assigned to grid row `i` and column `j`. */
interface Cell {
    d: IndexedPoint;
    i: number;
    j: number;
}

/**
 * Recursively bisects the point set and the `r × s` cell block together,
 * splitting along whichever axis has more cells to give.
 */
function dgrid(
    G: Cell[],
    P: IndexedPoint[],
    r: number,
    s: number,
    i = 0,
    j = 0,
): void {
    const N = P.length;
    if (N === 0) return;
    else if (N === 1) {
        G.push({ d: P[0], i, j });
    } else {
        if (r > s) {
            const r_half = Math.ceil(r / 2);
            const [P1, P2] = split(P, 1, r_half * s); // split by y
            dgrid(G, P1, r_half, s, i, j);
            dgrid(G, P2, r - r_half, s, i + r_half, j);
        } else {
            const s_half = Math.ceil(s / 2);
            const [P1, P2] = split(P, 0, s_half * r); // split by x
            dgrid(G, P1, r, s_half, i, j);
            dgrid(G, P2, r, s - s_half, i, j + s_half);
        }
    }
}

/** Sorts by `dim` and cuts the result at `pos`. */
function split(
    P: IndexedPoint[],
    dim: 0 | 1,
    pos: number,
): [IndexedPoint[], IndexedPoint[]] {
    const sorted = P.sort((a, b) => a[dim] - b[dim]);
    return [sorted.slice(0, pos), sorted.slice(pos)];
}

/**
 * Parameters for {@link gridify_dgrid}. Supply either an explicit `rows` and
 * `cols`, or an `aspect_ratio` from which both are derived.
 */
export type DGridParameters =
    | { rows: number; cols: number; aspect_ratio?: number }
    | { aspect_ratio: number; rows?: undefined; cols?: undefined };

/**
 * Removes overlaps by recursively bisecting the data into a dense grid, in the
 * manner of Hilasaca et al.'s Distance-preserving Grid. Fills the whole
 * visualization area, keeping the requested cell aspect ratio.
 *
 * @param D - Input points.
 * @returns One `[column, row]` cell per input point, in input order.
 */
export function gridify_dgrid(
    D: readonly Point[],
    parameters: DGridParameters,
): Point[] {
    const N = D.length;
    let rows: number;
    let cols: number;
    if (
        parameters != null &&
        "rows" in parameters &&
        "cols" in parameters &&
        parameters.rows !== undefined &&
        parameters.cols !== undefined
    ) {
        rows = parameters.rows;
        cols = parameters.cols;
    } else if (parameters != null && "aspect_ratio" in parameters) {
        rows = Math.floor(Math.sqrt(N * (parameters.aspect_ratio as number)));
        cols = Math.ceil(N / rows);
    } else {
        throw new Error("gridify_dgrid: expected either {rows, cols} or {aspect_ratio}");
    }
    const G: Cell[] = [];
    const added_index: IndexedPoint[] = D.map(([x, y], i) => [x, y, i]);
    dgrid(G, added_index, rows, cols);
    return G.sort((a, b) => a.d[2] - b.d[2]).map((g) => [g.j, g.i] as Point);
}
