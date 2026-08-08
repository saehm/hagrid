/**
 * Gridfit, after
 *
 *   Keim & Herrmann, "The Gridfit Algorithm: An Efficient and Effective
 *   Approach to Visualizing Large Amounts of Spatial Data",
 *   IEEE Visualization '98.
 *
 * The data space is partitioned hierarchically. At each level a node is cut
 * into four subregions, and the dividing lines are shifted off centre by the
 * smallest integer amount that gives every subregion at least as many grid
 * cells as it has points — the invariant `P(B) <= A(B)` of §3.2. Points are
 * assigned to the *equally sized* quadrants; only the cells they get are
 * adjusted.
 *
 * Partitioning and divider placement happen together, level by level — doing
 * them as separate passes leaves each node's children positioned against an
 * extent their parent no longer has.
 */
import type { Point } from "./types.js";

/**
 * A region as `[x_min, x_max, y_min, y_max]`, in grid cells and **half-open**:
 * it covers `x_min <= x < x_max`, so it holds exactly
 * `(x_max - x_min) * (y_max - y_min)` cells. Area and point containment must
 * use the same convention or regions come up a cell short.
 */
type Extent = [number, number, number, number];

/** `A(B)` of §3.2: how many grid cells a region holds. */
function area(x_min: number, x_max: number, y_min: number, y_max: number): number {
    return Math.max(0, x_max - x_min) * Math.max(0, y_max - y_min);
}

/**
 * The smallest shift of a dividing line that satisfies `fits`, searched
 * outward from the centre — `min{x ∈ ℤ | ...}` of §3.2.
 *
 * Both directions are tried, nearest first, so the line moves as little as
 * possible. Returns `null` if no position in `[lo, hi]` works.
 */
function shift(centre: number, lo: number, hi: number, fits: (m: number) => boolean): number | null {
    if (lo > hi) return null;
    const start = Math.min(Math.max(centre, lo), hi);
    if (fits(start)) return start;
    const reach = Math.max(start - lo, hi - start);
    for (let d = 1; d <= reach; ++d) {
        const up = start + d;
        if (up <= hi && fits(up)) return up;
        const down = start - d;
        if (down >= lo && fits(down)) return down;
    }
    return null;
}

/** Fills a region row by row; anything that will not fit is returned. */
function fill(
    indices: readonly number[],
    [x_min, x_max, y_min, y_max]: Extent,
    out: (Point | undefined)[],
): number[] {
    const leftovers: number[] = [];
    let k = 0;
    for (let x = x_min; x < x_max && k < indices.length; ++x) {
        for (let y = y_min; y < y_max && k < indices.length; ++y) {
            out[indices[k++]] = [x, y];
        }
    }
    for (; k < indices.length; ++k) leftovers.push(indices[k]);
    return leftovers;
}

/**
 * Recursively partitions `indices` into the region `E`, per §3.2.
 *
 * @returns Indices that could not be placed.
 */
function partition(
    data: readonly Point[],
    indices: readonly number[],
    E: Extent,
    out: (Point | undefined)[],
): number[] {
    const [x_min, x_max, y_min, y_max] = E;
    const count = indices.length;
    if (count === 0) return [];

    const w = x_max - x_min;
    const h = y_max - y_min;
    // Nothing left to divide, or nothing to gain by dividing.
    if (w <= 1 || h <= 1 || count === 1) return fill(indices, E, out);

    // Points belong to the *equally sized* quadrant, per the definition of P(B).
    const mx = x_min + Math.floor(w / 2);
    const my = y_min + Math.floor(h / 2);
    const q: number[][] = [[], [], [], []];
    for (const i of indices) {
        const [px, py] = data[i];
        const right = px >= mx ? 1 : 0;
        const top = py >= my ? 2 : 0;
        q[right + top].push(i);
    }
    const [P1, P2, P3, P4] = q.map((g) => g.length);

    // x^i: the left column must hold B1 and B3, the right column B2 and B4.
    const x = shift(
        mx,
        x_min + 1,
        x_max - 1,
        (m) => P1 + P3 <= area(x_min, m, y_min, y_max) && P2 + P4 <= area(m, x_max, y_min, y_max),
    );
    if (x === null) return fill(indices, E, out);

    // y1^i and y2^i: each column gets its own horizontal divider.
    const y1 = shift(
        my,
        y_min + 1,
        y_max - 1,
        (m) => P1 <= area(x_min, x, y_min, m) && P3 <= area(x_min, x, m, y_max),
    );
    const y2 = shift(
        my,
        y_min + 1,
        y_max - 1,
        (m) => P2 <= area(x, x_max, y_min, m) && P4 <= area(x, x_max, m, y_max),
    );
    if (y1 === null || y2 === null) return fill(indices, E, out);

    // E(B1..B4) of §3.2.
    const children: [number[], Extent][] = [
        [q[0], [x_min, x, y_min, y1]],
        [q[1], [x, x_max, y_min, y2]],
        [q[2], [x_min, x, y1, y_max]],
        [q[3], [x, x_max, y2, y_max]],
    ];

    const leftovers: number[] = [];
    for (const [group, extent] of children) {
        leftovers.push(...partition(data, group, extent, out));
    }
    return leftovers;
}

/**
 * Places `index` at the free cell nearest `[cx, cy]`, searching outward.
 *
 * §3.3: a straight-line partition is not always possible even when
 * `P(B) <= A(B)`, so the paper keeps the stragglers aside and places them at
 * the nearest unoccupied position in a second pass.
 */
function place_nearest(
    index: number,
    [cx, cy]: Point,
    E: Extent,
    taken: Set<number>,
    stride: number,
    out: (Point | undefined)[],
): boolean {
    const [x_min, x_max, y_min, y_max] = E;
    const reach = Math.max(x_max - x_min, y_max - y_min);
    for (let r = 0; r <= reach; ++r) {
        for (let dx = -r; dx <= r; ++dx) {
            for (let dy = -r; dy <= r; ++dy) {
                if (r > 0 && Math.abs(dx) !== r && Math.abs(dy) !== r) continue;
                const x = cx + dx;
                const y = cy + dy;
                if (x < x_min || x >= x_max || y < y_min || y >= y_max) continue;
                const key = (x - x_min) * stride + (y - y_min);
                if (taken.has(key)) continue;
                taken.add(key);
                out[index] = [x, y];
                return true;
            }
        }
    }
    return false;
}

/**
 * Removes overlaps by hierarchically partitioning the data space, shifting each
 * dividing line by the smallest amount that leaves every subregion enough grid
 * cells for its points.
 *
 * @param D - Input points, in grid-cell units.
 * @returns One grid cell per input point, in input order.
 */
export function gridify_gridfit(D: readonly Point[]): Point[] {
    const N = D.length;
    const out: (Point | undefined)[] = new Array(N);
    if (N === 0) return [];

    let x_min = Infinity;
    let x_max = -Infinity;
    let y_min = Infinity;
    let y_max = -Infinity;
    for (const [x, y] of D) {
        if (x < x_min) x_min = x;
        if (x > x_max) x_max = x;
        if (y < y_min) y_min = y;
        if (y > y_max) y_max = y;
    }
    x_min = Math.floor(x_min);
    y_min = Math.floor(y_min);
    x_max = Math.floor(x_max) + 1;
    y_max = Math.floor(y_max) + 1;

    // A(B0) >= P(B0) is the precondition for the whole recursion: the display
    // must have at least as many cells as there are points. Grow the region
    // outward, keeping its aspect, until it does.
    while (area(x_min, x_max, y_min, y_max) < N) {
        if (x_max - x_min <= y_max - y_min) x_max += 1;
        else y_max += 1;
    }

    const root: Extent = [x_min, x_max, y_min, y_max];
    const indices = Array.from({ length: N }, (_, i) => i);
    const leftovers = partition(D, indices, root, out);

    // Second pass for anything the partition could not seat (§3.3).
    const stride = y_max - y_min;
    const taken = new Set<number>();
    for (const cell of out) {
        if (cell) taken.add((cell[0] - x_min) * stride + (cell[1] - y_min));
    }
    for (const i of leftovers) {
        const [px, py] = D[i];
        const cx = Math.min(Math.max(Math.floor(px), x_min), x_max - 1);
        const cy = Math.min(Math.max(Math.floor(py), y_min), y_max - 1);
        place_nearest(i, [cx, cy], root, taken, stride, out);
    }

    return out.map((cell) => cell ?? [x_min, y_min]);
}
