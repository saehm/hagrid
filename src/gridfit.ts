import type { Point } from "./types.js";

/**
 * A region of the plane as `[x_min, x_max, y_min, y_max]`, measured in grid
 * cells rather than data units.
 */
export type Extent = [number, number, number, number];

/** A node of the quadtree: the indices it contains and the region it covers. */
export interface QuadNode {
    indices: number[];
    E: Extent;
    children?: QuadNode[];
    parent?: QuadNode;
}

/** A point placed at a concrete grid cell. */
export interface Placement {
    /** Index into the original input. */
    i: number;
    /** Assigned cell. */
    cell: Point;
}

/** The outcome of a {@link gridify_gridfit} run. */
export interface GridFitResult {
    /** Every non-empty quadtree node produced during subdivision. */
    nodes: QuadNode[];
    /** Points that were successfully placed. */
    final: Placement[];
    /** Points the fit could not place. */
    leftovers: number[];
    /** Set when the fit bailed out early, with the reason. */
    error: string | null;
}

function area([x_min, x_max, y_min, y_max]: Extent): number {
    return (x_max - x_min) * (y_max - y_min);
}

function extent_of(indices: readonly number[], data: readonly Point[]): Extent {
    let x_min = Infinity;
    let x_max = -Infinity;
    let y_min = Infinity;
    let y_max = -Infinity;
    for (const i of indices) {
        const [x, y] = data[i];
        if (x < x_min) x_min = x;
        if (x > x_max) x_max = x;
        if (y < y_min) y_min = y;
        if (y > y_max) y_max = y;
    }
    return [x_min, x_max, y_min, y_max];
}

/** Recursively subdivides `node` into four quadrants until cells are unit-sized. */
function quadtree(
    data: readonly Point[],
    node: QuadNode,
    [x_min, x_max, y_min, y_max]: Extent,
    nodes: QuadNode[],
): void {
    if (Math.abs(x_max - x_min) <= 1 || Math.abs(y_max - y_min) <= 1) {
        return;
    }

    const x_mid = Math.min(Math.max(x_min, Math.floor((x_max + x_min) / 2)), x_max - 1);
    const y_mid = Math.min(Math.max(y_min, Math.floor((y_max + y_min) / 2)), y_max - 1);

    const quadrants: Extent[] = [
        [x_min, x_mid, y_min, y_mid], // B1
        [x_mid, x_max, y_min, y_mid], // B2
        [x_min, x_mid, y_mid, y_max], // B3
        [x_mid, x_max, y_mid, y_max], // B4
    ];

    node.children = quadrants.map(([qx_min, qx_max, qy_min, qy_max], q) => {
        // The upper half of each axis is exclusive at the midpoint, so a point
        // on the boundary lands in exactly one quadrant.
        const low_x = q === 0 || q === 2;
        const low_y = q === 0 || q === 1;
        const indices = node.indices.filter((i) => {
            const [x, y] = data[i];
            const in_x = low_x ? x >= qx_min && x <= qx_max : x > qx_min && x <= qx_max;
            const in_y = low_y ? y >= qy_min && y <= qy_max : y > qy_min && y <= qy_max;
            return in_x && in_y;
        });
        return { indices, E: [qx_min, qx_max, qy_min, qy_max] as Extent, parent: node };
    });

    nodes.push(...node.children.filter((child) => child.indices.length > 0));

    node.children.forEach((child) => {
        if (child.indices.length > 0) quadtree(data, child, child.E, nodes);
    });
}

/**
 * Slides a dividing line until both sides have enough cells for their points.
 *
 * @returns The settled position, or `null` if no position satisfies both sides.
 */
function balance(
    start: number,
    lo: number,
    hi: number,
    P_low: number,
    P_high: number,
    area_low: (mid: number) => number,
    area_high: (mid: number) => number,
): number | null {
    let mid = start;
    const fits = (m: number) => P_low <= area_low(m) && P_high <= area_high(m);
    if (fits(mid)) return mid;

    const direction = P_low <= area_low(mid) ? -1 : P_high <= area_high(mid) ? 1 : 0;
    // direction 0 means neither side fits at any offset, so sliding is futile.
    if (direction === 0) return null;

    // The line can only occupy positions in [lo, hi], which bounds the walk and
    // rules out the runaway loop the original could fall into.
    for (let step = 0; step <= hi - lo; ++step) {
        mid += direction;
        if (mid < lo || mid > hi) return null;
        if (fits(mid)) return mid;
    }
    return null;
}

/** Assigns concrete cells to the points of a leaf node. */
function place_leaf(node: QuadNode, result: GridFitResult): void {
    const [x_min, x_max, y_min, y_max] = node.E;
    if (node.indices.length === 0) return;
    const remaining = node.indices.slice();
    const x_span = x_max - x_min;
    const y_span = y_max - y_min;

    for (let dx = 0; dx < x_span; ++dx) {
        for (let dy = 0; dy < y_span; ++dy) {
            if (remaining.length === 0) break;
            const i = remaining.pop() as number;
            result.final.push({ i, cell: [x_min + dx, y_min + dy] });
        }
    }
    if (remaining.length > 0) {
        result.leftovers.push(...remaining);
    }
}

/** Recursively fits each subtree into its share of the grid. */
function fit(node: QuadNode, result: GridFitResult): void {
    const [x_min, x_max, y_min, y_max] = node.E;
    if (!node.children) {
        place_leaf(node, result);
        return;
    }

    const y_mid = Math.floor((y_max + y_min) / 2);
    const [B1, B2, B3, B4] = node.children;
    const [P1, P2, P3, P4] = node.children.map((child) => child.indices.length);

    // Vertical divider: columns left of it must hold B1 and B3.
    const x_mid = balance(
        Math.floor((x_max + x_min) / 2),
        x_min,
        x_max,
        P1 + P3,
        P2 + P4,
        (m) => area([x_min, m, y_min, y_max]),
        (m) => area([m, x_max, y_min, y_max]),
    );
    if (x_mid === null) {
        result.error ??= "no vertical split leaves enough space for both halves";
        result.leftovers.push(...node.indices);
        return;
    }

    // Horizontal dividers, one per column.
    const y_1 = balance(
        y_mid, y_min, y_max, P1, P3,
        (m) => area([x_min, x_mid, y_min, m]),
        (m) => area([x_min, x_mid, m, y_max]),
    );
    const y_2 = balance(
        y_mid, y_min, y_max, P2, P4,
        (m) => area([x_mid, x_max, y_min, m]),
        (m) => area([x_mid, x_max, m, y_max]),
    );
    if (y_1 === null || y_2 === null) {
        result.error ??= "no horizontal split leaves enough space for both halves";
        result.leftovers.push(...node.indices);
        return;
    }

    // set new borders
    B1.E[1] = x_mid;
    B3.E[1] = x_mid;
    B2.E[0] = x_mid;
    B4.E[0] = x_mid;
    B1.E[3] = y_1;
    B3.E[2] = y_1;
    B2.E[3] = y_2;
    B4.E[2] = y_2;
    node.children.forEach((child) => fit(child, result));
}

/**
 * Removes overlaps by recursively fitting a quadtree of the data onto the grid,
 * sliding each subdivision line until both sides have room for their points.
 *
 * The original implementation referenced an unimported `d3` and could never
 * run; `d3.extent`/`d3.range` are now computed locally. The algorithm itself is
 * unchanged and still needs validation against the paper.
 *
 * @param D - Input points, in grid-cell units.
 */
export function gridify_gridfit(D: readonly Point[]): GridFitResult {
    const N = D.length;
    const result: GridFitResult = { nodes: [], final: [], leftovers: [], error: null };

    const root: QuadNode = { indices: Array.from({ length: N }, (_, i) => i), E: [0, 0, 0, 0] };
    root.E = extent_of(root.indices, D);
    quadtree(D, root, root.E, result.nodes);
    fit(root, result);

    return result;
}
