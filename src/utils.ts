import type { Accessor, Bounds, Point } from "./types.js";

/** Computes the axis-aligned bounding box of `A`. */
export function get_bounds(A: readonly Point[]): Bounds {
    let min_x = Infinity;
    let max_x = -Infinity;
    let min_y = Infinity;
    let max_y = -Infinity;

    A.forEach(([x, y]) => {
        min_x = Math.min(min_x, x);
        max_x = Math.max(max_x, x);
        min_y = Math.min(min_y, y);
        max_y = Math.max(max_y, y);
    });

    return {
        x: min_x,
        y: min_y,
        width: max_x - min_x,
        height: max_y - min_y,
    };
}

/**
 * Linearly remaps a point from one bounding box to another. Boxes are given as
 * `[x, y, width, height]`.
 */
export function remap(
    [x, y]: Point,
    [x1, y1, w1, h1]: [number, number, number, number],
    [x2, y2, w2, h2]: [number, number, number, number],
): Point {
    return [((x - x1) / w1) * w2 + x2, ((y - y1) / h1) * h2 + y2];
}

/** Options for {@link get_scales}. */
export interface ScaleOptions<T = Point> {
    /** Pads the shorter axis so the data's aspect ratio survives the mapping. */
    keep_aspect_ratio?: boolean;
    /** Rounds outputs to integers, as grid-cell indices require. */
    round?: boolean;
    /** Reads the x component of a datum. */
    x?: Accessor<T>;
    /** Reads the y component of a datum. */
    y?: Accessor<T>;
}

/**
 * Builds `[x_scale, y_scale]` mapping `data` onto the target ranges
 * `[[X_min, X_max], [Y_min, Y_max]]`.
 */
export function get_scales<T = Point>(
    data: readonly T[],
    [[X_min, X_max], [Y_min, Y_max]]: [Point, Point],
    {
        keep_aspect_ratio = false,
        round = true,
        x = (d) => (d as unknown as Point)[0],
        y = (d) => (d as unknown as Point)[1],
    }: ScaleOptions<T> = {},
): [Accessor<T>, Accessor<T>] {
    const X_span = X_max - X_min;
    const Y_span = Y_max - Y_min;
    let [x_min, x_max] = extent(data, x);
    let [y_min, y_max] = extent(data, y);
    const x_span = x_max - x_min;
    const y_span = y_max - y_min;

    if (keep_aspect_ratio) {
        let o: number;
        if (x_span > y_span) {
            o = (x_span - y_span) / 2;
            y_min -= o;
            y_max += o;
        } else {
            o = (y_span - x_span) / 2;
            x_min -= o;
            x_max += o;
        }
    }

    const f = round ? Math.round : (d: number) => d;
    const x_scale: Accessor<T> = (d) => f(((x(d) - x_min) / x_span) * X_span + X_min);
    const y_scale: Accessor<T> = (d) => f(((y(d) - y_min) / y_span) * Y_span + Y_min);
    return [x_scale, y_scale];
}

/** Returns `[min, max]` of `accessor` over `data`. */
export function extent<T>(data: Iterable<T>, accessor: Accessor<T>): Point {
    let min = Infinity;
    let max = -Infinity;
    for (const d of data) {
        const v = accessor(d);
        min = Math.min(min, v);
        max = Math.max(max, v);
    }
    return [min, max];
}

/** Euclidean distance between two points. */
export function distance([ax, ay]: Point, [bx, by]: Point): number {
    return Math.hypot(bx - ax, by - ay);
}
