import type { Accessor, Bounds, Point, Range } from "./types.js";

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
    [[X_min, X_max], [Y_min, Y_max]]: [Range, Range],
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
    let x_span = x_max - x_min;
    let y_span = y_max - y_min;

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
        // scale by the padded extent, or the padding only shifts the data
        x_span = x_max - x_min;
        y_span = y_max - y_min;
    }

    const f = round ? Math.round : (d: number) => d;
    // A zero span (all data on one line, or a single point) would divide 0 by
    // 0; centre that axis in its target range instead.
    const axis = (
        v: Accessor<T>,
        min: number,
        span: number,
        Min: number,
        Span: number,
    ): Accessor<T> =>
        span === 0
            ? () => f(Min + Span / 2)
            : (d) => f(((v(d) - min) / span) * Span + Min);
    return [axis(x, x_min, x_span, X_min, X_span), axis(y, y_min, y_span, Y_min, Y_span)];
}

/** Returns `[min, max]` of `accessor` over `data`. */
export function extent<T>(data: Iterable<T>, accessor: Accessor<T>): Range {
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
