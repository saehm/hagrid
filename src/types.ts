/** A 2D point, as `[x, y]`. */
export type Point = [number, number];

/**
 * A 2D point carrying its index in the original input, as `[x, y, i]`.
 * Several methods reorder points internally and use the trailing index to
 * restore the caller's ordering before returning.
 */
export type IndexedPoint = [number, number, number];

/** An axis-aligned bounding box. */
export interface Bounds {
    x: number;
    y: number;
    width: number;
    height: number;
}

/** A rectangular region given by opposing corners, used to bound `cmds`. */
export interface Region {
    x1: number;
    x2: number;
    y1: number;
    y2: number;
}

/** Maps a datum to a scalar, used to read out one dimension of a point. */
export type Accessor<T = Point> = (d: T) => number;

/**
 * Which way up the hexagons sit.
 *
 * `"pointy"` puts a vertex at the top and neighbours at multiples of 60°;
 * `"flat"` puts an edge at the top and neighbours at 30° off those. The two
 * lattices are the same lattice turned by 30°, so switching rotates the whole
 * grid and leaves the curve's cell ordering identical.
 */
export type HexOrientation = "pointy" | "flat";
