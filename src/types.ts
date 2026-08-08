/** A 2D point, as `[x, y]`. */
export type Point = [number, number];

/** A point carrying its index in the input, as `[x, y, i]`. */
export type IndexedPoint = [number, number, number];

/** An inclusive `[min, max]` interval along one axis. */
export type Range = [number, number];

/** An axis-aligned bounding box. */
export interface Bounds {
    x: number;
    y: number;
    width: number;
    height: number;
}

/** Maps a datum to a scalar, used to read out one dimension of a point. */
export type Accessor<T = Point> = (d: T) => number;

/**
 * Which way up the hexagons sit: `"pointy"` has a vertex at the top,
 * `"flat"` an edge. The two are the same lattice turned by 30°, so switching
 * rotates the grid and leaves the cell ordering unchanged.
 */
export type HexOrientation = "pointy" | "flat";
