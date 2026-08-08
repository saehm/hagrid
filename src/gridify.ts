import { get_bounds, remap } from "./utils.js";
import { gridify_hilbert } from "./hilbert.js";
import { gridify_gilbert } from "./gilbert.js";
import { gridify_gosper } from "./gosper.js";
import { gridify_dgrid } from "./dgrid.js";
import { gridify_nmap } from "./nmap.js";
import { gridify_cmds, type CMDSStep } from "./cmds.js";
import type { Point, Region } from "./types.js";

/** Overlap-removal methods reachable through {@link gridify}. */
export type GridifyMethod =
    | "hilbert"
    | "gilbert"
    | "gosper"
    | "dgrid"
    | "cmds"
    | "nmap";

/** Result of a method that produces a final layout directly. */
export interface GridifyResult {
    /** One grid position per input point, in input order. */
    positions: Point[];
    /** Wall-clock time of the layout, in milliseconds. */
    runtime: number;
}

/**
 * Result of `cmds`, which optimises iteratively rather than in one pass.
 *
 * `runtime` is not reported here: the work happens as `steps` is consumed, so
 * there is nothing meaningful to measure up front.
 */
export interface GridifyStepsResult {
    /** Successive layouts; the last one is the grid-snapped result. */
    steps: AsyncGenerator<CMDSStep, void, undefined>;
}

/** Parameters accepted by {@link gridify}, by method. */
export interface GridifyParameters {
    /** `hilbert`, `gosper`: extra curve levels beyond the fitted minimum. */
    pluslevel?: number;
    /** `hilbert`, `gilbert`: inflates the grid so cells outnumber points. */
    whitespace?: number;
    /** `gilbert`: explicit grid size. Derived from the data's aspect ratio if omitted. */
    width?: number;
    height?: number;
    /** `hilbert`: preserves the data's aspect ratio. */
    keep_aspect_ratio?: boolean;
    /** `gosper`: which way up the hexagons sit. Defaults to `"flat"`. */
    orientation?: import("./types.js").HexOrientation;
    /** `dgrid`: cell aspect ratio. Defaults to `1`. */
    aspect_ratio?: number;
    /** `dgrid`: explicit grid dimensions, overriding `aspect_ratio`. */
    rows?: number;
    cols?: number;
    /** `cmds`: boundary the layout is confined to. Defaults to the data's square bounds. */
    Gamma?: Region;
    /** `cmds`: starting node size. Defaults to a cell of the fitted grid, and shrinks if needed. */
    size?: number;
    /** `cmds`: weight of the boundary constraint. Defaults to `0.1`. */
    alpha?: number;
    /** `cmds`: maximum outer iterations. Defaults to `50`. */
    max_iter?: number;
    /** `nmap`: region to fill. Defaults to the data's bounding box. */
    BB?: import("./types.js").Bounds | null;
    /** `nmap`: pads the input so it fills a square grid. */
    squared?: boolean;
}

function timed(run: () => Point[]): GridifyResult {
    const start = performance.now();
    const positions = run();
    const end = performance.now();
    return { positions, runtime: end - start };
}

/**
 * Removes overlaps in a scatterplot by snapping points onto a grid.
 *
 * @param data - Input points as `[x, y]`.
 * @param method - Which technique to apply. `hilbert`, `gilbert` and `gosper`
 *   use space filling curves; `dgrid`, `cmds` and `nmap` are reimplementations
 *   of other published techniques.
 * @param parameters - Method-specific options; see {@link GridifyParameters}.
 *
 * @example
 * ```ts
 * const { positions } = gridify(data, "hilbert", { whitespace: 2 });
 * ```
 */
export function gridify(
    data: readonly Point[],
    method: "cmds",
    parameters?: GridifyParameters,
): GridifyStepsResult;
export function gridify(
    data: readonly Point[],
    method?: Exclude<GridifyMethod, "cmds">,
    parameters?: GridifyParameters,
): GridifyResult;
export function gridify(
    data: readonly Point[],
    method: GridifyMethod = "hilbert",
    parameters: GridifyParameters = {},
): GridifyResult | GridifyStepsResult {
    switch (method) {
        case "hilbert": {
            const pluslevel = parameters.pluslevel ?? 0;
            const whitespace = parameters.whitespace ?? 1;
            const level =
                Math.ceil(Math.log2(data.length * whitespace) / Math.log2(4)) + pluslevel;
            const size = Math.pow(2, level);
            const { x: ux, y: uy, width: w, height: h } = get_bounds(data);
            const D = data.map((d) => remap(d, [ux, uy, w, h], [0, 0, size, size]));
            return timed(() =>
                gridify_hilbert(D, {
                    level,
                    keep_aspect_ratio: parameters.keep_aspect_ratio,
                }),
            );
        }
        case "gilbert":
            // Unlike `hilbert` this keeps the grid proportional to the data,
            // so a wide dataset is not padded into a mostly-empty square.
            return timed(() =>
                gridify_gilbert(data, {
                    width: parameters.width,
                    height: parameters.height,
                    whitespace: parameters.whitespace,
                }),
            );
        case "gosper": {
            // Defaulted like hilbert's; previously this read an absent
            // `pluslevel` and produced a NaN level.
            const pluslevel = parameters.pluslevel ?? 0;
            const level = Math.ceil(Math.log2(data.length) / Math.log2(7) + 1) + pluslevel;
            return timed(() =>
                gridify_gosper(data, { level, orientation: parameters.orientation }),
            );
        }
        case "dgrid": {
            const p =
                parameters.rows !== undefined && parameters.cols !== undefined
                    ? { rows: parameters.rows, cols: parameters.cols }
                    : { aspect_ratio: parameters.aspect_ratio ?? 1 };
            return timed(() => gridify_dgrid(data, p));
        }
        case "cmds": {
            const Gamma =
                parameters.Gamma ??
                (() => {
                    const { x: min_x, y: min_y, width: wx, height: hy } = get_bounds(data);
                    const w = Math.max(wx, hy);
                    return { x1: min_x, x2: min_x + w, y1: min_y, y2: min_y + w };
                })();
            // `size` is now optional: left out, cmds derives it from the m x n
            // grid of the paper's §4.2 and shrinks it until the layout converges.
            return {
                steps: gridify_cmds(data, {
                    Gamma,
                    size: parameters.size,
                    alpha: parameters.alpha,
                    max_iter: parameters.max_iter,
                }),
            };
        }
        case "nmap":
            return timed(() =>
                gridify_nmap(data, { BB: parameters.BB, squared: parameters.squared }),
            );
        default:
            throw new Error(`gridify: unknown method ${JSON.stringify(method)}`);
    }
}
