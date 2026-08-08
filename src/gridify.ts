import { get_bounds, remap } from "./utils.js";
import { gridify_hilbert } from "./hilbert.js";
import { gridify_gilbert } from "./gilbert.js";
import { gridify_gridfit } from "./gridfit.js";
import { gridify_gosper } from "./gosper.js";
import { gridify_dgrid } from "./dgrid.js";
import { gridify_nmap } from "./nmap.js";
import { gridify_cmds, type CMDSStep } from "./cmds.js";
import type { Bounds, HexOrientation, Point } from "./types.js";

/** Overlap-removal methods reachable through {@link gridify}. */
export type GridifyMethod =
    | "hilbert"
    | "gilbert"
    | "gosper"
    | "gridfit"
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
    steps: Generator<CMDSStep, void, undefined>;
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
    orientation?: HexOrientation;
    /** `dgrid`: cell aspect ratio. Defaults to `1`. */
    aspect_ratio?: number;
    /** `dgrid`: explicit grid dimensions, overriding `aspect_ratio`. */
    rows?: number;
    cols?: number;
    /** `cmds`: boundary the layout is confined to. Defaults to the data's square bounds. */
    Gamma?: Bounds;
    /** `cmds`: starting node size. Defaults to a cell of the fitted grid, and shrinks if needed. */
    size?: number;
    /** `cmds`: weight of the boundary constraint. Defaults to `0.1`. */
    alpha?: number;
    /** `cmds`: rounds of Algorithm 1 before the grid is grown. */
    max_iter?: number;
    /** `cmds`: how many times the grid may grow before giving up. */
    max_rescales?: number;
    /** `cmds`: inner iterations per round, `N` in Algorithm 1. Defaults to `20`. */
    inner_iterations?: number;
    /** `nmap`: region to fill. Defaults to the data's bounding box. */
    BB?: Bounds | null;
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
            const pluslevel = parameters.pluslevel ?? 0;
            const whitespace = parameters.whitespace ?? 1;
            // A level-L curve has 7^L cells, so log7 of the point count is the
            // smallest level that fits. Each extra level costs 7x the cells,
            // so ask for slack with `whitespace` rather than `pluslevel`.
            const level =
                Math.ceil(Math.log2(data.length * whitespace) / Math.log2(7)) + pluslevel;
            return timed(() =>
                gridify_gosper(data, { level, orientation: parameters.orientation }),
            );
        }
        case "gridfit":
            // Expects grid-cell units; the region grows until it has at least
            // one cell per point.
            return timed(() => gridify_gridfit(data));
        case "dgrid":
            return timed(() =>
                gridify_dgrid(data, {
                    rows: parameters.rows,
                    cols: parameters.cols,
                    aspect_ratio: parameters.aspect_ratio,
                }),
            );
        case "cmds": {
            const Gamma =
                parameters.Gamma ??
                (() => {
                    const { x, y, width, height } = get_bounds(data);
                    const w = Math.max(width, height);
                    return { x, y, width: w, height: w };
                })();
            // Left out, `size` is derived from the m x n grid of §4.2 and
            // shrunk until the layout converges.
            return {
                steps: gridify_cmds(data, {
                    Gamma,
                    size: parameters.size,
                    alpha: parameters.alpha,
                    max_iter: parameters.max_iter,
                    max_rescales: parameters.max_rescales,
                    inner_iterations: parameters.inner_iterations,
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
