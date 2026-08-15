/**
 * Runs every layout off the main thread.
 *
 * All seven methods are slow enough to be felt at the point counts the demos
 * allow — `nmap` at 3 000 points and `cmds` at any size block for whole seconds,
 * and even the fast ones stutter a slider drag. None of them can be chunked:
 * `gridify` is one synchronous call, and `cmds` yields from a generator that
 * cannot be paused mid-round. A worker is the only place they can run without
 * the page going unresponsive.
 *
 * The worker rebuilds the dataset from its seed rather than receiving it, so a
 * request is a few numbers rather than a point array.
 */
import { build } from "./datasets.js";
import { layout, run_cmds } from "./layouts.js";
import { pack_frames, pack_points, transfers_of } from "./serialise.js";

self.onmessage = (event) => {
    const { id, method, dataset, n, seed, params } = event.data;
    const data = build(dataset, n, seed);

    const result =
        method === "cmds"
            ? run_cmds(data, params, (frames) =>
                  self.postMessage({ type: "progress", id, frames }),
              )
            : layout(method, data, params);

    if (result.error) {
        self.postMessage({ type: "done", id, error: result.error });
        return;
    }

    const packed = pack_frames(result.frames);
    const payload = {
        type: "done",
        id,
        error: null,
        source: pack_points(result.source),
        positions: packed.positions,
        stride: packed.stride,
        grids: packed.grids,
        count: result.frames.length,
        runtime: result.runtime,
        facts: result.facts,
        approximate_source: result.approximate_source ?? false,
    };
    self.postMessage(payload, transfers_of(payload));
};
