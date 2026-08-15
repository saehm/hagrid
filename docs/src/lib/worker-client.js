import { unpack_grid, unpack_points } from "./serialise.js";

/**
 * Main-thread side of the layout worker.
 *
 * A worker cannot be interrupted once a synchronous run has started, so a
 * superseded request terminates the worker outright and starts a fresh one.
 * Letting the old run finish would burn a core for seconds on a result nobody
 * is going to look at — which on a dragged slider is every result but the last.
 */
export function create_layout_runner() {
    let worker = null;
    let next_id = 0;
    let active = null;

    const spawn = () => {
        worker = new Worker(new URL("./layout-worker.js", import.meta.url), {
            type: "module",
        });
        worker.onmessage = (event) => {
            const message = event.data;
            if (!active || message.id !== active.id) return;
            if (message.type === "progress") {
                active.on_progress(message.frames);
                return;
            }
            const job = active;
            active = null;
            job.resolve(message.error ? { error: message.error } : hydrate(message));
        };
        worker.onerror = (event) => {
            const job = active;
            active = null;
            job?.resolve({ error: event.message || "the layout worker failed" });
        };
    };

    return {
        /** Queues a layout, cancelling whatever was running. */
        run(spec, on_progress = () => {}) {
            if (active) {
                // Mid-run and uninterruptible: throw the worker away.
                worker.terminate();
                worker = null;
                active.resolve({ cancelled: true });
                active = null;
            }
            if (!worker) spawn();
            const id = ++next_id;
            return new Promise((resolve) => {
                active = { id, resolve, on_progress };
                worker.postMessage({ id, ...spec });
            });
        },
        dispose() {
            worker?.terminate();
            worker = null;
            active = null;
        },
    };
}

/**
 * Rebuilds a result from its flat buffers.
 *
 * Frames are left packed and materialised one at a time by `frame(i)`: a `cmds`
 * run holds hundreds of them and only ever draws one, so rebuilding them all up
 * front would allocate megabytes of pairs to throw away.
 */
function hydrate(message) {
    const { positions, stride, grids, count } = message;
    const grid_cache = new Array(count);
    return {
        error: null,
        source: unpack_points(message.source),
        runtime: message.runtime,
        facts: message.facts,
        approximate_source: message.approximate_source,
        count,
        frame(i) {
            const k = Math.min(Math.max(i, 0), count - 1);
            if (!grid_cache[k]) grid_cache[k] = unpack_grid(grids[k]);
            return {
                positions: unpack_points(positions, k * stride, stride / 2),
                grid: grid_cache[k],
            };
        },
    };
}
