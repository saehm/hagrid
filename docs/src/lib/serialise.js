/**
 * Packing for the trip across the worker boundary.
 *
 * Layouts are arrays of `[x, y]` pairs, which structured-clone handles but
 * slowly: every pair is a separate object to walk and rebuild. A `cmds` run of
 * 800 points over 400 rounds is 320 000 of them. Flattening to one Float64Array
 * per result makes it a single transferable buffer with no copy at all.
 */

/** Flattens `[[x, y], ...]` into `[x0, y0, x1, y1, ...]`. */
export function pack_points(points) {
    const out = new Float64Array(points.length * 2);
    for (let i = 0; i < points.length; ++i) {
        out[2 * i] = points[i][0];
        out[2 * i + 1] = points[i][1];
    }
    return out;
}

/** Rebuilds `count` pairs starting at `offset` floats into the buffer. */
export function unpack_points(buffer, offset = 0, count = buffer.length / 2) {
    const out = new Array(count);
    for (let i = 0; i < count; ++i) {
        out[i] = [buffer[offset + 2 * i], buffer[offset + 2 * i + 1]];
    }
    return out;
}

/** Flattens every frame's positions into one buffer, laid end to end. */
export function pack_frames(frames) {
    if (!frames.length) return { positions: new Float64Array(0), stride: 0, grids: [] };
    const stride = frames[0].positions.length * 2;
    const out = new Float64Array(stride * frames.length);
    frames.forEach((frame, f) => {
        const p = frame.positions;
        for (let i = 0; i < p.length; ++i) {
            out[f * stride + 2 * i] = p[i][0];
            out[f * stride + 2 * i + 1] = p[i][1];
        }
    });
    return { positions: out, stride, grids: frames.map((frame) => pack_grid(frame.grid)) };
}

/** Hexagon centres are the only grid field big enough to be worth packing. */
export function pack_grid(grid) {
    if (grid?.type === "hex") return { ...grid, centres: pack_points(grid.centres) };
    return grid;
}

export function unpack_grid(grid) {
    if (grid?.type === "hex") return { ...grid, centres: unpack_points(grid.centres) };
    return grid;
}

/** Every transferable buffer in a packed payload, for the second postMessage argument. */
export function transfers_of(payload) {
    const out = [];
    const take = (v) => {
        if (v instanceof Float64Array) out.push(v.buffer);
    };
    take(payload.source);
    take(payload.positions);
    for (const g of payload.grids ?? []) take(g?.centres);
    return out;
}
