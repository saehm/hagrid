/**
 * The two ways a layout can be coloured, and the legends that make them
 * readable.
 *
 * They answer different questions. *Origin* asks whether the layout preserved
 * the arrangement — smooth colour means neighbours stayed neighbours, scrambled
 * colour means the structure was lost. *Movement* asks who paid for it, which
 * is where a method's failures actually live: a respectable mean hides the
 * handful of points flung across the plot, and this puts them in plain sight.
 */

/** Extent of a point set, with one scale for both axes. */
function unit_frame(data) {
    let x0 = Infinity;
    let x1 = -Infinity;
    let y0 = Infinity;
    let y1 = -Infinity;
    for (const [x, y] of data) {
        if (x < x0) x0 = x;
        if (x > x1) x1 = x;
        if (y < y0) y0 = y;
        if (y > y1) y1 = y;
    }
    const w = x1 - x0 || 1;
    const h = y1 - y0 || 1;
    // A per-axis fit would turn the colour wheel into an ellipse and make a wide
    // dataset look as though it had directional structure it does not have.
    const s = 1 / Math.max(w, h);
    return { x0, y0, s, dx: (1 - w * s) / 2, dy: (1 - h * s) / 2 };
}

/** Hue by direction from the centre of the data, depth by distance from it. */
export function origin_colour(px, py, dark) {
    const hue = ((Math.atan2(py, px) * 180) / Math.PI + 360) % 360;
    const rad = Math.min(1, Math.hypot(px, py) / 0.7071);
    const sat = 25 + 55 * rad;
    const lig = dark ? 72 - 22 * rad : 42 + 14 * (1 - rad);
    return `hsl(${hue.toFixed(0)} ${sat.toFixed(0)}% ${lig.toFixed(0)}%)`;
}

export function by_origin(data, dark) {
    const f = unit_frame(data);
    return data.map(([x, y]) =>
        origin_colour((x - f.x0) * f.s + f.dx - 0.5, (y - f.y0) * f.s + f.dy - 0.5, dark),
    );
}

/*
 * A sequential ramp held at roughly constant lightness rather than running dark
 * to bright. Most points barely move, so the low end carries most of the marks —
 * a ramp that starts near black would make the typical case invisible on the
 * dark theme, and one that starts near white would do the same on the light one.
 * Holding the middle keeps every stop legible against both.
 *
 * Blue at the cool end, not green. The ramp used to run teal → sage → yellow →
 * orange → red, which put its two ends in the one channel that red-green colour
 * blindness cannot read: the old green stop and the red stop were ΔE 5.6 apart
 * under deuteranopia against 19.7 for normal vision, so "barely moved" and
 * "flung across the layout" were the same colour to roughly one reader in
 * twelve men. Blue against red differs in the blue-yellow channel, which both
 * protanopia and deuteranopia preserve. Every well-separated pair here clears
 * ΔE 13.5 simulated, and the tightest of them is still ΔE 15.9 to normal vision.
 *
 * The route from blue to red is the fiddly part, and it is why the blues lean
 * violet rather than cyan. A ramp is a path through colour space, and the
 * obvious path — cyan-blue straight to amber — passes through the achromatic
 * axis, so it shows a muddy grey-green band across its middle: chroma fell to
 * 0.025 with a green cast, and a tenth of the ramp read as grey. Walking the
 * hue circle instead, blue → violet → magenta → red, keeps the chroma up but
 * fails the thing this ramp exists for: blue and magenta collapse to ΔE 2.4
 * under protanopia, worse than the green it replaced.
 *
 * So the crossing is unavoidable — what is negotiable is where it lands. Both
 * ends are held slightly on the magenta side of neutral, and the path crosses
 * the blue-yellow axis at a dusty rose rather than at grey: minimum chroma
 * 0.059, no sample below the 0.045 that reads as washed out. The old ramp's
 * smoothness, without its blind spot.
 *
 * Deliberately independent of the site's brown and pink: this encodes a
 * quantity, and a scale that moved whenever the branding did would stop being
 * comparable between one figure and the next.
 */
const RAMP = [
    [74, 108, 196],
    [138, 134, 207],
    [199, 154, 180],
    [229, 169, 122],
    [238, 148, 85],
    [226, 96, 74],
];

/** Samples the movement ramp at `t` in [0, 1]. */
export function movement_colour(t) {
    const u = Math.max(0, Math.min(1, t)) * (RAMP.length - 1);
    const i = Math.min(RAMP.length - 2, Math.floor(u));
    const k = u - i;
    const a = RAMP[i];
    const b = RAMP[i + 1];
    const c = [0, 1, 2].map((j) => Math.round(a[j] + (b[j] - a[j]) * k));
    return `rgb(${c[0]} ${c[1]} ${c[2]})`;
}

/**
 * Colours by how far each point moved.
 *
 * The scale tops out at the 95th percentile rather than the maximum: a single
 * outlier ten times the typical move would otherwise flatten every other point
 * onto the first stop of the ramp. Anything beyond simply saturates, which is
 * also what the histogram's tail is telling you.
 */
export function by_movement(displacements, top) {
    const t = top > 0 ? top : 1;
    return displacements.map((d) => movement_colour(d / t));
}

/** Ramp stops for drawing the movement scale as a gradient. */
export function movement_stops(n = 24) {
    return Array.from({ length: n }, (_, i) => movement_colour(i / (n - 1)));
}

/**
 * Gradient stops for an axis that runs to `max` while the ramp saturates at
 * `cap` — the 95th percentile, as the point colouring does.
 *
 * The two have to agree or a bar and the points it counts would be different
 * colours. Stops are positioned against `cap` and the last one is held out to
 * the end of the axis, which is what saturation looks like on a gradient.
 */
export function ramp_gradient(cap, max, n = 9) {
    const top = max > 0 ? max : 1;
    const span = Math.min(1, (cap > 0 ? cap : top) / top);
    const stops = movement_stops(n);
    const out = stops.map((colour, i) => ({
        colour,
        offset: (i / (stops.length - 1)) * span,
    }));
    if (span < 1) out.push({ colour: stops[stops.length - 1], offset: 1 });
    return out;
}

/** Hue wheel samples for the origin legend, as `{x, y, colour}` in [0,1]². */
export function origin_swatches(steps = 9, dark = true) {
    const out = [];
    for (let i = 0; i < steps; ++i) {
        for (let j = 0; j < steps; ++j) {
            const px = (i + 0.5) / steps - 0.5;
            const py = (j + 0.5) / steps - 0.5;
            if (Math.hypot(px, py) > 0.5) continue;
            out.push({ x: px + 0.5, y: py + 0.5, colour: origin_colour(px, py, dark) });
        }
    }
    return out;
}
