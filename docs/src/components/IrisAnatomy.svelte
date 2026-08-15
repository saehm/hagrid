<script>
    /**
     * A single petal or sepal with its two measurements marked, drawn from the
     * same geometry the glyph uses so the legend cannot drift from the plot.
     *
     * The blade runs from the end of the stem to the tip — that span is the
     * *length* — and the *width* is across it at its widest. Both dimensions
     * are set outside the shape with extension lines, so nothing is written
     * over the thing being measured.
     */
    import {
        PETAL_FILL,
        SEPAL_FILL,
        petal_path,
        sepal_path,
        arm_geometry,
        anatomy_frame,
    } from "../lib/iris-glyph.js";

    let {
        kind = "petal",
        /** A real row, so the proportions are honest. */
        row = [6.3, 3.3, 6.0, 2.5],
        species = "Iris-virginica",
        size = 5,
        /**
         * Shared with the sibling diagram so both draw to one scale and their
         * annotations line up. Computed per shape when omitted.
         */
        frame = null,
    } = $props();

    let petal = $derived(kind === "petal");
    let { len, wid, stem } = $derived(arm_geometry(kind, row, size));
    let d = $derived(petal ? petal_path(size, len, wid) : sepal_path(size, len, wid));
    let fill = $derived(petal ? PETAL_FILL[species] : SEPAL_FILL[species]);

    let box = $derived(frame ?? anatomy_frame(row, size));
    let dim_x = $derived(box.dim_x);      // length arrow, to the right
    let dim_y = $derived(box.dim_y);      // width arrow, below the tip
    let view = $derived(box.view);
    let waist = $derived(stem + len / 2);   // the blade's widest point

    // Markers are per-instance: two of these render side by side.
    let id = $derived(`arrow-${kind}`);
</script>

<svg viewBox={view} role="img" aria-label="{kind}, with its length and width marked">
    <defs>
        <!-- auto-start-reverse is what makes the tail arrowhead face outward -->
        <marker
            id={id} markerWidth="5" markerHeight="5" refX="5" refY="2.5"
            orient="auto-start-reverse" markerUnits="userSpaceOnUse"
        >
            <path d="M 0 0 L 5 2.5 L 0 5 z" fill="currentColor" />
        </marker>
    </defs>

    <path {d} {fill} stroke="rgba(0,0,0,.35)" stroke-width={size / 10} />
    <circle r={size / 2.5} fill="var(--sl-color-gray-4)" />

    <g class="dim">
        <!-- length: end of the stem to the tip -->
        <line class="ext" x1={0} y1={stem} x2={dim_x + 3} y2={stem} />
        <line class="ext" x1={0} y1={stem + len} x2={dim_x + 3} y2={stem + len} />
        <line
            x1={dim_x} y1={stem} x2={dim_x} y2={stem + len}
            marker-start="url(#{id})" marker-end="url(#{id})"
        />
        <text x={dim_x + 5} y={waist} dominant-baseline="middle">length</text>

        <!-- width: across the blade at its widest -->
        <line class="ext" x1={-wid / 2} y1={waist} x2={-wid / 2} y2={dim_y + 3} />
        <line class="ext" x1={wid / 2} y1={waist} x2={wid / 2} y2={dim_y + 3} />
        <line
            x1={-wid / 2} y1={dim_y} x2={wid / 2} y2={dim_y}
            marker-start="url(#{id})" marker-end="url(#{id})"
        />
        <text x={0} y={dim_y + 8} text-anchor="middle">width</text>
    </g>
</svg>

<style>
    svg { display: block; width: 100%; height: auto; }
    .dim { color: var(--sl-color-gray-3); }
    .dim line { stroke: currentColor; stroke-width: 0.9; }
    .dim line.ext { stroke-dasharray: 1.5 1.5; opacity: 0.5; }
    .dim text {
        fill: var(--sl-color-gray-2);
        font-size: 6.5px;
        font-family: var(--sl-font-system, sans-serif);
    }
</style>
