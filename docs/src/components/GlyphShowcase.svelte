<script>
    /**
     * Draws a 2D layout with a glyph at every point, and lets you move those
     * glyphs onto a grid with any of hagrid's techniques.
     *
     * Dataset-agnostic: supply the projected points, how far each glyph reaches
     * per unit of size, and a snippet that draws one. See `IrisShowcase` for a
     * worked example.
     *
     * The glyph size is computed, not searched for. A glyph's reach is
     * proportional to its size, so the closest pair of positions fixes the
     * largest glyph that cannot overlap its neighbour — and on a lattice that
     * distance is the cell pitch, a property of the technique rather than of
     * where the data happened to fall.
     */
    import * as hagrid from "@saehrimnir/hagrid";

    const AS_PROJECTED = "as projected";

    let {
        /** `[[x, y], …]`, already projected to two dimensions. */
        points,
        /**
         * Per point, how far its glyph reaches for `size = 1` — either an
         * array, or a function of the current layout for a glyph that changes
         * shape with the lattice. An image is a square on a square grid and a
         * hexagon on Gosper's, and the two reach differently.
         */
        reach,
        /** `(index, size, cell)` — draws one glyph, centred on the origin. */
        glyph,
        /** Optional legend, rendered under the readout. */
        legend = null,
        /** `(index)` — contents of the hover tooltip. Omit for no tooltip. */
        details = null,
        method = AS_PROJECTED,
        /** A readable glyph, which a projection typically cannot accommodate. */
        legible = 2,
        /** Pan and zoom the plot, with a minimap. Worth it when glyphs are dense. */
        zoom = false,
        width = 720,
        height = 460,
    } = $props();

    /**
     * `[value, label]`: the value is the library's own method key and must stay
     * lowercase, the label is the technique's name as it is written in prose.
     */
    const METHODS = [
        [AS_PROJECTED, "As projected"],
        ["hilbert", "Hilbert"],
        ["gilbert", "Gilbert"],
        ["gosper", "Gosper"],
        ["dgrid", "DGrid"],
        ["gridfit", "GridFit"],
        ["nmap", "NMap"],
    ];

    /**
     * The controls each technique exposes, mirroring `GridifyParameters`.
     * `pass` decides whether a value is forwarded — `0` means "unset" for the
     * explicit-grid options, which the library reads as a pinned grid if given.
     */
    const CONTROLS = {
        hilbert: [
            { key: "pluslevel", label: "Extra levels", min: 0, max: 3, step: 1, value: 0 },
            { key: "whitespace", label: "Whitespace (cells per point)", min: 1, max: 8, step: 0.1, value: 1 },
            { key: "keep_aspect_ratio", label: "Keep aspect ratio", type: "check", value: false },
        ],
        gilbert: [
            // Opens at slack rather than 1: Gilbert has to seat every point in
            // its own cell, and a grid with barely more cells than points sends
            // displaced glyphs a long way to find one. See the Gilbert page —
            // 2.4 captures most of the reduction available at any slack.
            { key: "whitespace", label: "Whitespace (cells per point)", min: 1, max: 8, step: 0.1, value: 2.4 },
            { key: "cols", label: "Columns (0 = derive)", min: 0, max: 40, step: 1, value: 0, pass: (v) => v > 0 },
            { key: "rows", label: "Rows (0 = derive)", min: 0, max: 40, step: 1, value: 0, pass: (v) => v > 0 },
        ],
        gosper: [
            { key: "pluslevel", label: "Extra levels", min: 0, max: 2, step: 1, value: 0 },
            { key: "whitespace", label: "Whitespace (cells per point)", min: 1, max: 8, step: 0.1, value: 1 },
            {
                key: "orientation", label: "Hexagon orientation", type: "select", value: "flat",
                options: [["flat", "flat-topped"], ["pointy", "pointy-topped"]],
            },
        ],
        dgrid: [
            { key: "aspect_ratio", label: "Cell aspect ratio", min: 0.25, max: 4, step: 0.25, value: 1 },
            { key: "cols", label: "Columns (0 = derive)", min: 0, max: 40, step: 1, value: 0, pass: (v) => v > 0 },
            { key: "rows", label: "Rows (0 = derive)", min: 0, max: 40, step: 1, value: 0, pass: (v) => v > 0 },
        ],
        nmap: [{ key: "squared", label: "Squared (uniform cells)", type: "check", value: true }],
        gridfit: [{ key: "__cells", label: "Cells per point", min: 1, max: 6, step: 0.5, value: 2 }],
        [AS_PROJECTED]: [],
    };

    let layout = $state(method);
    let params = $state(defaults(method));
    /** Off lets you overrule the computed size and watch it break. */
    let fit = $state(method !== AS_PROJECTED);
    let manual_size = $state(legible);

    function defaults(m) {
        return Object.fromEntries((CONTROLS[m] ?? []).map((c) => [c.key, c.value]));
    }
    function switch_to(m) {
        layout = m;
        params = defaults(m);
        // A projection has no cell to fit, so it opens at a size worth reading
        // — which is exactly where it falls apart.
        fit = m !== AS_PROJECTED;
        if (!fit) manual_size = legible;
    }

    let gridify_params = $derived.by(() => {
        const out = {};
        for (const c of CONTROLS[layout] ?? []) {
            if (c.key.startsWith("__")) continue;
            const v = params[c.key];
            if (c.pass && !c.pass(v)) continue;
            out[c.key] = v;
        }
        return out;
    });

    /** gridfit wants its input already in grid-cell units. */
    function to_grid(pts, per_point) {
        const side = Math.max(2, Math.ceil(Math.sqrt(pts.length * per_point)));
        let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
        for (const [x, y] of pts) {
            if (x < x0) x0 = x; if (x > x1) x1 = x;
            if (y < y0) y0 = y; if (y > y1) y1 = y;
        }
        return pts.map(([x, y]) => [
            Math.round(((x - x0) / (x1 - x0 || 1)) * (side - 1)),
            Math.round(((y - y0) / (y1 - y0 || 1)) * (side - 1)),
        ]);
    }

    let placed = $derived.by(() => {
        if (layout === AS_PROJECTED) return { points, error: null, ms: 0 };
        try {
            const t0 = performance.now();
            const input = layout === "gridfit" ? to_grid(points, params.__cells ?? 2) : points;
            const { positions } = hagrid.gridify(input, layout, gridify_params);
            return { points: positions, error: null, ms: performance.now() - t0 };
        } catch (e) {
            return { points, error: e.message, ms: 0 };
        }
    });

    const PAD = 26;

    /** Fits a point set into the plot area, keeping its aspect. */
    function fit_frame(pts) {
        let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
        for (const [x, y] of pts) {
            if (x < x0) x0 = x; if (x > x1) x1 = x;
            if (y < y0) y0 = y; if (y > y1) y1 = y;
        }
        const s = Math.min((width - 2 * PAD) / (x1 - x0 || 1), (height - 2 * PAD) / (y1 - y0 || 1));
        return {
            s,
            ox: PAD + (width - 2 * PAD - (x1 - x0) * s) / 2 - x0 * s,
            oy: PAD + (height - 2 * PAD - (y1 - y0) * s) / 2 - y0 * s,
        };
    }
    const apply = (f, pts) => pts.map(([x, y]) => [x * f.s + f.ox, y * f.s + f.oy]);

    let frame = $derived(fit_frame(placed.points));
    let screen = $derived(apply(frame, placed.points));

    /**
     * The input drawn into the same rectangle as the layout, so the two can be
     * compared. The lattice methods answer in grid indices, which share no
     * units with the projection, so both are fitted to the plot area rather
     * than one being expressed in the other's coordinates.
     */
    let source_screen = $derived(apply(fit_frame(points), points));

    const MAX_REACH = $derived(Math.max(...reach_now));

    /**
     * The lattice the current layout puts glyphs on, handed to the glyph
     * snippet so it can match its shape to the cell. Gosper is the only
     * hexagonal one, and it carries the orientation its cells were built with.
     */
    const cell = $derived({
        kind: layout === "gosper" ? "hex" : "square",
        orientation: params.orientation ?? "flat",
    });

    /** Resolved per layout, so a glyph may reach differently on each lattice. */
    const reach_now = $derived(typeof reach === "function" ? reach(layout, cell) : reach);

    /** Closest two positions — the cell pitch, on a lattice. */
    let pitch = $derived.by(() => {
        let best = Infinity;
        for (let i = 0; i < screen.length; ++i) {
            for (let j = i + 1; j < screen.length; ++j) {
                const d = Math.hypot(screen[i][0] - screen[j][0], screen[i][1] - screen[j][1]);
                if (d > 1e-9 && d < best) best = d;
            }
        }
        return Number.isFinite(best) ? best : 0;
    });

    // Every glyph inside a disc of half the pitch cannot reach its neighbour's,
    // so this is overlap-free by construction — nothing to search, nothing to
    // check afterwards.
    let cell_size = $derived(pitch / (2 * MAX_REACH));
    let glyph_size = $derived(fit ? cell_size : manual_size);

    /* ------------------------------------------------------------- zoom ---- */

    /*
     * A viewport over the plot, as scale plus translation. Glyph geometry is
     * left alone and the whole layer is transformed instead, so zooming costs
     * one attribute rather than a recomputation of every glyph — and the
     * computed glyph size, which is the point of this component, keeps meaning
     * what it says at zoom 1.
     */
    let view = $state({ k: 1, x: 0, y: 0 });
    const zoomed = $derived(view.k !== 1 || view.x !== 0 || view.y !== 0);
    const reset_view = () => (view = { k: 1, x: 0, y: 0 });

    // Layout changes rearrange everything under the viewport, so a pan held
    // across the change would leave you looking at empty space.
    $effect(() => {
        void layout;
        reset_view();
    });

    /** Plot coordinates from client coordinates, undoing the viewport. */
    function to_plot(event) {
        const box = event.currentTarget.getBoundingClientRect();
        const vx = ((event.clientX - box.left) / box.width) * width;
        const vy = ((event.clientY - box.top) / box.height) * height;
        return { vx: (vx - view.x) / view.k, vy: (vy - view.y) / view.k, box };
    }

    function on_wheel(event) {
        if (!zoom) return;
        event.preventDefault();
        const box = event.currentTarget.getBoundingClientRect();
        const px = ((event.clientX - box.left) / box.width) * width;
        const py = ((event.clientY - box.top) / box.height) * height;
        const k = Math.min(12, Math.max(1, view.k * (event.deltaY < 0 ? 1.15 : 1 / 1.15)));
        // Anchor on the pointer: the point under it stays under it.
        view = { k, x: px - ((px - view.x) / view.k) * k, y: py - ((py - view.y) / view.k) * k };
        clamp_view();
    }

    let dragging = $state(false);
    let drag_from = { x: 0, y: 0, vx: 0, vy: 0 };

    function on_down(event) {
        if (!zoom || view.k === 1) return;
        dragging = true;
        drag_from = { x: event.clientX, y: event.clientY, vx: view.x, vy: view.y };
        event.currentTarget.setPointerCapture?.(event.pointerId);
    }
    function on_up(event) {
        dragging = false;
        event.currentTarget.releasePointerCapture?.(event.pointerId);
    }
    function on_drag(event) {
        const box = event.currentTarget.getBoundingClientRect();
        const sx = width / box.width;
        view = {
            k: view.k,
            x: drag_from.vx + (event.clientX - drag_from.x) * sx,
            y: drag_from.vy + (event.clientY - drag_from.y) * sx,
        };
        clamp_view();
    }

    /** Keeps the plot from being dragged off its own frame. */
    function clamp_view() {
        const min_x = width - width * view.k;
        const min_y = height - height * view.k;
        view = {
            k: view.k,
            x: Math.min(0, Math.max(min_x, view.x)),
            y: Math.min(0, Math.max(min_y, view.y)),
        };
    }

    /** The minimap's frame, and the slice of it currently on screen. */
    const MAP = 96;
    const map_box = $derived({
        w: MAP,
        h: MAP * (height / width),
        vx: (-view.x / (width * view.k)) * MAP,
        vy: (-view.y / (height * view.k)) * (MAP * (height / width)),
        vw: MAP / view.k,
        vh: (MAP * (height / width)) / view.k,
    });

    /* ------------------------------------------------------------ hover ---- */

    /*
     * Hover is resolved by finding the nearest glyph to the pointer rather than
     * by putting a listener on each one: at a few hundred glyphs that is one
     * handler instead of hundreds, and it still works when the glyphs are only
     * a couple of pixels across — which is exactly the state the ungridded
     * layout is in, and the state worth being able to inspect.
     */
    let hover = $state(-1);
    let pointer = $state({ x: 0, y: 0 });

    function track(event) {
        if (dragging) return on_drag(event);
        const { vx, vy, box } = to_plot(event);

        let best = -1;
        let bd = Infinity;
        for (let i = 0; i < screen.length; ++i) {
            const d = Math.hypot(screen[i][0] - vx, screen[i][1] - vy);
            if (d < bd) { bd = d; best = i; }
        }
        // Generous, so tiny glyphs stay reachable, but not unbounded. The floor
        // is divided by the zoom so it stays about six *screen* pixels rather
        // than growing into a grab bag as you magnify.
        const slack = Math.max(reach_now[best] * glyph_size, 6 / view.k);
        hover = bd <= slack ? best : -1;
        pointer = { x: event.clientX - box.left, y: event.clientY - box.top };
    }

    let collisions = $derived.by(() => {
        let n = 0;
        for (let i = 0; i < screen.length; ++i) {
            const ri = reach_now[i] * glyph_size;
            for (let j = i + 1; j < screen.length; ++j) {
                const d = Math.hypot(screen[i][0] - screen[j][0], screen[i][1] - screen[j][1]);
                if (d < ri + reach_now[j] * glyph_size) ++n;
            }
        }
        return n;
    });

    /**
     * How far the average glyph had to move, as a share of the plot's
     * diagonal — the other side of the trade. Asking for more cells shrinks the
     * glyph but keeps points nearer where they started.
     *
     * Deliberately not measured in cells: that divides by the pitch, so it
     * climbs when cells shrink even as points move less on screen, which
     * inverts the very trade-off it is meant to show.
     */
    const DIAGONAL = Math.hypot(width - 2 * PAD, height - 2 * PAD);
    let displacement = $derived.by(() => {
        if (layout === AS_PROJECTED) return 0;
        let total = 0;
        for (let i = 0; i < screen.length; ++i) {
            total += Math.hypot(
                screen[i][0] - source_screen[i][0],
                screen[i][1] - source_screen[i][1],
            );
        }
        return total / screen.length / DIAGONAL;
    });
</script>

<!-- `hagrid-demo` opts into the markdown reset in custom.css; without it every
     element inside picks up Starlight's `* + *` top margin. -->
<figure class="showcase hagrid-demo">
    <div class="controls">
        <label>
            Layout
            <select value={layout} onchange={(e) => switch_to(e.currentTarget.value)}>
                {#each METHODS as [value, label] (value)}<option {value}>{label}</option>{/each}
            </select>
        </label>

        {#each CONTROLS[layout] ?? [] as c (c.key)}
            {#if c.type === "check"}
                <label class="check">
                    <input type="checkbox" bind:checked={params[c.key]} />
                    {c.label}
                </label>
            {:else if c.type === "select"}
                <label>
                    {c.label}
                    <select bind:value={params[c.key]}>
                        {#each c.options as [v, text]}<option value={v}>{text}</option>{/each}
                    </select>
                </label>
            {:else}
                <label class="range">
                    {c.label} <span class="num">{params[c.key]}</span>
                    <input type="range" min={c.min} max={c.max} step={c.step} bind:value={params[c.key]} />
                </label>
            {/if}
        {/each}
    </div>

    <div class="stage">
        <svg
            class="plot" class:grabbable={zoom && view.k > 1} class:grabbing={dragging}
            viewBox="0 0 {width} {height}" role="img"
            aria-label="Each data point drawn as a glyph"
            onpointermove={details || zoom ? track : undefined}
            onpointerleave={() => { hover = -1; dragging = false; }}
            onwheel={zoom ? on_wheel : undefined}
            onpointerdown={zoom ? on_down : undefined}
            onpointerup={zoom ? on_up : undefined}
        >
            <!-- One transform over the whole layer: the glyphs keep the
                 geometry the layout gave them, and the viewport moves. -->
            <g transform="translate({view.x} {view.y}) scale({view.k})">
                {#each screen as [x, y], i}
                    <g transform="translate({x} {y})">{@render glyph(i, glyph_size, cell)}</g>
                {/each}
                {#if details && hover >= 0}
                    <circle
                        cx={screen[hover][0]} cy={screen[hover][1]}
                        r={Math.max(reach_now[hover] * glyph_size, 4)}
                        class="halo" style="stroke-width: {1.5 / view.k}"
                    />
                {/if}
            </g>

            {#if zoom && zoomed}
                <!-- Minimap: the whole plot, with the slice you are looking at.
                     Only while zoomed — at 1:1 it would say nothing. -->
                <g class="minimap" transform="translate({width - MAP - 12} 12)" pointer-events="none">
                    <rect class="map-frame" width={map_box.w} height={map_box.h} rx="2" />
                    {#each screen as [x, y], i}
                        <circle
                            cx={(x / width) * map_box.w}
                            cy={(y / height) * map_box.h}
                            r="0.7"
                            class="map-dot"
                        />
                    {/each}
                    <rect
                        class="map-view"
                        x={map_box.vx} y={map_box.vy}
                        width={map_box.vw} height={map_box.vh}
                    />
                </g>
            {/if}
        </svg>

        {#if zoom && zoomed}
            <button type="button" class="reset" onclick={reset_view}>Reset view</button>
        {/if}

        {#if details && hover >= 0}
            <div
                class="tip"
                style="left: {pointer.x}px; top: {pointer.y}px;
                       --flip-x: {pointer.x > 0.62 * width ? -1 : 1};
                       --flip-y: {pointer.y > 0.62 * height ? -1 : 1}"
            >
                {@render details(hover)}
            </div>
        {/if}
    </div>

    <div class="controls under">
        <label class="range">
            Glyph size <span class="num">{glyph_size.toFixed(2)}</span>
            <input
                type="range" min="0.05" max="12" step="0.05"
                value={glyph_size}
                oninput={(e) => { fit = false; manual_size = +e.currentTarget.value; }}
            />
        </label>
        {#if layout !== AS_PROJECTED}
            <label class="check">
                <input type="checkbox" bind:checked={fit} />
                fit to cell
            </label>
        {/if}
    </div>

    <figcaption>
        {#if placed.error}
            <strong class="bad">{layout} refused this input:</strong> {placed.error}
        {:else}
            <strong>{points.length} glyphs</strong> ·
            cell {pitch.toFixed(1)} px → glyph <strong>{cell_size.toFixed(2)}</strong> ·
            <span class:bad={collisions > 0} class:good={collisions === 0}>
                {collisions} overlapping {collisions === 1 ? "pair" : "pairs"}
            </span>
            {#if layout !== AS_PROJECTED}· moved <strong>{(displacement * 100).toFixed(1)}%</strong> of the plot on average{/if}
            {#if placed.ms > 0}· {placed.ms.toFixed(1)} ms{/if}
        {/if}
    </figcaption>

    {#if legend}
        <div class="legends">{@render legend()}</div>
    {/if}
</figure>

<style>
    .showcase {
        margin: 1.75rem 0;
        border: 1px solid var(--sl-color-gray-5);
        border-radius: 0.5rem;
        overflow: hidden;
    }
    .controls {
        display: flex;
        flex-wrap: wrap;
        gap: 0.75rem 1rem;
        align-items: end;
        padding: 0.7rem 0.9rem;
        border-bottom: 1px solid var(--sl-color-gray-5);
        background: var(--sl-color-gray-6);
    }
    /* The size control sits under the plot it governs; the technique controls
       stay above, where they decide what is being drawn at all. */
    .controls.under {
        border-bottom: 0;
        border-top: 1px solid var(--sl-color-gray-5);
        background: transparent;
    }
    label {
        display: flex;
        flex-direction: column;
        gap: 0.2rem;
        font-size: var(--sl-text-xs);
        color: var(--sl-color-gray-2);
    }
    .range { min-width: 10rem; flex: 1 1 10rem; }
    .check { flex-direction: row; align-items: center; gap: 0.35rem; }
    .check input { width: auto; }
    .num { font-variant-numeric: tabular-nums; color: var(--sl-color-white); }
    select, input[type="range"] { width: 100%; }
    /* Only the plot gets the dark ground — legends sit on the panel. */
    .stage { position: relative; }
    .plot { display: block; width: 100%; height: auto; background: var(--sl-color-black); }
    /* Only once there is somewhere to drag to; at 1:1 a grab cursor promises
       something the plot will not do. */
    .plot.grabbable { cursor: grab; touch-action: none; }
    .plot.grabbing { cursor: grabbing; }

    .map-frame {
        fill: color-mix(in srgb, var(--sl-color-black) 82%, transparent);
        stroke: var(--sl-color-gray-5);
    }
    .map-dot { fill: var(--sl-color-gray-4); }
    /* The one saturated mark in the chrome, because its whole job is to say
       "you are here" against a field of grey dots. */
    .map-view {
        fill: color-mix(in srgb, var(--demo-accent, #ff3fa8) 18%, transparent);
        stroke: var(--demo-accent, #ff3fa8);
        stroke-width: 1;
    }

    .reset {
        position: absolute;
        right: 0.75rem;
        top: calc(12px + 96px * 0.64 + 0.5rem);
        font: inherit;
        font-size: var(--sl-text-xs);
        padding: 0.15rem 0.45rem;
        color: var(--sl-color-white);
        background: color-mix(in srgb, var(--sl-color-black) 82%, transparent);
        border: 1px solid var(--sl-color-gray-5);
        border-radius: 0.25rem;
        cursor: pointer;
    }
    .reset:hover { border-color: var(--demo-accent, #ff3fa8); }
    .halo {
        fill: none;
        stroke: var(--sl-color-white);
        stroke-width: 1.5;
        opacity: 0.8;
        pointer-events: none;
    }
    .tip {
        position: absolute;
        z-index: 3;
        pointer-events: none;
        /* Flips to whichever side of the pointer has room, so it never leaves
           the plot on the right or bottom edges. */
        transform: translate(
            calc(var(--flip-x) * 0.75rem + (var(--flip-x) - 1) * 50%),
            calc(var(--flip-y) * 0.75rem + (var(--flip-y) - 1) * 50%)
        );
        max-width: 15rem;
        padding: 0.4rem 0.55rem;
        border: 1px solid var(--sl-color-gray-5);
        border-radius: 0.35rem;
        background: var(--sl-color-gray-6);
        box-shadow: 0 4px 14px rgb(0 0 0 / 45%);
        font-size: var(--sl-text-xs);
        line-height: 1.4;
        color: var(--sl-color-gray-2);
    }
    /* Label takes the slack and wraps; the number stays compact and aligned.
       The other way round, one long label sets the column width and pushes
       every value off to the right. */
    .tip :global(dl) {
        margin: 0;
        display: grid;
        grid-template-columns: 1fr auto;
        gap: 0.1rem 0.6rem;
        align-items: baseline;
    }
    .tip :global(dt) { color: var(--sl-color-gray-3); min-width: 0; }
    .tip :global(dd) {
        margin: 0;
        text-align: right;
        white-space: nowrap;
        font-variant-numeric: tabular-nums;
        color: var(--sl-color-white);
    }
    .tip :global(.tip-title) { margin: 0 0 0.25rem; font-weight: 600; color: var(--sl-color-white); }
    figcaption {
        padding: 0.6rem 0.9rem;
        border-top: 1px solid var(--sl-color-gray-5);
        font-size: var(--sl-text-xs);
        color: var(--sl-color-gray-2);
    }
    .bad { color: var(--sl-color-red); }
    .good { color: var(--sl-color-green); }
    .legends {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(17rem, 1fr));
        gap: 1rem;
        padding: 0.9rem;
        border-top: 1px solid var(--sl-color-gray-5);
        background: var(--sl-color-gray-6);
    }
</style>
