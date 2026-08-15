<script>
    /**
     * Gridding as an *interaction* rather than a layout.
     *
     * Every other showcase grids the whole plot once and lives with the
     * displacement that costs. This one leaves the projection exactly as it is
     * — overlaps and all — and grids only what is under the cursor, into a
     * panel beside it. The layout you are reading never moves; the thing you
     * are inspecting becomes readable on demand.
     *
     * It is the original HaGrid demo's lens, rebuilt on this site's components.
     * The trade is different from the static case: displacement inside the lens
     * does not matter much, because the lens is a magnifier and not the map.
     */
    import * as hagrid from "@saehrimnir/hagrid";
    import ImageGlyph from "./ImageGlyph.svelte";

    let {
        /** The JSON from `build-images.mjs`. */
        data,
        width = 720,
        height = 520,
    } = $props();

    // BASE_URL arrives without a trailing slash; see ImagesShowcase.
    const BASE = import.meta.env.BASE_URL.replace(/\/+$/, "");
    const src = (name) => `${BASE}/${data.dir}/${name}.${data.extension}`;

    const PALETTE = [
        "#7aa2f7", "#e0af68", "#9ece6a", "#f7768e", "#bb9af7",
        "#7dcfff", "#ff9e64", "#73daca", "#c0caf5", "#e06c9f", "#a1b5d8",
    ];
    const COLOUR = Object.fromEntries(
        data.classes.map((c, i) => [c, PALETTE[i % PALETTE.length]]),
    );

    const PAD = 30;
    const THUMB = 11;
    const LENS = 260;

    let radius = $state(60);
    let method = $state("hilbert");
    let inside = $state(false);
    let pointer = $state({ x: 0, y: 0 });

    const METHODS = [
        ["hilbert", "Hilbert"],
        ["gilbert", "Gilbert"],
        ["gosper", "Gosper"],
        ["dgrid", "DGrid"],
    ];

    /** The projection, fitted to the plot box once. */
    const screen = $derived.by(() => {
        let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
        for (const [x, y] of data.points) {
            if (x < x0) x0 = x; if (x > x1) x1 = x;
            if (y < y0) y0 = y; if (y > y1) y1 = y;
        }
        const s = Math.min((width - 2 * PAD) / (x1 - x0 || 1), (height - 2 * PAD) / (y1 - y0 || 1));
        const ox = PAD + (width - 2 * PAD - (x1 - x0) * s) / 2 - x0 * s;
        const oy = PAD + (height - 2 * PAD - (y1 - y0) * s) / 2 - y0 * s;
        return data.points.map(([x, y]) => [x * s + ox, y * s + oy]);
    });

    /** Indices under the lens circle. */
    const picked = $derived.by(() => {
        if (!inside) return [];
        const out = [];
        for (let i = 0; i < screen.length; ++i) {
            const dx = screen[i][0] - pointer.x;
            const dy = screen[i][1] - pointer.y;
            if (dx * dx + dy * dy <= radius * radius) out.push(i);
        }
        return out;
    });

    /**
     * The selection, gridded. This is the whole idea: a handful of points, laid
     * out on demand, at interactive rates — the layout call is on the order of a
     * millisecond at these counts, so it can run on every pointer move.
     */
    const lens = $derived.by(() => {
        if (picked.length === 0) return null;
        // One point has nothing to lay out, and a degenerate set makes some
        // methods throw; show it as-is rather than failing under the cursor.
        if (picked.length === 1) return { cells: [[0, 0]], cols: 1, rows: 1, ms: 0 };
        const points = picked.map((i) => data.points[i]);
        const t0 = performance.now();
        let cells;
        try {
            cells = hagrid.gridify(points, method, { keep_aspect_ratio: true }).positions;
        } catch {
            return null;
        }
        const ms = performance.now() - t0;
        let cols = 0, rows = 0;
        for (const [cx, cy] of cells) {
            cols = Math.max(cols, cx + 1);
            rows = Math.max(rows, cy + 1);
        }
        // Gosper answers in the input's own space, not in cell indices, so the
        // cells are normalised to a 0..1 box and scaled by the panel instead.
        const lattice = method !== "gosper";
        return { cells, cols, rows, ms, lattice };
    });

    /** Places the panel in whichever corner the cursor is not heading into. */
    const panel = $derived({
        x: pointer.x < width / 2 ? pointer.x + radius + 12 : pointer.x - radius - 12 - LENS,
        y: pointer.y < height / 2 ? pointer.y + radius + 12 : pointer.y - radius - 12 - LENS,
    });

    /**
     * The lens shows things at a *readable* size, not at whatever size fills the
     * panel. Scaling the layout to the frame looks right at forty images and
     * absurd at three, where two thumbnails would each be half the panel wide.
     *
     * So the cell size is chosen first — as large as fits, capped — and the
     * layout is scaled to that pitch and centred. The block then grows with the
     * selection until it fills the panel, and stops.
     */
    const MAX_CELL = 44;

    const geometry = $derived.by(() => {
        if (!lens) return null;

        // Cells come back in cell indices for the lattice methods and in the
        // input's own coordinates for Gosper; dividing by the layout's own
        // pitch puts both in units of one cell.
        let pitch_data = Infinity;
        for (let i = 0; i < lens.cells.length; ++i) {
            for (let j = i + 1; j < lens.cells.length; ++j) {
                const d = Math.hypot(
                    lens.cells[i][0] - lens.cells[j][0],
                    lens.cells[i][1] - lens.cells[j][1],
                );
                if (d > 1e-9 && d < pitch_data) pitch_data = d;
            }
        }
        if (!Number.isFinite(pitch_data) || pitch_data <= 0) pitch_data = 1;

        let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
        for (const [x, y] of lens.cells) {
            if (x < x0) x0 = x; if (x > x1) x1 = x;
            if (y < y0) y0 = y; if (y > y1) y1 = y;
        }
        // How many cells the block spans, plus the one it starts on.
        const across = Math.max((x1 - x0) / pitch_data, (y1 - y0) / pitch_data) + 1;
        const step = Math.min(MAX_CELL, LENS / Math.max(1, across));
        const s = step / pitch_data;

        const ox = (LENS - (x1 - x0) * s) / 2;
        const oy = (LENS - (y1 - y0) * s) / 2;
        return {
            placed: lens.cells.map(([x, y]) => [(x - x0) * s + ox, (y - y0) * s + oy]),
            size: step / 2,
            shape: lens.lattice ? "square" : "hex",
        };
    });

    function move(event) {
        const rect = event.currentTarget.getBoundingClientRect();
        pointer = {
            x: ((event.clientX - rect.left) / rect.width) * width,
            y: ((event.clientY - rect.top) / rect.height) * height,
        };
        inside = true;
    }

    function wheel(event) {
        event.preventDefault();
        radius = Math.min(240, Math.max(20, radius - Math.sign(event.deltaY) * 8));
    }

    /** How many of the picked thumbnails overlap where they sit, unaided. */
    const collisions = $derived.by(() => {
        let n = 0;
        for (let a = 0; a < picked.length; ++a) {
            for (let b = a + 1; b < picked.length; ++b) {
                const dx = screen[picked[a]][0] - screen[picked[b]][0];
                const dy = screen[picked[a]][1] - screen[picked[b]][1];
                if (Math.abs(dx) < THUMB && Math.abs(dy) < THUMB) ++n;
            }
        }
        return n;
    });
</script>

<figure class="showcase hagrid-demo">
    <div class="controls">
        <label>
            Technique
            <select bind:value={method}>
                {#each METHODS as [value, label] (value)}<option {value}>{label}</option>{/each}
            </select>
        </label>
        <label class="range">
            Lens radius <span class="num">{radius}px</span>
            <input type="range" min="20" max="240" step="4" bind:value={radius} />
        </label>
        <span class="hint">Move the pointer over the plot; scroll to resize the lens.</span>
    </div>

    <div class="stage">
        <svg
            class="plot"
            viewBox="0 0 {width} {height}"
            role="img"
            aria-label="A projection of {data.points.length} images; hovering grids the neighbourhood under the pointer"
            onpointermove={move}
            onpointerleave={() => (inside = false)}
            onwheel={wheel}
        >
            {#each screen as [x, y], i}
                <g transform="translate({x} {y})">
                    <rect
                        x={-THUMB / 2 - 1} y={-THUMB / 2 - 1}
                        width={THUMB + 2} height={THUMB + 2}
                        fill={COLOUR[data.labels[i]]}
                    />
                    <image
                        href={src(data.names[i])}
                        x={-THUMB / 2} y={-THUMB / 2}
                        width={THUMB} height={THUMB}
                        preserveAspectRatio="xMidYMid slice"
                    />
                </g>
            {/each}

            {#if inside}
                <g pointer-events="none">
                    <circle
                        cx={pointer.x} cy={pointer.y} r={radius}
                        class="ring"
                    />
                    {#if lens && geometry}
                        <g transform="translate({panel.x} {panel.y})">
                            <rect
                                class="panel"
                                x="-8" y="-8"
                                width={LENS + 16} height={LENS + 16}
                                rx="4"
                            />
                            {#each picked as index, k}
                                {@const [gx, gy] = geometry.placed[k]}
                                <g transform="translate({gx} {gy})">
                                    <ImageGlyph
                                        src={src(data.names[index])}
                                        size={geometry.size}
                                        fill={COLOUR[data.labels[index]]}
                                        shape={geometry.shape}
                                        alt={data.labels[index]}
                                    />
                                </g>
                            {/each}
                        </g>
                    {/if}
                </g>
            {/if}
        </svg>
    </div>

    <figcaption>
        {#if lens}
            <strong>{picked.length}</strong> images under the lens ·
            <span class:bad={collisions > 0}>{collisions} overlapping in the plot</span> ·
            gridded in <strong>{lens.ms.toFixed(2)} ms</strong>
        {:else}
            <strong>{data.points.length} images</strong>, drawn where the projection put them —
            hover to grid a neighbourhood.
        {/if}
    </figcaption>
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
        align-items: end;
        gap: 0.5rem 1.1rem;
        padding: 0.7rem 0.8rem;
        border-bottom: 1px solid var(--sl-color-gray-5);
        background: var(--sl-color-gray-7, var(--sl-color-gray-6));
        font-size: var(--sl-text-xs);
    }
    .controls label {
        display: flex;
        flex-direction: column;
        gap: 0.2rem;
        color: var(--sl-color-gray-3);
    }
    .controls select {
        font: inherit;
        font-size: 0.78rem;
        padding: 0.15rem 0.3rem;
        color: var(--sl-color-white);
        background: var(--sl-color-black);
        border: 1px solid var(--sl-color-gray-5);
        border-radius: 0.25rem;
    }
    .range { min-width: 12rem; }
    .range input { width: 100%; accent-color: var(--sl-color-accent); }
    .num { color: var(--sl-color-white); font-variant-numeric: tabular-nums; }
    .hint { margin-left: auto; color: var(--sl-color-gray-4); }

    .stage { background: var(--sl-color-black); }
    .plot { display: block; width: 100%; height: auto; touch-action: none; }

    /* The ring is the only chrome on the plot, so it has to read against both
       a dense field of thumbnails and the empty background. */
    .ring {
        fill: none;
        stroke: var(--demo-accent, #ff3fa8);
        stroke-width: 2;
    }
    .panel {
        fill: color-mix(in srgb, var(--sl-color-black) 92%, transparent);
        stroke: var(--demo-accent, #ff3fa8);
        stroke-width: 1.5;
    }

    figcaption {
        padding: 0.55rem 0.8rem;
        border-top: 1px solid var(--sl-color-gray-5);
        background: var(--sl-color-gray-7, var(--sl-color-gray-6));
        font-size: var(--sl-text-xs);
        color: var(--sl-color-gray-3);
        font-variant-numeric: tabular-nums;
    }
    figcaption strong { color: var(--sl-color-white); }
    .bad { color: var(--sl-color-red); }
</style>
