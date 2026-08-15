<script>
    import Scatterplot from "./Scatterplot.svelte";
    import Histogram from "./Histogram.svelte";
    import ColourLegend from "./ColourLegend.svelte";
    import { build, datasets, dataset_keys } from "../lib/datasets.js";
    import { by_movement, by_origin } from "../lib/colour.js";
    import {
        auto_grid,
        count_distinct,
        displacements,
        pitch_of,
        summarise,
    } from "../lib/layouts.js";
    import { create_layout_runner } from "../lib/worker-client.js";

    let {
        method: fixed_method = "hilbert",
        dataset = "clusters",
        datasets: allowed = dataset_keys,
        n = 600,
        // cmds is orders of magnitude slower than the rest, so its pages cap
        // this rather than letting the slider reach a count that hangs the tab.
        max_points = 3000,
        pick_method = false,
    } = $props();

    const METHODS = [
        ["hilbert", "Hilbert — square grid"],
        ["gilbert", "Gilbert — rectangular grid"],
        ["gosper", "Gosper — hexagonal"],
        ["gridfit", "GridFit — hierarchical"],
        ["dgrid", "DGrid — dense grid"],
        ["nmap", "NMap-squared — treemap"],
        ["cmds", "CMDS — constrained MDS"],
    ];

    /**
     * The slack the demos open on, for the methods that can actually use it.
     *
     * Both GridFit and Gilbert have to seat every point in its own cell, so a
     * grid with barely more cells than points leaves the placement no freedom
     * and points travel a long way to find somewhere to sit. Swept over the
     * datasets on this site at N = 400, mean displacement falls steeply up to
     * about 2 cells per point and then flattens: Gilbert on clusters goes 4.73
     * cells at 1.0 to 1.76 at 2.4, and only to 0.63 by 6.0. Across the datasets
     * 2.4 captures 72–97% of the reduction available at any slack at all, which
     * is the knee — past it you are spending empty cells for very little.
     *
     * Hilbert and Gosper are left at 1 deliberately. They only exist at whole
     * curve levels, so this dial does nothing at all until it crosses a power
     * of four (or seven) — measured identical from 1.0 to 2.4 — and a default
     * of 2.4 there would promise a trade the method cannot make.
     */
    const SLACK = 2.4;

    /**
     * A range input stepping by 0.1 lands on values like 1.7000000000000002,
     * which is true but not worth showing. Formats to the precision the step
     * implies, so an integer step still reads as an integer.
     */
    const show = (value, step = 1) => {
        const places = Math.max(0, Math.ceil(-Math.log10(step)));
        return Number(value).toFixed(places);
    };

    const PARAMS = {
        hilbert: [
            { key: "whitespace", label: "Whitespace (cells per point)", type: "range", min: 1, max: 8, step: 0.1, value: 1 },
            { key: "pluslevel", label: "Extra levels", type: "range", min: 0, max: 3, step: 1, value: 0 },
            { key: "keep_aspect_ratio", label: "Keep aspect", type: "check", value: false },
        ],
        gilbert: [
            { key: "whitespace", label: "Whitespace (cells per point)", type: "range", min: 1, max: 8, step: 0.1, value: SLACK },
            { key: "cols", label: "Columns", type: "number", min: 1, max: 20000, step: 1, value: 0 },
            { key: "rows", label: "Rows", type: "number", min: 1, max: 20000, step: 1, value: 0 },
            {
                type: "note",
                value: "Whitespace proposes a grid and columns and rows show what it chose; edit either to pin your own, and move whitespace to propose again. A grid with fewer cells than there are points fails.",
            },
        ],
        gosper: [
            { key: "whitespace", label: "Whitespace (cells per point)", type: "range", min: 1, max: 8, step: 0.1, value: 1 },
            { key: "pluslevel", label: "Extra levels", type: "range", min: 0, max: 3, step: 1, value: 0 },
            {
                key: "orientation",
                label: "Hexagons",
                type: "select",
                value: "flat",
                options: [
                    ["flat", "flat-topped"],
                    ["pointy", "pointy-topped"],
                ],
            },
        ],
        gridfit: [
            // Defaults to slack rather than 1, for the reason given at SLACK:
            // on uniform data mean displacement measured 10.6 cells at 1.0
            // against 1.2 at 2.4. Drag it down to 1 and watch it jump.
            { key: "whitespace", label: "Whitespace (cells per point)", type: "range", min: 1, max: 6, step: 0.1, value: SLACK },
        ],
        dgrid: [
            { key: "aspect_ratio", label: "Cell aspect", type: "range", min: 0.25, max: 4, step: 0.25, value: 1 },
        ],
        nmap: [
            {
                type: "note",
                value: "Always run with squared: true — the NMap-squared method, which pads the point set to a square number of points so the treemap has a full grid to subdivide. Empty cells are the fillers. It equalises cell area, not cell shape: the aspect follows the region, so cells are square only when the data is.",
            },
        ],
        cmds: [
            { key: "max_iter", label: "Rounds", type: "range", min: 5, max: 400, step: 5, value: 60 },
            { key: "alpha", label: "Boundary α", type: "range", min: 0.01, max: 1, step: 0.01, value: 0.1 },
            { key: "max_rescales", label: "Grid growths", type: "range", min: 0, max: 30, step: 1, value: 8 },
            {
                type: "note",
                value: "Rounds are the quality dial: growing the grid early shrinks the nodes to a few percent of the plot area, while letting it run reaches about 70%. Cost is linear in rounds, so keep the point count low.",
            },
        ],
    };

    const COLOUR_MODES = [
        ["origin", "Origin", "Colour each point by where it started"],
        ["movement", "Movement", "Colour each point by how far it moved"],
    ];

    const VIEWS = [
        ["origins", "Origins", "Show where each point started, joined to where it ended up"],
        ["grid", "Grid", "Draw the cell boundaries"],
    ];

    let method = $state(fixed_method);
    let current_dataset = $state(dataset);
    let count = $state(n);
    let seed = $state(1);
    // Origins on by default: the layout alone shows where points ended up, and
    // the question a reader arrives with is what it cost to put them there.
    let views = $state({ grid: true, origins: true });
    let colour_mode = $state("origin");

    const schema = $derived(PARAMS[method] ?? []);
    let params = $state(defaults(fixed_method));

    function defaults(m) {
        return Object.fromEntries(
            (PARAMS[m] ?? []).filter((p) => p.key).map((p) => [p.key, p.value]),
        );
    }

    // Switching method must reset the parameters with it, or a stale `whitespace`
    // from the previous method leaks into the next one's defaults.
    $effect(() => {
        void method;
        params = defaults(method);
    });

    const data = $derived(build(current_dataset, count, seed));

    /*
     * gilbert's columns and rows are shown filled in rather than blank, so the
     * grid the library picked is visible before you touch anything.
     *
     * The dependencies are only the data and the whitespace dial, so editing a
     * field does not re-run this and overwrite what was just typed — but moving
     * whitespace, or changing the dataset, proposes a fresh grid. That is the
     * intended escape hatch from a pinned grid that no longer fits.
     */
    $effect(() => {
        if (method !== "gilbert") return;
        const g = auto_grid(data, params.whitespace ?? 1);
        params.cols = g.cols;
        params.rows = g.rows;
    });

    /*
     * Every method runs in the worker, not just the slow ones. Keeping one path
     * means the fast methods cannot quietly regress into blocking the page when
     * someone raises a point-count cap, and a dragged slider stays smooth
     * because each new request cancels the last.
     */
    const runner = create_layout_runner();
    let result = $state({ pending: true });
    let busy = $state(false);
    let progress = $state(0);
    let frame_index = $state(0);

    $effect(() => {
        const spec = {
            method,
            dataset: current_dataset,
            n: count,
            seed,
            params: { ...params },
        };
        busy = true;
        progress = 0;
        stop_playing();
        runner
            .run(spec, (frames) => (progress = frames))
            .then((r) => {
                // A cancelled run has already been superseded; its successor owns
                // the state now, so touching anything here would flicker the
                // canvas back to the layout the user just moved away from.
                if (r.cancelled) return;
                result = r;
                frame_index = r.count ? r.count - 1 : 0;
                busy = false;
            });
    });

    $effect(() => runner.dispose);

    const frame_count = $derived(result.count ?? 0);
    const frame = $derived(
        result.frame ? result.frame(Math.min(frame_index, frame_count - 1)) : null,
    );

    const pitch = $derived(frame ? pitch_of(frame.grid) : 1);
    const disp = $derived(
        frame && result.source ? displacements(result.source, frame.positions, pitch) : [],
    );
    const stats = $derived(summarise(disp));
    const distinct = $derived(frame ? count_distinct(frame.positions) : 0);
    const clean = $derived(distinct === data.length);

    let dark = $state(true);
    // Origin colours depend only on the input, so they survive a whole cmds
    // replay untouched; movement colours are per frame by definition.
    const point_colours = $derived(
        colour_mode === "origin"
            ? by_origin(data, dark)
            : by_movement(disp, stats.p95 || stats.max),
    );

    $effect(() => {
        const root = document.documentElement;
        const read = () => (dark = (root.dataset.theme ?? "dark") !== "light");
        read();
        const observer = new MutationObserver(read);
        observer.observe(root, { attributes: true, attributeFilter: ["data-theme"] });
        return () => observer.disconnect();
    });

    let playing = $state(null);
    function stop_playing() {
        if (playing) clearInterval(playing);
        playing = null;
    }
    function toggle_play() {
        if (playing) return stop_playing();
        if (frame_index >= frame_count - 1) frame_index = 0;
        playing = setInterval(() => {
            if (frame_index >= frame_count - 1) return stop_playing();
            ++frame_index;
        }, 40);
    }
    $effect(() => stop_playing);

    const notes = $derived(schema.filter((p) => p.type === "note"));
    const controls = $derived(schema.filter((p) => p.type !== "note"));
    const fmt = (v) => (typeof v === "number" ? v.toLocaleString() : v);
</script>

<figure class="demo hagrid-demo">
    <!-- Three registers, in the order the questions come: what am I looking at,
         how is it laid out, how is it drawn. Grouping matters most on the
         technique pages, where a method's own parameters sit next to sliders
         that belong to the data and read as one undifferentiated wall. -->
    <div class="bar">
        <section class="grp" aria-label="Data">
            <p class="grp-title">Data</p>
            <div class="ctls">
                <label class="ctl">
                    <span class="lbl">Dataset</span>
                    <select bind:value={current_dataset}>
                        {#each allowed as key (key)}
                            <option value={key}>{datasets[key].label}</option>
                        {/each}
                    </select>
                </label>

                <label class="ctl">
                    <span class="lbl">Points <b>{count.toLocaleString()}</b></span>
                    <input type="range" min="20" max={max_points} step="20" bind:value={count} />
                </label>

                <label class="ctl">
                    <span class="lbl">Seed <b>{seed}</b></span>
                    <input type="range" min="1" max="40" step="1" bind:value={seed} />
                </label>
            </div>
        </section>

        {#if pick_method || controls.length}
            <section class="grp" aria-label="Technique">
                <p class="grp-title">Technique</p>
                <div class="ctls">
                    {#if pick_method}
                        <label class="ctl wide">
                            <span class="lbl">Method</span>
                            <select bind:value={method}>
                                {#each METHODS as [key, label] (key)}
                                    <option value={key}>{label}</option>
                                {/each}
                            </select>
                        </label>
                    {/if}

                    {#each controls as p (p.key)}
                        {#if p.type === "check"}
                            <label class="ctl">
                                <span class="lbl">{p.label}</span>
                                <span class="boxwrap"><input type="checkbox" bind:checked={params[p.key]} /></span>
                            </label>
                        {:else if p.type === "number"}
                            <label class="ctl">
                                <span class="lbl">{p.label}</span>
                                <input type="number" min={p.min} max={p.max} step={p.step} bind:value={params[p.key]} />
                            </label>
                        {:else if p.type === "select"}
                            <label class="ctl">
                                <span class="lbl">{p.label}</span>
                                <select bind:value={params[p.key]}>
                                    {#each p.options as [v, l] (v)}
                                        <option value={v}>{l}</option>
                                    {/each}
                                </select>
                            </label>
                        {:else}
                            <label class="ctl">
                                <span class="lbl">{p.label} <b>{show(params[p.key], p.step)}</b></span>
                                <input type="range" min={p.min} max={p.max} step={p.step} bind:value={params[p.key]} />
                            </label>
                        {/if}
                    {/each}
                </div>
            </section>
        {/if}

        <section class="grp" aria-label="Vis">
            <p class="grp-title">Vis</p>
            <div class="ctls">
                <div class="toggles">
                    <div class="group">
                        <span class="lbl">Colour</span>
                        <div class="seg" role="group" aria-label="Colour points by">
                            {#each COLOUR_MODES as [key, label, title] (key)}
                                <button
                                    type="button"
                                    {title}
                                    aria-pressed={colour_mode === key}
                                    onclick={() => (colour_mode = key)}>{label}</button
                                >
                            {/each}
                        </div>
                    </div>
                    <div class="group">
                        <span class="lbl">Show</span>
                        <div class="seg loose" role="group" aria-label="Display options">
                            {#each VIEWS as [key, label, title] (key)}
                                <button
                                    type="button"
                                    {title}
                                    aria-pressed={views[key]}
                                    onclick={() => (views[key] = !views[key])}>{label}</button
                                >
                            {/each}
                        </div>
                    </div>
                </div>
            </div>
        </section>
    </div>

    <div class="stage">
        {#if frame}
            <Scatterplot
                source={result.source}
                positions={frame.positions}
                grid={frame.grid}
                colours={point_colours}
                show_grid={views.grid}
                show_links={views.origins}
                show_source={views.origins}
            />
        {:else}
            <p class="empty">
                {#if busy}
                    Computing{progress ? ` — round ${progress}` : ""}…
                {:else if result.error}
                    No layout
                {:else}
                    …
                {/if}
            </p>
        {/if}

        {#if busy && frame}
            <!-- The previous layout stays on screen while the next one computes,
                 so the canvas never blanks mid-drag; this says why it is stale. -->
            <p class="working">Computing{progress ? ` — round ${progress}` : ""}…</p>
        {/if}

        {#if result.error}
            <p class="banner">{result.error}</p>
        {/if}

        {#if frame_count > 1}
            <div class="transport">
                <button type="button" class="play" onclick={toggle_play}
                    >{playing ? "Pause" : "Play"}</button
                >
                <input
                    type="range"
                    min="0"
                    max={frame_count - 1}
                    step="1"
                    bind:value={frame_index}
                    oninput={stop_playing}
                    aria-label="Optimisation round"
                />
                <span class="round">{frame_index + 1} / {frame_count}</span>
            </div>
        {/if}

    </div>

    <div class="footer">
        <!-- A layout that left cells empty seated two points in one cell, which
             is the one failure this whole widget exists to make visible. Red
             text alone was too easy to scan past, so the band itself turns. -->
        <div class="metrics" class:danger={!clean}>
            <span><b>{data.length.toLocaleString()}</b> points</span>
            {#each result.facts ?? [] as [k, v] (k)}
                <span><b>{fmt(v)}</b> {k}</span>
            {/each}
            <span class={clean ? "ok" : "bad"}>
                <b>{distinct.toLocaleString()}</b> of {data.length.toLocaleString()} cells filled
            </span>
            <span>
                <b>{stats.mean.toFixed(2)}</b> cells moved{result.approximate_source ? "*" : ""}
            </span>
            <span><b>{result.runtime != null ? result.runtime.toFixed(1) : "—"}</b> ms</span>
        </div>

        {#if frame}
            <ColourLegend {dark} />
        {/if}

        {#if disp.length}
            <Histogram values={disp} summary={stats} />
        {/if}

        <figcaption>
            <p>{datasets[current_dataset].note}</p>
            {#each notes as note (note.value)}
                <p>{note.value}</p>
            {/each}
            {#if result.approximate_source}
                <p>
                    * <code>dgrid</code> assigns cells by rank rather than by scaling the
                    coordinates, so the originals sit in a fair reference frame rather than an exact
                    preimage. Read the movement figure as indicative.
                </p>
            {/if}
        </figcaption>
    </div>
</figure>

<style>
    /* One box, three registers: controls, plot, readout. The border lives here
       rather than on each band so the corners round once and the dividers can
       be hairlines instead of doubled-up edges. */
    .demo {
        margin: 1.75rem 0;
        display: flex;
        flex-direction: column;
        border: 1px solid var(--sl-color-gray-5);
        border-radius: 0.5rem;
        overflow: hidden;
    }

    /* --- toolbar ---
       Stacked groups, each an auto-fit grid of its own. Every control is one
       cell of a shared column track, and every control body is exactly --ctl-h
       tall, so labels line up across the row and controls line up down the
       column no matter which method's parameters are on screen or how the row
       wraps. Groups own the horizontal rules, so a group that renders nothing
       leaves no stray divider behind. */
    .bar {
        --ctl-h: 1.9rem;
        display: flex;
        flex-direction: column;
        padding: 0.7rem 0.8rem;
        border-bottom: 1px solid var(--sl-color-gray-5);
        background: var(--sl-color-gray-7, var(--sl-color-gray-6));
    }
    .grp + .grp {
        margin-top: 0.6rem;
        padding-top: 0.6rem;
        border-top: 1px solid var(--sl-color-gray-6);
    }
    /* Quieter and smaller than a control's own label, so the group reads as a
       heading over the controls rather than as another control beside them. */
    .grp-title {
        margin: 0 0 0.4rem;
        font-size: 0.6rem;
        line-height: 1;
        font-weight: 700;
        letter-spacing: 0.1em;
        text-transform: uppercase;
        color: var(--sl-color-gray-4);
    }
    .ctls {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(8.5rem, 1fr));
        align-items: start;
        gap: 0.6rem 0.85rem;
    }
    .ctl {
        display: grid;
        grid-template-rows: 1rem var(--ctl-h);
        align-items: center;
        min-width: 0;
    }
    .ctl.wide {
        grid-column: span 2;
    }
    .lbl {
        font-size: 0.66rem;
        line-height: 1rem;
        letter-spacing: 0.04em;
        text-transform: uppercase;
        color: var(--sl-color-gray-3);
        white-space: nowrap;
        overflow: hidden;
        text-overflow: ellipsis;
    }
    .lbl b {
        color: var(--sl-color-white);
        font-variant-numeric: tabular-nums;
        text-transform: none;
        letter-spacing: 0;
    }
    select,
    input[type="number"] {
        width: 100%;
        min-width: 0;
        height: var(--ctl-h);
        box-sizing: border-box;
        font: inherit;
        font-size: 0.78rem;
        padding: 0 0.35rem;
        color: var(--sl-color-white);
        background: var(--sl-color-black);
        border: 1px solid var(--sl-color-gray-5);
        border-radius: 0.25rem;
    }
    input[type="range"] {
        display: block;
        width: 100%;
        min-width: 0;
        height: var(--ctl-h);
        box-sizing: border-box;
        margin: 0;
        padding: 0;
        accent-color: var(--sl-color-accent);
    }
    .boxwrap {
        display: flex;
        align-items: center;
        height: var(--ctl-h);
    }
    .boxwrap input {
        margin: 0;
        accent-color: var(--sl-color-accent);
    }

    /* --- view toggles ---
       These change how the layout is drawn, not what is computed, which is why
       they live in their own group. */
    .toggles {
        grid-column: 1 / -1;
        display: flex;
        flex-wrap: wrap;
        align-items: center;
        gap: 0.5rem 1.4rem;
    }
    .group {
        display: flex;
        align-items: center;
        gap: 0.45rem;
    }
    .seg {
        display: flex;
        align-items: center;
    }
    .seg button {
        display: inline-flex;
        align-items: center;
        box-sizing: border-box;
        font: inherit;
        font-size: 0.7rem;
        line-height: 1;
        height: 1.65rem;
        padding: 0 0.6rem;
        border: 1px solid var(--sl-color-gray-5);
        background: transparent;
        color: var(--sl-color-gray-2);
        cursor: pointer;
    }
    /* Joined: the colour modes are one choice, so they read as one control. */
    .seg button:first-child {
        border-radius: 0.25rem 0 0 0.25rem;
    }
    .seg button:last-child {
        border-radius: 0 0.25rem 0.25rem 0;
    }
    .seg button:not(:first-child) {
        border-left: 0;
    }
    /* Loose: these are independent switches, so they must not look like a
       one-of-many choice. */
    .seg.loose {
        gap: 0.3rem;
    }
    .seg.loose button {
        border-radius: 0.25rem;
        border-left: 1px solid var(--sl-color-gray-5);
    }
    /* Hotpink is reserved for what is currently on. It is the only saturated
       colour in the chrome, so a glance answers "what am I looking at". */
    .seg button[aria-pressed="true"] {
        background: var(--demo-accent);
        border-color: var(--demo-accent);
        /* Flips with the theme: white on the bright pink is only 3.2:1, but the
           near-black token gives 5.8:1, and in light mode the pink is dark
           enough that the same token resolves to near-white at 5.4:1. */
        color: var(--sl-color-black);
        font-weight: 600;
    }
    .seg button:hover:not([aria-pressed="true"]) {
        border-color: var(--demo-accent);
        color: var(--sl-color-white);
    }

    /* --- stage --- */
    .stage {
        position: relative;
    }
    .empty {
        display: grid;
        place-items: center;
        height: 320px;
        margin: 0;
        background: var(--sl-color-black);
        color: var(--sl-color-gray-3);
        font-size: 0.85rem;
    }
    .working {
        position: absolute;
        top: 0.5rem;
        right: 0.5rem;
        margin: 0;
        padding: 0.2rem 0.5rem;
        border-radius: 0.25rem;
        font-size: 0.68rem;
        font-variant-numeric: tabular-nums;
        color: var(--sl-color-gray-2);
        background: color-mix(in srgb, var(--sl-color-gray-6) 88%, transparent);
        border: 1px solid var(--sl-color-gray-5);
    }
    .banner {
        margin: 0;
        padding: 0.5rem 0.7rem;
        font-size: 0.74rem;
        line-height: 1.5;
        color: var(--sl-color-red);
        background: color-mix(in srgb, var(--sl-color-red) 12%, var(--sl-color-black));
        border-top: 1px solid color-mix(in srgb, var(--sl-color-red) 40%, transparent);
    }
    .transport {
        display: flex;
        align-items: center;
        gap: 0.6rem;
        padding: 0.4rem 0.6rem;
        background: var(--sl-color-black);
        border-top: 1px solid var(--sl-color-gray-5);
    }
    .transport input[type="range"] {
        flex: 1;
    }
    .play {
        font: inherit;
        font-size: 0.72rem;
        padding: 0.24rem 0.6rem;
        min-width: 3.6rem;
        border: 1px solid var(--sl-color-gray-5);
        border-radius: 0.25rem;
        background: transparent;
        color: var(--sl-color-white);
        cursor: pointer;
    }
    .play:hover {
        border-color: var(--sl-color-accent);
    }
    .round {
        font-size: 0.72rem;
        font-variant-numeric: tabular-nums;
        color: var(--sl-color-gray-3);
        min-width: 4rem;
        text-align: right;
    }

    /* --- readout --- */
    .footer {
        padding: 0.7rem 0.8rem;
        border-top: 1px solid var(--sl-color-gray-5);
        background: var(--sl-color-gray-7, var(--sl-color-gray-6));
    }
    .metrics {
        display: flex;
        flex-wrap: wrap;
        align-items: baseline;
        gap: 0.3rem 1rem;
        padding-bottom: 0.6rem;
        border-bottom: 1px solid var(--sl-color-gray-6);
        font-size: 0.74rem;
        color: var(--sl-color-gray-3);
    }
    .metrics b {
        font-weight: 600;
        font-variant-numeric: tabular-nums;
        color: var(--sl-color-white);
    }
    .metrics .ok b {
        color: var(--sl-color-green, #3fa06a);
    }
    .metrics .bad,
    .metrics .bad b {
        color: var(--sl-color-red);
    }
    /* Bled out to the footer's edges with negative margins so it reads as a
       band across the widget rather than a tinted paragraph inside it. */
    .metrics.danger {
        margin: -0.7rem -0.8rem 0;
        padding: 0.7rem 0.8rem;
        background: color-mix(in srgb, var(--sl-color-red) 15%, transparent);
        border-bottom-color: color-mix(in srgb, var(--sl-color-red) 45%, transparent);
        /* The muted label grey drops to about 3:1 on the wash, so the band
           steps up to body ink; the failing figure keeps the red. */
        color: var(--sl-color-gray-2);
    }
    figcaption {
        margin-top: 0.6rem;
        padding-top: 0.6rem;
        border-top: 1px solid var(--sl-color-gray-6);
        display: flex;
        flex-direction: column;
        gap: 0.35rem;
    }
    figcaption p {
        margin: 0;
        font-size: 0.75rem;
        line-height: 1.55;
        color: var(--sl-color-gray-3);
    }

    :is(button, select, input):focus-visible {
        outline: 2px solid var(--sl-color-accent);
        outline-offset: 1px;
    }

    @media (max-width: 560px) {
        .ctls {
            grid-template-columns: repeat(auto-fit, minmax(7rem, 1fr));
        }
    }
</style>
