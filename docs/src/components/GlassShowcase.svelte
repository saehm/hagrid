<script>
    /**
     * Glass samples drawn as their composition: eight oxides stacked to a fixed
     * height, because they are percentages of a whole rather than eight
     * independent numbers.
     */
    import GlyphShowcase from "./GlyphShowcase.svelte";
    import CompositionGlyph from "./CompositionGlyph.svelte";
    import { box_reach } from "../lib/glyph-reach.js";
    import glass from "../data/glass-2d.json";

    let { method = "as projected", ...rest } = $props();

    /*
     * Two encodings share this glyph, so they must not share a palette.
     *
     * The bands are an *ordered* stack — the same eight oxides in the same
     * order every time — so they take a single-hue ramp, light at the bottom to
     * dark at the top. Order is the thing to read, and a ramp shows order.
     *
     * The frame is a *nominal* class, so it takes categorical hues. Because the
     * bands are all one hue, any coloured frame is unambiguous against them —
     * which was the problem before, when both used the same categorical set and
     * a blue band sat inside a blue frame.
     */
    const OXIDES = glass.oxides;
    const OXIDE_NAMES = OXIDES.map((j) => glass.columns[j]);
    const COLOURS = [
        "#dbe7f3", "#bcd0e6", "#9db9d8", "#7ea2cb",
        "#5f8bbd", "#4874a4", "#365a80", "#25405c",
    ];

    const CLASS_FILL = {
        build_wind_float: "#f7768e",
        "build_wind_non-float": "#ff9e64",
        vehic_wind_float: "#e0af68",
        containers: "#9ece6a",
        tableware: "#bb9af7",
        headlamps: "#ff4dd2",
    };
    const CLASSES = [...new Set(glass.labels)];

    // A stacked column is 1.24 x size wide and 2 x size tall, so its corner is
    // the half-diagonal out.
    // A stacked column is 0.62 wide against its full height, so its corner sits
    // further out than its edges — but only the projection can put a neighbour
    // there.
    const REACH = box_reach(glass.points.length, 0.62, 1);

    const parts_of = (i) => OXIDES.map((j) => glass.values[i][j]);

    let show_frames = $state(true);

    /** The sample nearest each class's average composition. */
    const EXEMPLAR = CLASSES.map((name) => {
        const rows = glass.values.filter((_, i) => glass.labels[i] === name);
        const mean = OXIDES.map((j) => rows.reduce((s, r) => s + r[j], 0) / rows.length);
        let best = -1, bd = Infinity;
        glass.values.forEach((r, i) => {
            if (glass.labels[i] !== name) return;
            const d = OXIDES.reduce((s, j, k) => s + (r[j] - mean[k]) ** 2, 0);
            if (d < bd) { bd = d; best = i; }
        });
        return { name, index: best };
    });

    const pretty = (n) => n.replace(/'/g, "");
</script>

<GlyphShowcase points={glass.points} reach={REACH} {method} {...rest}>
    {#snippet glyph(i, size)}
        <CompositionGlyph
            parts={parts_of(i)}
            colours={COLOURS}
            {size}
            outline={show_frames ? CLASS_FILL[glass.labels[i]] : null}
        />
    {/snippet}

    {#snippet details(i)}
        <p class="tip-title" style="color: {CLASS_FILL[glass.labels[i]]}">
            {glass.labels[i].replace(/_/g, " ")}
        </p>
        <!-- Three explicit columns — swatch, name, amount — so the rows line up
             whether or not a row has a swatch, and long names cannot push the
             numbers around. -->
        <div class="parts">
            {#each OXIDES as j, k}
                {#if glass.values[i][j] > 0}
                    <span class="swatch" style="background: {COLOURS[k]}"></span>
                    <span class="part">{pretty(glass.columns[j])}</span>
                    <span class="amount">{glass.values[i][j].toFixed(2)}%</span>
                {/if}
            {/each}
            <span class="rule"></span>
            <span class="part">RI</span>
            <span class="amount">{glass.values[i][0].toFixed(4)}</span>
        </div>
    {/snippet}

    {#snippet legend()}
        <div class="legend">
            <p class="legend-title">Reading a sample</p>
            <div class="reading">
                <svg viewBox="-20 -22 40 44" aria-hidden="true">
                    <CompositionGlyph
                        parts={parts_of(EXEMPLAR[0].index)}
                        colours={COLOURS}
                        size={18}
                        outline={CLASS_FILL[EXEMPLAR[0].name]}
                    />
                </svg>
                <p>
                    The eight oxides, stacked bottom to top and always filling the same height.
                    They are percentages of one sample, so the total says nothing and only the
                    split does — which is why these are stacked rather than drawn as bars.
                </p>
            </div>
            <ul class="oxides">
                {#each OXIDE_NAMES as n, k}
                    <li><span class="swatch" style="background: {COLOURS[k]}"></span>{pretty(n)}</li>
                {/each}
            </ul>
        </div>

        <div class="legend">
            <p class="legend-title">Use</p>
            <ul class="classes">
                {#each EXEMPLAR as e}
                    <li>
                        <svg viewBox="-20 -22 40 44" aria-hidden="true">
                            <CompositionGlyph
                                parts={parts_of(e.index)}
                                colours={COLOURS}
                                size={18}
                                outline={CLASS_FILL[e.name]}
                            />
                        </svg>
                        <span style="color: {CLASS_FILL[e.name]}">{e.name.replace(/_/g, " ")}</span>
                    </li>
                {/each}
            </ul>
            <label class="check">
                <input type="checkbox" bind:checked={show_frames} />
                frame each glyph by use
            </label>
            <p class="foot">
                Turn the frames off and the classes are still visible — headlamps carry barium,
                containers carry more potassium. The composition is doing the work, not the label.
            </p>
        </div>
    {/snippet}
</GlyphShowcase>

<style>
    .legend-title {
        margin: 0 0 0.5rem;
        font-size: var(--sl-text-xs);
        text-transform: uppercase;
        letter-spacing: 0.08em;
        color: var(--sl-color-gray-3);
    }
    .reading { display: flex; gap: 0.9rem; align-items: flex-start; }
    .reading svg { width: 3.4rem; flex: none; }
    .reading p {
        margin: 0;
        font-size: var(--sl-text-xs);
        color: var(--sl-color-gray-2);
        line-height: 1.45;
    }
    .oxides {
        display: grid;
        grid-template-columns: repeat(auto-fill, minmax(4rem, 1fr));
        gap: 0.15rem 0.6rem;
        margin: 0.7rem 0 0;
        padding: 0;
        list-style: none;
        font-size: var(--sl-text-xs);
        color: var(--sl-color-gray-3);
    }
    .oxides li { display: flex; align-items: center; gap: 0.35rem; }
    .swatch {
        display: inline-block;
        width: 0.6rem;
        height: 0.6rem;
        border-radius: 2px;
        flex: none;
    }

    /* The tooltip readout: swatch, name, amount, each in its own column so the
       rows align whether or not a row carries a swatch. */
    .parts {
        display: grid;
        grid-template-columns: auto 1fr auto;
        align-items: center;
        gap: 0.15rem 0.45rem;
    }
    .parts .part { color: var(--sl-color-gray-3); }
    .parts .amount {
        text-align: right;
        white-space: nowrap;
        font-variant-numeric: tabular-nums;
        color: var(--sl-color-white);
    }
    .parts .rule {
        width: 0.6rem;
        height: 1px;
        background: var(--sl-color-gray-5);
    }
    .classes {
        display: grid;
        grid-template-columns: repeat(auto-fill, minmax(4.6rem, 1fr));
        gap: 0.5rem;
        margin: 0;
        padding: 0;
        list-style: none;
    }
    .classes li {
        display: grid;
        grid-template-rows: auto auto;
        justify-items: center;
        gap: 0.2rem;
        font-size: var(--sl-text-xs);
        text-align: center;
    }
    .classes svg { display: block; width: 100%; max-width: 2.6rem; height: auto; }
    .check {
        display: flex;
        align-items: center;
        gap: 0.35rem;
        margin-top: 0.7rem;
        font-size: var(--sl-text-xs);
        color: var(--sl-color-gray-2);
    }
    .check input { width: auto; }
    .foot { margin: 0.5rem 0 0; font-size: var(--sl-text-xs); color: var(--sl-color-gray-3); }
</style>
