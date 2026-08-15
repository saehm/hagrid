<script>
    /**
     * Any dataset produced by `scripts/build-dataset.mjs`, drawn with a
     * general-purpose glyph and wired into {@link GlyphShowcase}.
     *
     * Adding a showcase is then: generate the JSON, import it, pick a glyph.
     * Bespoke glyphs — like the iris flower — get their own wrapper instead.
     */
    import GlyphShowcase from "./GlyphShowcase.svelte";
    import StarGlyph from "./StarGlyph.svelte";
    import BarGlyph from "./BarGlyph.svelte";
    import { box_reach } from "../lib/glyph-reach.js";

    let {
        /** The JSON from `build-dataset.mjs`. */
        data,
        /** `"star"` for a profile silhouette, `"bars"` for a small multiple. */
        glyph_kind = "star",
        method = "as projected",
        ...rest
    } = $props();

    /**
     * Colours by class. Chosen for hue separation at small sizes rather than
     * for prettiness — a glyph a few pixels across is mostly one flat colour,
     * so that colour has to carry the class on its own.
     */
    const PALETTE = [
        "#7aa2f7", "#e0af68", "#9ece6a", "#f7768e", "#bb9af7",
        "#7dcfff", "#ff9e64", "#73daca", "#c0caf5", "#e06c9f", "#a1b5d8",
    ];
    const CLASSES = [...new Set(data.labels)];
    const COLOUR = Object.fromEntries(CLASSES.map((c, i) => [c, PALETTE[i % PALETTE.length]]));

    // A star fills a circle of radius `size`, so it reaches the same in every
    // direction. Bars fill a square, which reaches its corner only where a
    // neighbour can be diagonal — see `box_reach`.
    const REACH =
        glyph_kind === "bars"
            ? box_reach(data.points.length)
            : data.points.map(() => 1);

    const Glyph = $derived(glyph_kind === "bars" ? BarGlyph : StarGlyph);

    // One representative item per class for the legend: the one closest to that
    // class's average profile, so it looks like the class rather than like an
    // outlier.
    const EXEMPLAR = CLASSES.map((c) => {
        const rows = data.values.filter((_, i) => data.labels[i] === c);
        const mean = data.columns.map((_, j) => rows.reduce((s, r) => s + r[j], 0) / rows.length);
        let best = -1, bd = Infinity;
        data.values.forEach((r, i) => {
            if (data.labels[i] !== c) return;
            const d = r.reduce((s, v, j) => s + (v - mean[j]) ** 2, 0);
            if (d < bd) { bd = d; best = i; }
        });
        return { label: c, index: best };
    });

    const LEGEND_SIZE = 17;
</script>

<GlyphShowcase points={data.points} reach={REACH} {method} {...rest}>
    {#snippet glyph(i, size)}
        <Glyph
            values={data.values[i]}
            ranges={data.ranges}
            {size}
            fill={COLOUR[data.labels[i]]}
        />
    {/snippet}

    {#snippet details(i)}
        <p class="tip-title" style="color: {COLOUR[data.labels[i]]}">{data.labels[i]}</p>
        <dl>
            {#each data.columns as c, j}
                <dt>{c}</dt>
                <dd>{data.values[i][j].toLocaleString(undefined, { maximumFractionDigits: 2 })}</dd>
            {/each}
        </dl>
    {/snippet}

    {#snippet legend()}
        <div class="legend">
            <p class="legend-title">Reading a glyph</p>
            <div class="reading">
                <svg viewBox="-22 -22 44 44" aria-hidden="true">
                    <Glyph values={data.values[EXEMPLAR[0].index]} ranges={data.ranges} size={LEGEND_SIZE} fill={COLOUR[EXEMPLAR[0].label]} guides />
                </svg>
                <p>
                    {#if glyph_kind === "bars"}
                        One bar per measurement, left to right in the order below, sharing a
                        baseline and a scale.
                    {:else}
                        One spoke per measurement, clockwise from the top in the order below.
                        The outline is the item's profile.
                    {/if}
                    Each is scaled against that column's range across the whole dataset, so two
                    glyphs are comparable.
                </p>
            </div>
            <ol class="columns">
                {#each data.columns as c, j}
                    <li><span class="idx">{j + 1}</span><span class="name">{c}</span></li>
                {/each}
            </ol>
        </div>

        <div class="legend">
            <p class="legend-title">Classes</p>
            <ul class="classes">
                {#each EXEMPLAR as e}
                    <li>
                        <svg viewBox="-22 -22 44 44" aria-hidden="true">
                            <Glyph values={data.values[e.index]} ranges={data.ranges} size={LEGEND_SIZE} fill={COLOUR[e.label]} guides />
                        </svg>
                        <span>{e.label}</span>
                    </li>
                {/each}
            </ul>
            <p class="foot">
                The item closest to each class's average profile — so it looks like the class
                rather than like an outlier.
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
    .reading { display: flex; gap: 0.8rem; align-items: flex-start; }
    .reading svg { width: 4.6rem; flex: none; color: var(--sl-color-gray-2); }
    .reading p {
        margin: 0;
        font-size: var(--sl-text-xs);
        color: var(--sl-color-gray-2);
        line-height: 1.45;
    }
    /* A grid rather than multi-column text: the numbers stay aligned and a
       long name wraps inside its own cell instead of shifting the rest. */
    .columns {
        display: grid;
        grid-template-columns: repeat(auto-fill, minmax(8.5rem, 1fr));
        gap: 0.1rem 0.7rem;
        margin: 0.7rem 0 0;
        padding: 0;
        list-style: none;
        font-size: var(--sl-text-xs);
        color: var(--sl-color-gray-3);
    }
    .columns li { display: flex; gap: 0.35rem; align-items: baseline; }
    .idx {
        min-width: 1.1rem;
        text-align: right;
        flex: none;
        color: var(--sl-color-gray-4);
        font-variant-numeric: tabular-nums;
    }
    .name { min-width: 0; }

    /* Also a grid: equal tracks put every exemplar on one baseline, which a
       flex row does not when the labels differ in width. */
    .classes {
        display: grid;
        grid-template-columns: repeat(auto-fill, minmax(4rem, 1fr));
        gap: 0.6rem 0.5rem;
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
        color: var(--sl-color-gray-2);
        text-align: center;
    }
    .classes svg {
        display: block;
        width: 100%;
        max-width: 3.2rem;
        height: auto;
        color: var(--sl-color-gray-2);
    }
    .foot {
        margin: 0.6rem 0 0;
        font-size: var(--sl-text-xs);
        color: var(--sl-color-gray-3);
    }
</style>
