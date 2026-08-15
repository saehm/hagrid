<script>
    /**
     * An image collection on the {@link GlyphShowcase} harness.
     *
     * Same controls, same overlap counting, same displacement readout as the
     * glyph showcases — the only difference is that the thing in each cell is a
     * photograph. That is the point of putting it here: the case for gridding
     * does not depend on the glyph being clever.
     */
    import GlyphShowcase from "./GlyphShowcase.svelte";
    import ImageGlyph from "./ImageGlyph.svelte";
    import { cell_reach } from "../lib/glyph-reach.js";

    let {
        /** The JSON from `build-images.mjs`. */
        data,
        method = "as projected",
        ...rest
    } = $props();

    /*
     * Files live under `public/`, so they need the site base in front. Astro
     * hands `BASE_URL` over without a trailing slash here, so the separator is
     * added rather than assumed — concatenating gave `/hagridshowcase/...`.
     */
    const BASE = import.meta.env.BASE_URL.replace(/\/+$/, "");
    const src = (name) => `${BASE}/${data.dir}/${name}.${data.extension}`;

    const PALETTE = [
        "#7aa2f7", "#e0af68", "#9ece6a", "#f7768e", "#bb9af7",
        "#7dcfff", "#ff9e64", "#73daca", "#c0caf5", "#e06c9f", "#a1b5d8",
    ];
    const COLOUR = Object.fromEntries(
        data.classes.map((c, i) => [c, PALETTE[i % PALETTE.length]]),
    );

    // Square on a square lattice, hexagon on Gosper's — either way the glyph is
    // the shape of its cell, so it tiles and reaches one half-pitch.
    const REACH = cell_reach(data.points.length);

    const counts = $derived(
        Object.fromEntries(
            data.classes.map((c) => [c, data.labels.filter((l) => l === c).length]),
        ),
    );

    /** The first item of each class, as its legend swatch. */
    const EXEMPLAR = data.classes.map((c) => ({
        label: c,
        index: data.labels.indexOf(c),
    }));

    const pretty = (s) => s.replace(/_/g, " ");
</script>

<GlyphShowcase points={data.points} reach={REACH} {method} {...rest}>
    {#snippet glyph(i, size, cell)}
        <ImageGlyph
            src={src(data.names[i])}
            {size}
            fill={COLOUR[data.labels[i]]}
            shape={cell?.kind ?? "square"}
            orientation={cell?.orientation ?? "flat"}
            alt={pretty(data.labels[i])}
        />
    {/snippet}

    {#snippet details(i)}
        <p class="tip-title" style="color: {COLOUR[data.labels[i]]}">{pretty(data.labels[i])}</p>
        <img class="tip-image" src={src(data.names[i])} alt={pretty(data.labels[i])} />
        <dl>
            <dt>file</dt><dd>{data.names[i]}</dd>
        </dl>
    {/snippet}

    {#snippet legend()}
        <div class="legend">
            <p class="legend-title">Classes</p>
            <ul class="classes">
                {#each EXEMPLAR as e}
                    <li>
                        <svg viewBox="-22 -22 44 44" aria-hidden="true">
                            <ImageGlyph
                                src={src(data.names[e.index])}
                                size={19}
                                fill={COLOUR[e.label]}
                                guides
                            />
                        </svg>
                        <span class="name">{pretty(e.label)}</span>
                        <span class="count">{counts[e.label]}</span>
                    </li>
                {/each}
            </ul>
            <p class="foot">
                The mount around each thumbnail is its class. Nothing else is encoded — the
                picture is the data, and the colour is only there so you can see whether the
                projection grouped the classes or not.
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
    .classes {
        display: grid;
        grid-template-columns: repeat(auto-fill, minmax(4.5rem, 1fr));
        gap: 0.6rem 0.5rem;
        margin: 0;
        padding: 0;
        list-style: none;
    }
    .classes li {
        display: grid;
        justify-items: center;
        gap: 0.15rem;
        font-size: var(--sl-text-xs);
        text-align: center;
    }
    .classes svg { display: block; width: 100%; max-width: 3rem; height: auto; }
    .classes .name { color: var(--sl-color-white); font-weight: 600; }
    .classes .count { color: var(--sl-color-gray-4); font-variant-numeric: tabular-nums; }
    .tip-image {
        display: block;
        width: 6rem;
        height: 6rem;
        object-fit: cover;
        border-radius: 3px;
        margin: 0.1rem 0 0.4rem;
    }
    .foot { margin: 0.6rem 0 0; font-size: var(--sl-text-xs); color: var(--sl-color-gray-3); }
</style>
