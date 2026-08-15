<script>
    /** Wheat kernels drawn as kernels. See `lib/seed-glyph.js`. */
    import GlyphShowcase from "./GlyphShowcase.svelte";
    import SeedGlyph from "./SeedGlyph.svelte";
    import { seed_reach, VARIETY_FILL } from "../lib/seed-glyph.js";
    import seeds from "../data/seeds-2d.json";

    let { method = "as projected", ...rest } = $props();

    const REACH = seeds.values.map((v) => seed_reach(v, seeds.ranges));
    const VARIETIES = [...new Set(seeds.labels)];

    /** The kernel nearest each variety's average, so it looks typical. */
    const EXEMPLAR = VARIETIES.map((name) => {
        const rows = seeds.values.filter((_, i) => seeds.labels[i] === name);
        const mean = seeds.columns.map((_, j) => rows.reduce((s, r) => s + r[j], 0) / rows.length);
        let best = -1, bd = Infinity;
        seeds.values.forEach((r, i) => {
            if (seeds.labels[i] !== name) return;
            const d = r.reduce((s, v, j) => s + (v - mean[j]) ** 2, 0);
            if (d < bd) { bd = d; best = i; }
        });
        return { name, index: best };
    });

    const CHANNELS = [
        ["area", "overall size"],
        ["length & width", "how slender or plump"],
        ["compactness", "how full the outline is"],
        ["asymmetry", "how far the tip leans"],
        ["groove length", "the crease down the middle"],
        ["perimeter", "outline weight — a ragged edge"],
    ];
</script>

<GlyphShowcase points={seeds.points} reach={REACH} {method} {...rest}>
    {#snippet glyph(i, size)}
        <SeedGlyph
            values={seeds.values[i]}
            ranges={seeds.ranges}
            {size}
            variety={seeds.labels[i]}
        />
    {/snippet}

    {#snippet details(i)}
        <p class="tip-title" style="color: {VARIETY_FILL[seeds.labels[i]]}">{seeds.labels[i]}</p>
        <dl>
            {#each seeds.columns as c, j}
                <dt>{c}</dt>
                <dd>{seeds.values[i][j].toFixed(2)}</dd>
            {/each}
        </dl>
    {/snippet}

    {#snippet legend()}
        <div class="legend">
            <p class="legend-title">Reading a kernel</p>
            <div class="reading">
                <svg viewBox="-26 -26 52 52" aria-hidden="true">
                    <SeedGlyph
                        values={seeds.values[EXEMPLAR[0].index]}
                        ranges={seeds.ranges}
                        size={17}
                        variety={EXEMPLAR[0].name}
                    />
                </svg>
                <dl>
                    {#each CHANNELS as [what, how]}
                        <div><dt>{what}</dt><dd>{how}</dd></div>
                    {/each}
                </dl>
            </div>
            <p class="foot">
                Each measurement drives the property of the kernel it actually describes, so the
                shape is legible without a key — no spoke-counting.
            </p>
        </div>

        <div class="legend">
            <p class="legend-title">Varieties</p>
            <ul class="varieties">
                {#each EXEMPLAR as e}
                    <li>
                        <svg viewBox="-26 -26 52 52" aria-hidden="true">
                            <SeedGlyph
                                values={seeds.values[e.index]}
                                ranges={seeds.ranges}
                                size={17}
                                variety={e.name}
                            />
                        </svg>
                        <span style="color: {VARIETY_FILL[e.name]}">{e.name}</span>
                    </li>
                {/each}
            </ul>
            <p class="foot">
                Rosa kernels are simply bigger; Canadian ones are small and pale. The varieties
                differ more in size than in shape, which the glyphs show and the projection does
                not.
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
    .reading svg { width: 4.4rem; flex: none; }
    dl { margin: 0; font-size: var(--sl-text-xs); color: var(--sl-color-gray-2); }
    dl div { display: flex; gap: 0.4rem; margin: 0 0 0.15rem; }
    dt { color: var(--sl-color-white); white-space: nowrap; }
    dd { margin: 0; color: var(--sl-color-gray-3); }
    /* Equal grid tracks, so every exemplar sits on one baseline whatever the
       label width. */
    .varieties {
        display: grid;
        grid-template-columns: repeat(auto-fill, minmax(4.4rem, 1fr));
        gap: 0.5rem;
        margin: 0;
        padding: 0;
        list-style: none;
    }
    .varieties li {
        display: grid;
        grid-template-rows: auto auto;
        justify-items: center;
        gap: 0.2rem;
        font-size: var(--sl-text-xs);
        text-align: center;
    }
    .varieties svg { display: block; width: 100%; max-width: 3.6rem; height: auto; }
    .foot { margin: 0.6rem 0 0; font-size: var(--sl-text-xs); color: var(--sl-color-gray-3); }
</style>
