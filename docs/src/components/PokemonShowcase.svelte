<script>
    /**
     * Pokémon base stats, drawn as the hexagonal radar the games use.
     *
     * The interesting part is the mismatch: the projection is computed from
     * stats, and the colour shows typing — which turns out to be nearly
     * unrelated. The glyphs are what let you see that rather than guess it.
     */
    import GlyphShowcase from "./GlyphShowcase.svelte";
    import PokemonGlyph from "./PokemonGlyph.svelte";
    import { STATS, STAT_SIDE, TYPE_COLOUR, types_of } from "../lib/pokemon.js";
    import pokemon from "../data/pokemon-2d.json";

    let { method = "as projected", ...rest } = $props();

    // The polygon is inscribed in a circle of radius `size`.
    const REACH = pokemon.points.map(() => 1);

    const primary = pokemon.labels.map((l) => types_of(l)[0]);
    const TYPES = [...new Set(primary)].sort(
        (a, b) =>
            primary.filter((t) => t === b).length - primary.filter((t) => t === a).length,
    );

    /** A recognisable example for each stat archetype. */
    const find = (name) => pokemon.names.findIndex((n) => n === name);
    const ARCHETYPES = [
        ["Shuckle", "all defence, no offence"],
        ["Chansey", "enormous HP, paper defence"],
        ["Mewtwo", "strong everywhere"],
        ["Magikarp", "nothing at all"],
    ]
        .map(([name, note]) => ({ name, note, index: find(name) }))
        .filter((a) => a.index >= 0);

    const total = (i) => pokemon.values[i].reduce((s, v) => s + v, 0) / 6;
</script>

<GlyphShowcase points={pokemon.points} reach={REACH} {method} {...rest}>
    {#snippet glyph(i, size)}
        <PokemonGlyph stats={pokemon.values[i]} label={pokemon.labels[i]} {size} />
    {/snippet}

    {#snippet details(i)}
        <p class="tip-title">
            #{pokemon.numbers[i]} {pokemon.names[i]}
        </p>
        <p class="typing">
            {#each types_of(pokemon.labels[i]) as t}
                <span class="chip" style="background: {TYPE_COLOUR[t]}">{t}</span>
            {/each}
        </p>
        <svg viewBox="-24 -24 48 48" class="tip-glyph" aria-hidden="true">
            <PokemonGlyph stats={pokemon.values[i]} label={pokemon.labels[i]} size={19} guides />
        </svg>
        <dl>
            {#each STATS as s, j}
                <dt>{s}</dt>
                <dd>{(pokemon.values[i][j] * 100).toFixed(0)}</dd>
            {/each}
        </dl>
    {/snippet}

    {#snippet legend()}
        <div class="legend">
            <p class="legend-title">Reading a Pokémon</p>
            <div class="reading">
                <svg viewBox="-26 -26 52 52" aria-hidden="true">
                    <PokemonGlyph
                        stats={pokemon.values[ARCHETYPES[0]?.index ?? 0]}
                        label={pokemon.labels[ARCHETYPES[0]?.index ?? 0]}
                        size={20}
                        guides
                    />
                </svg>
                <ul class="sides">
                    {#each STATS as s}
                        <li><b>{s}</b><span>{STAT_SIDE[s]}</span></li>
                    {/each}
                </ul>
            </div>
            <p class="foot">
                The games' own layout: physical stats down the right, special down the left, speed
                at the bottom. So a physically bulky Pokémon leans right and a special attacker
                leans left. A dual-typed Pokémon is <b>split down the middle</b> — primary type
                on the left, secondary on the right; a single type fills the whole shape.
            </p>
        </div>

        <!-- Both of these are shorter than the reading guide, so they share a
             column rather than leaving a void under a lone short block. -->
        <div class="stack">
            <div class="legend">
                <p class="legend-title">Shapes worth knowing</p>
                <ul class="archetypes">
                    {#each ARCHETYPES as a}
                        <li>
                            <svg viewBox="-24 -24 48 48" aria-hidden="true">
                                <PokemonGlyph
                                    stats={pokemon.values[a.index]}
                                    label={pokemon.labels[a.index]}
                                    size={19}
                                    guides
                                />
                            </svg>
                            <span class="who">{a.name}</span>
                            <span class="note">{a.note}</span>
                        </li>
                    {/each}
                </ul>
            </div>

            <div class="legend">
                <p class="legend-title">Primary type</p>
                <ul class="types">
                    {#each TYPES as t}
                        <li>
                            <span class="dot" style="background: {TYPE_COLOUR[t]}"></span>{t}
                            <small>{primary.filter((p) => p === t).length}</small>
                        </li>
                    {/each}
                </ul>
                <p class="foot">
                    Scattered, not clustered — which is the point. Typing decides what a Pokémon is
                    strong against, not how its stats are distributed.
                </p>
            </div>
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
    .stack { display: grid; gap: 1rem; align-content: start; }
    .reading { display: flex; gap: 0.9rem; align-items: flex-start; }
    .reading svg { width: 5rem; flex: none; color: var(--sl-color-gray-2); }
    .sides { margin: 0; padding: 0; list-style: none; font-size: var(--sl-text-xs); }
    .sides li { display: grid; grid-template-columns: 4.2rem 1fr; gap: 0.4rem; }
    .sides b { color: var(--sl-color-white); font-weight: 600; }
    .sides span { color: var(--sl-color-gray-3); }

    /* Two fixed columns, not auto-fill: four items wrapping 3 + 1 leaves a
       ragged row, and a `1fr` track that an unbreakable note can widen makes
       the glyphs come out different sizes. */
    .archetypes {
        display: grid;
        grid-template-columns: repeat(2, minmax(0, 1fr));
        gap: 0.6rem 0.9rem;
        margin: 0;
        padding: 0;
        list-style: none;
    }
    /* Glyph beside its label, the same way the reading guide above reads. */
    .archetypes li {
        display: grid;
        grid-template-columns: 3rem minmax(0, 1fr);
        align-items: center;
        column-gap: 0.5rem;
        font-size: var(--sl-text-xs);
    }
    .archetypes svg {
        grid-row: span 2;
        display: block;
        width: 3rem;
        height: auto;
        color: var(--sl-color-gray-2);
    }
    .archetypes .who { align-self: end; color: var(--sl-color-white); font-weight: 600; }
    .archetypes .note { align-self: start; color: var(--sl-color-gray-3); line-height: 1.3; }

    .types {
        display: grid;
        grid-template-columns: repeat(auto-fill, minmax(6.2rem, 1fr));
        gap: 0.1rem 0.6rem;
        margin: 0;
        padding: 0;
        list-style: none;
        font-size: var(--sl-text-xs);
        color: var(--sl-color-gray-2);
    }
    .types li { display: flex; align-items: center; gap: 0.35rem; }
    .types small { margin-left: auto; color: var(--sl-color-gray-4); font-variant-numeric: tabular-nums; }
    .dot { width: 0.6rem; height: 0.6rem; border-radius: 50%; flex: none; }

    .typing { margin: 0.15rem 0 0.35rem; display: flex; gap: 0.25rem; }
    .chip {
        padding: 0.05rem 0.35rem;
        border-radius: 999px;
        font-size: 0.68em;
        font-weight: 600;
        color: #12161c;
    }
    .tip-glyph { display: block; width: 6rem; height: auto; margin: 0 auto 0.4rem; }
    .foot { margin: 0.6rem 0 0; font-size: var(--sl-text-xs); color: var(--sl-color-gray-3); }
</style>
