<script>
    /**
     * The IRIS example from the VINCI'21 talk, wired into the generic
     * {@link GlyphShowcase}: 150 flowers projected to two dimensions, each drawn
     * from its own four measurements.
     *
     * This file supplies only what is specific to iris — the data, how to draw
     * one flower, and how to read it. Everything about layouts, sizing and
     * overlap counting lives in the harness.
     */
    import { IRIS } from "@saehrimnir/mistle";
    import GlyphShowcase from "./GlyphShowcase.svelte";
    import IrisGlyph from "./IrisGlyph.svelte";
    import IrisAnatomy from "./IrisAnatomy.svelte";
    import { reach, anatomy_frame } from "../lib/iris-glyph.js";
    import iris from "../data/iris-2d.json";

    let { method = "as projected", ...rest } = $props();

    const REACH = IRIS.values.map(reach);

    // One representative flower per species. They share a viewBox sized to the
    // largest of the three, so the species are compared at one scale rather
    // than each filling its own box.
    const SPECIES = [...new Set(IRIS.labels)].map((name) => ({
        name,
        index: IRIS.labels.indexOf(name),
    }));
    const SPECIES_SIZE = 4.5;
    const SPECIES_HALF = Math.ceil(
        Math.max(...SPECIES.map((s) => reach(IRIS.values[s.index]))) * SPECIES_SIZE + 3,
    );

    // Petal and sepal share a frame, so they draw to one scale, sit on a common
    // origin, and their dimension lines align across the pair.
    const ANATOMY_ROW = [6.3, 3.3, 6.0, 2.5];
    const ANATOMY_SIZE = 5;
    const ANATOMY_FRAME = anatomy_frame(ANATOMY_ROW, ANATOMY_SIZE);
</script>

<GlyphShowcase points={iris.points} reach={REACH} {method} {...rest}>
    {#snippet glyph(i, size)}
        <IrisGlyph index={i} {size} />
    {/snippet}

    {#snippet details(i)}
        <p class="tip-title">{IRIS.labels[i].replace("Iris-", "")}</p>
        <svg viewBox="-46 -46 92 92" class="tip-glyph" aria-hidden="true">
            <IrisGlyph index={i} size={7} />
        </svg>
        <dl>
            {#each IRIS.columns as c, j}
                <dt>{c.replace("length", " length").replace("width", " width")}</dt>
                <dd>{IRIS.values[i][j].toFixed(1)} cm</dd>
            {/each}
        </dl>
    {/snippet}

    {#snippet legend()}
        <div class="legend">
            <p class="legend-title">Reading a flower</p>
            <div class="anatomy">
                <figure>
                    <IrisAnatomy kind="petal" row={ANATOMY_ROW} size={ANATOMY_SIZE} frame={ANATOMY_FRAME} />
                    <figcaption><b>Petal</b> — the inner, darker arm. Three of them.</figcaption>
                </figure>
                <figure>
                    <IrisAnatomy kind="sepal" row={ANATOMY_ROW} size={ANATOMY_SIZE} frame={ANATOMY_FRAME} />
                    <figcaption><b>Sepal</b> — the outer, lighter arm. Three of them.</figcaption>
                </figure>
            </div>
            <p class="foot">
                A flower carries four numbers: the length and width of its petals and of its
                sepals. Each is scaled against that column's maximum across the dataset, so a
                bigger flower on screen really is a bigger flower.
            </p>
        </div>

        <div class="legend">
            <p class="legend-title">Species</p>
            <ul class="species">
                {#each SPECIES as s}
                    <li>
                        <svg
                            viewBox="{-SPECIES_HALF} {-SPECIES_HALF} {2 * SPECIES_HALF} {2 * SPECIES_HALF}"
                            aria-hidden="true"
                        >
                            <IrisGlyph index={s.index} size={SPECIES_SIZE} />
                        </svg>
                        <span>{s.name.replace("Iris-", "")}</span>
                    </li>
                {/each}
            </ul>
            <p class="foot">
                Petals are the darker shade, sepals the lighter. One exemplar each, so the shapes
                differ as the measurements do.
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
    .anatomy { display: flex; gap: 1.1rem; align-items: flex-start; }
    .anatomy figure { margin: 0; flex: 1 1 0; min-width: 0; }
    .anatomy figcaption {
        padding: 0.35rem 0 0;
        font-size: var(--sl-text-xs);
        color: var(--sl-color-gray-2);
        text-align: center;
        line-height: 1.35;
    }
    ul { margin: 0; padding: 0; list-style: none; font-size: var(--sl-text-xs); }
    .species { display: flex; gap: 0.6rem; }
    .species li {
        flex: 1 1 0;
        min-width: 0;
        display: flex;
        flex-direction: column;
        align-items: center;
        gap: 0.25rem;
        text-align: center;
        color: var(--sl-color-gray-2);
    }
    .species svg { display: block; width: 100%; max-width: 4.4rem; height: auto; }
    .foot {
        margin: 0.6rem 0 0;
        font-size: var(--sl-text-xs);
        color: var(--sl-color-gray-3);
    }
    .tip-glyph { display: block; width: 5rem; height: auto; margin: 0.1rem auto 0.4rem; }
</style>
