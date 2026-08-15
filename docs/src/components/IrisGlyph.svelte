<script>
    /**
     * One iris flower, drawn from its four measurements.
     *
     * Three petals carry petal length and width, three sepals carry sepal
     * length and width, and the two sets interleave so each measurement gets
     * its own visual channel. The disc in the middle is decoration, sized by
     * the `size` prop.
     *
     * Drawn around the origin — place it with a transform.
     */
    import { IRIS } from "@saehrimnir/mistle";
    import {
        PETAL_FILL,
        SEPAL_FILL,
        PETAL_ANGLES,
        SEPAL_ANGLES,
        EXAMPLE_ROW,
        arms,
        petal_path,
        sepal_path,
    } from "../lib/iris-glyph.js";

    let {
        /** Row of IRIS to draw, or -1 for a neutral example flower. */
        index = -1,
        /** Radius of the centre disc; every other length scales with it. */
        size = 2,
        stroke = "none",
        /** Overrides the species colour, e.g. to grey out a glyph. */
        petal_fill = null,
        sepal_fill = null,
    } = $props();

    let row = $derived(index >= 0 ? IRIS.values[index] : EXAMPLE_ROW);
    let species = $derived(index >= 0 ? IRIS.labels[index] : "default");
    let [sepal_length, sepal_width, petal_length, petal_width] = $derived(arms(row, size));
</script>

<g>
    {#each PETAL_ANGLES as angle}
        <g transform="rotate({angle})">
            <path
                d={petal_path(size, petal_length, petal_width)}
                {stroke}
                stroke-width={size / 10}
                fill={petal_fill ?? PETAL_FILL[species] ?? PETAL_FILL.default}
            />
        </g>
    {/each}
    {#each SEPAL_ANGLES as angle}
        <g transform="rotate({angle})">
            <path
                d={sepal_path(size, sepal_length, sepal_width)}
                {stroke}
                stroke-width={size / 10}
                fill={sepal_fill ?? SEPAL_FILL[species] ?? SEPAL_FILL.default}
            />
        </g>
    {/each}
    <circle r={size} fill="yellow" />
</g>
