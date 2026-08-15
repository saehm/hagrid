<script>
    /**
     * A star (radar) glyph: one spoke per measurement, at equal angles, each as
     * long as that measurement is large. Joining the spoke ends gives an
     * outline whose *shape* is the item's profile — which is what makes a field
     * of them readable at a glance, long before any single value is.
     *
     * Every column is normalised against its own range across the dataset, so
     * two glyphs can be compared. Scaling each row to itself would make every
     * item look average.
     *
     * Drawn around the origin, inside a circle of radius `size`.
     */
    let {
        /** The item's measurements, one per column. */
        values,
        /** `[[min, max], …]` per column, from the whole dataset. */
        ranges,
        size = 6,
        fill = "#7aa2f7",
        stroke = "none",
        /** Faint rings and spokes behind the outline. Off by default: in the plot
         * there are hundreds of them and they read as a broken grid. Legends turn
         * them on, where one large glyph benefits. */
        guides = false,
    } = $props();

    // A floor so a minimum value is still a visible nub rather than nothing.
    const FLOOR = 0.12;

    let unit = $derived(
        values.map((v, j) => {
            const [lo, hi] = ranges[j];
            const t = hi > lo ? (v - lo) / (hi - lo) : 0.5;
            return FLOOR + (1 - FLOOR) * Math.min(1, Math.max(0, t));
        }),
    );

    let spokes = $derived(
        unit.map((t, j) => {
            const a = (j / unit.length) * 2 * Math.PI - Math.PI / 2;
            return [Math.cos(a) * t * size, Math.sin(a) * t * size];
        }),
    );

    let outline = $derived(spokes.map(([x, y]) => `${x},${y}`).join(" "));
    let show_guides = $derived(guides && size >= 7);
</script>

<g>
    {#if show_guides}
        <circle r={size} class="ring" stroke-width={size / 40} />
        <circle r={size / 2} class="ring" stroke-width={size / 40} />
        {#each spokes as _, j}
            {@const a = (j / spokes.length) * 2 * Math.PI - Math.PI / 2}
            <line
                x1="0" y1="0" x2={Math.cos(a) * size} y2={Math.sin(a) * size}
                class="spoke" stroke-width={size / 40}
            />
        {/each}
    {/if}
    <polygon points={outline} {fill} {stroke} stroke-width={size / 12} fill-opacity="0.85" />
</g>

<style>
    .ring, .spoke { fill: none; stroke: currentColor; opacity: 0.28; }
</style>
