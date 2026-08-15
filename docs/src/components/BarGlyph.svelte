<script>
    /**
     * A bar glyph — a small multiple, one bar per measurement, sharing a
     * baseline and a scale.
     *
     * Where a star glyph reads as a silhouette, bars read as a sequence: it is
     * much easier to say *which* column is the tall one, and to compare the
     * same column across neighbours, because every bar sits at a fixed x. The
     * cost is that the profile is no longer one closed shape, so it takes
     * slightly longer to scan.
     *
     * Drawn around the origin, filling a square of side `2 × size`.
     */
    let {
        /** The item's measurements, one per column. */
        values,
        /** `[[min, max], …]` per column, from the whole dataset. */
        ranges,
        size = 6,
        fill = "#7aa2f7",
        /** Baseline behind the bars. Off by default — see StarGlyph. */
        guides = false,
    } = $props();

    const FLOOR = 0.06;

    let unit = $derived(
        values.map((v, j) => {
            const [lo, hi] = ranges[j];
            const t = hi > lo ? (v - lo) / (hi - lo) : 0.5;
            return FLOOR + (1 - FLOOR) * Math.min(1, Math.max(0, t));
        }),
    );

    let step = $derived((2 * size) / unit.length);
    let gap = $derived(Math.min(step * 0.22, size / 8));
    let show_guides = $derived(guides && size >= 7);
</script>

<g>
    {#if show_guides}
        <line
            x1={-size} y1={size} x2={size} y2={size}
            class="base" stroke-width={size / 30}
        />
    {/if}
    {#each unit as t, j}
        <rect
            x={-size + j * step + gap / 2}
            y={size - t * 2 * size}
            width={Math.max(step - gap, 0.1)}
            height={t * 2 * size}
            {fill}
            fill-opacity="0.9"
        />
    {/each}
</g>

<style>
    .base { stroke: currentColor; opacity: 0.35; }
</style>
