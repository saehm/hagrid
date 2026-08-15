<script>
    /**
     * A stacked column: parts of a whole, in a fixed order, filling a fixed
     * height.
     *
     * This is the encoding compositional data asks for. A star or a bar chart
     * treats the columns as independent magnitudes, but these sum to 100% — so
     * the total carries no information and only the split does. A stack shows
     * exactly that: every glyph is the same height, and the reading is where
     * the boundaries fall.
     *
     * The parts keep fixed colours rather than fixed positions in a legend, so
     * a band can be followed across neighbouring glyphs.
     *
     * Drawn around the origin, filling a square of side `2 × size`.
     */
    let {
        /** The parts, in display order. Rescaled to sum to 1. */
        parts,
        /** One colour per part. */
        colours,
        size = 6,
        /** A frame in the class colour; omitted when null. */
        outline = null,
    } = $props();

    let total = $derived(parts.reduce((a, b) => a + b, 0) || 1);

    /** Cumulative offsets, so each band knows where it starts. */
    let bands = $derived.by(() => {
        const out = [];
        let y = size; // bottom
        for (let i = 0; i < parts.length; ++i) {
            const h = (parts[i] / total) * 2 * size;
            out.push({ y: y - h, h, fill: colours[i % colours.length] });
            y -= h;
        }
        return out;
    });
</script>

<g>
    {#each bands as b}
        {#if b.h > 0}
            <rect x={-size * 0.62} y={b.y} width={size * 1.24} height={b.h} fill={b.fill} />
        {/if}
    {/each}
    {#if outline}
        <rect
            x={-size * 0.62}
            y={-size}
            width={size * 1.24}
            height={2 * size}
            fill="none"
            stroke={outline}
            stroke-width={Math.max(0.4, size / 6)}
        />
    {/if}
</g>
