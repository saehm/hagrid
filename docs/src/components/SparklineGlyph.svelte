<script>
    /**
     * The signal itself, drawn small.
     *
     * No encoding decision is taken here: the item *is* a curve, so the glyph is
     * that curve. Every other glyph in this section maps numbers onto a shape
     * someone chose; this one just plots them, which makes it the purest case
     * for the argument — the projection is computed from these curves, and the
     * plot shows you the curves.
     *
     * The vertical range is shared across the dataset rather than fitted per
     * item, so a tall beat really is a tall beat.
     *
     * Drawn around the origin, filling a square of side `2 × size`.
     */
    let {
        series,
        /** `[min, max]` across every series, so the glyphs are comparable. */
        range,
        size = 6,
        stroke = "#7dcfff",
        /** Baseline behind the trace. Off by default — hundreds of them tile into
         * a false grid; legends turn it on. */
        guides = false,
    } = $props();

    let path = $derived.by(() => {
        const [lo, hi] = range;
        const span = hi - lo || 1;
        const n = series.length;
        const step = (2 * size) / (n - 1);
        let d = "";
        for (let i = 0; i < n; ++i) {
            const x = -size + i * step;
            const t = (series[i] - lo) / span;
            const y = size - Math.min(1, Math.max(0, t)) * 2 * size;
            d += `${i === 0 ? "M" : "L"} ${x.toFixed(2)} ${y.toFixed(2)} `;
        }
        return d;
    });

    let mid = $derived.by(() => {
        const [lo, hi] = range;
        return size - ((0 - lo) / (hi - lo || 1)) * 2 * size;
    });
    let show_guides = $derived(guides && size >= 8);
</script>

<g>
    {#if show_guides}
        <line x1={-size} y1={mid} x2={size} y2={mid} class="base" stroke-width={size / 40} />
    {/if}
    <path
        d={path}
        fill="none"
        {stroke}
        stroke-width={Math.max(0.35, size / 7)}
        stroke-linejoin="round"
        stroke-linecap="round"
    />
</g>

<style>
    .base { stroke: currentColor; opacity: 0.3; }
</style>
