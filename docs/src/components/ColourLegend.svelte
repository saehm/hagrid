<script>
    /**
     * Legend for the plot's colour. Without one the colour is decoration — a
     * reader can see that two points differ but not what the difference means,
     * which is the whole payload of the encoding.
     *
     * Only origin, and it stays put when the plot switches to movement colours.
     * The movement ramp has its own key drawn on the histogram's axis below,
     * where a colour and a number share the same position — a better legend
     * than a detached swatch strip, and one that is on screen at all times. So
     * swapping this panel to a second copy of it bought nothing and cost the
     * reader the explanation of the encoding they started on.
     */
    import { origin_swatches } from "../lib/colour.js";

    let { dark = true } = $props();

    const swatches = $derived(origin_swatches(11, dark));
</script>

<div class="legend">
    <svg viewBox="0 0 1 1" class="disc" aria-hidden="true">
        {#each swatches as sw, i (i)}
            <rect x={sw.x - 0.05} y={1 - sw.y - 0.05} width="0.1" height="0.1" fill={sw.colour} />
        {/each}
    </svg>
    <p>
        <b>Where each point started.</b> Hue is the direction from the centre of the data, depth is
        the distance from it. Smooth colour means neighbours stayed neighbours; scrambled colour
        means the arrangement was lost.
    </p>
</div>

<style>
    .legend {
        display: flex;
        align-items: center;
        gap: 0.7rem;
        padding-top: 0.6rem;
    }
    p {
        margin: 0;
        font-size: 0.72rem;
        line-height: 1.55;
        color: var(--sl-color-gray-3);
    }
    b {
        color: var(--sl-color-gray-2);
        font-weight: 600;
    }
    .disc {
        flex: none;
        width: 2.6rem;
        height: 2.6rem;
    }
    @media (max-width: 560px) {
        .legend {
            align-items: start;
        }
    }
</style>
