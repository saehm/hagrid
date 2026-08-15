<script>
    /**
     * One recording as its spectrogram: time left to right, frequency bottom to
     * top, brightness for energy.
     *
     * The image is a 16x16 greyscale PNG built at generation time, so a glyph
     * is a single `<image>` rather than 256 rectangles — 300 of those would be
     * 76,800 nodes. Colour carries the class instead, as a frame, which also
     * means one image serves whichever labelling is selected.
     *
     * Drawn around the origin, filling a square of side `2 × size`.
     */
    let {
        /** A `data:image/png;base64,…` spectrogram. */
        src,
        size = 6,
        /** Class colour, drawn as the frame. */
        colour = "#7dcfff",
    } = $props();
</script>

<g>
    <image
        href={src}
        x={-size}
        y={-size}
        width={2 * size}
        height={2 * size}
        preserveAspectRatio="none"
    />
    <rect
        x={-size}
        y={-size}
        width={2 * size}
        height={2 * size}
        fill="none"
        stroke={colour}
        stroke-width={Math.max(0.5, size / 5)}
    />
</g>

<style>
    /* The source is 16x16; smoothing it turns a spectrogram into a smudge. */
    image { image-rendering: pixelated; }
</style>
