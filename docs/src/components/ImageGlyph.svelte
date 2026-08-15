<script>
    /**
     * A photograph, framed in its class colour.
     *
     * The one "glyph" in this section that encodes nothing at all: the item is
     * a picture, so the picture is the glyph, and every argument about mapping
     * columns onto shapes falls away. What remains is the argument about room —
     * a thumbnail needs a cell far more obviously than a star does.
     *
     * The frame carries the label rather than an outline around the image,
     * because a stroke on a photograph competes with the photograph. Drawn as a
     * filled square behind it, the class reads as a mount around the picture.
     */
    let {
        /** Already resolved against the site base. */
        src,
        /**
         * The cell's inradius: half the side of a square cell, and half the
         * across-flats distance of a hexagonal one. Either way it is half the
         * pitch, which is what makes the two shapes comparable.
         */
        size,
        /** The class colour, used for the mount. */
        fill = "#888",
        /** `"square"` or `"hex"` — matched to the lattice the layout produced. */
        shape = "square",
        /** Which way up the hexagons sit, when hexagonal. */
        orientation = "flat",
        /** Enlarge the mount so the colour survives at legend sizes. */
        guides = false,
        alt = "",
    } = $props();

    // Thin at plot sizes, generous in a legend where there is room for it.
    const inset = $derived(guides ? size * 0.22 : Math.max(1, size * 0.12));

    /*
     * The largest hexagon that fits the cell is the cell: with neighbours one
     * pitch away across the flats, the inradius is `size` and the circumradius
     * is `size · 2/√3`. Using the circumradius as the half-width — as a square
     * glyph would — draws a hexagon 15% too small and leaves a visible gap
     * along every edge.
     */
    const R = $derived((size * 2) / Math.sqrt(3));

    /** Vertices at 0° for a flat top, 30° for a pointy one. */
    const hexagon = (r, deg0) =>
        Array.from({ length: 6 }, (_, i) => {
            const a = ((deg0 + i * 60) * Math.PI) / 180;
            return `${(Math.cos(a) * r).toFixed(2)},${(Math.sin(a) * r).toFixed(2)}`;
        }).join(" ");

    const deg0 = $derived(orientation === "pointy" ? 30 : 0);
    const uid = $props.id();
</script>

{#if shape === "hex"}
    <defs>
        <clipPath id="hexclip-{uid}">
            <polygon points={hexagon(R - inset, deg0)} />
        </clipPath>
    </defs>
    <polygon points={hexagon(R, deg0)} {fill} />
    <!-- Square source, hexagonal hole: the image covers the hexagon's bounding
         box and the clip takes the middle of it. -->
    <image
        href={src}
        x={-(R - inset)}
        y={-(R - inset)}
        width={2 * (R - inset)}
        height={2 * (R - inset)}
        preserveAspectRatio="xMidYMid slice"
        clip-path="url(#hexclip-{uid})"
        aria-label={alt}
    />
{:else}
    <rect x={-size} y={-size} width={2 * size} height={2 * size} rx={size * 0.08} {fill} />
    <image
        href={src}
        x={-size + inset}
        y={-size + inset}
        width={2 * (size - inset)}
        height={2 * (size - inset)}
        preserveAspectRatio="xMidYMid slice"
        aria-label={alt}
    />
{/if}
