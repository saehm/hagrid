<script>
    /**
     * One wheat kernel, drawn from its seven measurements.
     *
     * See `lib/seed-glyph.js` for which measurement drives what. The point of
     * drawing the thing the data describes, rather than a star, is that the
     * profile is legible without a key: a plump pale kernel and a long ruddy
     * one are obviously different grains, and that reads at sizes where
     * counting spokes would not.
     *
     * Drawn around the origin — place it with a transform.
     */
    import { seed_geometry, seed_path, groove_path, VARIETY_FILL } from "../lib/seed-glyph.js";

    let {
        values,
        ranges,
        size = 6,
        /** Variety name, or a colour directly. */
        variety = "Kama",
        fill = null,
    } = $props();

    let g = $derived(seed_geometry(values, ranges, size));
    let colour = $derived(fill ?? VARIETY_FILL[variety] ?? "#d8a657");
</script>

<g>
    <path
        d={seed_path(g)}
        fill={colour}
        stroke="rgba(0,0,0,.45)"
        stroke-width={g.stroke}
        stroke-linejoin="round"
    />
    <path
        d={groove_path(g)}
        fill="none"
        stroke="rgba(0,0,0,.4)"
        stroke-width={g.stroke * 1.1}
        stroke-linecap="round"
    />
</g>
