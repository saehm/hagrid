<script>
    /**
     * Base stats as the hexagonal radar the games use, with the two types shown
     * as the two halves of the fill.
     *
     * Two encodings, two channels: the *shape* is the stat spread, the *colours*
     * are the typing. An outline was the first attempt at the second type and it
     * failed badly — to be visible at glyph size the stroke has to be thick, and
     * a few hundred thick rings turn the plot into a mesh where no shape can be
     * read. Splitting the fill down the middle costs no extra ink and leaves the
     * silhouette intact.
     *
     * Single-typed Pokémon are simply one colour, so a split reads as "two
     * types" without any key.
     *
     * Drawn around the origin, inside a circle of radius `size`.
     */
    import { corners, half, to_points, TYPE_COLOUR, types_of } from "../lib/pokemon.js";

    let {
        /** The six normalised stats, in dataset order. */
        stats,
        /** e.g. `"Grass Poison"`. */
        label,
        size = 6,
        /** The hexagon behind the shape. Off by default: 250 of them overlapping
         * read as a broken grid rather than as guides. */
        guides = false,
    } = $props();

    let types = $derived(types_of(label));
    let primary = $derived(TYPE_COLOUR[types[0]] ?? "#a8a77a");
    let secondary = $derived(types[1] ? TYPE_COLOUR[types[1]] : null);

    let shape = $derived(corners(stats, size));
    let points = $derived(to_points(shape));
    let right = $derived(secondary ? to_points(half(shape, true)) : "");
    let frame = $derived(to_points(corners([1, 1, 1, 1, 1, 1], size)));
    let show_guides = $derived(guides && size >= 9);
</script>

<g>
    {#if show_guides}
        <polygon points={frame} class="frame" stroke-width={size / 40} />
    {/if}
    <polygon {points} fill={primary} />
    {#if secondary && right}
        <polygon points={right} fill={secondary} />
    {/if}
    <!-- A hairline only, to separate touching glyphs. Any thicker and it
         competes with the fill it is supposed to be delimiting. -->
    <polygon {points} fill="none" class="edge" stroke-width={Math.max(0.25, size / 22)} />
</g>

<style>
    .frame { fill: none; stroke: currentColor; opacity: 0.25; }
    .edge { stroke: rgb(0 0 0 / 45%); }
</style>
