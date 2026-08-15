/**
 * How far a rectangular glyph reaches, which depends on where its neighbours
 * can be.
 *
 * `GlyphShowcase` sizes glyphs by treating each as a disc of radius
 * `reach × size`: two glyphs cannot touch if their discs do not. For a round
 * glyph that is exact. For a rectangular one it depends on the layout, and
 * getting it wrong is expensive in a way that is easy to miss.
 *
 * On the raw projection a neighbour can be at any angle, including the 45°
 * that brings two corners together, so the reach is to the corner:
 * `√(w² + h²)` in units of the half-extent.
 *
 * On a lattice it is not. Neighbours sit on the axes, one pitch away, and
 * along an axis the glyph only extends by its half-width — so squares of half
 * side `pitch/2` tile the plane exactly. The diagonal neighbour is `pitch·√2`
 * away and its corner merely touches. Charging the corner reach there shrinks
 * every glyph by √2, about 41%, for a clearance that nothing needs: a gap of
 * `pitch(1 - 1/√2)` around every cell.
 *
 * The same holds for a hexagon matched to a hexagonal cell, whose neighbours
 * meet it across the flats at one inradius.
 */
export function box_reach(count, half_width = 1, half_height = half_width) {
    const corner = Math.hypot(half_width, half_height);
    const edge = Math.max(half_width, half_height);

    return (layout, cell) => {
        let reach;
        if (layout === "as projected") reach = corner;
        else if (cell?.kind === "hex") reach = hex_reach(half_width, half_height, cell.orientation);
        else reach = edge;
        return new Array(count).fill(reach);
    };
}

/**
 * A glyph that takes the *shape of its cell* — square on a square lattice,
 * hexagon on a hexagonal one — and so tiles whatever it is given.
 *
 * It reaches one half-pitch on every lattice. Only the raw projection, where
 * the glyph is a square and a neighbour may be diagonal, charges the corner.
 */
export function cell_reach(count) {
    return (layout) => new Array(count).fill(layout === "as projected" ? Math.SQRT2 : 1);
}

/**
 * A box on a *hexagonal* lattice, which is the case that looks obvious and is
 * not. Only two of the six neighbours sit on an axis; the other four are 30°
 * off, one pitch away, which places them `√3/2 · pitch` across and half a pitch
 * up. A square sized to clear the axial neighbour — half-side `pitch/2`, the
 * square-lattice answer — still overlaps those four, because both of their
 * offsets are under its full width.
 *
 * Two axis-aligned boxes miss each other if they are clear on *either* axis, so
 * the largest size is the tightest of the six directions:
 *
 *     s ≤ min over neighbours of  max( |dx| / 2w , |dy| / 2h )
 *
 * For a square that gives `2/√3 ≈ 1.1547` — the glyph is *smaller* on a
 * hexagonal grid than on a square one, the opposite of what a round glyph does.
 * A glyph shaped like the cell has no such problem, which is why the image
 * showcase draws a hexagon on Gosper and asks for a reach of 1.
 */
function hex_reach(half_width, half_height, orientation = "flat") {
    // Flat-topped cells meet their neighbours at 30° off the x-axis, pointy
    // ones along it. Three directions suffice; the other three mirror them.
    const base = orientation === "pointy" ? 0 : 30;
    let u = Infinity;
    for (let k = 0; k < 6; ++k) {
        const a = ((base + k * 60) * Math.PI) / 180;
        const dx = Math.abs(Math.cos(a));
        const dy = Math.abs(Math.sin(a));
        u = Math.min(u, Math.max(dx / (2 * half_width), dy / (2 * half_height)));
    }
    return 1 / (2 * u);
}
