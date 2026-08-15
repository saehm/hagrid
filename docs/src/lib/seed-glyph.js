/**
 * Geometry for the wheat-kernel glyph, shared by the glyph and its legend.
 *
 * The seeds dataset measures seven things about a kernel, and a kernel is a
 * shape — so rather than abstracting them onto spokes, each measurement drives
 * the property of the drawn kernel that it actually describes:
 *
 *   area          overall size
 *   length/width  the aspect: how slender or plump
 *   compactness   how full the outline is between base and tip
 *   asymmetry     how far the tip leans off the axis
 *   groove length the crease down the middle
 *   perimeter     outline weight — a long perimeter for a given area means a
 *                 more ragged edge
 *
 * Column order is `[area, perimeter, compactness, length, width, asymmetry,
 * groove]`, as `build-dataset.mjs` writes it.
 */

const clamp01 = (t) => Math.min(1, Math.max(0, t));

/** Each measurement as 0..1 against its range across the dataset. */
function unit(values, ranges) {
    return values.map((v, j) => {
        const [lo, hi] = ranges[j];
        return hi > lo ? clamp01((v - lo) / (hi - lo)) : 0.5;
    });
}

export function seed_geometry(values, ranges, size) {
    const [area, perimeter, compactness, length, width, asymmetry, groove] = unit(values, ranges);

    // Area sets the scale; length and width set the shape at that scale.
    const scale = 0.78 + 0.44 * area;
    const half_length = size * scale * (0.62 + 0.38 * length);
    const half_width = size * scale * (0.30 + 0.34 * width);

    return {
        half_length,
        half_width,
        // How far the side control points push out: a compact kernel is full
        // between base and tip, a loose one waisted.
        bulge: 0.72 + 0.55 * compactness,
        // The tip leans off-axis; every kernel leans a little, the asymmetric
        // ones a lot.
        lean: (0.12 + 0.75 * asymmetry) * half_width,
        // The crease runs up from the base.
        groove: (0.35 + 0.5 * groove) * 2 * half_length,
        stroke: size * (0.035 + 0.05 * perimeter),
    };
}

/** How far the glyph reaches from its centre, per unit of `size`. */
export function seed_reach(values, ranges) {
    const g = seed_geometry(values, ranges, 1);
    return Math.max(g.half_length, g.half_width + g.lean) + g.stroke;
}

/**
 * A kernel pointing up, drawn about the origin: two cubics from base to tip,
 * the tip displaced sideways by `lean`.
 */
export function seed_path({ half_length: L, half_width: W, bulge, lean }) {
    const tip_x = lean;
    const tip_y = -L;
    const cy1 = L * 0.45;
    const cy2 = -L * 0.35;
    return (
        `M 0 ${L} ` +
        `C ${W * bulge} ${cy1} ${W * bulge} ${cy2} ${tip_x} ${tip_y} ` +
        `C ${-W * bulge + lean * 0.6} ${cy2} ${-W * bulge} ${cy1} 0 ${L} Z`
    );
}

/** The crease, following the same lean as the outline. */
export function groove_path({ half_length: L, groove, lean }) {
    const top = L - groove;
    return `M 0 ${L * 0.86} Q ${lean * 0.35} ${(L * 0.86 + top) / 2} ${lean * 0.7} ${top}`;
}

export const VARIETY_FILL = {
    Kama: "#d8a657",
    Rosa: "#c2704a",
    Canadian: "#e3c88a",
};
