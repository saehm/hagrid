/**
 * The Pokémon glyph, which is the one case in this section where the right
 * shape is not a design decision — the games already draw base stats as a
 * hexagonal radar, so anyone who has played one can read it without a legend.
 *
 * Six stats, and the arrangement is theirs too: HP on top, the physical stats
 * down the right, the special stats down the left, speed at the bottom. That
 * puts Attack opposite Sp. Def and Defense opposite Sp. Atk, so a physically
 * bulky Pokémon leans right and a special attacker leans left.
 */

export const STATS = ["HP", "Attack", "Defense", "Sp. Atk", "Sp. Def", "Speed"];

/** Degrees clockwise from the top, indexed as the data is. */
export const STAT_ANGLE = [0, 60, 120, -60, -120, 180];

/** Where each stat sits, for the legend. */
export const STAT_SIDE = {
    HP: "top",
    Attack: "upper right",
    Defense: "lower right",
    "Sp. Atk": "upper left",
    "Sp. Def": "lower left",
    Speed: "bottom",
};

/** The series' own type palette — familiar, and already distinct by design. */
export const TYPE_COLOUR = {
    Normal: "#a8a77a",
    Fire: "#ee8130",
    Water: "#6390f0",
    Electric: "#f7d02c",
    Grass: "#7ac74c",
    Ice: "#96d9d6",
    Fighting: "#c22e28",
    Poison: "#a33ea1",
    Ground: "#e2bf65",
    Flying: "#a98ff3",
    Psychic: "#f95587",
    Bug: "#a6b91a",
    Rock: "#b6a136",
    Ghost: "#735797",
    Dragon: "#6f35fc",
    Dark: "#705746",
    Steel: "#b7b7ce",
    Fairy: "#d685ad",
};

/** `"Grass Poison"` to `["Grass", "Poison"]`. */
export const types_of = (label) => label.trim().split(/\s+/);

/**
 * The six corners of the stat hexagon, for a glyph of radius `size`.
 *
 * A floor keeps a minimum stat as a visible nub rather than collapsing the
 * polygon onto its own centre, where the shape would stop being readable.
 */
/**
 * The order the corners must be *drawn* in, which is not the order the stats
 * are stored in.
 *
 * The stats sit at 0°, 60°, 120°, 300°, 240°, 180° so that physical runs down
 * the right and special down the left. Emitting the polygon in that order walks
 * back and forth across the circle and produces a self-intersecting star rather
 * than a hexagon. Sorting by angle once, here, keeps every consumer honest.
 */
const DRAW_ORDER = STAT_ANGLE.map((angle, i) => ({ angle, i }))
    .sort((a, b) => a.angle - b.angle)
    .map(({ i }) => i);

export function corners(values, size, floor = 0.14) {
    return DRAW_ORDER.map((i) => {
        const v = values[i];
        const t = floor + (1 - floor) * Math.min(1, Math.max(0, v));
        const a = ((STAT_ANGLE[i] - 90) * Math.PI) / 180;
        return [Math.cos(a) * t * size, Math.sin(a) * t * size];
    });
}

/**
 * The part of a polygon on one side of the vertical axis.
 *
 * Used to fill the two halves of the hexagon in the two type colours, which is
 * how a dual typing gets shown without an outline. An outline was the obvious
 * answer and the wrong one: at glyph sizes it has to be thick enough to see,
 * and hundreds of thick rings turn the plot into a mesh. A split fill costs no
 * extra ink at all.
 *
 * Sutherland–Hodgman against a single edge — six vertices, so this is cheaper
 * than the clip path it replaces.
 */
export function half(points, keep_right = true) {
    const inside = ([x]) => (keep_right ? x >= 0 : x <= 0);
    const cross = (a, b) => {
        // Where the segment crosses x = 0.
        const t = a[0] / (a[0] - b[0]);
        return [0, a[1] + t * (b[1] - a[1])];
    };

    const out = [];
    for (let i = 0; i < points.length; ++i) {
        const a = points[i];
        const b = points[(i + 1) % points.length];
        const a_in = inside(a);
        const b_in = inside(b);
        if (a_in) out.push(a);
        if (a_in !== b_in) out.push(cross(a, b));
    }
    return out;
}

export const to_points = (pts) =>
    pts.map(([x, y]) => `${x.toFixed(2)},${y.toFixed(2)}`).join(" ");
