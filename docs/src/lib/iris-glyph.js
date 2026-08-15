/**
 * Geometry and palette for the iris flower glyph, shared by the glyph itself
 * and by the legend that explains it — so the legend cannot drift from what is
 * actually drawn.
 *
 * Transcribed from the VINCI'21 talk (renecutura.eu/vinci21_hagrid).
 */
import { IRIS } from "@saehrimnir/mistle";

export const PETAL_FILL = {
    "Iris-setosa": "#5d5faa",
    "Iris-versicolor": "#67479c",
    "Iris-virginica": "#aa80b9",
    default: "purple",
};

export const SEPAL_FILL = {
    "Iris-setosa": "#b8bcdf",
    "Iris-versicolor": "#baaed5",
    "Iris-virginica": "#a893bd",
    default: "mediumpurple",
};

/** A neutral flower, for when no row is named. */
export const EXAMPLE_ROW = [5, 3, 8, 3];

const COLUMN_MAX = IRIS.columns.map((c) => IRIS.statistics[c].max);

/**
 * Each measurement against its own column's maximum across the dataset, mapped
 * to [0, 5]. Scaling per column rather than per row is what makes two glyphs
 * comparable: a bigger flower on screen is a bigger flower.
 */
export const scale = (value, column) => (value / COLUMN_MAX[column]) * 5;

/** `[sepal_length, sepal_width, petal_length, petal_width]`, in glyph units. */
export function arms(row, size) {
    return [
        scale(row[0], 0) * size,
        scale(row[1], 1) * size,
        scale(row[2], 2) * size,
        scale(row[3], 3) * size,
    ];
}

/** Petals sit on a 4-unit stem, sepals on a 2-unit one. */
export const PETAL_STEM = 4;
export const SEPAL_STEM = 2;

export const PETAL_ANGLES = [0, 120, 240];
export const SEPAL_ANGLES = [60, 180, 300];

/**
 * One arm, pointing down the +y axis: a stem of `stem × size`, then a blade
 * `len` long and `wid` across at its widest.
 *
 * Petals close with quadratics and sepals with cubics, which is what gives the
 * sepals their narrower tip.
 */
export function petal_path(size, len, wid) {
    const s = PETAL_STEM * size;
    return `M ${size} 0
         v ${s}
         Q ${wid / 2} ${s} ${wid / 2} ${s + len / 2}
         Q ${wid / 2} ${s + len} 0 ${s + len}
         Q ${-wid / 2} ${s + len} ${-wid / 2} ${s + len / 2}
         Q ${-wid / 2} ${s} ${-size} ${s}
         v ${-s} Z`;
}

export function sepal_path(size, len, wid) {
    const s = SEPAL_STEM * size;
    return `M ${size} 0
         v ${s}
         Q ${wid / 2} ${s} ${wid / 2} ${s + len / 2}
         C ${wid / 2} ${s + len / 2 + len / 5} ${wid / 5} ${s + len} 0 ${s + len}
         C ${-wid / 5} ${s + len} ${-wid / 2} ${s + len / 2 + len / 5} ${-wid / 2} ${s + len / 2}
         Q ${-wid / 2} ${s} ${-size} ${s}
         v ${-s} Z`;
}

/**
 * How far a glyph reaches from its centre, per unit of `size`. Multiply by the
 * size to get the radius; two glyphs whose centres are further apart than the
 * sum of their radii cannot overlap.
 */
export function reach(row) {
    return Math.max(PETAL_STEM + scale(row[2], 2), SEPAL_STEM + scale(row[0], 0));
}

/** Blade length, blade width and stem length of one arm, in glyph units. */
export function arm_geometry(kind, row, size) {
    const petal = kind === "petal";
    return {
        len: petal ? scale(row[2], 2) * size : scale(row[0], 0) * size,
        wid: petal ? scale(row[3], 3) * size : scale(row[1], 1) * size,
        stem: (petal ? PETAL_STEM : SEPAL_STEM) * size,
    };
}

/**
 * A viewBox and dimension-line positions that suit *both* arms at once.
 *
 * Petals and sepals differ in every measurement, so left to themselves the two
 * diagrams get different viewBoxes, render at different scales once stretched
 * to the same column width, and end at different heights — which puts their
 * captions on different baselines. Sharing a frame anchors both at the flower's
 * centre, draws them to one scale, and lines the annotations up across the pair.
 */
export function anatomy_frame(row, size, { label_width = 24, gap = 8, pad = 6 } = {}) {
    const p = arm_geometry("petal", row, size);
    const s = arm_geometry("sepal", row, size);
    const wid = Math.max(p.wid, s.wid);
    const bottom = Math.max(p.stem + p.len, s.stem + s.len);

    const dim_x = wid / 2 + gap;
    const dim_y = bottom + gap + 2;
    const left = -wid / 2 - pad;
    const right = dim_x + 5 + label_width;
    const top = -pad;
    const height = dim_y + 10 + pad;

    return { dim_x, dim_y, view: `${left} ${top} ${right - left} ${height - top}` };
}
