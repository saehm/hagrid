/**
 * Choosing how to measure similarity, before choosing how to project.
 *
 * Euclidean distance on raw columns is a modelling assumption, not a default,
 * and it is wrong in different ways for different data. The approach here is to
 * make the metric a *transform*: if Euclidean distance in the transformed space
 * is the distance you meant, then everything downstream — the projection and
 * the quality metrics that score it — works unchanged and honestly.
 *
 * A preservation metric cannot tell you which transform is right, since each
 * one measures faithfulness to its own notion of similarity. What it can do is
 * pick the best projection *given* a choice. To compare the choices themselves
 * we need something external, so the class labels are used for that and for
 * nothing else.
 */

export const zscore = (rows) => {
    const d = rows[0].length;
    const n = rows.length;
    const mean = Array.from({ length: d }, (_, j) => rows.reduce((s, r) => s + r[j], 0) / n);
    const std = Array.from({ length: d }, (_, j) =>
        Math.sqrt(rows.reduce((s, r) => s + (r[j] - mean[j]) ** 2, 0) / n) || 1,
    );
    return rows.map((r) => r.map((v, j) => (v - mean[j]) / std[j]));
};

/**
 * Centred log-ratio, the standard geometry for parts of a whole.
 *
 * Compositional data carries no information in its total — only in the ratios
 * between parts — so a plain euclidean distance on percentages counts the
 * closure as if it were signal. CLR maps the simplex into a real space where
 * euclidean distance *is* the Aitchison distance.
 *
 * Every part must be strictly positive, which is the catch: a structural zero
 * is not a small number, it is an absence, and no replacement value is honest
 * about that. Only pass parts that are always present.
 */
export function clr(rows) {
    return rows.map((r) => {
        const logs = r.map((v) => Math.log(v));
        const mean = logs.reduce((a, b) => a + b, 0) / logs.length;
        return logs.map((v) => v - mean);
    });
}

/** Splits a matrix into the named column groups. */
export const take = (rows, indices) => rows.map((r) => indices.map((j) => r[j]));

/** Side-by-side, so differently-treated blocks can be recombined. */
export const concat = (...blocks) =>
    blocks[0].map((_, i) => blocks.flatMap((b) => b[i]));

/**
 * How often a point's k nearest neighbours share its label.
 *
 * External to every transform, so it is the one number that can be compared
 * across representations. It says whether a notion of similarity agrees with
 * what the data is actually about — not whether a projection is faithful.
 */
export function label_agreement(rows, labels, k = 10) {
    const n = rows.length;
    let total = 0;
    for (let i = 0; i < n; ++i) {
        const near = rows
            .map((r, j) => {
                let d = 0;
                for (let c = 0; c < r.length; ++c) d += (r[c] - rows[i][c]) ** 2;
                return [d, j];
            })
            .sort((a, b) => a[0] - b[0])
            .slice(1, k + 1);
        total += near.filter(([, j]) => labels[j] === labels[i]).length / k;
    }
    return total / n;
}
