/**
 * Scoring projections with `@saehrimnir/sickle`.
 *
 * No single metric decides this. They disagree by construction, and each one
 * can be gamed by a projection that is bad in the way it does not measure:
 *
 *   auc log R_NX   rank preservation across every neighbourhood size at once
 *   trustworthiness   are the neighbours it shows real ones?      (false positives)
 *   continuity        did it keep the neighbours that existed?    (false negatives)
 *   spearman ρ        do distances stay in the same order overall? (global shape)
 *   neighbourhood hit do same-label points end up together?        (class structure)
 *
 * Trustworthiness and continuity trade against each other; a projection that
 * shatters a cluster scores well on the first and badly on the second. Spearman
 * catches layouts that are locally tidy but globally scrambled. The hit rate is
 * the only label-aware one, so it is included but never allowed to decide alone
 * — optimising for it would just be asking which projection best separates
 * classes we already know.
 *
 * Ranking is therefore by mean rank across the metrics rather than by any one.
 */
import * as sickle from "@saehrimnir/sickle";

const K = 12;

/** The grid to search. */
export function candidates() {
    const out = [];
    for (const seed of [1, 2]) {
        for (const n_neighbors of [10, 25]) {
            for (const MN_ratio of [0.5, 1.0]) {
                out.push({
                    method: "PaCMAP",
                    params: { seed, n_neighbors, MN_ratio, FP_ratio: 2.0 },
                });
            }
        }
        for (const n_neighbors of [10, 25]) {
            for (const min_dist of [0.05, 0.3]) {
                out.push({ method: "UMAP", params: { seed, n_neighbors, min_dist } });
            }
        }
        for (const perplexity of [20, 40]) {
            out.push({ method: "TSNE", params: { seed, perplexity, epsilon: 5 } });
        }
    }
    return out;
}

/** Some metrics answer with `{ value, … }` and some with a bare number. */
const scalar = (r) => (typeof r === "number" ? r : r.value);

/**
 * Every metric for one projection. Higher is better throughout.
 *
 * `labels` is optional, and should be left out when the labels have little to
 * do with the space being projected — ranking on a hit rate that is barely
 * above chance is ranking on noise, and it drags the genuinely informative
 * metrics around with it.
 */
export function evaluate(high, low, labels = null) {
    const cr = sickle.coRanking(high, low);
    const out = {
        auc: scalar(sickle.aucLogRnx(cr)),
        trust: scalar(sickle.trustworthiness(cr, K)),
        cont: scalar(sickle.continuity(cr, K)),
        rho: Math.abs(scalar(sickle.spearmanRho(high, low))),
    };
    for (const [name, values] of Object.entries(labels ?? {})) {
        out[`hit_${name}`] = scalar(sickle.neighborhoodHit(low, values));
    }
    return out;
}

/**
 * Ranks every candidate on each metric, then sorts by the mean rank. Prints the
 * table and returns the winner.
 */
export function report(scored) {
    const METRICS = Object.keys(scored[0].scores);
    const ranks = new Map(scored.map((s) => [s, []]));
    for (const m of METRICS) {
        [...scored]
            .sort((a, b) => b.scores[m] - a.scores[m])
            .forEach((s, i) => ranks.get(s).push(i + 1));
    }
    for (const s of scored) {
        const r = ranks.get(s);
        s.mean_rank = r.reduce((a, b) => a + b, 0) / r.length;
    }
    scored.sort((a, b) => a.mean_rank - b.mean_rank);

    const head = ["method".padEnd(8), "parameters".padEnd(46)];
    for (const m of METRICS) head.push(m.padStart(11));
    head.push("mean rank".padStart(10));
    console.log(`\n${head.join("  ")}`);

    for (const s of scored.slice(0, 12)) {
        const row = [s.method.padEnd(8), JSON.stringify(s.params).padEnd(46)];
        for (const m of METRICS) row.push(s.scores[m].toFixed(3).padStart(11));
        row.push(s.mean_rank.toFixed(1).padStart(10));
        console.log(row.join("  "));
    }
    console.log(`… ${scored.length} candidates evaluated, best ${scored.length > 12 ? "12 " : ""}shown\n`);

    const best = scored[0];
    console.log(`chosen: ${best.method} ${JSON.stringify(best.params)}  (mean rank ${best.mean_rank.toFixed(1)})`);
    return best;
}
