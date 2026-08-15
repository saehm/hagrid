<script>
    /** Heartbeats drawn as heartbeats. See `scripts/build-ecg.mjs`. */
    import GlyphShowcase from "./GlyphShowcase.svelte";
    import SparklineGlyph from "./SparklineGlyph.svelte";
    import { box_reach } from "../lib/glyph-reach.js";
    import ecg from "../data/ecg-2d.json";

    let { method = "as projected", ...rest } = $props();

    const CLASS_FILL = {
        Normal: "#7dcfff",
        "R-on-T PVC": "#f7768e",
        PVC: "#e0af68",
        Supraventricular: "#9ece6a",
        Unclassifiable: "#bb9af7",
    };

    // A shared vertical range, so a tall beat is a tall beat. Clipped at the
    // 0.2nd/99.8th percentile: a handful of extreme samples would otherwise
    // flatten every other trace to a line.
    const RANGE = (() => {
        const all = ecg.series.flat().sort((a, b) => a - b);
        const at = (q) => all[Math.floor(q * (all.length - 1))];
        return [at(0.002), at(0.998)];
    })();

    // A trace fills its box, so it reaches the corner only off-lattice.
    const REACH = box_reach(ecg.points.length);

    /** The beat nearest each class's average shape. */
    const EXEMPLAR = ecg.classes.map((name) => {
        const rows = ecg.series.filter((_, i) => ecg.labels[i] === name);
        const mean = ecg.series[0].map((_, j) => rows.reduce((s, r) => s + r[j], 0) / rows.length);
        let best = -1, bd = Infinity;
        ecg.series.forEach((r, i) => {
            if (ecg.labels[i] !== name) return;
            const d = r.reduce((s, v, j) => s + (v - mean[j]) ** 2, 0);
            if (d < bd) { bd = d; best = i; }
        });
        return { name, index: best, count: rows.length };
    });
</script>

<GlyphShowcase points={ecg.points} reach={REACH} {method} {...rest}>
    {#snippet glyph(i, size)}
        <SparklineGlyph
            series={ecg.series[i]}
            range={RANGE}
            {size}
            stroke={CLASS_FILL[ecg.labels[i]] ?? "#a1b5d8"}
        />
    {/snippet}

    {#snippet details(i)}
        <p class="tip-title" style="color: {CLASS_FILL[ecg.labels[i]]}">{ecg.labels[i]}</p>
        <svg viewBox="-26 -18 52 36" class="tip-trace" aria-hidden="true">
            <SparklineGlyph
                series={ecg.series[i]}
                range={RANGE}
                size={16}
                stroke={CLASS_FILL[ecg.labels[i]] ?? "#a1b5d8"}
                guides
            />
        </svg>
        <dl>
            <dt>beat</dt><dd>#{i + 1}</dd>
            <dt>samples</dt><dd>{ecg.series[i].length}</dd>
            <dt>peak</dt><dd>{Math.max(...ecg.series[i]).toFixed(2)}</dd>
            <dt>trough</dt><dd>{Math.min(...ecg.series[i]).toFixed(2)}</dd>
        </dl>
    {/snippet}

    {#snippet legend()}
        <div class="legend">
            <p class="legend-title">Reading a beat</p>
            <div class="reading">
                <svg viewBox="-22 -22 44 44" aria-hidden="true">
                    <SparklineGlyph
                        series={ecg.series[EXEMPLAR[0].index]}
                        range={RANGE}
                        size={18}
                        stroke={CLASS_FILL[EXEMPLAR[0].name]}
                    />
                </svg>
                <p>
                    140 samples of one heartbeat, left to right. Nothing is encoded — this is the
                    signal the projection was computed from, drawn small. Every trace shares one
                    vertical scale, so heights are comparable between glyphs.
                </p>
            </div>
        </div>

        <div class="legend">
            <p class="legend-title">Beat types</p>
            <ul class="classes">
                {#each EXEMPLAR as e}
                    <li>
                        <svg viewBox="-26 -16 52 32" aria-hidden="true">
                            <SparklineGlyph
                                series={ecg.series[e.index]}
                                range={RANGE}
                                size={15}
                                stroke={CLASS_FILL[e.name]}
                                guides
                            />
                        </svg>
                        <span class="name" style="color: {CLASS_FILL[e.name]}">{e.name}</span>
                        <span class="count">{e.count} beats</span>
                    </li>
                {/each}
            </ul>
            <p class="foot">
                The beat nearest each type's average shape. A normal beat has a narrow spike; the
                ventricular ones are broader and displaced, which is visible in the glyph before
                you know any of the labels.
            </p>
        </div>
    {/snippet}
</GlyphShowcase>

<style>
    .legend-title {
        margin: 0 0 0.5rem;
        font-size: var(--sl-text-xs);
        text-transform: uppercase;
        letter-spacing: 0.08em;
        color: var(--sl-color-gray-3);
    }
    .reading { display: flex; gap: 0.9rem; align-items: flex-start; }
    .reading svg { width: 4.6rem; flex: none; color: var(--sl-color-gray-2); }
    .reading p {
        margin: 0;
        font-size: var(--sl-text-xs);
        color: var(--sl-color-gray-2);
        line-height: 1.45;
    }
    /* One row per type: the names are long and uneven ("Supraventricular" next
       to "PVC"), so side-by-side columns either wrap mid-word or leave a ragged
       gap. A list reads cleanly at any width. */
    .classes {
        display: flex;
        flex-direction: column;
        gap: 0.15rem;
        margin: 0;
        padding: 0;
        list-style: none;
    }
    .classes li {
        display: grid;
        grid-template-columns: 3.4rem 1fr auto;
        align-items: center;
        gap: 0.5rem;
        padding: 0.1rem 0;
        font-size: var(--sl-text-xs);
    }
    .classes svg { display: block; width: 100%; height: auto; color: var(--sl-color-gray-2); }
    .classes .name { font-weight: 600; }
    .classes .count { color: var(--sl-color-gray-4); font-variant-numeric: tabular-nums; }
    .tip-trace { display: block; width: 8rem; height: auto; margin: 0.1rem 0 0.35rem; }
    .foot { margin: 0.6rem 0 0; font-size: var(--sl-text-xs); color: var(--sl-color-gray-3); }
</style>
