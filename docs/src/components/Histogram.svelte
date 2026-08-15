<script module>
    let instances = 0;
</script>

<script>
    /**
     * Distribution of how far points moved: a histogram over a box plot, both on
     * one shared scale, so a bump in the tail lines up vertically with the
     * outlier dots that explain it.
     *
     * The mean lives on the box plot rather than on the histogram. Drawn over the
     * bars it was a rule cutting through the data it was meant to summarise; on
     * the box it sits beside the median, which is the comparison actually worth
     * making — the gap between the two is what tells you the tail is dragging.
     */
    import { ramp_gradient } from "../lib/colour.js";

    let { values = [], summary = null, unit = "cells" } = $props();

    // Unique per instance so two charts on one page cannot share a gradient def.
    const gid = `hagrid-ramp-${(instances += 1)}`;

    const BINS = 44;
    const ML = 34;
    const MR = 10;
    const MT = 6;
    const HIST = 74;
    const BOX = 30;
    const RAMP_H = 8;
    const MB = 22;
    const HEIGHT = MT + HIST + BOX + RAMP_H + MB;
    const ramp_y = MT + HIST + BOX + 2;

    let body = $state(null);
    let width = $state(720);
    let hover = $state(null);

    const s = $derived(
        summary ?? {
            mean: 0, median: 0, p95: 0, max: 0, min: 0,
            q1: 0, q3: 0, lo_fence: 0, hi_fence: 0, outliers: [],
        },
    );
    const top = $derived(Math.max(s.max, 1e-9));
    const counts = $derived.by(() => {
        const c = new Array(BINS).fill(0);
        for (const v of values) c[Math.min(BINS - 1, Math.floor((v / top) * BINS))]++;
        return c;
    });
    const peak = $derived(Math.max(1, ...counts));

    const pw = $derived(Math.max(80, width - ML - MR));
    const bw = $derived(pw / BINS);
    const X = $derived((v) => ML + Math.min(1, v / top) * pw);
    const box_y = MT + HIST + BOX / 2;

    function on_move(event) {
        const rect = body.getBoundingClientRect();
        const px = ((event.clientX - rect.left) / rect.width) * width;
        const bin = Math.floor((px - ML) / bw);
        if (bin < 0 || bin >= BINS || !counts[bin]) return (hover = null);
        hover = {
            bin,
            x: ((ML + (bin + 0.5) * bw) / width) * rect.width,
            lo: (bin / BINS) * top,
            hi: ((bin + 1) / BINS) * top,
            count: counts[bin],
        };
    }

    $effect(() => {
        if (!body) return;
        const observer = new ResizeObserver((entries) => {
            width = Math.max(240, Math.round(entries[0].contentRect.width));
        });
        observer.observe(body);
        return () => observer.disconnect();
    });

    const diamond = $derived((x, y, r) => `${x},${y - r} ${x + r},${y} ${x},${y + r} ${x - r},${y}`);

    /*
     * The bars are filled from the same ramp the plot uses to colour points by
     * movement, laid across the axis in user space so each bar picks up the
     * colour at its own displacement. A bar spans a range of values, so a flat
     * fill per bin would be a small lie; one gradient across the whole plot
     * gets the ends of every bar right for free.
     *
     * The ramp saturates at the 95th percentile, exactly as the point colouring
     * does, so the stops are positioned against p95 and the last one is held out
     * to the end of the axis. Without that the chart and the plot would disagree
     * about what a colour means.
     */
    const ramp = $derived(ramp_gradient(s.p95, top));
</script>

<figure class="hist">
    <div
        class="body"
        bind:this={body}
        onpointermove={on_move}
        onpointerleave={() => (hover = null)}
        role="img"
        aria-label="Distribution of how far each point moved, in {unit}. Mean {s.mean.toFixed(
            2,
        )}, median {s.median.toFixed(2)}, quartiles {s.q1.toFixed(2)} to {s.q3.toFixed(
            2,
        )}, longest move {s.max.toFixed(2)}."
    >
        <svg viewBox="0 0 {width} {HEIGHT}" style="height: {HEIGHT}px">
            <defs>
                <linearGradient id={gid} gradientUnits="userSpaceOnUse" x1={ML} x2={ML + pw}>
                    {#each ramp as stop, i (i)}
                        <stop offset="{(stop.offset * 100).toFixed(2)}%" stop-color={stop.colour} />
                    {/each}
                </linearGradient>
            </defs>

            {#each [peak, Math.round(peak / 2)] as v (v)}
                {#if v}
                    <line class="axis" x1={ML} x2={ML + pw} y1={MT + HIST - (v / peak) * HIST} y2={MT + HIST - (v / peak) * HIST} />
                    <text class="tick" x={ML - 6} y={MT + HIST - (v / peak) * HIST + 3} text-anchor="end">{v}</text>
                {/if}
            {/each}

            {#each counts as c, b (b)}
                {#if c}
                    <rect
                        class="bar"
                        class:on={hover?.bin === b}
                        fill="url(#{gid})"
                        x={ML + b * bw + 0.75}
                        y={MT + HIST - (c / peak) * HIST}
                        width={Math.max(0.75, bw - 1.5)}
                        height={(c / peak) * HIST}
                    />
                {/if}
            {/each}
            <line class="axis" x1={ML} x2={ML + pw} y1={MT + HIST} y2={MT + HIST} />

            <!-- box plot, on the histogram's own scale -->
            <line class="whisker" x1={X(s.lo_fence)} x2={X(s.q1)} y1={box_y} y2={box_y} />
            <line class="whisker" x1={X(s.q3)} x2={X(s.hi_fence)} y1={box_y} y2={box_y} />
            <line class="cap" x1={X(s.lo_fence)} x2={X(s.lo_fence)} y1={box_y - 4} y2={box_y + 4} />
            <line class="cap" x1={X(s.hi_fence)} x2={X(s.hi_fence)} y1={box_y - 4} y2={box_y + 4} />
            <rect class="box" x={X(s.q1)} y={box_y - 7} width={Math.max(1, X(s.q3) - X(s.q1))} height={14} />
            <line class="median" x1={X(s.median)} x2={X(s.median)} y1={box_y - 7} y2={box_y + 7} />
            {#each s.outliers as o, i (i)}
                <circle class="outlier" fill="url(#{gid})" cx={X(o)} cy={box_y} r="2" />
            {/each}
            <polygon class="mean" points={diamond(X(s.mean), box_y, 4.5)}>
                <title>mean {s.mean.toFixed(2)} {unit}</title>
            </polygon>

            <!-- The colour scale, drawn on the axis it encodes. Sharing the
                 geometry is the legend: a bar's colour and its position carry
                 the same number, so there is nothing to cross-reference. -->
            <rect
                class="ramp-key"
                fill="url(#{gid})"
                x={ML}
                y={ramp_y}
                width={pw}
                height={RAMP_H}
            />
            {#if s.p95 > 0 && s.p95 < top}
                <line
                    class="cap-rule"
                    x1={X(s.p95)}
                    x2={X(s.p95)}
                    y1={ramp_y - 2}
                    y2={ramp_y + RAMP_H + 2}
                />
            {/if}

            {#each [0, 0.25, 0.5, 0.75, 1] as t (t)}
                <text
                    class="tick"
                    x={ML + t * pw}
                    y={HEIGHT - 6}
                    text-anchor={t === 0 ? "start" : t === 1 ? "end" : "middle"}>{(t * top).toFixed(2)}</text
                >
            {/each}
        </svg>

        {#if hover}
            <div class="tip" style="left: {hover.x}px">
                {hover.lo.toFixed(2)}–{hover.hi.toFixed(2)}
                {unit} · {hover.count}
                {hover.count === 1 ? "point" : "points"}
            </div>
        {/if}
    </div>

    <figcaption>
        <b>How far points moved</b>, in {unit}. The bar under the axis is the movement colour
        scale, which the bars above are filled from and which colours the plot whenever Colour is
        set to Movement{#if s.p95 > 0 && s.p95 < top}; it saturates at the 95th percentile, marked
        on the scale at {s.p95.toFixed(2)}{/if}. The box spans the middle half, {s.q1.toFixed(2)} to
        {s.q3.toFixed(2)}, with the median at {s.median.toFixed(2)} and the mean
        <span class="key-mean">◆</span> at {s.mean.toFixed(2)}. Whiskers reach 1.5 × IQR; the
        {s.outliers.length}
        {s.outliers.length === 1 ? "dot" : "dots"} beyond
        {s.outliers.length === 1 ? "is a point that" : "are the points that"} travelled furthest, out
        to {s.max.toFixed(2)}.
    </figcaption>
</figure>

<style>
    .hist {
        margin: 0.55rem 0 0;
    }
    .body {
        position: relative;
    }
    svg {
        display: block;
        width: 100%;
        touch-action: none;
    }
    /* No `fill` here: an SVG presentation attribute loses to any CSS rule, so a
       fill declared on .bar would quietly win over the gradient. Hover is a
       stroke rather than a fill swap, which would throw away the encoding at
       exactly the moment the reader is inspecting it. */
    .bar {
        stroke: transparent;
        stroke-width: 1;
    }
    .bar.on {
        stroke: var(--sl-color-white);
    }
    .axis {
        stroke: var(--sl-color-gray-5);
    }
    .tick {
        fill: var(--sl-color-gray-3);
        font-size: 10px;
        font-variant-numeric: tabular-nums;
    }
    .box {
        fill: color-mix(in srgb, var(--sl-color-white) 10%, transparent);
        stroke: var(--sl-color-gray-4);
    }
    .ramp-key {
        rx: 1;
    }
    .cap-rule {
        stroke: var(--sl-color-white);
        stroke-width: 1.5;
    }
    .median {
        stroke: var(--sl-color-white);
        stroke-width: 2;
    }
    .whisker,
    .cap {
        stroke: var(--sl-color-gray-4);
    }
    .outlier {
        stroke: var(--sl-color-black);
        stroke-width: 0.75;
    }
    /* Mean and median are structural marks, not data, so they stay achromatic —
       a warm diamond would read as another point at the far end of the ramp. */
    .mean {
        fill: var(--sl-color-white);
        stroke: var(--sl-color-black);
        stroke-width: 1;
    }
    .key-mean {
        color: var(--sl-color-white);
    }
    .tip {
        position: absolute;
        top: -0.4rem;
        transform: translate(-50%, -100%);
        pointer-events: none;
        white-space: nowrap;
        background: var(--sl-color-gray-6);
        border: 1px solid var(--sl-color-gray-5);
        border-radius: 3px;
        padding: 0.18rem 0.4rem;
        font-size: 0.68rem;
        font-variant-numeric: tabular-nums;
    }
    figcaption {
        margin-top: 0.3rem;
        font-size: 0.72rem;
        line-height: 1.55;
        color: var(--sl-color-gray-3);
    }
    figcaption b {
        color: var(--sl-color-gray-2);
        font-weight: 600;
    }
</style>
