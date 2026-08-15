<script>
    /**
     * Spoken digits: every point is a recording, drawn as its spectrogram, and
     * clicking one plays the actual audio.
     *
     * Nothing is synthesised here. The clips stream from the dataset's own
     * repository on demand, so the page ships no audio and loads none until you
     * ask for it.
     *
     * The projection carries two structures at once — what was said and who
     * said it — so the colouring can be switched between them.
     */
    import GlyphShowcase from "./GlyphShowcase.svelte";
    import SpectrogramGlyph from "./SpectrogramGlyph.svelte";
    import { box_reach } from "../lib/glyph-reach.js";
    import fsdd from "../data/fsdd-2d.json";

    let { method = "as projected", ...rest } = $props();

    const DIGIT_FILL = {
        0: "#7aa2f7", 2: "#9ece6a", 4: "#e0af68", 6: "#f7768e", 8: "#bb9af7",
    };
    const SPEAKER_FILL = {
        george: "#7dcfff", jackson: "#ff9e64", lucas: "#73daca", nicolas: "#e06c9f",
    };

    // A spectrogram fills a square of side 2 x size, so its corner is sqrt(2) out.
    // A spectrogram fills its box, so it reaches the corner only off-lattice.
    const REACH = box_reach(fsdd.points.length);

    let colour_by = $state("digit");
    let playing = $state(-1);

    let colour = $derived(
        colour_by === "digit"
            ? (i) => DIGIT_FILL[fsdd.digit[i]]
            : (i) => SPEAKER_FILL[fsdd.speaker[i]],
    );

    /** One element reused, so a quick series of clicks does not pile up. */
    let element = null;
    function hear(i) {
        if (typeof Audio === "undefined") return;
        if (!element) element = new Audio();
        element.src = `${fsdd.audio_base}/${fsdd.files[i]}`;
        element.currentTime = 0;
        element.play().catch(() => { /* blocked or offline; nothing to do */ });
        playing = i;
    }

    /** One example per class, for the legend. */
    const first_of = (pred) => fsdd.points.findIndex((_, i) => pred(i));
    const DIGIT_EXAMPLES = fsdd.digits.map((d) => ({
        label: String(d),
        index: first_of((i) => fsdd.digit[i] === d),
        fill: DIGIT_FILL[d],
    }));
    const SPEAKER_EXAMPLES = fsdd.speakers.map((s) => ({
        label: s,
        index: first_of((i) => fsdd.speaker[i] === s),
        fill: SPEAKER_FILL[s],
    }));
    let examples = $derived(colour_by === "digit" ? DIGIT_EXAMPLES : SPEAKER_EXAMPLES);
</script>

<GlyphShowcase points={fsdd.points} reach={REACH} {method} {...rest}>
    {#snippet glyph(i, size)}
        <g
            role="button"
            tabindex="0"
            aria-label="Play {fsdd.digit[i]} spoken by {fsdd.speaker[i]}"
            class="hit"
            onclick={() => hear(i)}
            onkeydown={(e) => {
                if (e.key === "Enter" || e.key === " ") { e.preventDefault(); hear(i); }
            }}
        >
            <rect
                x={-size * 1.4} y={-size * 1.4}
                width={size * 2.8} height={size * 2.8}
                fill="transparent"
            />
            <SpectrogramGlyph src={fsdd.images[i]} {size} colour={colour(i)} />
        </g>
    {/snippet}

    {#snippet details(i)}
        <p class="tip-title" style="color: {colour(i)}">
            “{fsdd.digit[i]}” — {fsdd.speaker[i]}
        </p>
        <img class="tip-spectrogram" src={fsdd.images[i]} alt="" />
        <p class="tip-note">Click to hear it.</p>
    {/snippet}

    {#snippet legend()}
        <div class="legend">
            <p class="legend-title">Reading a recording</p>
            <div class="reading">
                <svg viewBox="-22 -22 44 44" aria-hidden="true">
                    <SpectrogramGlyph
                        src={fsdd.images[examples[0].index]}
                        size={18}
                        colour={examples[0].fill}
                    />
                </svg>
                <p>
                    A spectrogram: time left to right, frequency low at the bottom to high at the
                    top, brightness for energy. The frame is the class.
                </p>
            </div>
            <p class="foot">
                <strong>Click any glyph to hear it.</strong> These are the original recordings,
                streamed from the dataset's repository — nothing here is synthesised, and no
                audio is downloaded until you ask for it.
                {#if playing >= 0}
                    <br />Last played: <strong>{fsdd.digit[playing]}</strong>, spoken by
                    {fsdd.speaker[playing]}.
                {/if}
            </p>
        </div>

        <div class="legend">
            <p class="legend-title">Colour by</p>
            <div class="switch">
                <button
                    type="button" class:on={colour_by === "digit"}
                    onclick={() => (colour_by = "digit")}
                >digit</button>
                <button
                    type="button" class:on={colour_by === "speaker"}
                    onclick={() => (colour_by = "speaker")}
                >speaker</button>
            </div>
            <ul class="classes">
                {#each examples as e}
                    <li>
                        <button type="button" onclick={() => hear(e.index)}>
                            <svg viewBox="-22 -22 44 44" aria-hidden="true">
                                <SpectrogramGlyph src={fsdd.images[e.index]} size={18} colour={e.fill} />
                            </svg>
                            <span style="color: {e.fill}">{e.label}</span>
                        </button>
                    </li>
                {/each}
            </ul>
            <p class="foot">
                The layout does not change — only the colouring. Both structures are in the same
                projection, so switching shows which one a given clump is about.
            </p>
        </div>
    {/snippet}
</GlyphShowcase>

<style>
    .hit { cursor: pointer; }
    .hit:focus-visible { outline: 2px solid var(--sl-color-accent-high); outline-offset: 2px; }
    .legend-title {
        margin: 0 0 0.5rem;
        font-size: var(--sl-text-xs);
        text-transform: uppercase;
        letter-spacing: 0.08em;
        color: var(--sl-color-gray-3);
    }
    .reading { display: flex; gap: 0.9rem; align-items: flex-start; }
    .reading svg { width: 4.4rem; flex: none; }
    .reading p {
        margin: 0;
        font-size: var(--sl-text-xs);
        color: var(--sl-color-gray-2);
        line-height: 1.45;
    }
    .switch { display: flex; gap: 0.3rem; margin-bottom: 0.6rem; }
    .switch button {
        flex: 1 1 0;
        padding: 0.25rem 0.4rem;
        font-size: var(--sl-text-xs);
        color: var(--sl-color-gray-2);
        background: none;
        border: 1px solid var(--sl-color-gray-5);
        border-radius: 0.3rem;
        cursor: pointer;
    }
    .switch button.on {
        color: var(--sl-color-white);
        border-color: var(--sl-color-accent);
        background: color-mix(in srgb, var(--sl-color-accent) 18%, transparent);
    }
    .classes {
        display: grid;
        grid-template-columns: repeat(auto-fill, minmax(3.6rem, 1fr));
        gap: 0.4rem;
        margin: 0;
        padding: 0;
        list-style: none;
    }
    .classes button {
        display: grid;
        grid-template-rows: auto auto;
        justify-items: center;
        gap: 0.15rem;
        width: 100%;
        padding: 0.15rem;
        background: none;
        border: 1px solid transparent;
        border-radius: 0.3rem;
        cursor: pointer;
        font-size: var(--sl-text-xs);
    }
    .classes button:hover { border-color: var(--sl-color-gray-5); }
    .classes svg { display: block; width: 100%; max-width: 3rem; height: auto; }
    .foot { margin: 0.6rem 0 0; font-size: var(--sl-text-xs); color: var(--sl-color-gray-3); }
    .tip-spectrogram {
        display: block;
        width: 8rem;
        height: 8rem;
        image-rendering: pixelated;
        border-radius: 2px;
    }
    .tip-note { margin: 0.35rem 0 0; color: var(--sl-color-gray-4); }
</style>
