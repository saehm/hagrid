<script>
    /**
     * Draws a layout and the points it came from in one shared transform.
     *
     * The transform is isotropic on purpose. Fitting x and y independently would
     * make the canvas tidier and every cell a lie — square cells would render as
     * rectangles and hexagons as ellipsoids, which is exactly the property most
     * of these pages are about.
     */
    let {
        source = [],
        positions = [],
        grid = null,
        colours = null,
        show_grid = true,
        show_links = true,
        show_source = true,
        height = null,
        highlight = -1,
    } = $props();

    let canvas = $state(null);
    let container = $state(null);
    let width = $state(600);

    /*
     * The canvas takes its shape from the layout it has to hold, not a fixed
     * 3:2. The scale is isotropic — it has to be, or square cells would render
     * as rectangles — so a 100:1 dataset in a 3:2 box is a thin line of points
     * across an empty field. Matching the box to the content instead spends the
     * height on the data. Clamped at both ends so a degenerate aspect cannot
     * collapse the plot or run it off the screen.
     */
    const box_height = $derived.by(() => {
        if (height) return height;
        const { x0, y0, x1, y1 } = extent_of();
        const aspect = (x1 - x0) / (y1 - y0 || 1);
        const fitted = Number.isFinite(aspect) && aspect > 0 ? width / aspect : width * 0.62;
        return Math.round(Math.min(520, Math.max(200, fitted)));
    });

    const PAD = 14;

    function extent_of() {
        let x0 = Infinity;
        let y0 = Infinity;
        let x1 = -Infinity;
        let y1 = -Infinity;
        const consider = (x, y) => {
            if (x < x0) x0 = x;
            if (x > x1) x1 = x;
            if (y < y0) y0 = y;
            if (y > y1) y1 = y;
        };
        for (const [x, y] of positions) consider(x, y);
        if (show_source) for (const [x, y] of source) consider(x, y);
        if (grid?.type === "regions") {
            for (const c of grid.cells) {
                consider(c.x, c.y);
                consider(c.x + c.width, c.y + c.height);
            }
        }
        if (!Number.isFinite(x0)) return { x0: 0, y0: 0, x1: 1, y1: 1 };
        // A degenerate axis would divide by zero; give it something to occupy.
        if (x1 - x0 < 1e-9) ((x0 -= 0.5), (x1 += 0.5));
        if (y1 - y0 < 1e-9) ((y0 -= 0.5), (y1 += 0.5));
        return { x0, y0, x1, y1 };
    }

    function css(name, fallback) {
        if (!container) return fallback;
        const v = getComputedStyle(container).getPropertyValue(name).trim();
        return v || fallback;
    }

    function draw() {
        if (!canvas) return;
        const dpr = window.devicePixelRatio || 1;
        const w = width;
        const h = box_height;
        canvas.width = Math.round(w * dpr);
        canvas.height = Math.round(h * dpr);
        canvas.style.width = `${w}px`;
        canvas.style.height = `${h}px`;

        const ctx = canvas.getContext("2d");
        ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
        ctx.clearRect(0, 0, w, h);

        const { x0, y0, x1, y1 } = extent_of();
        const k = Math.min((w - 2 * PAD) / (x1 - x0), (h - 2 * PAD) / (y1 - y0));
        // Centre the drawing in whichever direction has slack, and flip y so the
        // data's up is up.
        const ox = (w - (x1 - x0) * k) / 2 - x0 * k;
        const oy = (h - (y1 - y0) * k) / 2 + y1 * k;
        const X = (x) => x * k + ox;
        const Y = (y) => -y * k + oy;

        const c_grid = css("--demo-grid", "#c8c8d0");
        const c_source = css("--demo-source", "#8b8b96");
        const c_target = css("--demo-target", "#2f6fed");
        const c_link = css("--demo-link", "#b9c6e8");
        const c_hit = css("--demo-highlight", "#ff3fa8");

        if (show_grid && grid) draw_grid(ctx, grid, X, Y, k, c_grid);

        if (show_links && show_source) {
            ctx.strokeStyle = c_link;
            ctx.lineWidth = 1;
            ctx.beginPath();
            const n = Math.min(source.length, positions.length);
            for (let i = 0; i < n; ++i) {
                ctx.moveTo(X(source[i][0]), Y(source[i][1]));
                ctx.lineTo(X(positions[i][0]), Y(positions[i][1]));
            }
            ctx.stroke();
        }

        // 2.5px is the floor rather than the size: a dot has to stay legible
        // against the grid lines and the link threads at any point count, and
        // below that it reads as noise. It still grows with the pitch when the
        // cells are large enough to carry it.
        const r = Math.max(2.5, Math.min(4, k * 0.22));

        if (show_source) {
            ctx.strokeStyle = c_source;
            ctx.lineWidth = 1;
            ctx.beginPath();
            for (const [x, y] of source) {
                ctx.moveTo(X(x) + r, Y(y));
                ctx.arc(X(x), Y(y), r, 0, Math.PI * 2);
            }
            ctx.stroke();
        }

        if (colours) {
            // One path per point rather than one batched path: the whole point
            // of the encoding is that every mark carries its own hue.
            for (let i = 0; i < positions.length; ++i) {
                ctx.fillStyle = colours[i];
                ctx.beginPath();
                ctx.arc(X(positions[i][0]), Y(positions[i][1]), r, 0, Math.PI * 2);
                ctx.fill();
            }
        } else {
            ctx.fillStyle = c_target;
            ctx.beginPath();
            for (const [x, y] of positions) {
                ctx.moveTo(X(x) + r, Y(y));
                ctx.arc(X(x), Y(y), r, 0, Math.PI * 2);
            }
            ctx.fill();
        }

        if (highlight >= 0 && highlight < positions.length) {
            ctx.strokeStyle = c_hit;
            ctx.lineWidth = 2;
            ctx.beginPath();
            const [px, py] = positions[highlight];
            ctx.arc(X(px), Y(py), r + 3, 0, Math.PI * 2);
            ctx.stroke();
        }
    }

    /**
     * Grid lines are anchored to a cell boundary and stepped by the pitch.
     * Stepping from the viewport edge instead leaves the lattice floating by an
     * arbitrary fraction of a cell, which makes seated points look mis-seated.
     */
    function draw_grid(ctx, g, X, Y, k, colour) {
        ctx.strokeStyle = colour;
        ctx.lineWidth = 1;

        if (g.type === "rect") {
            // Too dense to read, and it would just fill the panel with ink.
            if (k * g.pitch < 4) return;
            const top = g.y0 + g.rows * g.pitch;
            const right = g.x0 + g.cols * g.pitch;
            ctx.beginPath();
            for (let i = 0; i <= g.cols; ++i) {
                const gx = X(g.x0 + i * g.pitch);
                ctx.moveTo(gx, Y(g.y0));
                ctx.lineTo(gx, Y(top));
            }
            for (let j = 0; j <= g.rows; ++j) {
                const gy = Y(g.y0 + j * g.pitch);
                ctx.moveTo(X(g.x0), gy);
                ctx.lineTo(X(right), gy);
            }
            ctx.stroke();
            return;
        }

        if (g.type === "hex") {
            if (k * g.radius < 3) return;
            ctx.beginPath();
            for (const [cx, cy] of g.centres) {
                for (let v = 0; v < 6; ++v) {
                    const a = g.angle + (v * Math.PI) / 3;
                    const px = X(cx + Math.cos(a) * g.radius);
                    const py = Y(cy + Math.sin(a) * g.radius);
                    if (v === 0) ctx.moveTo(px, py);
                    else ctx.lineTo(px, py);
                }
                ctx.closePath();
            }
            ctx.stroke();
            return;
        }

        if (g.type === "regions") {
            ctx.beginPath();
            for (const c of g.cells) {
                const x = X(c.x);
                const y = Y(c.y + c.height);
                ctx.rect(x, y, c.width * k, c.height * k);
            }
            ctx.stroke();
        }
    }

    $effect(() => {
        // Touch every input so the effect re-runs when any of them changes.
        void source;
        void positions;
        void grid;
        void colours;
        void show_grid;
        void show_links;
        void show_source;
        void highlight;
        void width;
        void box_height;
        draw();
    });

    $effect(() => {
        if (!container) return;
        const observer = new ResizeObserver((entries) => {
            width = Math.max(240, Math.floor(entries[0].contentRect.width));
        });
        observer.observe(container);
        return () => observer.disconnect();
    });
</script>

<div class="plot" bind:this={container}>
    <canvas bind:this={canvas}></canvas>
</div>

<style>
    .plot {
        width: 100%;
        --demo-grid: color-mix(in srgb, var(--sl-color-gray-4) 55%, transparent);
        --demo-source: var(--sl-color-gray-3);
        --demo-target: var(--sl-color-accent);
        --demo-link: color-mix(in srgb, var(--sl-color-accent) 35%, transparent);
        --demo-highlight: var(--demo-accent, #ff3fa8);
    }
    canvas {
        display: block;
        background: var(--sl-color-black);
    }
</style>
