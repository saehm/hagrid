## Hagrid — Removes overlaps in scatterplots using space filling curves.

<a href="#"><img src="https://raw.githubusercontent.com/saehm/hagrid/refs/heads/master/hagrid.svg" width=80 align="left" hspace="10" vspace="6"></a>

Draw a scatterplot with glyphs or images instead of points and they overlap, which is where the plot stops being readable.
Hagrid moves each point to its own cell on a grid, keeping it as close to where it started as the grid allows — seven techniques behind one call, each documented with what it refuses and what it gets wrong.

<br>

[![CI](https://img.shields.io/github/actions/workflow/status/saehm/hagrid/ci.yml?branch=master&label=tests)](https://github.com/saehm/hagrid/actions/workflows/ci.yml)
[![Docs](https://img.shields.io/github/actions/workflow/status/saehm/hagrid/docs.yml?branch=master&label=docs)](https://saehm.github.io/hagrid/)
[![Coverage](https://img.shields.io/codecov/c/github/saehm/hagrid?label=coverage)](https://app.codecov.io/gh/saehm/hagrid)
![NPM Downloads](https://img.shields.io/npm/dw/%40saehrimnir%2Fhagrid)
[![License](https://img.shields.io/github/license/saehm/hagrid)](https://raw.githubusercontent.com/saehm/hagrid/refs/heads/master/LICENCE)

### Installation

If you use npm, install with `npm install @saehrimnir/hagrid`, and use it with

```js
import * as hagrid from "@saehrimnir/hagrid";
```

Otherwise download the files [here](https://github.com/saehm/hagrid/releases/latest), or use for instance [jsDelivr](https://cdn.jsdelivr.net/npm/@saehrimnir/hagrid) this way:

```html
<script src="https://cdn.jsdelivr.net/npm/@saehrimnir/hagrid"></script>
```

The package ships ESM, CommonJS and a browser build, with TypeScript types included. `d3-delaunay` is a runtime dependency of `cmds` and is left external to the ESM and CommonJS builds; the browser build bundles it.

### Quick start

```js
import * as hagrid from "@saehrimnir/hagrid";

const { positions, runtime } = hagrid.gridify(data, "hilbert");
```

`data` is an array of `[x, y]` pairs. `positions` holds one grid position per input point, **in input order**, so `positions[i]` is where `data[i]` ended up, and `runtime` is the layout time in milliseconds.

Or call a method directly, which skips the parameter defaulting `gridify` does for you:

```js
const positions = hagrid.gridify_hilbert(data, { level });
const positions = hagrid.gridify_gilbert(data);
const positions = hagrid.gridify_gosper(data, { level });
const positions = hagrid.gridify_dgrid(data, { aspect_ratio: 1 });
const positions = hagrid.gridify_nmap(data);
```

The methods do not all answer in the same coordinate space, which trips people up: `gosper`, `nmap` and `cmds` return positions in the input's own space, `hilbert`, `gilbert` and `dgrid` return integer cell indices, and `gridfit` takes *and* returns grid units. See [Which space do I get back?](https://saehm.github.io/hagrid/getting-started/#which-space-do-i-get-back).

### Techniques

**Space filling curves** —
[Hilbert](https://saehm.github.io/hagrid/techniques/hilbert/) ·
[Gilbert](https://saehm.github.io/hagrid/techniques/gilbert/) ·
[Gosper](https://saehm.github.io/hagrid/techniques/gosper/)

**Reimplementations of published methods** —
[DGrid](https://saehm.github.io/hagrid/techniques/dgrid/) ·
[CMDS](https://saehm.github.io/hagrid/techniques/cmds/) ·
[NMap](https://saehm.github.io/hagrid/techniques/nmap/) ·
[GridFit](https://saehm.github.io/hagrid/techniques/gridfit/)

| name | grid | parameters |
| ---- | ---- | ---- |
| `hilbert` | square, `2^level` a side | `pluslevel`, `whitespace`, `keep_aspect_ratio` |
| `gilbert` | arbitrary `cols × rows` | `whitespace`, or explicit `cols` and `rows` |
| `gosper` | hexagonal | `pluslevel`, `whitespace`, `orientation` |
| `dgrid` | dense, fills the area | `aspect_ratio`, or explicit `rows` and `cols` |
| `cmds` | constrained MDS, iterative | `alpha`, `Gamma`, `size`, `max_iter`, `max_rescales`, `inner_iterations` |
| `nmap` | treemap rectangles | `BB`, `squared` |
| `gridfit` | hierarchical partitioning | none; input is already in grid-cell units |

Elongated data is the usual reason to choose between them: `hilbert` must pad to a square, so a 1000×10 extent leaves about 1% of its cells usable, where `gilbert` matches the shape exactly. Round glyphs are the other: `gosper`'s hexagons hold a circle 7.5% wider than a square cell of the same area, which is 15.5% more glyph for free.

### `cmds` is iterative

Alone among the methods, `cmds` optimises over many steps rather than in one pass, so it returns a generator of intermediate layouts. The last yield is the grid-snapped result:

```js
const { steps } = hagrid.gridify(data, "cmds");
for (const [layout, proximity_graph, size] of steps) {
    positions = layout;   // draw here to animate the optimisation
    pitch = size;         // the node size in force, which changes as the grid grows
}
```

The generator is synchronous, so that loop runs the whole optimisation in one blocking pass. Drive it a step at a time to animate without freezing the page — [the technique page](https://saehm.github.io/hagrid/techniques/cmds/) has the pattern.

### Performance

Both curve implementations are built around not allocating, which is most of the difference:

| | before | after |
| ---- | ---- | ---- |
| `hilbert_encode` | 0.12 M/s | **10.1 M/s** |
| `gosper_encode` | — | 6.6 M/s |
| `gosper` end to end, N = 200 000 | 21 197 ms | **192 ms** |

`gosper` no longer builds its curve at all: the direct indexing of Uher et al. computes an index arithmetically in `O(level)` and never materialises anything. Hilbert indices are kept in floating point — exact to 2^53, so level 26 — and a non-power-of-two size or an out-of-range level is rejected rather than silently returning wrong coordinates.

### Resources

- [Documentation](https://saehm.github.io/hagrid/)
- [Getting started](https://saehm.github.io/hagrid/getting-started/) — install, the shape of the API, and which space each method answers in
- [Playground](https://saehm.github.io/hagrid/playground/) — every method over datasets chosen to show where each one struggles
- [Showcase](https://saehm.github.io/hagrid/showcase/glyphs/iris/) — why gridify at all: seven worked examples with glyphs
- [API reference](https://saehm.github.io/hagrid/api/readme/) — generated from the source by TypeDoc

The techniques and their evaluation are described in *Hagrid: using Hilbert and Gosper curves to gridify scatterplots*, Journal of Visualization 25(6), 2022 — [doi:10.1007/s12650-022-00854-7](https://doi.org/10.1007/s12650-022-00854-7).

### Licence

hagrid is licensed under the **GNU Lesser General Public License v3.0 or later**; see `LICENCE`. Third-party notices, including the BSD 2-Clause notice for the gilbert algorithm, are in `THIRD-PARTY.md`.

### Development

```sh
pnpm install
pnpm test        # behavioural tests
pnpm typecheck
pnpm build
pnpm docs:dev    # builds the library, then serves the docs site
```

The tests check what each method promises — a complete tiling, no overlaps, containment, curve continuity — rather than pinning recorded output. Layouts may therefore shift between versions; the guarantees listed on each technique page are what stays fixed.

> **Status:** all seven techniques are implemented and documented. `cmds` is the least trustworthy of them — it rarely satisfies its own end condition, and a run that stops with overlaps left can snap two nodes into one cell; prefer `hilbert`, `gilbert`, `gosper` or `gridfit`. Each technique page lists what its method refuses and what it gets wrong, with a demo you can make fail.
