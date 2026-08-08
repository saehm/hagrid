# HAGRID
<a href="#"><img src="img/hagrid.png" width=120 align="left" hspace="10" vspace="4"></a>
Hagrid is a JavaScript library for removing overlaps in scatterplots. When using glyphs or images instead of points in a scatterplot, overlaps happen which impedes the readability of the scatterplot. Hagrid uses space filling curves to remove the overlaps in such a scatterplot.

## Installation

```sh
npm install @saehrimnir/hagrid
```

```js
import * as hagrid from "@saehrimnir/hagrid";
```

The package ships ESM, CommonJS and a browser build, with TypeScript types included. For a script tag, use the UMD bundle and the `hagrid` global:

```html
<script src="https://cdn.jsdelivr.net/npm/@saehrimnir/hagrid"></script>
```

## Usage

Hagrid provides a method `gridify`:

```js
const { positions, runtime } = hagrid.gridify(data, method, parameters);
```

`data` is an array of `[x, y]` pairs. `positions` holds one grid position per input point, in input order, and `runtime` is the layout time in milliseconds.

Or use a method directly, which skips the parameter defaulting `gridify` does for you:

```js
const positions = hagrid.gridify_hilbert(data, { level });
const positions = hagrid.gridify_gosper(data, { level });
const positions = hagrid.gridify_dgrid(data, { aspect_ratio: 1 });
const positions = hagrid.gridify_nmap(data);
```

### Methods & Parameters

The methods `hilbert` and `gosper` are ours; `dgrid`, `cmds`, and `nmap` are reimplementations or translations of other techniques.

| name | description | parameters |
| ---- | ---- | ---- |
| `hilbert` | uses the Hilbert curve, produces rectangular grid cells | `pluslevel` (default `0`), `whitespace` (default `1`), `keep_aspect_ratio` |
| `gosper` | uses the Gosper curve, produces hexagonal grid cells | `pluslevel` (default `0`), `scale_factor` (default `0.8`) |
| `dgrid` | fills the whole visualization area, keeping a given aspect ratio of the grid cells | `aspect_ratio` (default `1`), or explicit `rows` and `cols` |
| `cmds` | similar to MDS but with constraints | `alpha`, `Gamma`, `size`, `max_iter` |
| `nmap` | a space filling treemap approach | `BB`, `squared` |

### `cmds` is iterative

`cmds` optimises over many steps rather than in one pass, so it returns a generator of intermediate layouts instead of a final position array. The last yielded layout is the grid-snapped result:

```js
const { steps } = hagrid.gridify(data, "cmds", { max_iter: 50 });
let positions;
for await (const [layout, proximity_graph] of steps) {
    positions = layout; // draw here to animate the optimisation
}
```

## Migrating from 0.1.x

`0.2.0` is the first release since the library was ported to TypeScript. The layouts produced by `hilbert`, `dgrid` and `nmap` are unchanged — they are pinned by fixtures captured from `0.1.1` — but the surrounding API changed:

- **`gridify` returns `{ positions, runtime }`** instead of an array with a `runtime` property bolted on. Replace `const result = gridify(...)` with `const { positions } = gridify(...)`.
- **`gridify(data, "cmds", …)` returns `{ steps }`**, an async generator. Previously it returned a generator that could not complete: the final proximity graph was constructed with the wrong arity, so the grid-snapped layout always threw before reaching the caller.
- **`gridify(data, "gosper", …)` works without an explicit `pluslevel`.** It previously read an absent `pluslevel` and computed a `NaN` curve level.
- **The two Gosper implementations no longer collide.** `gosper.js` and `gosper2.js` both exported `gosper_encode`, `gosper_decode` and `gridify_gosper`, and which one you got depended on whether you called the export directly or went through `gridify`. The direct-indexing variant of Uher et al. is now exported under its own names — `gridify_gosper_uher`, `gosper_uher_encode`, `gosper_uher_decode`, `gosper_uher_size` — and is **known to produce incorrect coordinates**; prefer `gridify_gosper`.
- **`gridify_gridfit` returns `{ nodes, final, leftovers, error }`** and actually runs. It previously referenced an unimported `d3` and threw a `ReferenceError` immediately.
- **`d3-delaunay` is no longer bundled** into the ESM and CommonJS builds; it is a peer of your own module graph, as its `dependencies` entry always implied. The browser build still bundles it.
- The package now has an `exports` map. Deep imports into `dist/` are no longer supported.

## Development

```sh
pnpm install
pnpm test        # golden + behavioural tests
pnpm typecheck
pnpm build
```

The fixtures in `test/__goldens__/` were captured from the `0.1.1` bundle. A test failing there means the port drifted; re-record only for a deliberate, explained behaviour change.

## Known issues

- `gridify_gosper_uher` does not yet reproduce `gridify_gosper`'s coordinates.
- `nmap`'s internal `sizes()` splits a region using the first half's point count for both halves, so subdivisions are not proportional. Preserved for now because fixing it changes every `nmap` output.
- `gridfit` runs but has not been validated against its paper.
