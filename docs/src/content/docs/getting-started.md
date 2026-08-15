---
title: Getting started
description: Install hagrid and remove the overlaps from a scatterplot.
---

## Install

```sh
npm install @saehrimnir/hagrid
```

```js
import * as hagrid from "@saehrimnir/hagrid";
```

The package ships ESM, CommonJS and a browser build, with TypeScript types included. For a
script tag, use the UMD bundle and the `hagrid` global:

```html
<script src="https://cdn.jsdelivr.net/npm/@saehrimnir/hagrid"></script>
```

`d3-delaunay` is a runtime dependency and is **not** bundled into the ESM and CommonJS
builds, so your own module graph resolves it. The browser build does bundle it, since a
`<script>` tag has no resolver.

## Remove the overlaps

```js
const { positions, runtime } = hagrid.gridify(data, "hilbert");
```

`data` is an array of `[x, y]` pairs. `positions` holds one grid position per input point,
**in input order**, so `positions[i]` is where `data[i]` ended up. `runtime` is the layout
time in milliseconds.

Or call a method directly, which skips the parameter defaulting `gridify` does for you:

```js
const positions = hagrid.gridify_hilbert(data, { level });
const positions = hagrid.gridify_gilbert(data);
const positions = hagrid.gridify_gosper(data, { level });
const positions = hagrid.gridify_dgrid(data, { aspect_ratio: 1 });
const positions = hagrid.gridify_nmap(data);
```

## Which space do I get back?

This trips people up, so it is worth being explicit. The methods do not all answer in the
same coordinate space:

| method | returns |
| ---- | ---- |
| Gosper, NMap, CMDS | positions in the input's own coordinate space |
| Hilbert, Gilbert, DGrid | integer cell indices on a `cols × rows` lattice |
| GridFit | grid units — and it must be *given* grid units |

So drawing the original points and the layout on the same axes needs care. Normalising each
by its own bounding box is the tempting move and it is wrong: the two bounding boxes differ,
so it invents displacement that is not in the layout. Project the originals into the
method's output space instead, and use one transform for both.

For the lattice methods that projection is the same affine the method applies internally —
the data's bounding box stretched onto the lattice, so the point at the left edge lands in
column $0$ and the one at the right edge in column $\text{cols} - 1$:

$$
x' = (x - x_{\min}) \cdot \frac{\text{cols} - 1}{w}, \qquad
y' = (y - y_{\min}) \cdot \frac{\text{rows} - 1}{h}
$$

```js
const bb = hagrid.utils.get_bounds(data);
const sx = bb.width > 0 ? (cols - 1) / bb.width : 0;
const sy = bb.height > 0 ? (rows - 1) / bb.height : 0;
const source = data.map(([x, y]) => [(x - bb.x) * sx, (y - bb.y) * sy]);
```

## Methods and parameters

Hilbert and Gosper are ours; DGrid, CMDS, NMap and GridFit are reimplementations
or translations of other published techniques.

| name | description | parameters |
| ---- | ---- | ---- |
| [Hilbert](/techniques/hilbert/) | Hilbert curve, rectangular cells on a square grid | `pluslevel`, `whitespace`, `keep_aspect_ratio` |
| [Gilbert](/techniques/gilbert/) | generalized Hilbert curve on a grid matching the data's aspect ratio | `whitespace`, or explicit `cols` and `rows` |
| [Gosper](/techniques/gosper/) | Gosper curve, hexagonal cells | `pluslevel`, `whitespace`, `orientation` |
| [DGrid](/techniques/dgrid/) | fills the whole area, keeping a given cell aspect ratio | `aspect_ratio`, or explicit `rows` and `cols` |
| [CMDS](/techniques/cmds/) | constrained MDS; iterative | `alpha`, `Gamma`, `size`, `max_iter`, `max_rescales`, `inner_iterations` |
| [NMap](/techniques/nmap/) | space filling treemap | `BB`, `squared` |
| [GridFit](/techniques/gridfit/) | hierarchical partitioning, keeps points near their cell | none |

## CMDS is iterative

CMDS optimises over many steps rather than in one pass, so it returns a generator of
intermediate layouts instead of a final position array:

```js
const { steps } = hagrid.gridify(data, "cmds");
let positions, pitch;
for (const [layout, proximity_graph, size] of steps) {
    positions = layout;
    pitch = size;
}
```

The generator is synchronous, so pulling every step runs the whole optimisation in one
blocking loop. To animate without freezing the page, drive it a step at a time:

```js
const step = () => {
    const { value, done } = steps.next();
    if (done) return;
    draw(value[0], value[2]);
    requestAnimationFrame(step);
};
step();
```

## Licence

hagrid is licensed under the **GNU Lesser General Public License v3.0 or later**. Third-party
notices, including the BSD 2-Clause notice for the Gilbert algorithm, ship with the package
in `THIRD-PARTY.md`.
