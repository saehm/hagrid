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
const positions = hagrid.gridify_gilbert(data);
const positions = hagrid.gridify_gosper(data, { level });
const positions = hagrid.gridify_dgrid(data, { aspect_ratio: 1 });
const positions = hagrid.gridify_nmap(data);
```

### Methods & Parameters

The methods `hilbert` and `gosper` are ours; `dgrid`, `cmds`, and `nmap` are reimplementations or translations of other techniques.

| name | description | parameters |
| ---- | ---- | ---- |
| `hilbert` | uses the Hilbert curve, produces rectangular grid cells on a square grid | `pluslevel` (default `0`), `whitespace` (default `1`), `keep_aspect_ratio` |
| `gilbert` | generalized Hilbert curve on a grid matching the data's aspect ratio | `whitespace` (default `1`), or explicit `width` and `height` |
| `gosper` | uses the Gosper curve, produces hexagonal grid cells | `pluslevel` (default `0`), `orientation` (default `"flat"`) |
| `dgrid` | fills the whole visualization area, keeping a given aspect ratio of the grid cells | `aspect_ratio` (default `1`), or explicit `rows` and `cols` |
| `cmds` | similar to MDS but with constraints | `alpha`, `Gamma`, `size`, `max_iter`, `max_rescales` |
| `nmap` | a space filling treemap approach | `BB`, `squared` |
| `gridfit` | hierarchical partitioning that keeps points near their original cell | none; input is in grid-cell units |

### `cmds` is iterative

`cmds` optimises over many steps rather than in one pass, so it returns a generator of intermediate
layouts instead of a final position array. Each step also carries the node size in force, which
changes as the grid grows. The last yielded layout is the grid-snapped result:

```js
const { steps } = hagrid.gridify(data, "cmds");
let positions, pitch;
for await (const [layout, proximity_graph, size] of steps) {
    positions = layout;   // draw here to animate the optimisation
    pitch = size;
}
```

`size` is a *starting* node size rather than a fixed one — see below.

## Migrating from 0.1.x

`0.3.0` is the first release since the library was ported to TypeScript. The layouts produced by `hilbert`, `dgrid` and `nmap` are unchanged — they are pinned by fixtures captured from `0.1.1` — but the surrounding API changed:

- **`gridify` returns `{ positions, runtime }`** instead of an array with a `runtime` property bolted on. Replace `const result = gridify(...)` with `const { positions } = gridify(...)`.
- **`gridify(data, "cmds", …)` returns `{ steps }`**, an async generator. Previously it returned a generator that could not complete: the final proximity graph was constructed with the wrong arity, so the grid-snapped layout always threw before reaching the caller.
- **`gridify(data, "gosper", …)` works without an explicit `pluslevel`.** It previously read an absent `pluslevel` and computed a `NaN` curve level.
- **There is now one Gosper implementation.** `gosper.js` and `gosper2.js` both exported `gosper_encode`, `gosper_decode` and `gridify_gosper`, and which one you got depended on whether you called the export directly or went through `gridify`. The curve-building one has been removed in favour of the direct indexing of Uher et al., which keeps the plain names — `gosper_encode`, `gosper_decode`, `gosper_size`, `gosper_curve`, `gridify_gosper`. Layouts differ from `0.1.x`. `scale_factor` is gone (the fit is isotropic now); `orientation` replaces it. `gosper_encode` returns `-1` outside the island, and `CurveIndex` is gone with the scan it existed to accelerate.
- **`gridify_gridfit` returns `{ nodes, final, leftovers, error }`** and actually runs. It previously referenced an unimported `d3` and threw a `ReferenceError` immediately.
- **`d3-delaunay` is no longer bundled** into the ESM and CommonJS builds; it is a peer of your own module graph, as its `dependencies` entry always implied. The browser build still bundles it.
- The package now has an `exports` map. Deep imports into `dist/` are no longer supported.

## Licence

hagrid is licensed under the **GNU Lesser General Public License v3.0 or later**; see `COPYING.LESSER`
and `COPYING`. Third-party notices, including the BSD 2-Clause notice for the gilbert algorithm, are
in `THIRD-PARTY.md`.

## Development

```sh
pnpm install
pnpm test        # golden + behavioural tests
pnpm typecheck
pnpm build
```

The fixtures in `test/__goldens__/` were captured from the `0.1.1` bundle. A test failing there means the port drifted; re-record only for a deliberate, explained behaviour change.

## Performance

Both curve implementations were rebuilt around not allocating. `rotate()` used to return a fresh
`[x, y]` array on every level of every call, which dominated everything else:

| | before | after |
| ---- | ---- | ---- |
| `hilbert_encode` | 0.12 M/s | **10.1 M/s** |
| `hilbert_decode` | 1.88 M/s | 2.02 M/s |
| `gosper_encode` | — | 6.6 M/s |
| `gosper_decode` | — | 8.8 M/s |

End to end, `gosper` no longer builds the curve at all:

| N | before (curve + index) | after (direct) |
| ---- | ---- | ---- |
| 5 000 | 80 ms | 7 ms |
| 20 000 | 510 ms | 14 ms |
| 50 000 | 1 047 ms | 33 ms |
| 200 000 | 21 197 ms | **192 ms** |

`hilbert_encode`/`hilbert_decode` also no longer corrupt silently above level 15. The index was
passed through `>>` and `&`, which coerce to 32 bits, so a level-16 grid returned wrong
coordinates with no error. The index is now kept in floating point — exact to 2^53, i.e. level 26
— and anything beyond that, or a non-power-of-two size, is rejected outright.

## The Gosper implementation

There used to be two, and they collided (see *Migrating* above). The survivor is the direct
indexing of Uher et al., which computes a curve index arithmetically in `O(level)` and never
materializes anything. The other built the whole curve as an L-system and matched points against
it; it was dropped because it lost on every axis measured — 110x slower end to end, `O(7^level)`
memory, it threw on collinear input, and it was not even faster at producing the full curve for
drawing (17.7 ms vs 10.7 ms at level 6).

The one thing it did differently was semantics for points outside the island: it snapped them to
the nearest cell, where `gosper_encode` returns `-1`. Inside the island the two rules agree
exactly — verified on thousands of sampled points at levels 2 to 4 — so this only matters when
indexing arbitrary points, not when gridding a dataset.

### Getting it working

The direct implementation never worked. `decode` was always correct, but only 1 index in `7^level`
survived a `decode → encode` round trip. Two independent faults in the localization step (paper
Algorithm 1), where fixing either alone changes nothing — the round trip stays at 14% or 0%:

1. **`xy2cube` was not the inverse of `cube2xy`.** It used a *pointy-topped* hexagon convention
   while `cube2xy`, and therefore all of `decode`, is *flat-topped*. The two disagreed about the
   lattice itself.
2. **`B()` inverted the sign test on all three axes** relative to paper eq. (16), which specifies
   `x` dominant and `x < 0 → 2`, `x > 0 → 5`, likewise `y` and `z`. Its threshold was also
   `10e-3` where the paper gives `1e-5`.

Structurally, the loop also walked xy space applying cumulative `1/√7^i` scalings and `-α·i`
rotations, where the paper applies one constant matrix `T = [[2/7, -1/7], [1/7, 3/7]]` (eq. 13) to
*axial* coordinates, repeatedly. The file already contained a commented-out `transform()` that was
exactly this matrix — the correct approach had been written, then abandoned for the broken one.

Algorithm 2 (index conversion) was already correct, including its `E(P2, ·) = [4,0,1,2,3,6,5]`
transfer table, which matches the paper's Table 1. `encode` now inverts `decode` exactly at every
level tested — all 16 807 indices at level 5.

### Points outside the island

The climb records only each hexagon's position *relative to its parent* and never checks which
root island it ended in, so an outside point used to fold onto an arbitrary interior cell, very
often one at the opposite end of the curve. A contained point always finishes on the root hexagon
`(0, 0)`, so `gosper_encode` returns `-1` for anything else.

That exposed a larger bug: **the island is centred on the origin, but `gridify_gosper` never moved
the data there.** `gosper_size` scales the island so the data's bounding circle fits inside its
inscribed circle, which only holds about a shared centre. For ordinary data like `[0,100]²`,
centred at `(50,50)`, most points sat off the tile — 26 112 of 48 000 measured across square,
wide, degenerate-line and clustered inputs. It now centres the data and shifts the result back,
and the same measurement gives zero.

The discarded `(q, l)` is in fact the *island address*, so the scheme extends over the whole
plane: islands tile by translation (forced, since `T` is linear), and

```
global_position(island, index) = decode(index) + rotate(cube2xy(island, size·√7^level), −level·α)
```

reconstructs any cell — verified across ~50 islands out to six island widths. That is an unbounded
addressing scheme rather than a single global ordering, and it is not currently exposed as API.

## Hexagon orientation

`"pointy"` puts neighbours at multiples of 60°, `"flat"` at 30° off those. They are the same
lattice turned by 30° — the basis sends axial `(1,0)` to 0° and 30° respectively, at equal length:

```js
hagrid.gosper_curve(level, "pointy");
hagrid.gridify_gosper(data, { level, orientation: "pointy" });
hagrid.gosper_encode(p, level, size, "pointy");
```

The default is `"flat"`, the paper's own lattice. Because the rotation is rigid, cell ordering and
step length are unchanged and indexing round-trips exactly in either orientation. Note that
switching turns the lattice *underneath* the data rather than carrying each point to a matching
cell, so points are reassigned — that is usually the point of the option.

## N-Map tiles the plot area

`nmap` partitions the plotting area, so its rectangles should cover the bounding box exactly. They
did not: `sizes()` gave the first half's share to *both* halves, and since the left half is pinned
to the region's left edge and the right half to its right edge, shares that do not add up to the
full size leave a gap down the middle.

The outer extent still looked right, which is why it went unnoticed — the interior simply did not
tile:

| N | covered before | covered after |
| ---- | ---- | ---- |
| 8 | 100% | 100% |
| 50 | 62.6% | **100%** |
| 101 | 61.6% | **100%** |
| 500 | 88.9% | **100%** |
| 1 001 | 89.0% | **100%** |

N = 8 was always correct because an even split has `pA === pB`, so the wrong share happened to be
the right one — the bug only ever fired on odd counts.

The rectangles are now available directly, which is what you want for drawing a treemap:

```js
const cells = hagrid.nmap_regions(data);   // [{ x, y, width, height, index }, ...]
```

### nmap needs spread on both axes

`nmap` cuts a region at the midpoint between the two halves of its points. If every point in that
region shares the split coordinate, the boundary lands on the region's own edge, one side gets zero
extent, and the affine scale becomes `share / 0`. That used to propagate silently, so the layout
came back full of `null`s with nothing to indicate why. It now throws with the offending axis named.

**Distinct points are not enough.** All of these are 100% unique and all of them fail:

| input | n | outcome |
| ---- | ---- | ---- |
| four points on a vertical line | 4 | throws |
| 10×10 integer lattice | 100 | throws |
| coordinates rounded to integers | 194 | throws |
| points along a diagonal | 16 | fine |
| continuous random | 200 | fine |

The property that matters is spread on *both* axes within every subdivision, which is why the
diagonal survives (every `x` and every `y` distinct) and the lattice does not. Quantised data —
counts, scores, rounded measurements — commonly trips it. Jitter the input slightly, or use
`gilbert` or `hilbert`, which have no such requirement.

A real fix would fall back to a proportional split when a region is degenerate, placing tied points
in adjacent cells. That is not done yet.

## CMDS convergence

`cmds` never terminated: it burned every iteration on every input and then threw. Beyond the arity
bug (see *Migrating*), the update rule itself was already a faithful implementation of the paper's
eq. (4), with its recommended `N = 20` inner iterations and `α = 0.1`. What was missing was
everything §4.2 says about making that rule converge:

- the **`m × n` grid** that determines the node size, rather than an ad-hoc `w / (√N + 1)`;
- **`Γ` inset by half a node**, so a node whose centre is inside is drawn completely;
- the input **uniformly scaled** into `Γ` — `gridify` sets `Γ` to the data's own bounding square, so
  without this the extreme points start outside and can never get in;
- **growing the grid** by a row or column and rescaling the nodes when a run does not converge.

The termination test was also wrong: it summed the overlap over every edge and compared it to a
constant `1`, which across ~570 edges demands an average overlap below `0.002`. The paper asks for
"no overlap", which is per-pair.

Without the rescale loop the layout oscillates forever. The overlap correction enters the ideal
distance as `d_ij = l_ij + δ_ij`, so once a pair is separated `δ_ij` drops back to zero, the spring
pulls them together, and the overlap returns — measured overlap alternated exactly `0 → 5.19 → 0`
with period two.

**This is improved rather than solved.** `cmds` now terminates and converges on most inputs, where
before it converged on none. But it typically has to shrink the nodes to around 1% of the plot area,
where the paper reports converged layouts at roughly 70%. Something in the implementation still
differs from the paper, and `cmds` remains the least trustworthy method here — prefer `hilbert`,
`gilbert` or `gosper`. When a run does not converge, the final grid snap can place two nodes in the
same cell.

## Elongated data

`hilbert` needs a square, power-of-two grid, so a wide dataset is either stretched out of shape or
padded into a mostly-empty square — for a 1000x10 extent, roughly 1% of the cells are usable, so
holding the same number of points needs ~100x more of them. `gosper` fits its island isotropically
and fares no better: a 1000x10 dataset comes back at aspect 2.8 rather than 100.

`gilbert` is a generalized Hilbert curve that fills an arbitrary `width x height` rectangle, so the
grid can simply match the data:

```js
const { positions } = hagrid.gridify(data, "gilbert");
const { positions } = hagrid.gridify(data, "gilbert", { width: 200, height: 5 });
```

| N, at 1000x10 aspect | `gilbert` | `hilbert` (padded square) |
| ---- | ---- | ---- |
| 1 000 | 1 000 cells (250x4) | 262 144 |
| 10 000 | 10 000 cells (1000x10) | 1 048 576 |
| 100 000 | 100 000 cells (3125x32) | 16 777 216 |

It costs about 6x more per query than `hilbert` — a recursive descent rather than bit twiddling —
which two orders of magnitude fewer cells more than repays. It also has no power-of-two constraint
and no level ceiling.

One wrinkle: the curve is discontinuous on 506 of the 2304 shapes up to 48x48 — consecutive indices
land two cells apart instead of one. The mapping is still a bijection, so every point gets its own
cell and nothing visibly breaks, but the point of a curve is that neighbours along it are neighbours
in the plane, and collision resolution walks the index assuming that. Every failing shape has an
even short side and an odd long side, and nothing outside that family fails, so `gilbert_dimensions`
steps off it by growing the long side by one — a single extra row or column. Passing `width` and
`height` explicitly can still land on a discontinuous shape.

## Gridfit

`gridfit` (Keim & Herrmann, IEEE Visualization '98) hierarchically partitions the data space,
shifting each dividing line by the smallest integer amount that leaves every subregion at least as
many grid cells as it has points — the invariant `P(B) <= A(B)` of §3.2.

It could never run: it referenced an unimported `d3` and threw a `ReferenceError` on the first line
of work. Once that was fixed it still dropped points, for two reasons:

- **The extent convention was inconsistent.** `area()` counted `x_max - x_min` cells while the point
  filter accepted `x_min <= x <= x_max`, one column more, so every region was allocated one fewer
  cell than it took points for. Harmless at the root, fatal a few levels down where regions are two
  or three cells wide.
- **Partitioning and divider placement were separate passes.** The whole quadtree was built against
  unshifted extents, and only then walked shifting dividers — so every node's children had been
  positioned against an extent their parent no longer had. The paper interleaves the two, level by
  level, and it is now a single recursive pass.

It also now searches both directions for the smallest shift, as `min{x ∈ ℤ | ...}` requires, and
implements the second pass of §3.3 that seats anything a straight-line partition could not.

Every point is placed in its own cell, and the mean displacement from its original position stays
under one cell — which is optimization goal (1) of §2.1:

| N | time | mean displacement |
| ---- | ---- | ---- |
| 1 000 | 2 ms | 0.73 cells |
| 5 000 | 5 ms | 0.73 cells |
| 20 000 | 19 ms | 0.75 cells |

`gridfit` takes coordinates already in grid-cell units, and is the most robust method here on awkward
data: identical points, single rows or columns, and heavy duplicates all place cleanly, where `nmap`
rejects them.

## Known issues

- `hilbert` can leave a few points sharing a cell when the input has exact duplicates: its collision
  walk resolves most but not all of them. On a 200-point set with 195 distinct positions it returns
  196 distinct cells. Pre-existing, and byte-identical to `0.1.1`. `gridfit` and `gilbert` handle the
  same input cleanly.

- `gosper` does not preserve the data's aspect ratio for elongated inputs; use `gilbert` for those.
- `nmap` cannot lay out data where a subdivision has no spread on the axis it is cut along, and now
  throws rather than returning `NaN`. See below.
