# Third-party notices

hagrid is licensed under the GNU Lesser General Public License v3.0 or later
(see `COPYING.LESSER` and `COPYING`). It incorporates the work below, whose
licence terms are reproduced as those licences require. Including these notices
does not change hagrid's own licence.

---

## Generalized Hilbert ("gilbert") curve — `src/gilbert.ts`

The algorithm implemented in `src/gilbert.ts` is Jakub Červený's, published at
<https://github.com/jakubcerveny/gilbert> under the BSD 2-Clause licence. Our
implementation follows it closely enough to be treated as derived from it.

> Copyright (c) 2018, Jakub Červený
>
> Redistribution and use in source and binary forms, with or without
> modification, are permitted provided that the following conditions are met:
>
> 1. Redistributions of source code must retain the above copyright notice, this
>    list of conditions and the following disclaimer.
>
> 2. Redistributions in binary form must reproduce the above copyright notice,
>    this list of conditions and the following disclaimer in the documentation
>    and/or other materials provided with the distribution.
>
> THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
> AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
> IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
> DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
> FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
> DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
> SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
> CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
> OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
> OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

---

## Algorithms implemented from published papers

These are reimplementations written from the papers. Algorithms are not subject
to copyright, so no licence attaches; they are listed as attribution.

- **`src/gosper.ts`** — Uher, Gajdoš, Snášel, Lai & Radecký, "Hierarchical
  Hexagonal Clustering and Indexing", *Symmetry* 11(6):731, 2019.
  <https://doi.org/10.3390/sym11060731>
- **`src/dgrid.ts`** — Hilasaca, Marcílio-Jr, Eler, Martins & Paulovich, a
  distance-preserving grid layout.
- **`src/nmap.ts`** — Duarte, Sikansi, Fatore, Fadel & Paulovich, "Nmap: A
  Novel Neighborhood Preservation Space-filling Algorithm", IEEE TVCG 2014. The
  `squared` padding follows
  <https://github.com/sebastian-meier/nmap-squared.js>.
- **`src/cmds.ts`** — constrained MDS for overlap removal.

## Runtime dependencies

- **d3-delaunay** — ISC licence. Used by `src/cmds.ts`. Bundled into the browser
  build only; the ESM and CommonJS builds leave it external.
- **delaunator** — ISC licence, via d3-delaunay.

Both are permissive and compatible with the LGPL.
