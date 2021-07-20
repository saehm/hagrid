
function dgrid(G, P, r, s, i=0, j=0) {
    const N = P.length;
    if (N == 0) return
    else if (N == 1) {
        G.push({
            d: P[0],
            i: i,
            j: j,
        });
    } else {
        if (r > s) {
            const r_half = Math.ceil(r / 2);
            const [P1, P2] = split(P, 1, r_half * s); // split by y
            dgrid(G, P1, (r_half), s, i, j);
            dgrid(G, P2, (r - r_half), s, (i + r_half), j);
        } else {
            const s_half = Math.ceil(s / 2);
            const [P1, P2] = split(P, 0, s_half * r); // split by x
            dgrid(G, P1, r, s_half, i, j);
            dgrid(G, P2, r, (s - s_half), i, (j + s_half));
        }
    }
}

function split(P, dim, pos) {
    const sorted = P.sort((a, b) => a[dim] - b[dim]);
    return [sorted.slice(0, pos), sorted.slice(pos)];
}

export function gridify_dgrid(D, parameters) {
    const N = D.length;
    let rows;
    let cols;
    if ("rows" in parameters && "cols" in parameters) {
        rows = parameters.rows;
        cols = parameters.cols;
    } else if (parameters == "square") {
        rows = Math.ceil(Math.sqrt(N));
        cols = Math.ceil(Math.sqrt(N));
    } else if ("aspect_ratio" in parameters) {
        rows = Math.floor(Math.sqrt(N * parameters.aspect_ratio));
        cols = Math.ceil(N / rows)
    } else {
        throw "wrong parameters!"
    }
    let G = [];
    let added_index = D.map(([x, y], i) => [x, y, i]);
    dgrid(G, added_index, rows, cols);
    return G.sort((a, b) => a.d[2] - b.d[2]).map(g => [g.j, g.i])
}
