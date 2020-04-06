import {Delaunay} from "d3-delaunay";

function gamma([px, py], G) {
    const lx = Math.abs(G.x1 - px); // left border
    const rx = Math.abs(G.x2 - px); // right border
    const ly = Math.abs(G.y1 - py); // top border
    const ry = Math.abs(G.y2 - py); // bottom border
    let min_value = Infinity;
    let min_index = null;
    [lx, rx, ly, ry].forEach((d, i) => {
        if (d < min_value) {
            min_value = d;
            min_index = i;
        }
    })
    if (min_index == 0) {
        return [G.x1, py];
    } else if (min_index == 1) {
        return [G.x2, py];
    } else if (min_index == 2) {
        return [px, G.y1];
    } else if (min_index == 3) {
        return [px, G.y2];
    }
}

function line_intersects_rectangle(line, [rx, ry, rw, rh]) {
    const l = line_crosses_line(line, [[rx, ry], [rx, ry+rh]]);
    const r = line_crosses_line(line, [[rx + rw, ry], [rx + rw, ry + rh]]);
    const t = line_crosses_line(line, [[rx, ry], [rx + rw, ry]]);
    const b = line_crosses_line(line, [[rx, ry + rh], [rx + rw, ry + rh]]);

    const intersections = [l, r, t, b].filter(d => !!d)
    if (intersections.length != 2) {
        return 0;
    } else if (intersections.length == 1) {
        const [px, py] = line[1];
        const [ix, iy] = intersections[0];
        return Math.sqrt(Math.pow(px - ix, 2) + Math.pow(py - iy, 2));
    } else {
        //console.log(intersections)
        const [[i0x, i0y], [i1x, i1y]] = intersections;
        return Math.sqrt(Math.pow(i1x - i0x, 2) + Math.pow(i1y - i0y, 2));
    }
}

function line_crosses_line([[ax0, ay0], [ax1, ay1]], [[bx0, by0], [bx1, by1]]) {
    const uA = ((bx1-bx0)*(ay0-by0) - (by1-by0)*(ax0-bx0)) / ((by1-by0)*(ax1-ax0) - (bx1-bx0)*(ay1-ay0));
    const uB = ((ax1-ax0)*(ay0-by0) - (ay1-ay0)*(ax0-bx0)) / ((by1-by0)*(ax1-ax0) - (bx1-bx0)*(ay1-ay0));
    if (uA >= 0 && uA <= 1 && uB >= 0 && uB <= 1) {
        return [ax1 + (uA * (ax1-ax0)), ay1 + (uA * (ay1-ay0))];
    } else {
        return null;
    }
}

export function gridify_cmds(D, Gamma, size, alpha=1) {
    const N = D.length;
    const delaunay = Delaunay.from(D);
    const proximity_graph = D.map((d, i) => {
        return [...delaunay.neighbors(i)].map((j) => {
            const neighbor = D[j];
            return {
                "source": i,
                "target": j,
                "distance": Math.sqrt(Math.pow(d[0] - neighbor[0], 2) + Math.pow(d[1] - neighbor[1], 2))
            }
        })
    }).flat();

    const [kw, kh] = size;

    let deltas = new Array(N).fill(0)
    deltas = deltas.map(() => new Array(N).fill(0));

    // copy array
    let G = D.map(([x, y]) => [x, y, 0, 0])

    proximity_graph.forEach(({source: i, target: j, distance}) => {
        const l = distance//Math.sqrt(Math.pow(pjx - pix, 2) + Math.pow(pjy - piy, 2));
        let d = 0;
        for (let k = 0; k < N; ++k) {
            if (k == i || k == j) continue;
            const kx = G[k][0] - kw / 2;
            const ky = G[k][1] - kh / 2;
            d += line_intersects_rectangle([G[i], G[j]], [kx, ky, kw, kh])
        }
        deltas[i][j] = d + l;
    }) 

    let nominator = new Array(N).fill(0)
    nominator = nominator.map(() => [0, 0]);
    let denominator = new Array(N).fill(0);

    proximity_graph.forEach(({source: i, target: j, distance}) => {
        let pi = G[i];
        let pj = G[j];
        let dij = deltas[i][j];
        const wij = (1/(dij**2));
        nominator[i][0] += wij * (pj[0] + dij * ((pi[0] - pj[0]) / distance))
        nominator[i][1] += wij * (pj[1] + dij * ((pi[1] - pj[1]) / distance))
        denominator[i] += wij
    })

    G = G.map(([px, py], i) => {
        const gamma_i = gamma([px, py], Gamma);
        return [
            nominator[i][0] + alpha * gamma_i[0] / denominator[i],
            nominator[i][1] + alpha * gamma_i[1] / denominator[i],
        ]
    })


    return G
}