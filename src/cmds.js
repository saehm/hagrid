import {Delaunay} from "d3-delaunay";

// gamma puts a point outside the boundary back inside,
// 
function gamma([px, py], {x1, x2, y1, y2}) {
    let [nx, ny] = [px, py];
    if (px < x1) nx = x1;
    if (px > x2) nx = x2;
    if (py < y1) ny = y1;
    if (py > y2) ny = y2;
    return [nx, ny]
}

function distance([ax, ay], [bx, by]) {
    return Math.sqrt(Math.pow(ax - bx, 2) + Math.pow(ay - by, 2));
}

function compute_overlap([px, py], [qx, qy], s) {
    const dx = Math.abs(px - qx);
    const dy = Math.abs(py - qy);
    if (dx < s || dy < s) {
        const ox = Math.abs(Math.min(dx - s, 0));
        const oy = Math.abs(Math.min(dy - s, 0));
        const dd = distance([px, py], [qx, qy]);
        if (ox < oy) {
            return (dd * ox) / dx;
        } else if (oy <= ox) {
            return (dd * oy) / dy;
        }
    } else {
        return 0
    }
}

// computes all the stuff needed for cmds
function create_proximity_graph(A, D, s) {
    const delaunay = Delaunay.from(A);
    return A.map((d, i) => {
        return [...delaunay.neighbors(i)].map((j) => {
            const pd = D[i];
            const qd = D[j];
            const pa = A[i];
            const qa = A[j];
            // maybe better to save
            const odistance = distance(pd, qd);
            const overlap = compute_overlap(pa, qa, s);
            return {
                "source": i,
                "target": j,
                "distance": odistance,
                "overlap": overlap,
                "delta": odistance + overlap
            }
        })
    })
    .flat()
}

// end condition is
// no overlaps anymore and all points are inside the boundaries
function end_condition(proximity_graph, G, {x1, x2, y1, y2}) {
    const sum_overlaps = proximity_graph
        .map(d => d.overlap)
        .reduce((a, b) => a + b);

    const inside_Gamma = G
        .map(([px, py]) => px >= x1 && px <= x2 && py >= y1 && py <= y2)
        .reduce((a, b) => a && b);

    return sum_overlaps == 0 && inside_Gamma
}

export async function gridify_cmds(D, parameters) {
    const N = D.length;
    const Gamma = parameters.Gamma;
    const [kw, kh] = parameters.size;
    const alpha = parameters.alpha ? parameters.alpha : .1;
    let max_iter = parameters.max_iter ? parameters.max_iter : 5000;
    
    // construct proximity graph and 
    // calculate the ideal distances according to (3)
    let proximity_graph = create_proximity_graph(D, D, kw);
    // copy D
    let G = D.map(([x, y]) => [x, y]);
    
    do {
        // update step according to (6)
        let nominator = new Array(N).fill(0)
        nominator = nominator.map(() => [0, 0]);
        let denominator = new Array(N).fill(0);

        proximity_graph.forEach(({source: i, target: j, distance, delta}) => {
            let pi = G[i];
            let pj = G[j];
            const wij = (1/Math.pow(delta, 2));
            nominator[i][0] += (wij * (pj[0] + delta * ((pi[0] - pj[0]) / distance)))
            nominator[i][1] += (wij * (pj[1] + delta * ((pi[1] - pj[1]) / distance)))
            denominator[i] += wij
        })

        G = G.map(([px, py], i) => {
            const gamma_i = gamma([px, py], Gamma);
            return [
                (nominator[i][0] + alpha * gamma_i[0]) / (denominator[i] + alpha),
                (nominator[i][1] + alpha * gamma_i[1]) / (denominator[i] + alpha),
            ]
        })
        // reconstruct the proximity graph
        // with the new layout

        // should augment G with edges from pairs of nodes that overlap? not done
        proximity_graph = create_proximity_graph(G, D, kw);
    } while (--max_iter >= 0 && !end_condition(proximity_graph, G, Gamma))
    
    // round to grid
    G = G.map(([px, py]) => {
        return [
            Math.floor((px - kw / 2) / kw) * kw + (kw / 2),
            Math.floor((py - kh / 2) / kh) * kh + (kh / 2),
        ];
    })
    return G
}