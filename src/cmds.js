import {Delaunay} from "d3-delaunay";

// gamma puts a point outside the boundary back inside,
// 
function gamma([px, py], {x1, x2, y1, y2}) {
    let [nx, ny] = [px, py];
    if (px < x1) nx = x1;
    else if (px > x2) nx = x2;
    if (py < y1) ny = y1;
    else if (py > y2) ny = y2;
    return [nx, ny]
}

function distance([ax, ay], [bx, by]) {
    return Math.hypot(bx - ax, by - ay);//Math.sqrt(((ax - bx) ** 2) + ((ay - by) ** 2));
}

function compute_overlap([px, py], [qx, qy], s) {
    const dx = Math.abs(px - qx);
    const dy = Math.abs(py - qy);
    if (dx < s || dy < s) {
        const ox = Math.abs(Math.min(dx - s, 0));
        const oy = Math.abs(Math.min(dy - s, 0));
        const dd = distance([px, py], [qx, qy]);
        if (ox < oy) {
            return dx == 0 ? (dd * oy) / dy : (dd * ox) / dx;
        } else if (oy <= ox) {
            return dy == 0 ? (dd * ox) / dx : (dd * oy) / dy;
        }
    } else {
        return 0
    }
}

// computes all the stuff needed for cmds
function create_proximity_graph(A, B, s, augment=false) {
    const delaunay = Delaunay.from(A);
    return A.map((_, i) => {
        //return [...delaunay.neighbors(i)].map((j) => {
        const neighbors = [...delaunay.neighbors(i)];
        return A.map((_, j) => {
            /* const pd = D[i];
            const qd = D[j]; */
            const pa = A[i];
            const qa = A[j];
            const po = B[i];
            const qo = B[j];
            // maybe better to save
            const odistance = distance(po, qo);
            const overlap = compute_overlap(pa, qa, s);
            return {
                "source": i,
                "target": j,
                "distance": odistance, //l_ij
                "overlap": overlap, // delta_ij
                "delta": odistance + overlap // d_ij
            }
        }).filter(({target, overlap}) => {
            const keep_overlap = overlap > 1e-2
            const keep_neighbor = neighbors.findIndex(j => target == j) >= 0;
            return augment ? keep_overlap || keep_neighbor : keep_neighbor;
        })
    }).flat()
}

// end condition is
// no overlaps anymore and all points are inside the boundaries
function end_condition(proximity_graph, G, {x1, x2, y1, y2}) {
    const tol = 1;
    const sum_overlaps = proximity_graph
        .reduce((a, b) => a + b.overlap, 0);

    //console.log("sum overlaps", sum_overlaps, proximity_graph.filter(d => d.overlap > tol).length)
    const inside_Gamma = G
        .map(([px, py]) => px >= x1 && px <= x2 && py >= y1 && py <= y2)
        //.map(([px, py]) => (px - x1 < tol) && (px - x2 < tol) && (py - y1 < tol) && (py - y2 < tol))
        .reduce((a, b) => a && b);

    /* console.log("inside Gamma", G
    .map(([px, py]) => px >= x1 - s/2 && px <= x2 + s/2 && py >= y1 - s/2 && py <= y2 + s/2)
        .filter(d => !d).length, 
        inside_Gamma) */
    //console.log(inside_Gamma, sum_overlaps)
    return sum_overlaps < tol && inside_Gamma
}

export async function* gridify_cmds(D, parameters) {
    const N = D.length;
    const Gamma = parameters.Gamma;
    const s = parameters.size;
    let alpha = parameters.alpha ? parameters.alpha : .1;
    let max_iter = parameters.max_iter ? parameters.max_iter : 50;
    
    // copy D
    let G = D.map(([x, y]) => [x, y]);
    // construct proximity graph and 
    // calculate the ideal distances according to (3)
    let proximity_graph = create_proximity_graph(G, D, s, false);
    
    do {
        for (let i = 0; i < 20; ++i) {
            // update step according to (4)
            const nominator = Array.from({length: N}, () => [0, 0]);
            const denominator = Array.from({length: N}, () => 0);

            proximity_graph.forEach(({source: i, target: j, delta}) => {
                let pi = G[i];
                let pj = G[j];
                //if (proximity_graph.findIndex(({source, target}) => source == j && target == i) < k) return;
                const wij = 1 / Math.pow(delta, 2);
                const dist = distance(pi, pj);
                nominator[i][0] += (wij * (pj[0] + delta * ((pi[0] - pj[0]) / dist)))
                nominator[i][1] += (wij * (pj[1] + delta * ((pi[1] - pj[1]) / dist)))
                denominator[i] += wij
            })

            G = G.map((p, i) => {
                const [gamma_px, gamma_py] = gamma(p, Gamma, s);
                return [
                    (nominator[i][0] + alpha * gamma_px) / (denominator[i] + alpha),
                    (nominator[i][1] + alpha * gamma_py) / (denominator[i] + alpha),
                ]
            })

        }
        // reconstruct the proximity graph
        // with the new layout
        proximity_graph = create_proximity_graph(G, D, s, true);

        yield [G, proximity_graph];

    } while (--max_iter > 0 && !end_condition(proximity_graph, G, Gamma, s))
    // round to grid
    G = G.map(d => {
        return d.map(v => {
            return Math.round((v - s / 2) / s) * s + (s / 2)
        })
    });

    proximity_graph = create_proximity_graph(G, s);
    //console.log("fin")
    yield [G, proximity_graph]
    //return G
}