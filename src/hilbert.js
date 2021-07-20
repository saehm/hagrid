import { get_scales, distance } from "./utils.js";

function rotate(size, px, py, rx, ry) {
    if (ry == 0) {
        if (rx == 1) {
            px = size - 1 - px;
            py = size - 1 - py;
        }
        return [py, px];
    }
    return [px, py];
}

/**
 * 
 * @param {*} param0 
 * @param {*} size 
 */
export function hilbert_encode([px, py], size) {
    let n = 0;
    for (let s = size / 2; s > 0; s /= 2) {
        const rx = (px & s) > 0;
        const ry = (py & s) > 0;
        n += (s ** 2) * ((3 * rx) ^ ry);
        [px, py] = rotate(size, px, py, rx, ry);
    }
    return n;
}

export function hilbert_decode(n, size) {
    let t = n;
    let [px, py] = [0, 0];
    for (let s = 1; s < size; s *= 2) {
        const rx = 1 & (t / 2);
        const ry = 1 & (t ^ rx);
        [px, py] = rotate(s, px, py, rx, ry);
        px += (s * rx);
        py += (s * ry);
        t = t >> 2;
    }
    return [px, py];
}

function hilbert_collision(P, p, d, i, size) {
    const e = size ** 2;
    const valid = (p) => p >= 0 && p <= e;
    let pl = p;
    let pr = p;
    while (true) {
        ++pl; --pr;
        const el = !P.has(pl);
        const er = !P.has(pr);
        const vl = valid(pl);
        const vr = valid(pr);
        if (vl && el && !er) {
            P.set(pl, i);
            return;
        } else if (!el && vr && er) {
            P.set(pr, i);
            return;
        } else if (el && er) {
            if (vl && vr) {
                let dl = distance(d, hilbert_decode(pl, size));
                let dr = distance(d, hilbert_decode(pr, size));
                P.set(dl < dr ? pl : pr, i);
                return;
            } else if (vl) {
                P.set(pl, i);
                return;
            } else if (vr) {
                P.set(pr, i);
                return;
            }
        }
    } 
}

export function gridify_hilbert(D, {level, keep_aspect_ratio=false}) {
    const size = 2 ** level;
    const curveLength = 4 ** level;
    const N = D.length;
    const P = new Map();
    const Y = new Array(N).fill(0);
    const scales = get_scales(D, [[0, size - 1], [0, size - 1]], {keep_aspect_ratio: keep_aspect_ratio, round: true});
    D.forEach((d, i) => {
        const [x, y] = scales.map(s => s(d));
        let p = hilbert_encode([x, y], size);
        if (P.has(p)) { // collision detected
            hilbert_collision(P, p, [x, y], i, curveLength);
        } else {
            P.set(p, i);
        }
    })

    for (const [p, i] of P.entries()) {
        Y[i] = hilbert_decode(p, size);
    }

    return Y;
}
