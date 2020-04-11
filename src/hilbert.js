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

function distance([ax, ay], [bx, by]) {
    return Math.sqrt(Math.pow(ax - bx, 2) + Math.pow(ay - by, 2));
}

function hilbert_collision_1(P, p, d, i, size) {
    const distance_increment = distance(d, hilbert_decode(p + 1, size));
    const distance_decrement = distance(d, hilbert_decode(p - 1, size));
    const direction = distance_increment < distance_decrement ? 1 : -1; 
    //const direction = 1
    let flying = i;    
    while (P.has(p)) {
        // swap elements
        const tmp = P.get(p);
        //P.delete(p);
        P.set(p, flying);
        flying = tmp;
        // hop to next place in direction
        p += direction;
    }
    P.set(p, flying);
}

function hilbert_collision_2(P, p, i) {
    while (P.has(p)) {
        p += 1
    }
    P.set(p, i)
}

export async function gridify_hilbert(D, {level, collision}) {
    const size = 2 ** level;
    const N = D.length;
    const P = new Map();
    const Y = new Array(N).fill(0);

    D.forEach((d, i) => {
        const [x, y] = [Math.round(d[0]), Math.round(d[1])];
        let p = hilbert_encode([x, y], size);
        if (P.has(p)) { // collision detected
            // distinguish direction
            if (collision == "new") {
                hilbert_collision_1(P, p, d, i, size);
            } else {
                hilbert_collision_2(P, p, i)
            }
        } else {
            P.set(p, i);
        }
    })

    for (const [p, i] of P) {
        Y[i] = hilbert_decode(p, size);
    }

    return Y;
}
