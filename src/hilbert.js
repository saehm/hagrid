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
    for (let s = size >> 1; s > 0; s >>= 1) {
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
    for (let s = 1; s < size; s <<= 1) {
        const rx = 1 & (t / 2);
        const ry = 1 & (t ^ rx);
        [px, py] = rotate(s, px, py, rx, ry);
        px += (s * rx);
        py += (s * ry);
        t = t >> 2;
    }
    return [px, py];
}