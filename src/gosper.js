function create_gosper_fractal(level) {
    const t1 = ["a", "b", "b", "a", "a", "a", "b"];
    const d1 = [0, 5, 3, 4, 0, 0, 1];

    const t2 = ["a", "b", "b", "b", "a", "a", "b"];
    const d2 = [1, 0, 0, 4, 3, 5, 0];

    const SQRT7 = Math.sqrt(7);

    const result = [{
        "s": SQRT7,
        "t": ["a"],
        "d": [0]
    }];

    for (let l = 1; l < level; ++l) {
        const previous_result = result[l - 1];
        result.push({
            "s": previous_result.s * (1 / SQRT7),
            "t": [],
            "d": [],
        });

        const t = previous_result.t;
        const d = previous_result.d;

        for (let i = 0; i < t.length; ++i) {
            if (t[i] == "a") {
                result[l].t.push(...t1);
                result[l].d.push(...add_mod_6(d[i], d1));
            } else {
                result[l].t.push(...t2);
                result[l].d.push(...add_mod_6(d[i], d2));
            }
        }
    }

    return result;
}

function add_mod_6(m, d) {
    return d.map(e => (m + e) % 6);
}

function generate_level(level) {
    const k1 = .5;
    const k2 = Math.sqrt(3) / 2

    const d_cos = [1, k1, -k1, -1, -k1, k1];
    const d_sin = [0, k2, k2, 0, -k2, -k2];
    scale = level.s;

    const n = level.d.length + 1;
    const V = [[0, 0]];
    for (let i = 1; i < n - 1; ++i) {
        const d = level.d[i];
        const [px, py] = V[i - 1];
        V.push([
            px + scale * d_cos[d],
            py + scale * d_sin[d],
        ]);
    }

    return V;
}


export function gosper_encode([x, y], level) {
    const size = level.s /// Math.sqrt(3);
    let cube = xy2cube([x, y], size);
    cube = cube_round(cube);
    let [nx, ny] = cube2xy(cube, size);
    return knn.search([nx, ny], 1).first.element.index;
}

export function gosper_decode(i) {
    return i >= 0 ? curve[i] : null;
}

function xy2cube([x, y], size) {
    return axial2cube(xy2axial([x, y], size));
}

function cube_round([x, y, z]) {
    let [l, m, n] = [x, y, z].map(d => Math.round(d));
    let [dx, dy, dz] = [l - x, m - y, n - z].map(d => Math.abs(d));
    if (dx > dy && dx > dz) {
        l = -m - n;
    } else if (dy > dz) {
        m = -l - n;
    } else {
        n = -l - m; 
    }
    return [l, m, n].map(d => parseInt(Math.round(d)));
}

function cube2xy([x, y, z], size) {
    return axial2xy(cube2axial([x, y, z]), size);
}

function axial2cube([q, r]) {
    return [q, -q -r, r];
}

function cube2axial([x, y, z]){
    return [x, z];
}

function axial2xy([q, r], size) {
    return [(3 / 2 * q) * size, (Math.sqrt(3) * (q / 2 + r) * size)];
}

function xy2axial([x, y], size){
    return [(x * 2 / 3) / size, (-x / 3 + y * Math.sqrt(3) / 3) / size];
}
