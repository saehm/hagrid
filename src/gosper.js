import { extent, get_scales, distance } from "./utils.js";

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
    const scale = level.s;
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

export function gosper_encode(p, level) {
    const gosper_fractal = create_gosper_fractal(level + 1);
    const V = generate_level(gosper_fractal[level]);
    return find_nearest(p, V)
}

export function gosper_decode(i, level) {
    const gosper_fractal = create_gosper_fractal(level + 1);
    const V = generate_level(gosper_fractal[level]);
    return V[i]
}

export function gosper_curve(level) {
    const gosper_fractal = create_gosper_fractal(level + 1);
    return generate_level(gosper_fractal[level]);
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

function gosper_collision(P, p, d, i, V) {
    const e = V.length - 1;
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
                const dl = distance(d, V[pl]);
                const dr = distance(d, V[pr]);
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

export function gridify_gosper(D, {level, scale_factor = 0.8}) {
    const N = D.length;
    const P = new Map();
    const Y = new Array(N).fill(0);
    const V = gosper_curve(level);
    // distance between grid cells (grid cell size)
    const r = distance(V[0], V[1]) // Math.sqrt(3);

    const grid_extent = [extent(V, d => d[0] * scale_factor), extent(V, d => d[1]* scale_factor)]; 
    const scales = get_scales(D, grid_extent, {round: false});

    D.forEach((d, i) => {
        let [x, y] = scales.map(s => s(d))
        let p = find_nearest([x, y], V);
        if (P.has(p)) {
            gosper_collision(P, p, [x, y], i, V);
        } else {
            P.set(p, i)
        }
    })

    for (const [p, i] of P.entries()) {
        Y[i] = V[p];
    }

    return Y;
}

function find_nearest(p, list) {
    let nearest_index = -1;
    let nearest_distance = Infinity;

    list.forEach((q, i) => {
        const dist = distance(p, q);
        if (dist < nearest_distance) {
            nearest_index = i;
            nearest_distance = dist;
        }
    })

    return nearest_index;
}