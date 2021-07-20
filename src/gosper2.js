/*
https://doi.org/10.3390/sym11060731

Hierarchical Hexagonal Clustering and Indexing
by Vojtěch Uher 1,Petr Gajdoš 1,Václav Snášel 1, Yu-Chi Lai 2 and Michal Radecký 1

1 Department of Computer Science, VŠB-Technical University of Ostrava, Ostrava-Poruba 708 00, Czech Republic
2 Department of Computer Science and Information Engineering, National Taiwan University of Science and Technology, 43, Sec.4, Keelung Rd., Taipei 106, Taiwan
*/
const REC_SQRT7 = 0.3779644730092272; // 1 / Math.sqrt(7)
const SQRT3_DIV_3 = 0.5773502691896257; // Math.sqrt(3) / 3
const ALPHA = 0.3334731722518321; // Math.asin(Math.sqrt(3) / (2 * Math.sqrt(7)))
const SQRT7 = 2.6457513110645907; // Math.sqrt(7)
// Arry for mapping points onto the center pattern index
const pattern_index = [[5, 1, 3], [2, 4, 6]];

// transform to node gosper pattern
const transform_index = [4, 0, 1, 2, 3, 6, 5];

// precomputed array of circles inscribed into Gosper islands of different levels
const inscribed_circles = [0.755928946000, 0.755928946000, 0.750121467308, 0.746782631146, 0.746782631146, 0.746577727521, 0.746348363909, 0.746348363909, 0.746344578768, 0.746327538283, 0.746327538283, 0.746327538283, 0.746326555879, 0.746326555879, 0.746326555879, 0.746326510616, 0.746326510616, 0.746326510616, 0.746326508597, 0.746326508597, 0.746326508597]

function xy2cube([x, y], size) {
    const cx = (x * SQRT3_DIV_3 - y * 1 / 3) / size;
    const cz = y * 2 / 3 / size;
    const cy = -cx - cz;
    return [cx, cy, cz];
}

function cube2xy([cx, cy, cz], size) {
    const x = (3 / 2 * cx) * size;
    const y = (Math.sqrt(3) * (cx / 2 + cz)) * size;
    return [x, y];
}

function cube_abs_diff([ix, iy, iz], [dx, dy, dz]) {
    let [cx, cy, cz] = [
        ix - dx,
        iy - dy,
        iz - dz,
    ].map(Math.abs);

    if (cx > cy && cx > cz) {
        ix = -iy - iz;
    } else if (cy > cz) {
        iy = -ix - iz;
    } else {
        iz = -ix - iy;
    }

    return [ix, iy, iz];
}

function B([x, y, z], thd = 10e-3) {
    const [ax, ay, az] = [x, y, z].map(v => Math.abs(v));
    if (ax < thd && ay < thd && az < thd) {
        return 0;
    } else if (ax > ay && ax > az) { // x dominant
        return x > 0 ? 2 : 5;
    } else if (ay > ax && ay > az) { // y dominant
        return y > 0 ? 4 : 1;
    } else if (az > ax && az > ay) { // z dominant
        return z > 0 ? 6 : 3;
    }
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

function transform([x, y, z]) {
    const nx = (2 * x - z) / 7;
    const nz = (x + 3 * z) / 7;
    return [nx, -nx - nz, nz];
}

function inverse_transform([x, y, z]) {
    const nx = (3 * x + z);
    const nz = (-x + 2 * z);
    return [nx, -nx - nz, nz];
}

export function get_size(D, level) {
    let min_x = Infinity;
    let max_x = -Infinity;
    let min_y = Infinity;
    let max_y = -Infinity;
    D.forEach(([x, y]) => {
        min_x = x < min_x ? x : min_x;
        max_x = x > max_x ? x : max_x;
        min_y = y < min_y ? y : min_y;
        max_y = y > max_y ? y : max_y;
    })

    const half_diagional = Math.sqrt(Math.pow(max_x - min_x, 2) + Math.pow(max_y - min_y, 2)) / 2;
    const s = half_diagional / inscribed_circles[level];
    return s * Math.pow(REC_SQRT7, level);
}

const f_yArr = [-2.0943951023931953, 0, 0, -2.0943951023931953, 0, 2.0943951023931953, 0];
function f_y(i) {
    //[-120, 0, 0, -120, 0, 120, 0].map(d => d / 180 * Math.PI)
    return f_yArr[i];
}

const gammaArr = [0, -1.0471975511965976, -2.0943951023931953, -3.141592653589793, 0, 1.0471975511965976, 2.0943951023931953];
function gamma(i) {
    //[0, 60, 120, 180, 0, -60, -120].map(d => -d / 180 * Math.PI)
    return gammaArr[i];
}

function rotate([x, y], alpha) {
    const sin_alpha = Math.sin(alpha);
    const cos_alpha = Math.cos(alpha);
    return [cos_alpha * x - sin_alpha * y,
    sin_alpha * x + cos_alpha * y]
}

function scale([x, y], s) {
    return [x * s, y * s];
}

//export function gosper_encode(b, level, size) {
export function gosper_encode([xp, yp], level, size) {
    //let size = Math.pow(1 / SQRT7, level - 1);
    const b = new Array(level);
    let [x, y, z] = xy2cube([xp, yp], size);
    let [rx, ry, rz] = cube_round([x, y, z]);
    let p = cube2xy([rx, ry, rz], size);
    for (let i = 1; i <= level; ++i) {
        p = scale(p, Math.pow(REC_SQRT7, i));
        p = rotate(p, -ALPHA * i);
        [x, y, z] = xy2cube(p, size);
        [rx, ry, rz] = cube_round([x, y, z]);
        b[level - i] = B([x - rx, y - ry, z - rz]);
        if (i == level) continue;
        p = cube2xy([rx, ry, rz], size);
        p = scale(p, Math.pow(SQRT7, i));
        p = rotate(p, ALPHA * i);
    }
    //let b = new Array(level); // final position array
    /* let sign; // sign of decimal residue;
    let mini; // hexagon index according to the center pattern

    let [dx, dy, dz] = xy2cube([x, y], size);
    let [ix, iy, iz] = [dx, dy, dz].map(Math.round);
    [ix, iy, iz] = cube_abs_diff([ix, iy, iz], [dx, dy, dz]);

    for (let l = 0; l < level; ++l) {
        //hexcode *= 6;
        [dx, dy, dz] = transform([ix, iy, iz]);
        [ix, iy, iz] = [dx, dy, dz].map(Math.round);

        dx -= ix;
        dy -= iy;
        dz -= iz;

        let dominant = (Math.abs(dz) > Math.abs(dx)) ? 1 : 0 + (Math.abs(dz) > Math.abs(dy)) ? 1 : 0;
        dominant = (dominant == 2 ? dominant : (Math.abs(dy) > Math.abs(dx) ? 1 : 0));
        sign = [dx, dy, dz][dominant] < 0 ? 1 : 0;

        if (Math.abs([dx, dy, dz][dominant]) < 0.00001) {
            mini = 0;
        } else {
            mini = pattern_index[sign][dominant];
        }

        b[level - l - 1] = mini;
    } */

    // next step
    let rotDir = 0;
    let order = true;
    const k = new Array(level);
    for (let i = 0; i < level; ++i) {
        let bi = b[i];
        if (bi != 0 && rotDir != 0) {
            bi += 2 * rotDir;
            if (bi < 1) {
                bi += 6;
            } else if (bi > 6) {
                bi -= 6;
            }
        }
        let ki = transform_index[bi];
        if (ki == 0 || ki == 3) {
            --rotDir;
            if (rotDir < -1) rotDir = 1;
        } else if (ki == 5) {
            ++rotDir;
            if (rotDir > 1) rotDir = -1;
        }
        if (!order) {
            k[i] = 6 - ki;
        } else {
            k[i] = ki;
        }
        if (ki == 0 || ki == 4 || ki == 5) {
            order = !order;
        }
    }
    return k//.reduce((a, b) => a * 7 + b);
}

/**
 * 
 * @param {Array<Number>} t - hierarchical index of gosper curve 
 * @param {Number} level - level of gosper curve
 * @param {Number} size - grid size of goster curve 
 * @returns {[Number, Number]} - 2d position of index {@link t} on a gospercurve with level {@link level}, and grid size {@link size}.
 */
export function gosper_decode(t, level, size) {
    //console.log("decode", t, level, size)
    /* const k = new Uint8Array(level - 1).fill(0);
    for (let i = level - 1; i >= 0; --i) {
        const mod_t = t % 7;
        k[i] = mod_t;
        t = (t - mod_t) / 7;
    } */
    const k = t
    let d = [-size, 0, 0];
    for (let l = 0; l < level - 1; ++l) {
        d = inverse_transform(d);
    }
    d = cube2xy(d, 1);

    //rotate([-size * Math.pow(Math.sqrt(7), level - 1), 0], (0) * ALPHA);
    let c = [0, 0];
    let ord = true;
    for (let i = 0; i < level; ++i) {
        let ki = (i > 0 && !ord) ? (6 - k[i]) : k[i];
        const dk = rotate(d, gamma(ki));
        if (ki != 4) {
            c = [c[0] + dk[0], c[1] + dk[1]];
        }
        d = rotate(d, ALPHA + f_y(ki));
        d = scale(d, REC_SQRT7);
        if (ki == 0 || ki == 4 || ki == 5) {
            ord = !ord;
        }
    }
    return c;
}

function distance([ax, ay], [bx, by]) {
    return Math.hypot(bx - ax, by - ay);
}

function collision(P, p, d, i, level, size) {
    const e = 7 ** level;
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
                let dl = distance(d, gosper_decode(pl, level, size));
                let dr = distance(d, gosper_decode(pr, level, size));
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

export function gridify_gosper(D, { level: level }) {
    const size = get_size(D, level);
    const N = D.length;
    const P = new Map();
    const Y = new Array(N).fill(0);

    D.forEach((d, i) => {
        const [x, y] = d.map(Math.round);
        let p = gosper_encode([x, y], level, size);
        if (P.get(p)) { // collision detected
            collision(P, p, d, i, level, size);
        } else {
            P.set(p, i)
        }
    })

    for (const [p, i] of P.entries()) {
        console.log("i", i, "p", p)
        Y[i] = gosper_decode(p, level, size);
    }
    console.log(P)
    return Y;
}

function addToLevel(p, dir) {
    let overhang = true;
    const q = p.slice();
    const n = p.length - 1;
    for (let i = n; i >= 0 && overhang; --i) {
        let s = q[i] + dir;
        if (s > 6) {
            s = s % 7;
            overhang = true;
        } else if (s < 0) {
            s = 6;
        } else {
            overhang = false;
        }
        q[i] = s;
    }
    return q;
}
