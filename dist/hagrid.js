(function (global, factory) {
    typeof exports === 'object' && typeof module !== 'undefined' ? factory(exports) :
    typeof define === 'function' && define.amd ? define(['exports'], factory) :
    (global = typeof globalThis !== 'undefined' ? globalThis : global || self, factory(global.hagrid = {}));
}(this, (function (exports) { 'use strict';

    var version = "0.1.0";

    function get_bounds(A) {
        let min_x = Infinity;
        let max_x = -Infinity;
        let min_y = Infinity;
        let max_y = -Infinity;

        A.forEach(([x, y, i]) => {
            min_x = Math.min(min_x, x);
            max_x = Math.max(max_x, x);
            min_y = Math.min(min_y, y);
            max_y = Math.max(max_y, y);
        });

        return {
            "x": min_x,
            "y": min_y,
            "width": max_x - min_x,
            "height": max_y - min_y,
        }
    }

    function remap([x, y], [x1, y1, w1, h1], [x2, y2, w2, h2]) {
        return [
            (x - x1) / w1 * w2 + x2, 
            (y - y1) / h1 * h2 + y2,
        ];
    }


    function get_scales(data, [[X_min, X_max], [Y_min, Y_max]], {
        keep_aspect_ratio = false, 
        round = true, 
        x = d => d[0], 
        y = d => d[1]
    }) {
        const X_span = X_max - X_min;
        const Y_span = Y_max - Y_min; 
        let [x_min, x_max] = extent$1(data, x);
        let [y_min, y_max] = extent$1(data, y);
        let x_span = x_max - x_min;
        let y_span = y_max - y_min;

        if (keep_aspect_ratio) {
            let o;
            if (x_span > y_span) {
                o = (x_span - y_span) / 2;
                y_min -= o;
                y_max += o;
            } else {
                o = (y_span - x_span) / 2;
                x_min -= o;
                x_max += o;
            }
        }

        const f = round ? Math.round : d => d;
        const x_scale = d => f((x(d) - x_min) / x_span * X_span + X_min);
        const y_scale = d => f((y(d) - y_min) / y_span * Y_span + Y_min);
        return [x_scale, y_scale];
    }

    function extent$1(data, accessor) {
        let min = Infinity;
        let max = -Infinity;
        for (const d of data) {
            const v = accessor(d);
            min = Math.min(min, v);
            max = Math.max(max, v);
        }
        return [min, max]    
    }

    function distance$2([ax, ay], [bx, by]) {
        return Math.hypot(bx - ax, by - ay);//Math.sqrt(Math.pow(ax - bx, 2) + Math.pow(ay - by, 2));
    }

    var utils = /*#__PURE__*/Object.freeze({
        __proto__: null,
        get_bounds: get_bounds,
        remap: remap,
        get_scales: get_scales,
        extent: extent$1,
        distance: distance$2
    });

    function rotate$1(size, px, py, rx, ry) {
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
    function hilbert_encode([px, py], size) {
        let n = 0;
        for (let s = size / 2; s > 0; s /= 2) {
            const rx = (px & s) > 0;
            const ry = (py & s) > 0;
            n += (s ** 2) * ((3 * rx) ^ ry);
            [px, py] = rotate$1(size, px, py, rx, ry);
        }
        return n;
    }

    function hilbert_decode(n, size) {
        let t = n;
        let [px, py] = [0, 0];
        for (let s = 1; s < size; s *= 2) {
            const rx = 1 & (t / 2);
            const ry = 1 & (t ^ rx);
            [px, py] = rotate$1(s, px, py, rx, ry);
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
                    let dl = distance$2(d, hilbert_decode(pl, size));
                    let dr = distance$2(d, hilbert_decode(pr, size));
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

    function gridify_hilbert(D, {level, keep_aspect_ratio=false}) {
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
        });

        for (const [p, i] of P.entries()) {
            Y[i] = hilbert_decode(p, size);
        }

        return Y;
    }

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
        const k2 = Math.sqrt(3) / 2;
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

    function gosper_encode$1(p, level) {
        const gosper_fractal = create_gosper_fractal(level + 1);
        const V = generate_level(gosper_fractal[level]);
        return find_nearest(p, V)
    }

    function gosper_decode$1(i, level) {
        const gosper_fractal = create_gosper_fractal(level + 1);
        const V = generate_level(gosper_fractal[level]);
        return V[i]
    }

    function gosper_curve(level) {
        const gosper_fractal = create_gosper_fractal(level + 1);
        return generate_level(gosper_fractal[level]);
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
                    const dl = distance$2(d, V[pl]);
                    const dr = distance$2(d, V[pr]);
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

    function gridify_gosper$1(D, {level, scale_factor = 0.8}) {
        const N = D.length;
        const P = new Map();
        const Y = new Array(N).fill(0);
        const V = gosper_curve(level);
        // distance between grid cells (grid cell size)
        distance$2(V[0], V[1]); // Math.sqrt(3);

        const grid_extent = [extent$1(V, d => d[0] * scale_factor), extent$1(V, d => d[1]* scale_factor)]; 
        const scales = get_scales(D, grid_extent, {round: false});

        D.forEach((d, i) => {
            let [x, y] = scales.map(s => s(d));
            let p = find_nearest([x, y], V);
            if (P.has(p)) {
                gosper_collision(P, p, [x, y], i, V);
            } else {
                P.set(p, i);
            }
        });

        for (const [p, i] of P.entries()) {
            Y[i] = V[p];
        }

        return Y;
    }

    function find_nearest(p, list) {
        let nearest_index = -1;
        let nearest_distance = Infinity;

        list.forEach((q, i) => {
            const dist = distance$2(p, q);
            if (dist < nearest_distance) {
                nearest_index = i;
                nearest_distance = dist;
            }
        });

        return nearest_index;
    }

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
                const [P1, P2] = split$1(P, 1, r_half * s); // split by y
                dgrid(G, P1, (r_half), s, i, j);
                dgrid(G, P2, (r - r_half), s, (i + r_half), j);
            } else {
                const s_half = Math.ceil(s / 2);
                const [P1, P2] = split$1(P, 0, s_half * r); // split by x
                dgrid(G, P1, r, s_half, i, j);
                dgrid(G, P2, r, (s - s_half), i, (j + s_half));
            }
        }
    }

    function split$1(P, dim, pos) {
        const sorted = P.sort((a, b) => a[dim] - b[dim]);
        return [sorted.slice(0, pos), sorted.slice(pos)];
    }

    function gridify_dgrid(D, parameters) {
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
            cols = Math.ceil(N / rows);
        } else {
            throw "wrong parameters!"
        }
        let G = [];
        let added_index = D.map(([x, y], i) => [x, y, i]);
        dgrid(G, added_index, rows, cols);
        return G.sort((a, b) => a.d[2] - b.d[2]).map(g => [g.j, g.i])
    }

    function nmap(G, P, R, direction=true) {
        const N = P.length;
        const N_half = Math.floor(N/2);
        if (N == 1) {
            R.d = P[0];
            G.push(R);
        } else {
            const dim = direction ? 0 : 1;
            P = P.sort((a, b) => b[dim] - a[dim]);
            let [Da, Db] = split(P, N_half);
            let Ra = Object.assign({}, R);
            let Rb = Object.assign({}, R);
            if (direction) {
                const [wRa, wRb] = sizes(Da, Db, R.width);
                const bh = (Da[Da.length - 1][0] + Db[0][0]) / 2;
                Ra.width = bh - R.x;
                Rb.x += Ra.width; 
                Rb.width = R.width - Ra.width;
                const HRa = [
                    [wRa / Ra.width, 0, R.x * (1 - wRa / Ra.width)],
                    [0, 1, 0]
                ];
                const HRb = [
                    [wRb / Rb.width, 0, (R.x + R.width) * (1 - wRb / Rb.width)],
                    [0, 1, 0]
                ];          
                nmap(G, affine_transform_points(Da, HRa), affine_transform_rect(Ra, HRa), !direction);
                nmap(G, affine_transform_points(Db, HRb), affine_transform_rect(Rb, HRb), !direction);
            } else {
                const [hRa, hRb] = sizes(Da, Db, R.height);
                const bv = (Da[Da.length - 1][1] + Db[0][1]) / 2;
                Ra.height = bv - R.y;
                Rb.y = R.y + Ra.height; 
                Rb.height = R.height - Ra.height;
                const VRa = [
                    [1, 0, 0],
                    [0, hRa / Ra.height, R.y * (1 - hRa / Ra.height)]
                ];
                const VRb = [
                    [1, 0, 0],
                    [0, hRb / Rb.height, (R.y + R.height) * (1 - hRb / Rb.height)]
                ];
                nmap(G, affine_transform_points(Da, VRa), affine_transform_rect(Ra, VRa), !direction);
                nmap(G, affine_transform_points(Db, VRb), affine_transform_rect(Rb, VRb), !direction);
            }
        }
    }

    function affine_transform_points(D, [[scale_x, shear_x, translate_x], [scale_y, shear_y, translate_y]]) {
        return D.map(([x, y, i]) => [
            x * scale_x + y * shear_x + translate_x,
            x * scale_y + y * shear_y + translate_y,
            i
        ]);
    }

    function affine_transform_rect({x, y, width, height}, [[scale_x, shear_x, translate_x], [scale_y, shear_y, translate_y]]) {
        const nx = x * scale_x + y * shear_x + translate_x;
        const ny = x * scale_y + y * shear_y + translate_y;
        const nwidth = (x + width) * scale_x + (y + height) * shear_x + translate_x - nx;
        const nheight = (x + width) * scale_y + (y + height) * shear_y + translate_y - ny;
        return {
            "x": nx, 
            "y": ny,
            "width": nwidth,
            "height": nheight,
        }
    }

    function split(P, pos) {
        const N = P.length;
        return [P.slice(0, pos), P.slice(pos, N)];
    }

    function sizes(Da, Db, size) {
        const pA = Da.length;
        const pB = Db.length;
        return [pA / (pA + pB) * size, pA / (pA + pB) * size]
    }

    // https://github.com/sebastian-meier/nmap-squared.js/blob/master/nmap-squared.js
    function squared(D, {x: x0, y: y0, width, height}) {
        const BB2 = get_bounds(D);
        let sx = width / BB2.width;
        let sy = height / BB2.height;
        if (sx < sy) {
            sy = sx;
        } else {
            sx = sy;
        }
        D = D.map(([x, y, i]) => [(x - x0) * sx + 1, (y - y0) * sy + 1, i]);
        const grid_size = Math.ceil(Math.ceil(Math.sqrt(D.length)) / 4) * 4;
        const sq_amount = Math.pow(grid_size, 2);
        let sq_missing = sq_amount - D.length;
        let ep_num = 1 / sq_missing;

        // border method
        let extra = [];
        // top and bottom border
        for (let x = 0; x <= grid_size - ep_num; x += ep_num) {
            extra.push([(width * sx) / grid_size * x + ep_num * Math.random(), ep_num * Math.random(), undefined]);
            extra.push([(width * sx) / grid_size * x + ep_num * Math.random(), (height * sy) + 1 + ep_num * Math.random(), undefined]);
        }

        // left and right border
        for (let y = 0; y <= grid_size - ep_num; y += ep_num) {
            extra.push([ep_num * Math.random(), ((height * sy) / grid_size * y) + ep_num * Math.random(), undefined]);
            extra.push([(width * sx) + 1 + ep_num * Math.random(), ((height * sy) / grid_size * y) + ep_num * Math.random(), undefined]);
        }

        while (sq_missing > 0) {
            extra = calcDist(extra, D);
            D.push(extra[extra.length - 1]);
            --sq_missing;
        }

        return D;
    }

    function calcDist(A, B) {
        for (let i = 0; i < A.length; ++i) {
            const a = A[i];
            for (let j = 0; j < B.length; ++j) {
                const b = B[j];
                let t_dist = Math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2);
                if (b[2] == undefined) {
                    t_dist *= 2;
                }
                if (t_dist < a.dist) {
                    a.dist = t_dist;
                }
            }
        }
        return A.sort((a, b) => a.dist - b.dist);
    }

    function gridify_nmap(D, parameters) {
        let G = [];
        console.log(parameters);
        if (!"BB" in parameters || parameters.BB == null) {
            console.log("no BB");
            parameters.BB = get_bounds(D);
        }
        let added_index = D.map(([x, y], i) => [x, y, i]);
        if ("squared" in parameters && parameters.squared == true) {
            added_index = squared(added_index, parameters.BB);
        }
        nmap(G, added_index, parameters.BB, false);
        return G
            .filter(d => d.d[2] != undefined)
            .sort((a, b) => a.d[2] - b.d[2])
            .map(({x, y, width, height}) => [x + width / 2, y + height / 2])
    }

    const epsilon$1 = 1.1102230246251565e-16;
    const splitter = 134217729;
    const resulterrbound = (3 + 8 * epsilon$1) * epsilon$1;

    // fast_expansion_sum_zeroelim routine from oritinal code
    function sum(elen, e, flen, f, h) {
        let Q, Qnew, hh, bvirt;
        let enow = e[0];
        let fnow = f[0];
        let eindex = 0;
        let findex = 0;
        if ((fnow > enow) === (fnow > -enow)) {
            Q = enow;
            enow = e[++eindex];
        } else {
            Q = fnow;
            fnow = f[++findex];
        }
        let hindex = 0;
        if (eindex < elen && findex < flen) {
            if ((fnow > enow) === (fnow > -enow)) {
                Qnew = enow + Q;
                hh = Q - (Qnew - enow);
                enow = e[++eindex];
            } else {
                Qnew = fnow + Q;
                hh = Q - (Qnew - fnow);
                fnow = f[++findex];
            }
            Q = Qnew;
            if (hh !== 0) {
                h[hindex++] = hh;
            }
            while (eindex < elen && findex < flen) {
                if ((fnow > enow) === (fnow > -enow)) {
                    Qnew = Q + enow;
                    bvirt = Qnew - Q;
                    hh = Q - (Qnew - bvirt) + (enow - bvirt);
                    enow = e[++eindex];
                } else {
                    Qnew = Q + fnow;
                    bvirt = Qnew - Q;
                    hh = Q - (Qnew - bvirt) + (fnow - bvirt);
                    fnow = f[++findex];
                }
                Q = Qnew;
                if (hh !== 0) {
                    h[hindex++] = hh;
                }
            }
        }
        while (eindex < elen) {
            Qnew = Q + enow;
            bvirt = Qnew - Q;
            hh = Q - (Qnew - bvirt) + (enow - bvirt);
            enow = e[++eindex];
            Q = Qnew;
            if (hh !== 0) {
                h[hindex++] = hh;
            }
        }
        while (findex < flen) {
            Qnew = Q + fnow;
            bvirt = Qnew - Q;
            hh = Q - (Qnew - bvirt) + (fnow - bvirt);
            fnow = f[++findex];
            Q = Qnew;
            if (hh !== 0) {
                h[hindex++] = hh;
            }
        }
        if (Q !== 0 || hindex === 0) {
            h[hindex++] = Q;
        }
        return hindex;
    }

    function estimate(elen, e) {
        let Q = e[0];
        for (let i = 1; i < elen; i++) Q += e[i];
        return Q;
    }

    function vec(n) {
        return new Float64Array(n);
    }

    const ccwerrboundA = (3 + 16 * epsilon$1) * epsilon$1;
    const ccwerrboundB = (2 + 12 * epsilon$1) * epsilon$1;
    const ccwerrboundC = (9 + 64 * epsilon$1) * epsilon$1 * epsilon$1;

    const B$1 = vec(4);
    const C1 = vec(8);
    const C2 = vec(12);
    const D = vec(16);
    const u = vec(4);

    function orient2dadapt(ax, ay, bx, by, cx, cy, detsum) {
        let acxtail, acytail, bcxtail, bcytail;
        let bvirt, c, ahi, alo, bhi, blo, _i, _j, _0, s1, s0, t1, t0, u3;

        const acx = ax - cx;
        const bcx = bx - cx;
        const acy = ay - cy;
        const bcy = by - cy;

        s1 = acx * bcy;
        c = splitter * acx;
        ahi = c - (c - acx);
        alo = acx - ahi;
        c = splitter * bcy;
        bhi = c - (c - bcy);
        blo = bcy - bhi;
        s0 = alo * blo - (s1 - ahi * bhi - alo * bhi - ahi * blo);
        t1 = acy * bcx;
        c = splitter * acy;
        ahi = c - (c - acy);
        alo = acy - ahi;
        c = splitter * bcx;
        bhi = c - (c - bcx);
        blo = bcx - bhi;
        t0 = alo * blo - (t1 - ahi * bhi - alo * bhi - ahi * blo);
        _i = s0 - t0;
        bvirt = s0 - _i;
        B$1[0] = s0 - (_i + bvirt) + (bvirt - t0);
        _j = s1 + _i;
        bvirt = _j - s1;
        _0 = s1 - (_j - bvirt) + (_i - bvirt);
        _i = _0 - t1;
        bvirt = _0 - _i;
        B$1[1] = _0 - (_i + bvirt) + (bvirt - t1);
        u3 = _j + _i;
        bvirt = u3 - _j;
        B$1[2] = _j - (u3 - bvirt) + (_i - bvirt);
        B$1[3] = u3;

        let det = estimate(4, B$1);
        let errbound = ccwerrboundB * detsum;
        if (det >= errbound || -det >= errbound) {
            return det;
        }

        bvirt = ax - acx;
        acxtail = ax - (acx + bvirt) + (bvirt - cx);
        bvirt = bx - bcx;
        bcxtail = bx - (bcx + bvirt) + (bvirt - cx);
        bvirt = ay - acy;
        acytail = ay - (acy + bvirt) + (bvirt - cy);
        bvirt = by - bcy;
        bcytail = by - (bcy + bvirt) + (bvirt - cy);

        if (acxtail === 0 && acytail === 0 && bcxtail === 0 && bcytail === 0) {
            return det;
        }

        errbound = ccwerrboundC * detsum + resulterrbound * Math.abs(det);
        det += (acx * bcytail + bcy * acxtail) - (acy * bcxtail + bcx * acytail);
        if (det >= errbound || -det >= errbound) return det;

        s1 = acxtail * bcy;
        c = splitter * acxtail;
        ahi = c - (c - acxtail);
        alo = acxtail - ahi;
        c = splitter * bcy;
        bhi = c - (c - bcy);
        blo = bcy - bhi;
        s0 = alo * blo - (s1 - ahi * bhi - alo * bhi - ahi * blo);
        t1 = acytail * bcx;
        c = splitter * acytail;
        ahi = c - (c - acytail);
        alo = acytail - ahi;
        c = splitter * bcx;
        bhi = c - (c - bcx);
        blo = bcx - bhi;
        t0 = alo * blo - (t1 - ahi * bhi - alo * bhi - ahi * blo);
        _i = s0 - t0;
        bvirt = s0 - _i;
        u[0] = s0 - (_i + bvirt) + (bvirt - t0);
        _j = s1 + _i;
        bvirt = _j - s1;
        _0 = s1 - (_j - bvirt) + (_i - bvirt);
        _i = _0 - t1;
        bvirt = _0 - _i;
        u[1] = _0 - (_i + bvirt) + (bvirt - t1);
        u3 = _j + _i;
        bvirt = u3 - _j;
        u[2] = _j - (u3 - bvirt) + (_i - bvirt);
        u[3] = u3;
        const C1len = sum(4, B$1, 4, u, C1);

        s1 = acx * bcytail;
        c = splitter * acx;
        ahi = c - (c - acx);
        alo = acx - ahi;
        c = splitter * bcytail;
        bhi = c - (c - bcytail);
        blo = bcytail - bhi;
        s0 = alo * blo - (s1 - ahi * bhi - alo * bhi - ahi * blo);
        t1 = acy * bcxtail;
        c = splitter * acy;
        ahi = c - (c - acy);
        alo = acy - ahi;
        c = splitter * bcxtail;
        bhi = c - (c - bcxtail);
        blo = bcxtail - bhi;
        t0 = alo * blo - (t1 - ahi * bhi - alo * bhi - ahi * blo);
        _i = s0 - t0;
        bvirt = s0 - _i;
        u[0] = s0 - (_i + bvirt) + (bvirt - t0);
        _j = s1 + _i;
        bvirt = _j - s1;
        _0 = s1 - (_j - bvirt) + (_i - bvirt);
        _i = _0 - t1;
        bvirt = _0 - _i;
        u[1] = _0 - (_i + bvirt) + (bvirt - t1);
        u3 = _j + _i;
        bvirt = u3 - _j;
        u[2] = _j - (u3 - bvirt) + (_i - bvirt);
        u[3] = u3;
        const C2len = sum(C1len, C1, 4, u, C2);

        s1 = acxtail * bcytail;
        c = splitter * acxtail;
        ahi = c - (c - acxtail);
        alo = acxtail - ahi;
        c = splitter * bcytail;
        bhi = c - (c - bcytail);
        blo = bcytail - bhi;
        s0 = alo * blo - (s1 - ahi * bhi - alo * bhi - ahi * blo);
        t1 = acytail * bcxtail;
        c = splitter * acytail;
        ahi = c - (c - acytail);
        alo = acytail - ahi;
        c = splitter * bcxtail;
        bhi = c - (c - bcxtail);
        blo = bcxtail - bhi;
        t0 = alo * blo - (t1 - ahi * bhi - alo * bhi - ahi * blo);
        _i = s0 - t0;
        bvirt = s0 - _i;
        u[0] = s0 - (_i + bvirt) + (bvirt - t0);
        _j = s1 + _i;
        bvirt = _j - s1;
        _0 = s1 - (_j - bvirt) + (_i - bvirt);
        _i = _0 - t1;
        bvirt = _0 - _i;
        u[1] = _0 - (_i + bvirt) + (bvirt - t1);
        u3 = _j + _i;
        bvirt = u3 - _j;
        u[2] = _j - (u3 - bvirt) + (_i - bvirt);
        u[3] = u3;
        const Dlen = sum(C2len, C2, 4, u, D);

        return D[Dlen - 1];
    }

    function orient2d(ax, ay, bx, by, cx, cy) {
        const detleft = (ay - cy) * (bx - cx);
        const detright = (ax - cx) * (by - cy);
        const det = detleft - detright;

        if (detleft === 0 || detright === 0 || (detleft > 0) !== (detright > 0)) return det;

        const detsum = Math.abs(detleft + detright);
        if (Math.abs(det) >= ccwerrboundA * detsum) return det;

        return -orient2dadapt(ax, ay, bx, by, cx, cy, detsum);
    }

    const EPSILON = Math.pow(2, -52);
    const EDGE_STACK = new Uint32Array(512);

    class Delaunator {

        static from(points, getX = defaultGetX, getY = defaultGetY) {
            const n = points.length;
            const coords = new Float64Array(n * 2);

            for (let i = 0; i < n; i++) {
                const p = points[i];
                coords[2 * i] = getX(p);
                coords[2 * i + 1] = getY(p);
            }

            return new Delaunator(coords);
        }

        constructor(coords) {
            const n = coords.length >> 1;
            if (n > 0 && typeof coords[0] !== 'number') throw new Error('Expected coords to contain numbers.');

            this.coords = coords;

            // arrays that will store the triangulation graph
            const maxTriangles = Math.max(2 * n - 5, 0);
            this._triangles = new Uint32Array(maxTriangles * 3);
            this._halfedges = new Int32Array(maxTriangles * 3);

            // temporary arrays for tracking the edges of the advancing convex hull
            this._hashSize = Math.ceil(Math.sqrt(n));
            this._hullPrev = new Uint32Array(n); // edge to prev edge
            this._hullNext = new Uint32Array(n); // edge to next edge
            this._hullTri = new Uint32Array(n); // edge to adjacent triangle
            this._hullHash = new Int32Array(this._hashSize).fill(-1); // angular edge hash

            // temporary arrays for sorting points
            this._ids = new Uint32Array(n);
            this._dists = new Float64Array(n);

            this.update();
        }

        update() {
            const {coords, _hullPrev: hullPrev, _hullNext: hullNext, _hullTri: hullTri, _hullHash: hullHash} =  this;
            const n = coords.length >> 1;

            // populate an array of point indices; calculate input data bbox
            let minX = Infinity;
            let minY = Infinity;
            let maxX = -Infinity;
            let maxY = -Infinity;

            for (let i = 0; i < n; i++) {
                const x = coords[2 * i];
                const y = coords[2 * i + 1];
                if (x < minX) minX = x;
                if (y < minY) minY = y;
                if (x > maxX) maxX = x;
                if (y > maxY) maxY = y;
                this._ids[i] = i;
            }
            const cx = (minX + maxX) / 2;
            const cy = (minY + maxY) / 2;

            let minDist = Infinity;
            let i0, i1, i2;

            // pick a seed point close to the center
            for (let i = 0; i < n; i++) {
                const d = dist(cx, cy, coords[2 * i], coords[2 * i + 1]);
                if (d < minDist) {
                    i0 = i;
                    minDist = d;
                }
            }
            const i0x = coords[2 * i0];
            const i0y = coords[2 * i0 + 1];

            minDist = Infinity;

            // find the point closest to the seed
            for (let i = 0; i < n; i++) {
                if (i === i0) continue;
                const d = dist(i0x, i0y, coords[2 * i], coords[2 * i + 1]);
                if (d < minDist && d > 0) {
                    i1 = i;
                    minDist = d;
                }
            }
            let i1x = coords[2 * i1];
            let i1y = coords[2 * i1 + 1];

            let minRadius = Infinity;

            // find the third point which forms the smallest circumcircle with the first two
            for (let i = 0; i < n; i++) {
                if (i === i0 || i === i1) continue;
                const r = circumradius(i0x, i0y, i1x, i1y, coords[2 * i], coords[2 * i + 1]);
                if (r < minRadius) {
                    i2 = i;
                    minRadius = r;
                }
            }
            let i2x = coords[2 * i2];
            let i2y = coords[2 * i2 + 1];

            if (minRadius === Infinity) {
                // order collinear points by dx (or dy if all x are identical)
                // and return the list as a hull
                for (let i = 0; i < n; i++) {
                    this._dists[i] = (coords[2 * i] - coords[0]) || (coords[2 * i + 1] - coords[1]);
                }
                quicksort(this._ids, this._dists, 0, n - 1);
                const hull = new Uint32Array(n);
                let j = 0;
                for (let i = 0, d0 = -Infinity; i < n; i++) {
                    const id = this._ids[i];
                    if (this._dists[id] > d0) {
                        hull[j++] = id;
                        d0 = this._dists[id];
                    }
                }
                this.hull = hull.subarray(0, j);
                this.triangles = new Uint32Array(0);
                this.halfedges = new Uint32Array(0);
                return;
            }

            // swap the order of the seed points for counter-clockwise orientation
            if (orient2d(i0x, i0y, i1x, i1y, i2x, i2y) < 0) {
                const i = i1;
                const x = i1x;
                const y = i1y;
                i1 = i2;
                i1x = i2x;
                i1y = i2y;
                i2 = i;
                i2x = x;
                i2y = y;
            }

            const center = circumcenter(i0x, i0y, i1x, i1y, i2x, i2y);
            this._cx = center.x;
            this._cy = center.y;

            for (let i = 0; i < n; i++) {
                this._dists[i] = dist(coords[2 * i], coords[2 * i + 1], center.x, center.y);
            }

            // sort the points by distance from the seed triangle circumcenter
            quicksort(this._ids, this._dists, 0, n - 1);

            // set up the seed triangle as the starting hull
            this._hullStart = i0;
            let hullSize = 3;

            hullNext[i0] = hullPrev[i2] = i1;
            hullNext[i1] = hullPrev[i0] = i2;
            hullNext[i2] = hullPrev[i1] = i0;

            hullTri[i0] = 0;
            hullTri[i1] = 1;
            hullTri[i2] = 2;

            hullHash.fill(-1);
            hullHash[this._hashKey(i0x, i0y)] = i0;
            hullHash[this._hashKey(i1x, i1y)] = i1;
            hullHash[this._hashKey(i2x, i2y)] = i2;

            this.trianglesLen = 0;
            this._addTriangle(i0, i1, i2, -1, -1, -1);

            for (let k = 0, xp, yp; k < this._ids.length; k++) {
                const i = this._ids[k];
                const x = coords[2 * i];
                const y = coords[2 * i + 1];

                // skip near-duplicate points
                if (k > 0 && Math.abs(x - xp) <= EPSILON && Math.abs(y - yp) <= EPSILON) continue;
                xp = x;
                yp = y;

                // skip seed triangle points
                if (i === i0 || i === i1 || i === i2) continue;

                // find a visible edge on the convex hull using edge hash
                let start = 0;
                for (let j = 0, key = this._hashKey(x, y); j < this._hashSize; j++) {
                    start = hullHash[(key + j) % this._hashSize];
                    if (start !== -1 && start !== hullNext[start]) break;
                }

                start = hullPrev[start];
                let e = start, q;
                while (q = hullNext[e], orient2d(x, y, coords[2 * e], coords[2 * e + 1], coords[2 * q], coords[2 * q + 1]) >= 0) {
                    e = q;
                    if (e === start) {
                        e = -1;
                        break;
                    }
                }
                if (e === -1) continue; // likely a near-duplicate point; skip it

                // add the first triangle from the point
                let t = this._addTriangle(e, i, hullNext[e], -1, -1, hullTri[e]);

                // recursively flip triangles from the point until they satisfy the Delaunay condition
                hullTri[i] = this._legalize(t + 2);
                hullTri[e] = t; // keep track of boundary triangles on the hull
                hullSize++;

                // walk forward through the hull, adding more triangles and flipping recursively
                let n = hullNext[e];
                while (q = hullNext[n], orient2d(x, y, coords[2 * n], coords[2 * n + 1], coords[2 * q], coords[2 * q + 1]) < 0) {
                    t = this._addTriangle(n, i, q, hullTri[i], -1, hullTri[n]);
                    hullTri[i] = this._legalize(t + 2);
                    hullNext[n] = n; // mark as removed
                    hullSize--;
                    n = q;
                }

                // walk backward from the other side, adding more triangles and flipping
                if (e === start) {
                    while (q = hullPrev[e], orient2d(x, y, coords[2 * q], coords[2 * q + 1], coords[2 * e], coords[2 * e + 1]) < 0) {
                        t = this._addTriangle(q, i, e, -1, hullTri[e], hullTri[q]);
                        this._legalize(t + 2);
                        hullTri[q] = t;
                        hullNext[e] = e; // mark as removed
                        hullSize--;
                        e = q;
                    }
                }

                // update the hull indices
                this._hullStart = hullPrev[i] = e;
                hullNext[e] = hullPrev[n] = i;
                hullNext[i] = n;

                // save the two new edges in the hash table
                hullHash[this._hashKey(x, y)] = i;
                hullHash[this._hashKey(coords[2 * e], coords[2 * e + 1])] = e;
            }

            this.hull = new Uint32Array(hullSize);
            for (let i = 0, e = this._hullStart; i < hullSize; i++) {
                this.hull[i] = e;
                e = hullNext[e];
            }

            // trim typed triangle mesh arrays
            this.triangles = this._triangles.subarray(0, this.trianglesLen);
            this.halfedges = this._halfedges.subarray(0, this.trianglesLen);
        }

        _hashKey(x, y) {
            return Math.floor(pseudoAngle(x - this._cx, y - this._cy) * this._hashSize) % this._hashSize;
        }

        _legalize(a) {
            const {_triangles: triangles, _halfedges: halfedges, coords} = this;

            let i = 0;
            let ar = 0;

            // recursion eliminated with a fixed-size stack
            while (true) {
                const b = halfedges[a];

                /* if the pair of triangles doesn't satisfy the Delaunay condition
                 * (p1 is inside the circumcircle of [p0, pl, pr]), flip them,
                 * then do the same check/flip recursively for the new pair of triangles
                 *
                 *           pl                    pl
                 *          /||\                  /  \
                 *       al/ || \bl            al/    \a
                 *        /  ||  \              /      \
                 *       /  a||b  \    flip    /___ar___\
                 *     p0\   ||   /p1   =>   p0\---bl---/p1
                 *        \  ||  /              \      /
                 *       ar\ || /br             b\    /br
                 *          \||/                  \  /
                 *           pr                    pr
                 */
                const a0 = a - a % 3;
                ar = a0 + (a + 2) % 3;

                if (b === -1) { // convex hull edge
                    if (i === 0) break;
                    a = EDGE_STACK[--i];
                    continue;
                }

                const b0 = b - b % 3;
                const al = a0 + (a + 1) % 3;
                const bl = b0 + (b + 2) % 3;

                const p0 = triangles[ar];
                const pr = triangles[a];
                const pl = triangles[al];
                const p1 = triangles[bl];

                const illegal = inCircle(
                    coords[2 * p0], coords[2 * p0 + 1],
                    coords[2 * pr], coords[2 * pr + 1],
                    coords[2 * pl], coords[2 * pl + 1],
                    coords[2 * p1], coords[2 * p1 + 1]);

                if (illegal) {
                    triangles[a] = p1;
                    triangles[b] = p0;

                    const hbl = halfedges[bl];

                    // edge swapped on the other side of the hull (rare); fix the halfedge reference
                    if (hbl === -1) {
                        let e = this._hullStart;
                        do {
                            if (this._hullTri[e] === bl) {
                                this._hullTri[e] = a;
                                break;
                            }
                            e = this._hullPrev[e];
                        } while (e !== this._hullStart);
                    }
                    this._link(a, hbl);
                    this._link(b, halfedges[ar]);
                    this._link(ar, bl);

                    const br = b0 + (b + 1) % 3;

                    // don't worry about hitting the cap: it can only happen on extremely degenerate input
                    if (i < EDGE_STACK.length) {
                        EDGE_STACK[i++] = br;
                    }
                } else {
                    if (i === 0) break;
                    a = EDGE_STACK[--i];
                }
            }

            return ar;
        }

        _link(a, b) {
            this._halfedges[a] = b;
            if (b !== -1) this._halfedges[b] = a;
        }

        // add a new triangle given vertex indices and adjacent half-edge ids
        _addTriangle(i0, i1, i2, a, b, c) {
            const t = this.trianglesLen;

            this._triangles[t] = i0;
            this._triangles[t + 1] = i1;
            this._triangles[t + 2] = i2;

            this._link(t, a);
            this._link(t + 1, b);
            this._link(t + 2, c);

            this.trianglesLen += 3;

            return t;
        }
    }

    // monotonically increases with real angle, but doesn't need expensive trigonometry
    function pseudoAngle(dx, dy) {
        const p = dx / (Math.abs(dx) + Math.abs(dy));
        return (dy > 0 ? 3 - p : 1 + p) / 4; // [0..1]
    }

    function dist(ax, ay, bx, by) {
        const dx = ax - bx;
        const dy = ay - by;
        return dx * dx + dy * dy;
    }

    function inCircle(ax, ay, bx, by, cx, cy, px, py) {
        const dx = ax - px;
        const dy = ay - py;
        const ex = bx - px;
        const ey = by - py;
        const fx = cx - px;
        const fy = cy - py;

        const ap = dx * dx + dy * dy;
        const bp = ex * ex + ey * ey;
        const cp = fx * fx + fy * fy;

        return dx * (ey * cp - bp * fy) -
               dy * (ex * cp - bp * fx) +
               ap * (ex * fy - ey * fx) < 0;
    }

    function circumradius(ax, ay, bx, by, cx, cy) {
        const dx = bx - ax;
        const dy = by - ay;
        const ex = cx - ax;
        const ey = cy - ay;

        const bl = dx * dx + dy * dy;
        const cl = ex * ex + ey * ey;
        const d = 0.5 / (dx * ey - dy * ex);

        const x = (ey * bl - dy * cl) * d;
        const y = (dx * cl - ex * bl) * d;

        return x * x + y * y;
    }

    function circumcenter(ax, ay, bx, by, cx, cy) {
        const dx = bx - ax;
        const dy = by - ay;
        const ex = cx - ax;
        const ey = cy - ay;

        const bl = dx * dx + dy * dy;
        const cl = ex * ex + ey * ey;
        const d = 0.5 / (dx * ey - dy * ex);

        const x = ax + (ey * bl - dy * cl) * d;
        const y = ay + (dx * cl - ex * bl) * d;

        return {x, y};
    }

    function quicksort(ids, dists, left, right) {
        if (right - left <= 20) {
            for (let i = left + 1; i <= right; i++) {
                const temp = ids[i];
                const tempDist = dists[temp];
                let j = i - 1;
                while (j >= left && dists[ids[j]] > tempDist) ids[j + 1] = ids[j--];
                ids[j + 1] = temp;
            }
        } else {
            const median = (left + right) >> 1;
            let i = left + 1;
            let j = right;
            swap(ids, median, i);
            if (dists[ids[left]] > dists[ids[right]]) swap(ids, left, right);
            if (dists[ids[i]] > dists[ids[right]]) swap(ids, i, right);
            if (dists[ids[left]] > dists[ids[i]]) swap(ids, left, i);

            const temp = ids[i];
            const tempDist = dists[temp];
            while (true) {
                do i++; while (dists[ids[i]] < tempDist);
                do j--; while (dists[ids[j]] > tempDist);
                if (j < i) break;
                swap(ids, i, j);
            }
            ids[left + 1] = ids[j];
            ids[j] = temp;

            if (right - i + 1 >= j - left) {
                quicksort(ids, dists, i, right);
                quicksort(ids, dists, left, j - 1);
            } else {
                quicksort(ids, dists, left, j - 1);
                quicksort(ids, dists, i, right);
            }
        }
    }

    function swap(arr, i, j) {
        const tmp = arr[i];
        arr[i] = arr[j];
        arr[j] = tmp;
    }

    function defaultGetX(p) {
        return p[0];
    }
    function defaultGetY(p) {
        return p[1];
    }

    const epsilon = 1e-6;

    class Path {
      constructor() {
        this._x0 = this._y0 = // start of current subpath
        this._x1 = this._y1 = null; // end of current subpath
        this._ = "";
      }
      moveTo(x, y) {
        this._ += `M${this._x0 = this._x1 = +x},${this._y0 = this._y1 = +y}`;
      }
      closePath() {
        if (this._x1 !== null) {
          this._x1 = this._x0, this._y1 = this._y0;
          this._ += "Z";
        }
      }
      lineTo(x, y) {
        this._ += `L${this._x1 = +x},${this._y1 = +y}`;
      }
      arc(x, y, r) {
        x = +x, y = +y, r = +r;
        const x0 = x + r;
        const y0 = y;
        if (r < 0) throw new Error("negative radius");
        if (this._x1 === null) this._ += `M${x0},${y0}`;
        else if (Math.abs(this._x1 - x0) > epsilon || Math.abs(this._y1 - y0) > epsilon) this._ += "L" + x0 + "," + y0;
        if (!r) return;
        this._ += `A${r},${r},0,1,1,${x - r},${y}A${r},${r},0,1,1,${this._x1 = x0},${this._y1 = y0}`;
      }
      rect(x, y, w, h) {
        this._ += `M${this._x0 = this._x1 = +x},${this._y0 = this._y1 = +y}h${+w}v${+h}h${-w}Z`;
      }
      value() {
        return this._ || null;
      }
    }

    class Polygon {
      constructor() {
        this._ = [];
      }
      moveTo(x, y) {
        this._.push([x, y]);
      }
      closePath() {
        this._.push(this._[0].slice());
      }
      lineTo(x, y) {
        this._.push([x, y]);
      }
      value() {
        return this._.length ? this._ : null;
      }
    }

    class Voronoi {
      constructor(delaunay, [xmin, ymin, xmax, ymax] = [0, 0, 960, 500]) {
        if (!((xmax = +xmax) >= (xmin = +xmin)) || !((ymax = +ymax) >= (ymin = +ymin))) throw new Error("invalid bounds");
        this.delaunay = delaunay;
        this._circumcenters = new Float64Array(delaunay.points.length * 2);
        this.vectors = new Float64Array(delaunay.points.length * 2);
        this.xmax = xmax, this.xmin = xmin;
        this.ymax = ymax, this.ymin = ymin;
        this._init();
      }
      update() {
        this.delaunay.update();
        this._init();
        return this;
      }
      _init() {
        const {delaunay: {points, hull, triangles}, vectors} = this;

        // Compute circumcenters.
        const circumcenters = this.circumcenters = this._circumcenters.subarray(0, triangles.length / 3 * 2);
        for (let i = 0, j = 0, n = triangles.length, x, y; i < n; i += 3, j += 2) {
          const t1 = triangles[i] * 2;
          const t2 = triangles[i + 1] * 2;
          const t3 = triangles[i + 2] * 2;
          const x1 = points[t1];
          const y1 = points[t1 + 1];
          const x2 = points[t2];
          const y2 = points[t2 + 1];
          const x3 = points[t3];
          const y3 = points[t3 + 1];

          const dx = x2 - x1;
          const dy = y2 - y1;
          const ex = x3 - x1;
          const ey = y3 - y1;
          const ab = (dx * ey - dy * ex) * 2;

          if (Math.abs(ab) < 1e-9) {
            // degenerate case (collinear diagram)
            // almost equal points (degenerate triangle)
            // the circumcenter is at the infinity, in a
            // direction that is:
            // 1. orthogonal to the halfedge.
            let a = 1e9;
            // 2. points away from the center; since the list of triangles starts
            // in the center, the first point of the first triangle
            // will be our reference
            const r = triangles[0] * 2;
            a *= Math.sign((points[r] - x1) * ey - (points[r + 1] - y1) * ex);
            x = (x1 + x3) / 2 - a * ey;
            y = (y1 + y3) / 2 + a * ex;
          } else {
            const d = 1 / ab;
            const bl = dx * dx + dy * dy;
            const cl = ex * ex + ey * ey;
            x = x1 + (ey * bl - dy * cl) * d;
            y = y1 + (dx * cl - ex * bl) * d;
          }
          circumcenters[j] = x;
          circumcenters[j + 1] = y;
        }

        // Compute exterior cell rays.
        let h = hull[hull.length - 1];
        let p0, p1 = h * 4;
        let x0, x1 = points[2 * h];
        let y0, y1 = points[2 * h + 1];
        vectors.fill(0);
        for (let i = 0; i < hull.length; ++i) {
          h = hull[i];
          p0 = p1, x0 = x1, y0 = y1;
          p1 = h * 4, x1 = points[2 * h], y1 = points[2 * h + 1];
          vectors[p0 + 2] = vectors[p1] = y0 - y1;
          vectors[p0 + 3] = vectors[p1 + 1] = x1 - x0;
        }
      }
      render(context) {
        const buffer = context == null ? context = new Path : undefined;
        const {delaunay: {halfedges, inedges, hull}, circumcenters, vectors} = this;
        if (hull.length <= 1) return null;
        for (let i = 0, n = halfedges.length; i < n; ++i) {
          const j = halfedges[i];
          if (j < i) continue;
          const ti = Math.floor(i / 3) * 2;
          const tj = Math.floor(j / 3) * 2;
          const xi = circumcenters[ti];
          const yi = circumcenters[ti + 1];
          const xj = circumcenters[tj];
          const yj = circumcenters[tj + 1];
          this._renderSegment(xi, yi, xj, yj, context);
        }
        let h0, h1 = hull[hull.length - 1];
        for (let i = 0; i < hull.length; ++i) {
          h0 = h1, h1 = hull[i];
          const t = Math.floor(inedges[h1] / 3) * 2;
          const x = circumcenters[t];
          const y = circumcenters[t + 1];
          const v = h0 * 4;
          const p = this._project(x, y, vectors[v + 2], vectors[v + 3]);
          if (p) this._renderSegment(x, y, p[0], p[1], context);
        }
        return buffer && buffer.value();
      }
      renderBounds(context) {
        const buffer = context == null ? context = new Path : undefined;
        context.rect(this.xmin, this.ymin, this.xmax - this.xmin, this.ymax - this.ymin);
        return buffer && buffer.value();
      }
      renderCell(i, context) {
        const buffer = context == null ? context = new Path : undefined;
        const points = this._clip(i);
        if (points === null || !points.length) return;
        context.moveTo(points[0], points[1]);
        let n = points.length;
        while (points[0] === points[n-2] && points[1] === points[n-1] && n > 1) n -= 2;
        for (let i = 2; i < n; i += 2) {
          if (points[i] !== points[i-2] || points[i+1] !== points[i-1])
            context.lineTo(points[i], points[i + 1]);
        }
        context.closePath();
        return buffer && buffer.value();
      }
      *cellPolygons() {
        const {delaunay: {points}} = this;
        for (let i = 0, n = points.length / 2; i < n; ++i) {
          const cell = this.cellPolygon(i);
          if (cell) cell.index = i, yield cell;
        }
      }
      cellPolygon(i) {
        const polygon = new Polygon;
        this.renderCell(i, polygon);
        return polygon.value();
      }
      _renderSegment(x0, y0, x1, y1, context) {
        let S;
        const c0 = this._regioncode(x0, y0);
        const c1 = this._regioncode(x1, y1);
        if (c0 === 0 && c1 === 0) {
          context.moveTo(x0, y0);
          context.lineTo(x1, y1);
        } else if (S = this._clipSegment(x0, y0, x1, y1, c0, c1)) {
          context.moveTo(S[0], S[1]);
          context.lineTo(S[2], S[3]);
        }
      }
      contains(i, x, y) {
        if ((x = +x, x !== x) || (y = +y, y !== y)) return false;
        return this.delaunay._step(i, x, y) === i;
      }
      *neighbors(i) {
        const ci = this._clip(i);
        if (ci) for (const j of this.delaunay.neighbors(i)) {
          const cj = this._clip(j);
          // find the common edge
          if (cj) loop: for (let ai = 0, li = ci.length; ai < li; ai += 2) {
            for (let aj = 0, lj = cj.length; aj < lj; aj += 2) {
              if (ci[ai] == cj[aj]
              && ci[ai + 1] == cj[aj + 1]
              && ci[(ai + 2) % li] == cj[(aj + lj - 2) % lj]
              && ci[(ai + 3) % li] == cj[(aj + lj - 1) % lj]
              ) {
                yield j;
                break loop;
              }
            }
          }
        }
      }
      _cell(i) {
        const {circumcenters, delaunay: {inedges, halfedges, triangles}} = this;
        const e0 = inedges[i];
        if (e0 === -1) return null; // coincident point
        const points = [];
        let e = e0;
        do {
          const t = Math.floor(e / 3);
          points.push(circumcenters[t * 2], circumcenters[t * 2 + 1]);
          e = e % 3 === 2 ? e - 2 : e + 1;
          if (triangles[e] !== i) break; // bad triangulation
          e = halfedges[e];
        } while (e !== e0 && e !== -1);
        return points;
      }
      _clip(i) {
        // degenerate case (1 valid point: return the box)
        if (i === 0 && this.delaunay.hull.length === 1) {
          return [this.xmax, this.ymin, this.xmax, this.ymax, this.xmin, this.ymax, this.xmin, this.ymin];
        }
        const points = this._cell(i);
        if (points === null) return null;
        const {vectors: V} = this;
        const v = i * 4;
        return V[v] || V[v + 1]
            ? this._clipInfinite(i, points, V[v], V[v + 1], V[v + 2], V[v + 3])
            : this._clipFinite(i, points);
      }
      _clipFinite(i, points) {
        const n = points.length;
        let P = null;
        let x0, y0, x1 = points[n - 2], y1 = points[n - 1];
        let c0, c1 = this._regioncode(x1, y1);
        let e0, e1 = 0;
        for (let j = 0; j < n; j += 2) {
          x0 = x1, y0 = y1, x1 = points[j], y1 = points[j + 1];
          c0 = c1, c1 = this._regioncode(x1, y1);
          if (c0 === 0 && c1 === 0) {
            e0 = e1, e1 = 0;
            if (P) P.push(x1, y1);
            else P = [x1, y1];
          } else {
            let S, sx0, sy0, sx1, sy1;
            if (c0 === 0) {
              if ((S = this._clipSegment(x0, y0, x1, y1, c0, c1)) === null) continue;
              [sx0, sy0, sx1, sy1] = S;
            } else {
              if ((S = this._clipSegment(x1, y1, x0, y0, c1, c0)) === null) continue;
              [sx1, sy1, sx0, sy0] = S;
              e0 = e1, e1 = this._edgecode(sx0, sy0);
              if (e0 && e1) this._edge(i, e0, e1, P, P.length);
              if (P) P.push(sx0, sy0);
              else P = [sx0, sy0];
            }
            e0 = e1, e1 = this._edgecode(sx1, sy1);
            if (e0 && e1) this._edge(i, e0, e1, P, P.length);
            if (P) P.push(sx1, sy1);
            else P = [sx1, sy1];
          }
        }
        if (P) {
          e0 = e1, e1 = this._edgecode(P[0], P[1]);
          if (e0 && e1) this._edge(i, e0, e1, P, P.length);
        } else if (this.contains(i, (this.xmin + this.xmax) / 2, (this.ymin + this.ymax) / 2)) {
          return [this.xmax, this.ymin, this.xmax, this.ymax, this.xmin, this.ymax, this.xmin, this.ymin];
        }
        return P;
      }
      _clipSegment(x0, y0, x1, y1, c0, c1) {
        while (true) {
          if (c0 === 0 && c1 === 0) return [x0, y0, x1, y1];
          if (c0 & c1) return null;
          let x, y, c = c0 || c1;
          if (c & 0b1000) x = x0 + (x1 - x0) * (this.ymax - y0) / (y1 - y0), y = this.ymax;
          else if (c & 0b0100) x = x0 + (x1 - x0) * (this.ymin - y0) / (y1 - y0), y = this.ymin;
          else if (c & 0b0010) y = y0 + (y1 - y0) * (this.xmax - x0) / (x1 - x0), x = this.xmax;
          else y = y0 + (y1 - y0) * (this.xmin - x0) / (x1 - x0), x = this.xmin;
          if (c0) x0 = x, y0 = y, c0 = this._regioncode(x0, y0);
          else x1 = x, y1 = y, c1 = this._regioncode(x1, y1);
        }
      }
      _clipInfinite(i, points, vx0, vy0, vxn, vyn) {
        let P = Array.from(points), p;
        if (p = this._project(P[0], P[1], vx0, vy0)) P.unshift(p[0], p[1]);
        if (p = this._project(P[P.length - 2], P[P.length - 1], vxn, vyn)) P.push(p[0], p[1]);
        if (P = this._clipFinite(i, P)) {
          for (let j = 0, n = P.length, c0, c1 = this._edgecode(P[n - 2], P[n - 1]); j < n; j += 2) {
            c0 = c1, c1 = this._edgecode(P[j], P[j + 1]);
            if (c0 && c1) j = this._edge(i, c0, c1, P, j), n = P.length;
          }
        } else if (this.contains(i, (this.xmin + this.xmax) / 2, (this.ymin + this.ymax) / 2)) {
          P = [this.xmin, this.ymin, this.xmax, this.ymin, this.xmax, this.ymax, this.xmin, this.ymax];
        }
        return P;
      }
      _edge(i, e0, e1, P, j) {
        while (e0 !== e1) {
          let x, y;
          switch (e0) {
            case 0b0101: e0 = 0b0100; continue; // top-left
            case 0b0100: e0 = 0b0110, x = this.xmax, y = this.ymin; break; // top
            case 0b0110: e0 = 0b0010; continue; // top-right
            case 0b0010: e0 = 0b1010, x = this.xmax, y = this.ymax; break; // right
            case 0b1010: e0 = 0b1000; continue; // bottom-right
            case 0b1000: e0 = 0b1001, x = this.xmin, y = this.ymax; break; // bottom
            case 0b1001: e0 = 0b0001; continue; // bottom-left
            case 0b0001: e0 = 0b0101, x = this.xmin, y = this.ymin; break; // left
          }
          // Note: this implicitly checks for out of bounds: if P[j] or P[j+1] are
          // undefined, the conditional statement will be executed.
          if ((P[j] !== x || P[j + 1] !== y) && this.contains(i, x, y)) {
            P.splice(j, 0, x, y), j += 2;
          }
        }
        if (P.length > 4) {
          for (let i = 0; i < P.length; i+= 2) {
            const j = (i + 2) % P.length, k = (i + 4) % P.length;
            if (P[i] === P[j] && P[j] === P[k]
            || P[i + 1] === P[j + 1] && P[j + 1] === P[k + 1])
              P.splice(j, 2), i -= 2;
          }
        }
        return j;
      }
      _project(x0, y0, vx, vy) {
        let t = Infinity, c, x, y;
        if (vy < 0) { // top
          if (y0 <= this.ymin) return null;
          if ((c = (this.ymin - y0) / vy) < t) y = this.ymin, x = x0 + (t = c) * vx;
        } else if (vy > 0) { // bottom
          if (y0 >= this.ymax) return null;
          if ((c = (this.ymax - y0) / vy) < t) y = this.ymax, x = x0 + (t = c) * vx;
        }
        if (vx > 0) { // right
          if (x0 >= this.xmax) return null;
          if ((c = (this.xmax - x0) / vx) < t) x = this.xmax, y = y0 + (t = c) * vy;
        } else if (vx < 0) { // left
          if (x0 <= this.xmin) return null;
          if ((c = (this.xmin - x0) / vx) < t) x = this.xmin, y = y0 + (t = c) * vy;
        }
        return [x, y];
      }
      _edgecode(x, y) {
        return (x === this.xmin ? 0b0001
            : x === this.xmax ? 0b0010 : 0b0000)
            | (y === this.ymin ? 0b0100
            : y === this.ymax ? 0b1000 : 0b0000);
      }
      _regioncode(x, y) {
        return (x < this.xmin ? 0b0001
            : x > this.xmax ? 0b0010 : 0b0000)
            | (y < this.ymin ? 0b0100
            : y > this.ymax ? 0b1000 : 0b0000);
      }
    }

    const tau = 2 * Math.PI, pow = Math.pow;

    function pointX(p) {
      return p[0];
    }

    function pointY(p) {
      return p[1];
    }

    // A triangulation is collinear if all its triangles have a non-null area
    function collinear(d) {
      const {triangles, coords} = d;
      for (let i = 0; i < triangles.length; i += 3) {
        const a = 2 * triangles[i],
              b = 2 * triangles[i + 1],
              c = 2 * triangles[i + 2],
              cross = (coords[c] - coords[a]) * (coords[b + 1] - coords[a + 1])
                    - (coords[b] - coords[a]) * (coords[c + 1] - coords[a + 1]);
        if (cross > 1e-10) return false;
      }
      return true;
    }

    function jitter(x, y, r) {
      return [x + Math.sin(x + y) * r, y + Math.cos(x - y) * r];
    }

    class Delaunay {
      static from(points, fx = pointX, fy = pointY, that) {
        return new Delaunay("length" in points
            ? flatArray(points, fx, fy, that)
            : Float64Array.from(flatIterable(points, fx, fy, that)));
      }
      constructor(points) {
        this._delaunator = new Delaunator(points);
        this.inedges = new Int32Array(points.length / 2);
        this._hullIndex = new Int32Array(points.length / 2);
        this.points = this._delaunator.coords;
        this._init();
      }
      update() {
        this._delaunator.update();
        this._init();
        return this;
      }
      _init() {
        const d = this._delaunator, points = this.points;

        // check for collinear
        if (d.hull && d.hull.length > 2 && collinear(d)) {
          this.collinear = Int32Array.from({length: points.length/2}, (_,i) => i)
            .sort((i, j) => points[2 * i] - points[2 * j] || points[2 * i + 1] - points[2 * j + 1]); // for exact neighbors
          const e = this.collinear[0], f = this.collinear[this.collinear.length - 1],
            bounds = [ points[2 * e], points[2 * e + 1], points[2 * f], points[2 * f + 1] ],
            r = 1e-8 * Math.hypot(bounds[3] - bounds[1], bounds[2] - bounds[0]);
          for (let i = 0, n = points.length / 2; i < n; ++i) {
            const p = jitter(points[2 * i], points[2 * i + 1], r);
            points[2 * i] = p[0];
            points[2 * i + 1] = p[1];
          }
          this._delaunator = new Delaunator(points);
        } else {
          delete this.collinear;
        }

        const halfedges = this.halfedges = this._delaunator.halfedges;
        const hull = this.hull = this._delaunator.hull;
        const triangles = this.triangles = this._delaunator.triangles;
        const inedges = this.inedges.fill(-1);
        const hullIndex = this._hullIndex.fill(-1);

        // Compute an index from each point to an (arbitrary) incoming halfedge
        // Used to give the first neighbor of each point; for this reason,
        // on the hull we give priority to exterior halfedges
        for (let e = 0, n = halfedges.length; e < n; ++e) {
          const p = triangles[e % 3 === 2 ? e - 2 : e + 1];
          if (halfedges[e] === -1 || inedges[p] === -1) inedges[p] = e;
        }
        for (let i = 0, n = hull.length; i < n; ++i) {
          hullIndex[hull[i]] = i;
        }

        // degenerate case: 1 or 2 (distinct) points
        if (hull.length <= 2 && hull.length > 0) {
          this.triangles = new Int32Array(3).fill(-1);
          this.halfedges = new Int32Array(3).fill(-1);
          this.triangles[0] = hull[0];
          inedges[hull[0]] = 1;
          if (hull.length === 2) {
            inedges[hull[1]] = 0;
            this.triangles[1] = hull[1];
            this.triangles[2] = hull[1];
          }
        }
      }
      voronoi(bounds) {
        return new Voronoi(this, bounds);
      }
      *neighbors(i) {
        const {inedges, hull, _hullIndex, halfedges, triangles, collinear} = this;

        // degenerate case with several collinear points
        if (collinear) {
          const l = collinear.indexOf(i);
          if (l > 0) yield collinear[l - 1];
          if (l < collinear.length - 1) yield collinear[l + 1];
          return;
        }

        const e0 = inedges[i];
        if (e0 === -1) return; // coincident point
        let e = e0, p0 = -1;
        do {
          yield p0 = triangles[e];
          e = e % 3 === 2 ? e - 2 : e + 1;
          if (triangles[e] !== i) return; // bad triangulation
          e = halfedges[e];
          if (e === -1) {
            const p = hull[(_hullIndex[i] + 1) % hull.length];
            if (p !== p0) yield p;
            return;
          }
        } while (e !== e0);
      }
      find(x, y, i = 0) {
        if ((x = +x, x !== x) || (y = +y, y !== y)) return -1;
        const i0 = i;
        let c;
        while ((c = this._step(i, x, y)) >= 0 && c !== i && c !== i0) i = c;
        return c;
      }
      _step(i, x, y) {
        const {inedges, hull, _hullIndex, halfedges, triangles, points} = this;
        if (inedges[i] === -1 || !points.length) return (i + 1) % (points.length >> 1);
        let c = i;
        let dc = pow(x - points[i * 2], 2) + pow(y - points[i * 2 + 1], 2);
        const e0 = inedges[i];
        let e = e0;
        do {
          let t = triangles[e];
          const dt = pow(x - points[t * 2], 2) + pow(y - points[t * 2 + 1], 2);
          if (dt < dc) dc = dt, c = t;
          e = e % 3 === 2 ? e - 2 : e + 1;
          if (triangles[e] !== i) break; // bad triangulation
          e = halfedges[e];
          if (e === -1) {
            e = hull[(_hullIndex[i] + 1) % hull.length];
            if (e !== t) {
              if (pow(x - points[e * 2], 2) + pow(y - points[e * 2 + 1], 2) < dc) return e;
            }
            break;
          }
        } while (e !== e0);
        return c;
      }
      render(context) {
        const buffer = context == null ? context = new Path : undefined;
        const {points, halfedges, triangles} = this;
        for (let i = 0, n = halfedges.length; i < n; ++i) {
          const j = halfedges[i];
          if (j < i) continue;
          const ti = triangles[i] * 2;
          const tj = triangles[j] * 2;
          context.moveTo(points[ti], points[ti + 1]);
          context.lineTo(points[tj], points[tj + 1]);
        }
        this.renderHull(context);
        return buffer && buffer.value();
      }
      renderPoints(context, r) {
        if (r === undefined && (!context || typeof context.moveTo !== "function")) r = context, context = null;
        r = r == undefined ? 2 : +r;
        const buffer = context == null ? context = new Path : undefined;
        const {points} = this;
        for (let i = 0, n = points.length; i < n; i += 2) {
          const x = points[i], y = points[i + 1];
          context.moveTo(x + r, y);
          context.arc(x, y, r, 0, tau);
        }
        return buffer && buffer.value();
      }
      renderHull(context) {
        const buffer = context == null ? context = new Path : undefined;
        const {hull, points} = this;
        const h = hull[0] * 2, n = hull.length;
        context.moveTo(points[h], points[h + 1]);
        for (let i = 1; i < n; ++i) {
          const h = 2 * hull[i];
          context.lineTo(points[h], points[h + 1]);
        }
        context.closePath();
        return buffer && buffer.value();
      }
      hullPolygon() {
        const polygon = new Polygon;
        this.renderHull(polygon);
        return polygon.value();
      }
      renderTriangle(i, context) {
        const buffer = context == null ? context = new Path : undefined;
        const {points, triangles} = this;
        const t0 = triangles[i *= 3] * 2;
        const t1 = triangles[i + 1] * 2;
        const t2 = triangles[i + 2] * 2;
        context.moveTo(points[t0], points[t0 + 1]);
        context.lineTo(points[t1], points[t1 + 1]);
        context.lineTo(points[t2], points[t2 + 1]);
        context.closePath();
        return buffer && buffer.value();
      }
      *trianglePolygons() {
        const {triangles} = this;
        for (let i = 0, n = triangles.length / 3; i < n; ++i) {
          yield this.trianglePolygon(i);
        }
      }
      trianglePolygon(i) {
        const polygon = new Polygon;
        this.renderTriangle(i, polygon);
        return polygon.value();
      }
    }

    function flatArray(points, fx, fy, that) {
      const n = points.length;
      const array = new Float64Array(n * 2);
      for (let i = 0; i < n; ++i) {
        const p = points[i];
        array[i * 2] = fx.call(that, p, i, points);
        array[i * 2 + 1] = fy.call(that, p, i, points);
      }
      return array;
    }

    function* flatIterable(points, fx, fy, that) {
      let i = 0;
      for (const p of points) {
        yield fx.call(that, p, i, points);
        yield fy.call(that, p, i, points);
        ++i;
      }
    }

    // gamma puts a point outside the boundary back inside,
    // 
    function gamma$1([px, py], {x1, x2, y1, y2}) {
        let [nx, ny] = [px, py];
        if (px < x1) nx = x1;
        else if (px > x2) nx = x2;
        if (py < y1) ny = y1;
        else if (py > y2) ny = y2;
        return [nx, ny]
    }

    function distance$1([ax, ay], [bx, by]) {
        return Math.hypot(bx - ax, by - ay);//Math.sqrt(((ax - bx) ** 2) + ((ay - by) ** 2));
    }

    function compute_overlap([px, py], [qx, qy], s) {
        const dx = Math.abs(px - qx);
        const dy = Math.abs(py - qy);
        if (dx < s || dy < s) {
            const ox = Math.abs(Math.min(dx - s, 0));
            const oy = Math.abs(Math.min(dy - s, 0));
            const dd = distance$1([px, py], [qx, qy]);
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
                const odistance = distance$1(po, qo);
                const overlap = compute_overlap(pa, qa, s);
                return {
                    "source": i,
                    "target": j,
                    "distance": odistance, //l_ij
                    "overlap": overlap, // delta_ij
                    "delta": odistance + overlap // d_ij
                }
            }).filter(({target, overlap}) => {
                const keep_overlap = overlap > 1e-2;
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

    async function* gridify_cmds(D, parameters) {
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
                    const dist = distance$1(pi, pj);
                    nominator[i][0] += (wij * (pj[0] + delta * ((pi[0] - pj[0]) / dist)));
                    nominator[i][1] += (wij * (pj[1] + delta * ((pi[1] - pj[1]) / dist)));
                    denominator[i] += wij;
                });

                G = G.map((p, i) => {
                    const [gamma_px, gamma_py] = gamma$1(p, Gamma);
                    return [
                        (nominator[i][0] + alpha * gamma_px) / (denominator[i] + alpha),
                        (nominator[i][1] + alpha * gamma_py) / (denominator[i] + alpha),
                    ]
                });

            }
            // reconstruct the proximity graph
            // with the new layout
            proximity_graph = create_proximity_graph(G, D, s, true);

            yield [G, proximity_graph];

        } while (--max_iter > 0 && !end_condition(proximity_graph, G, Gamma))
        // round to grid
        G = G.map(d => {
            return d.map(v => {
                return Math.round((v - s / 2) / s) * s + (s / 2)
            })
        });

        proximity_graph = create_proximity_graph(G, s);
        //console.log("fin")
        yield [G, proximity_graph];
        //return G
    }

    function area([x_min, x_max, y_min, y_max]) {
        return (x_max - x_min) * (y_max - y_min);
    }

    function extent(indices, data) {
        return [
            ...d3.extent(indices, (i) => data[i][0]),
            ...d3.extent(indices, (i) => data[i][1]),
        ];
    }

    function quadtree(data, node, [x_min, x_max, y_min, y_max], nodes) {
        if (Math.abs(x_max - x_min) <= 1 || Math.abs(y_max - y_min) <= 1) {
            return;
        }

        let x_mid = Math.min(
            Math.max(x_min, Math.floor((x_max + x_min) / 2)),
            x_max - 1
        );
        let y_mid = Math.min(
            Math.max(y_min, Math.floor((y_max + y_min) / 2)),
            y_max - 1
        );

        node.children = [];
        //B1
        let child;
        child = node.filter(
            i =>
                data[i][0] >= x_min &&
                data[i][0] <= x_mid &&
                data[i][1] >= y_min &&
                data[i][1] <= y_mid
        );
        child.E = [x_min, x_mid, y_min, y_mid];
        node.children.push(child);

        // B2
        child = node.filter(
            i =>
                data[i][0] > x_mid &&
                data[i][0] <= x_max &&
                data[i][1] >= y_min &&
                data[i][1] <= y_mid
        );
        child.E = [x_mid, x_max, y_min, y_mid];
        node.children.push(child);

        //B3
        child = node.filter(
            i =>
                data[i][0] >= x_min &&
                data[i][0] <= x_mid &&
                data[i][1] > y_mid &&
                data[i][1] <= y_max
        );
        child.E = [x_min, x_mid, y_mid, y_max];
        node.children.push(child);

        //B4
        child = node.filter(
            i =>
                data[i][0] > x_mid &&
                data[i][0] <= x_max &&
                data[i][1] > y_mid &&
                data[i][1] <= y_max
        );
        child.E = [x_mid, x_max, y_mid, y_max];
        node.children.push(child);

        nodes.push(...node.children.filter(child => child.length > 0));

        node.children.forEach((child, i) => {
            child.parent = node;
            //console.log("quad", node.E, i, child.E)
            if (child.length > 0) quadtree(data, child, child.E, nodes);
        });
    }

    function fit(node, count, data) {
        const [x_min, x_max, y_min, y_max] = node.E;
        // if node is leaf, then place the points on the grid
        if (!node.children) {
            if (node.length == 0) return;
            const x_span = x_max - x_min;
            const y_span = y_max - y_min;

            for (let dx = 0; dx < x_span; ++dx) {
                for (let dy = 0; dy < y_span; ++dy) {
                    if (node.length > 0) {
                        const i = node.pop();
                        const el = [x_min + dx, y_min + dy];
                        el.i = i;
                        nodes.final.push(el);
                    }
                }
            }
            if (node.length > 0) {
                nodes.leftovers.push(...node);
            }
            return;
        }
        let x_mid = Math.floor((x_max + x_min) / 2);
        let y_mid = Math.floor((y_max + y_min) / 2);

        const [B1, B2, B3, B4] = node.children;
        const [P1, P2, P3, P4] = node.children.map((child) => child.length);

        // Search for x
        const P13 = P1 + P3;
        const P24 = P2 + P4;
        let A13 = area([x_min, x_mid, y_min, y_max]);
        let A24 = area([x_mid, x_max, y_min, y_max]);
        if (P13 > A13 && P24 > A24) {
            console.log(P13, A13, P24, A24);
            throw "Not enough space, should not happen, X";
        }
        let direction = P13 <= A13 ? -1 : P24 <= A24 ? 1 : 0;
        while (!(P13 <= A13 && P24 <= A24)) {
            ++count;
            x_mid += direction;
            if (x_min > x_mid || x_mid > x_max) {
                console.log(x_min, x_mid, x_max);
                return;
            }
            A13 = area([x_min, x_mid, y_min, y_max]);
            A24 = area([x_mid, x_max, y_min, y_max]);
            let new_direction = P13 <= A13 ? -1 : P24 <= A24 ? 1 : 0;
            if (!(P13 <= A13 && P24 <= A24) && new_direction == -direction) {
                console.log("split the line x!", P13, A13, P24, A24);
            }
            console.log("xs", x_min, x_mid, x_max, P13, A13, P24, A24);
        }

        // search fro y_1
        let y_1 = y_mid;
        let A1 = area([x_min, x_mid, y_min, y_mid]);
        let A3 = area([x_min, x_mid, y_mid, y_max]);
        if (P1 + P3 > A1 + A3) {
            console.log(P1, A1, P3, A3);
            console.log("xs", x_min, x_mid, x_max);
            console.log("ys", y_min, y_mid, y_max, B1.E[3], B2.E[3], B3.E[3], B4.E[3]);
            console.log(node.map(i => data[i]));
            console.log(node, B1.E, B3.E, y_max);
            throw "Not enough space, should not happen, Y_1";
        }
        direction = P1 <= A1 ? -1 : P3 <= A3 ? 1 : 0;
        while (!(P1 <= A1 && P3 <= A3)) {
            //} && y_1 < y_max && y_1 > y_min) {
            //if (direction == 0) throw "y1 dir 0";
            //++count;
            y_1 += direction;
            if (y_min > y_1 || y_1 > y_max) {
                console.log(y_min, y_1, y_max);
                //nodes.leftovers.push(...node);
                return;
            }
            A1 = area([x_min, x_mid, y_min, y_1]);
            A3 = area([x_min, x_mid, y_1, y_max]);
            let new_direction = P1 <= A1 ? -1 : P3 <= A3 ? 1 : 0;
            if (!(P1 <= A1 && P3 <= A3) && new_direction == -direction) {
                console.log("split the line! y_1", P1, A1, P3, A3);
            }
        }

        // search for y_2
        let y_2 = y_mid;
        let A2 = area([x_mid, x_max, y_min, y_mid]);
        let A4 = area([x_mid, x_max, y_mid, y_max]);
        if (P2 > A2 && P4 > A4) {
            console.log(P2, A2, P4, A4);
            throw "Not enough space, should not happen, Y_2";
        }
        direction = P2 <= A2 ? -1 : P4 <= A4 ? 1 : 0;
        while (!(P2 <= A2 && P4 <= A4)) {
            //} && y_2 < y_max && y_2 > y_min) {
            //if (direction == 0) throw "y2 dir 0";
            ++count;
            y_2 += direction;
            if (y_2 < y_min || y_2 > y_max) {
                console.log(y_min, y_2, y_max);
                //nodes.leftovers.push(...node);
                return;
            }
            A2 = area([x_mid, x_max, y_min, y_2]);
            A4 = area([x_mid, x_max, y_2, y_max]);
            let new_direction = P2 <= A2 ? -1 : P4 <= A4 ? 1 : 0;
            if (!(P2 <= A2 && P4 <= A4) && new_direction == -direction) {
                console.log("split the line! y_2", P2, A2, P4, A4);
            }
        }

        // set new borders
        B1.E[1] = x_mid;
        B3.E[1] = x_mid;
        B2.E[0] = x_mid;
        B4.E[0] = x_mid;
        B1.E[3] = y_1;
        B3.E[2] = y_1;
        B2.E[3] = y_2;
        B4.E[2] = y_2;
        node.children.forEach((child) => fit(child, count, data));
    }

    let nodes = [];
    function gridify_gridfit(D) {
        const N = D.length;
        nodes = [];
        nodes.final = [];
        nodes.leftovers = [];
        let count = 0;
        // do the gridfit;

        // create quadtree;
        const root = d3.range(0, N);
        root.E = extent(root, D);
        root.P = N;
        quadtree(D, root, root.E, nodes);

        // fit to grid;
        try {
            fit(root, count, D);
        } catch (e) {
            console.error(e);
        } finally {
            console.log("count:", count);
            return nodes;
        }
    }

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

    // transform to node gosper pattern
    const transform_index = [4, 0, 1, 2, 3, 6, 5];

    // precomputed array of circles inscribed into Gosper islands of different levels
    const inscribed_circles = [0.755928946000, 0.755928946000, 0.750121467308, 0.746782631146, 0.746782631146, 0.746577727521, 0.746348363909, 0.746348363909, 0.746344578768, 0.746327538283, 0.746327538283, 0.746327538283, 0.746326555879, 0.746326555879, 0.746326555879, 0.746326510616, 0.746326510616, 0.746326510616, 0.746326508597, 0.746326508597, 0.746326508597];

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

    function inverse_transform([x, y, z]) {
        const nx = (3 * x + z);
        const nz = (-x + 2 * z);
        return [nx, -nx - nz, nz];
    }

    function get_size(D, level) {
        let min_x = Infinity;
        let max_x = -Infinity;
        let min_y = Infinity;
        let max_y = -Infinity;
        D.forEach(([x, y]) => {
            min_x = x < min_x ? x : min_x;
            max_x = x > max_x ? x : max_x;
            min_y = y < min_y ? y : min_y;
            max_y = y > max_y ? y : max_y;
        });

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
    function gosper_encode([xp, yp], level, size) {
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
    function gosper_decode(t, level, size) {
        //console.log("decode", t, level, size)
        /* const k = new Uint8Array(level - 1).fill(0);
        for (let i = level - 1; i >= 0; --i) {
            const mod_t = t % 7;
            k[i] = mod_t;
            t = (t - mod_t) / 7;
        } */
        const k = t;
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

    function gridify_gosper(D, { level: level }) {
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
                P.set(p, i);
            }
        });

        for (const [p, i] of P.entries()) {
            console.log("i", i, "p", p);
            Y[i] = gosper_decode(p, level, size);
        }
        console.log(P);
        return Y;
    }

    //import { gridify_gridfit } from "./gridfit.js";

    function gridify(data, method = "hilbert", parameters = {}) {
        if (method == "hilbert") {
            if (!Object.keys(parameters).includes("pluslevel")) {
                parameters.pluslevel = 0;
            }
            if (!("whitespace" in parameters)) {
                parameters.whitespace = 1;
            }
            const level = Math.ceil(Math.log2(data.length * parameters.whitespace) / Math.log2(4)) + parameters.pluslevel;
            const size = Math.pow(2, level);
            const { x: ux, y: uy, width: w, height: h } = get_bounds(data);
            const original_bounding_box = [ux, uy, w, h];
            const gridded_bounding_box = [0, 0, size, size];
            const D = data.map((d) => remap(d, original_bounding_box, gridded_bounding_box));
            const start = performance.now();
            const res = gridify_hilbert(D, { level: level });
            const end = performance.now();
            res.runtime = end - start;
            return res;
        } else if (method ===  "gosper") {
            const level = Math.ceil(Math.log2(data.length) / Math.log2(7) + 1) + parameters.pluslevel;
            const start = performance.now();
            const res = gridify_gosper(data, { level: level });
            const end = performance.now();
            res.runtime = end - start;
            return res;
        } else if (method === "dgrid") {
            if (parameters && Object.keys(parameters).length == 0) {
                parameters = {};
                parameters.aspect_ratio = 1;
            }
            const start = performance.now();
            const res = gridify_dgrid(data, parameters);
            const end = performance.now();
            res.runtime = end - start;
            return res;
        } else if (method === "cmds") {
            parameters.alpha = parameters.alpha ? parameters.alpha : 0.1;
            parameters.Gamma = parameters.Gamma ? parameters.Gamma : (() => {
                const { x: min_x, y: min_y, width: wx, height: hy } = get_bounds(data);
                const w = Math.max(wx, hy);
                return {
                    x1: min_x,
                    x2: min_x + w,
                    y1: min_y,
                    y2: min_y + w,
                }
            })();
            parameters.size = parameters.size ? parameters.size : (() => {
                const { x1, x2, y1, y2 } = parameters.Gamma;
                const w = Math.max(x2 - x1, y2 - y1);
                const area = w ** 2 / ((Math.sqrt(data.length) + 1) ** 2);
                return Math.sqrt(area);
            })();
            const start = performance.now();
            const res = gridify_cmds(data, parameters);
            const end = performance.now();
            res.runtime = end - start;
            return res;
        } else if (method === "nmap") {
            const start = performance.now();
            const res =  gridify_nmap(data, parameters);
            const end = performance.now();
            res.runtime = end - start;
            return res;
        } else {
            throw "not a valid method!";
        }
    }

    exports.gosper_curve = gosper_curve;
    exports.gosper_decode = gosper_decode$1;
    exports.gosper_encode = gosper_encode$1;
    exports.gridify = gridify;
    exports.gridify_cmds = gridify_cmds;
    exports.gridify_dgrid = gridify_dgrid;
    exports.gridify_gosper = gridify_gosper$1;
    exports.gridify_gridfit = gridify_gridfit;
    exports.gridify_hilbert = gridify_hilbert;
    exports.gridify_nmap = gridify_nmap;
    exports.hilbert_decode = hilbert_decode;
    exports.hilbert_encode = hilbert_encode;
    exports.utils = utils;
    exports.version = version;

    Object.defineProperty(exports, '__esModule', { value: true });

})));
