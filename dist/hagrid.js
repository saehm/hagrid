(function (global, factory) {
    typeof exports === 'object' && typeof module !== 'undefined' ? factory(exports) :
    typeof define === 'function' && define.amd ? define(['exports'], factory) :
    (global = global || self, factory(global.hagrid = {}));
}(this, (function (exports) { 'use strict';

    var version = "0.1.0";

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

    function hilbert_encode([px, py], size) {
        let n = 0;
        for (let s = size / 2; s > 0; s /= 2) {
            const rx = (px & s) > 0;
            const ry = (py & s) > 0;
            n += (s ** 2) * ((3 * rx) ^ ry);
            [px, py] = rotate(size, px, py, rx, ry);
        }
        return n;
    }

    function hilbert_decode(n, size) {
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
            p += 1;
        }
        P.set(p, i);
    }

    async function gridify_hilbert(D, {level, collision}) {
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
                    hilbert_collision_2(P, p, i);
                }
            } else {
                P.set(p, i);
            }
        });

        for (const [p, i] of P) {
            Y[i] = hilbert_decode(p, size);
        }

        return Y;
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

    // Arry for mapping points onto the center pattern index
    const pattern_index = [[5, 1, 3],[2, 4, 6]];

    // transform to node gosper pattern
    const transform_index = [4, 0, 1, 2, 3, 6, 5];

    // precomputed array of circles inscribed into Gosper islands of different levels
    const inscribed_circles = [0.755928946000, 0.755928946000, 0.750121467308, 0.746782631146, 0.746782631146, 0.746577727521, 0.746348363909, 0.746348363909, 0.746344578768, 0.746327538283, 0.746327538283, 0.746327538283, 0.746326555879, 0.746326555879, 0.746326555879, 0.746326510616, 0.746326510616, 0.746326510616, 0.746326508597, 0.746326508597, 0.746326508597];

    function xy2cube([x, y], size) {
        const cx = (x * SQRT3_DIV_3 - y * 1/3) / size;
        const cz = y * 2/3 / size;
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
    function f_y(i) {
        //[-120, 0, 0, -120, 0, 120, 0].map(d => d / 180 * Math.PI)
        return [-2.0943951023931953, 0, 0, -2.0943951023931953, 0, 2.0943951023931953, 0][i];
    }

    function gamma(i) {
        //[0, 60, 120, 180, 0, -60, -120].map(d => -d / 180 * Math.PI)
        return [0, -1.0471975511965976, -2.0943951023931953, -3.141592653589793, 0, 1.0471975511965976, 2.0943951023931953][i]
    }

    function rotate$1([x, y], alpha) {
        const sin_alpha = Math.sin(alpha);
        const cos_alpha = Math.cos(alpha);
        return [cos_alpha * x - sin_alpha * y,
                sin_alpha * x + cos_alpha * y]
    }

    function scale([x, y], s) {
        return [x * s, y * s];
    }

    function gosper_encode([x, y], level, size) {
        let b = new Array(level); // final position array
        let sign; // sign of decimal residue;
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
        }

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
        return k;
    }

    function gosper_decode(k, level, size) {
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
            const dk = rotate$1(d, gamma(ki));
            if (ki != 4) {
                c = [c[0] + dk[0], c[1] + dk[1]];
            }
            d = rotate$1(d, ALPHA + f_y(ki));
            d = scale(d, REC_SQRT7);
            if (ki == 0 || ki == 4 || ki == 5) {
                ord = !ord;
            }
        }
        return c;
    }

    function gridify_gosper(D, {level: level}) {
        const size = get_size(D, level);
        const N = D.length;
        const P = new Map();
        const Y = new Array(N).fill(0);

        D.forEach((d, i) => {
            let p = gosper_encode(d, level, size);
            if (P.has(p)) {
                //let q = P.get(p);
                while (P.has(p)) {
                    p = p + 1;
                }
            }
            P.set(p, d);
            Y[i] = gosper_decode(p, level, size);
        });

        return Y;
    }

    let times = [];

    function dgrid(G, P, r, s, i=0, j=0) {
        const N = P.length;
        if (N == 0) return
        else if (N == 1) {
            times.push(performance.now());
            G.push({
                d: P[0],
                i: i,
                j: j,
            });
        } else {
            if (r > s) {
                const r_half = Math.ceil(r / 2);
                const [P1, P2] = split(P, 1, r_half * s); // split by y
                dgrid(G, P1, (r_half), s, i, j);
                dgrid(G, P2, (r - r_half), s, (i + r_half), j);
            } else {
                const s_half = Math.ceil(s / 2);
                const [P1, P2] = split(P, 0, s_half * r); // split by x
                dgrid(G, P1, r, s_half, i, j);
                dgrid(G, P2, r, (s - s_half), i, (j + s_half));
            }
        }
    }

    function split(P, dim, pos) {
        const sorted = P.sort((a, b) => a[dim] - b[dim]);
        return [sorted.slice(0, pos), sorted.slice(pos)];
    }

    async function gridify_dgrid(D, size) {
        console.log(size);
        const N = D.length;
        let rows;
        let cols;
        if ("rows" in size && "cols" in size) {
            rows = size.rows;
            cols = size.cols;
        } else if (size == "square") {
            rows = Math.ceil(Math.sqrt(N));
            cols = Math.ceil(Math.sqrt(N));
        } else if ("aspect_ratio" in size) {
            rows = Math.floor(Math.sqrt(N * size.aspect_ratio));
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
            let [Da, Db] = split$1(P, N_half);
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

    function getBounds(D) {
        let min_x = Infinity;
        let max_x = -Infinity;
        let min_y = Infinity;
        let max_y = -Infinity;

        D.forEach(([x, y, i]) => {
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

    function split$1(P, pos) {
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
        const BB2 = getBounds(D);
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

    async function gridify_nmap(D, parameters) {
        let G = [];
        console.log(parameters);
        if (!"BB" in parameters || parameters.BB == null) {
            console.log("no BB");
            parameters.BB = getBounds(D);
        }
        let added_index = D.map(([x, y], i) => [x, y, i]);
        if ("squared" in parameters && parameters.squared == true) {
            added_index = squared(added_index, parameters.BB);
        }
        nmap(G, added_index, parameters.BB, false);
        return G.sort((a, b) => a.d[2] - b.d[2]).filter(d => d.d[2] != undefined)//.map(g => [g.j, g.i])
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
            if (orient(i0x, i0y, i1x, i1y, i2x, i2y)) {
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
                while (q = hullNext[e], !orient(x, y, coords[2 * e], coords[2 * e + 1], coords[2 * q], coords[2 * q + 1])) {
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
                while (q = hullNext[n], orient(x, y, coords[2 * n], coords[2 * n + 1], coords[2 * q], coords[2 * q + 1])) {
                    t = this._addTriangle(n, i, q, hullTri[i], -1, hullTri[n]);
                    hullTri[i] = this._legalize(t + 2);
                    hullNext[n] = n; // mark as removed
                    hullSize--;
                    n = q;
                }

                // walk backward from the other side, adding more triangles and flipping
                if (e === start) {
                    while (q = hullPrev[e], orient(x, y, coords[2 * q], coords[2 * q + 1], coords[2 * e], coords[2 * e + 1])) {
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

    // return 2d orientation sign if we're confident in it through J. Shewchuk's error bound check
    function orientIfSure(px, py, rx, ry, qx, qy) {
        const l = (ry - py) * (qx - px);
        const r = (rx - px) * (qy - py);
        return Math.abs(l - r) >= 3.3306690738754716e-16 * Math.abs(l + r) ? l - r : 0;
    }

    // a more robust orientation test that's stable in a given triangle (to fix robustness issues)
    function orient(rx, ry, qx, qy, px, py) {
        const sign = orientIfSure(px, py, rx, ry, qx, qy) ||
        orientIfSure(rx, ry, qx, qy, px, py) ||
        orientIfSure(qx, qy, px, py, rx, ry);
        return sign < 0;
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
          const bl = dx * dx + dy * dy;
          const cl = ex * ex + ey * ey;
          const ab = (dx * ey - dy * ex) * 2;

          if (!ab) {
            // degenerate case (collinear diagram)
            x = (x1 + x3) / 2 - 1e8 * ey;
            y = (y1 + y3) / 2 + 1e8 * ex;
          }
          else if (Math.abs(ab) < 1e-8) {
            // almost equal points (degenerate triangle)
            x = (x1 + x3) / 2;
            y = (y1 + y3) / 2;
          } else {
            const d = 1 / ab;
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
        if (points === null) return;
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
          if (cell) yield cell;
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
        let e0, e1;
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

    const tau = 2 * Math.PI;

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
            r = 1e-8 * Math.sqrt((bounds[3] - bounds[1])**2 + (bounds[2] - bounds[0])**2);
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
          this.triangles[1] = hull[1];
          this.triangles[2] = hull[1];
          inedges[hull[0]] = 1;
          if (hull.length === 2) inedges[hull[1]] = 0;
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
        let dc = (x - points[i * 2]) ** 2 + (y - points[i * 2 + 1]) ** 2;
        const e0 = inedges[i];
        let e = e0;
        do {
          let t = triangles[e];
          const dt = (x - points[t * 2]) ** 2 + (y - points[t * 2 + 1]) ** 2;
          if (dt < dc) dc = dt, c = t;
          e = e % 3 === 2 ? e - 2 : e + 1;
          if (triangles[e] !== i) break; // bad triangulation
          e = halfedges[e];
          if (e === -1) {
            e = hull[(_hullIndex[i] + 1) % hull.length];
            if (e !== t) {
              if ((x - points[e * 2]) ** 2 + (y - points[e * 2 + 1]) ** 2 < dc) return e;
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
      renderPoints(context, r = 2) {
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
        if (px > x2) nx = x2;
        if (py < y1) ny = y1;
        if (py > y2) ny = y2;
        return [nx, ny]
    }

    function distance$1([ax, ay], [bx, by]) {
        return Math.sqrt(Math.pow(ax - bx, 2) + Math.pow(ay - by, 2));
    }

    function compute_overlap([px, py], [qx, qy], s) {
        const dx = Math.abs(px - qx);
        const dy = Math.abs(py - qy);
        if (dx < s || dy < s) {
            const ox = Math.abs(Math.min(dx - s, 0));
            const oy = Math.abs(Math.min(dy - s, 0));
            const dd = distance$1([px, py], [qx, qy]);
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
                const odistance = distance$1(pd, qd);
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

    async function gridify_cmds(D, parameters) {
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
            let nominator = new Array(N).fill(0);
            nominator = nominator.map(() => [0, 0]);
            let denominator = new Array(N).fill(0);

            proximity_graph.forEach(({source: i, target: j, distance, delta}) => {
                let pi = G[i];
                let pj = G[j];
                const wij = (1/Math.pow(delta, 2));
                nominator[i][0] += (wij * (pj[0] + delta * ((pi[0] - pj[0]) / distance)));
                nominator[i][1] += (wij * (pj[1] + delta * ((pi[1] - pj[1]) / distance)));
                denominator[i] += wij;
            });

            G = G.map(([px, py], i) => {
                const gamma_i = gamma$1([px, py], Gamma);
                return [
                    (nominator[i][0] + alpha * gamma_i[0]) / (denominator[i] + alpha),
                    (nominator[i][1] + alpha * gamma_i[1]) / (denominator[i] + alpha),
                ]
            });
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
        });
        return G
    }

    exports.get_size = get_size;
    exports.gosper_decode = gosper_decode;
    exports.gosper_encode = gosper_encode;
    exports.gridify_cmds = gridify_cmds;
    exports.gridify_dgrid = gridify_dgrid;
    exports.gridify_gosper = gridify_gosper;
    exports.gridify_hilbert = gridify_hilbert;
    exports.gridify_nmap = gridify_nmap;
    exports.hilbert_decode = hilbert_decode;
    exports.hilbert_encode = hilbert_encode;
    exports.version = version;

    Object.defineProperty(exports, '__esModule', { value: true });

})));
