(function (global, factory) {
    typeof exports === 'object' && typeof module !== 'undefined' ? factory(exports) :
    typeof define === 'function' && define.amd ? define(['exports'], factory) :
    (global = global || self, factory(global.hagrid = {}));
}(this, (function (exports) { 'use strict';

    var version = "0.0.2";

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
        for (let s = size >> 1; s > 0; s >>= 1) {
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

    function gosper_encode([x, y], level) {
        const size = level.s; /// Math.sqrt(3);
        let cube = xy2cube([x, y], size);
        cube = cube_round(cube);
        let [nx, ny] = cube2xy(cube, size);
        return knn.search([nx, ny], 1).first.element.index;
    }

    function gosper_decode(i) {
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

    exports.gosper_decode = gosper_decode;
    exports.gosper_encode = gosper_encode;
    exports.hilbert_decode = hilbert_decode;
    exports.hilbert_encode = hilbert_encode;
    exports.version = version;

    Object.defineProperty(exports, '__esModule', { value: true });

})));
