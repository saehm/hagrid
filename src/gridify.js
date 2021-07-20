import { remap, get_bounds } from "./utils.js";
import { gridify_hilbert } from "./hilbert.js";
import { gridify_gosper } from "./gosper2.js";
import { gridify_dgrid } from "./dgrid.js";
import { gridify_nmap } from "./nmap.js";
import { gridify_cmds } from "./cmds.js";
//import { gridify_gridfit } from "./gridfit.js";

export function gridify(data, method = "hilbert", parameters = {}) {
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
        const gridded_bounding_box = [0, 0, size, size]
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
            parameters = {}
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
            const { x1, x2, y1, y2 } = parameters.Gamma
            const w = Math.max(x2 - x1, y2 - y1);
            const area = w ** 2 / ((Math.sqrt(data.length) + 1) ** 2)
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