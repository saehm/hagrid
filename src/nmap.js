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

function getBounds(D) {
    let min_x = Infinity;
    let max_x = -Infinity;
    let min_y = Infinity;
    let max_y = -Infinity;

    D.forEach(([x, y, i]) => {
        min_x = Math.min(min_x, x)
        max_x = Math.max(max_x, x)
        min_y = Math.min(min_y, y)
        max_y = Math.max(max_y, y)
    })

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
    let extra = []
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
        D.push(extra[extra.length - 1])
        --sq_missing;
    }

    return D;
}

function calcDist(A, B) {
    for (let i = 0; i < A.length; ++i) {
        const a = A[i]
        for (let j = 0; j < B.length; ++j) {
            const b = B[j]
            let t_dist = Math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2)
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

export async function gridify_nmap(D, parameters) {
    let G = [];
    console.log(parameters)
    if (!"BB" in parameters || parameters.BB == null) {
        console.log("no BB")
        parameters.BB = getBounds(D)
    }
    let added_index = D.map(([x, y], i) => [x, y, i]);
    if ("squared" in parameters && parameters.squared == true) {
        added_index = squared(added_index, parameters.BB);
    }
    nmap(G, added_index, parameters.BB, false);
    return G.sort((a, b) => a.d[2] - b.d[2]).filter(d => d.d[2] != undefined)//.map(g => [g.j, g.i])
}