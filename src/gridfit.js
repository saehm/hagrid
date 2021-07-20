function area([x_min, x_max, y_min, y_max]) {
    return (x_max - x_min) * (y_max - y_min);
}

function extent(indices, data) {
    return [
        ...d3.extent(indices, (i) => data[i][0]),
        ...d3.extent(indices, (i) => data[i][1]),
    ];
}

function remove_from_tree(node, value) {
    if (!node.includes(value)) return;
    const i = node.findIndex(n => n == value);
    node.splice(i, 1);
    node.children.forEach(child => remove_from_tree(child, value))
}

function get_direction(Pa, Aa, Pb, Ab) {
    return Pa <= Aa ? -1 : Pb <= Ab ? 1 : 0
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
        console.log("xs", x_min, x_mid, x_max, P13, A13, P24, A24)
    }

    // search fro y_1
    let y_1 = y_mid;
    let A1 = area([x_min, x_mid, y_min, y_mid]);
    let A3 = area([x_min, x_mid, y_mid, y_max]);
    if (P1 + P3 > A1 + A3) {
        console.log(P1, A1, P3, A3);
        console.log("xs", x_min, x_mid, x_max);
        console.log("ys", y_min, y_mid, y_max, B1.E[3], B2.E[3], B3.E[3], B4.E[3]);
        console.log(node.map(i => data[i]))
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
export function gridify_gridfit(D) {
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
