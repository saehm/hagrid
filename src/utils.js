export function get_bounds(A) {
    let min_x = Infinity;
    let max_x = -Infinity;
    let min_y = Infinity;
    let max_y = -Infinity;

    A.forEach(([x, y, i]) => {
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

export function remap([x, y], [x1, y1, w1, h1], [x2, y2, w2, h2]) {
    return [
        (x - x1) / w1 * w2 + x2, 
        (y - y1) / h1 * h2 + y2,
    ];
}


export function get_scales(data, [[X_min, X_max], [Y_min, Y_max]], {
    keep_aspect_ratio = false, 
    round = true, 
    x = d => d[0], 
    y = d => d[1]
}) {
    const X_span = X_max - X_min;
    const Y_span = Y_max - Y_min; 
    let [x_min, x_max] = extent(data, x);
    let [y_min, y_max] = extent(data, y);
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

export function extent(data, accessor) {
    let min = Infinity;
    let max = -Infinity;
    for (const d of data) {
        const v = accessor(d);
        min = Math.min(min, v);
        max = Math.max(max, v);
    }
    return [min, max]    
}

export function distance([ax, ay], [bx, by]) {
    return Math.hypot(bx - ax, by - ay);//Math.sqrt(Math.pow(ax - bx, 2) + Math.pow(ay - by, 2));
}