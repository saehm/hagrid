import * as d3 from d3;
import * as druid from "@saehrimnir/druidjs";

function size_increase(X, Y) {
    const x = d3.scaleLinear()
        .domain(d3.extent(X, d => d[0]))
        .range(d3.extent(Y, d => d[0]))
    const y = d3.scaleLinear()
        .domain(d3.extent(X, d => d[1]))
        .range(d3.extent(Y, d => d[1]))
    const remap = ([x0, y0]) => [x(x0), y(y0)];

    const X_hull = d3.polygonHull(X.map(remap))
    const Y_hull = d3.polygonHull(Y)
    const X_size = d3.polygonArea(X_hull)
    const Y_size = d3.polygonArea(Y_hull)
    return {
        value: X_size / Y_size,
        measures: [X_hull, Y_hull],
    }    
}

function euclidean_distances(X, Y) {
    const x = d3.scaleLinear()
        .domain(d3.extent(X, d => d[0]))
        .range(d3.extent(Y, d => d[0]))
    const y = d3.scaleLinear()
        .domain(d3.extent(X, d => d[1]))
        .range(d3.extent(Y, d => d[1]))
    const remap = ([x0, y0]) => [x(x0), y(y0)];

    const distances = d3.zip(X.map(remap), Y).map(([a, b]) => druid.euclidean(a, b));
    return {
        value: druid.neumair_sum(distances) / X.length,
        measures: distances,
    }    
}

function knn_preservation(X, Y, k=15, DX, DY) {
    DX = !!DX ? distance_matrix(X) : DX;
    DY = !!DY ? distance_matrix(Y) : DY;
    let N = X.length;
    let K = new Array(N).fill(0)
    for (let i = 0; i < N; ++i) {
        let knnX = DX[i]
            .map((d,i) => [d, i])
            .sort((a, b) => a[0] - b[0])
            .slice(0, k)
            .map(d => d[1]);
        let knnY = DY[i]
            .map((d,i) => [d, i])
            .sort((a, b) => a[0] - b[0])
            .slice(0, k)
            .map(d => d[1]);
        for (let j = 0; j < k; ++j) {
            if (knnY.find(d => d == knnX[j])) K[i] += 1
        }
        K[i] /= k;
    }
    return {
        value: druid.neumair_sum(K) / K.length,
        measures: K
    }
}

function alt_cc(X, Y) {
    const N = X.length;
    let DX = distance_matrix(X);
    let DY = distance_matrix(Y);
    const DXf = DX.flat();
    const DYf = DY.flat();

    const mx = d3.mean(DXf);
    const my = d3.mean(DYf);
    const ox = d3.deviation(DXf);
    const oy = d3.deviation(DYf);

    const l = druid.neumair_sum(DX.map(dx => druid.neumair_sum(dx.map(d => d - mx))))
    const r = druid.neumair_sum(DY.map(dy => druid.neumair_sum(dy.map(d => d - my))))
    
    
    const o = ox * oy
    return (l * r / o)
}

function crosscorrelation_all(X, Y) {
    const N = X.length;
    const DX = distance_matrix(X);
    const DY = distance_matrix(Y);

    let DXf = DX.flat();
    let DYf = DY.flat();

    const DXmax = d3.max(DXf, Math.abs)
    DXf = DXf.map(x => x / DXmax)
    const DYmax = d3.max(DYf, Math.abs)
    DYf = DYf.map(y => y / DYmax)

    const ox = d3.deviation(DXf);
    const oy = d3.deviation(DYf);
    const mx = d3.mean(DXf);
    const my = d3.mean(DYf);
    console.log(ox, oy, mx, my)
    let CC = new Array(N).fill(0)
    //CC = CC.map(() => new Array(N).fill(0));
    for (let i = 0; i < N; ++i) {
        for (let j = 0; j < N; ++j) {
            CC[i] += crosscorrelation(DXf[i * N + j], DYf[i * N + j], mx, my, ox, oy);
        }
    }
    CC = CC.map(c => c/(ox*oy))
    return {
        value: druid.neumair_sum(CC),
        measures: CC,
    };
}

function crosscorrelation(dx, dy, mx, my, ox, oy) {
    const l = dx - mx;
    const r = dy - my;
    const o = ox * oy;
    console.log(dx, dy, mx, my, ox, oy, l, r, (l*r))
    
    return (l*r)///o;
}

function distance_matrix(A) {
    const N = A.length;
    let M = new Array(N).fill(0)
    M = M.map(() => new Array(N).fill(0));
    for (let i = 0; i < N; ++i) {
        let [x1, y1] = A[i];
        for (let j = i + 1; j < N; ++j) {
            let [x2, y2] = A[j];
            let d = Math.sqrt(Math.pow(x1-x2, 2) + Math.pow(y1-y2, 2));
            M[i][j] = d;
            M[j][i] = d;
        }
    }
    return M;
}



function distance_correlation(X, Y) {
    const N = X.length;
    const C = new druid.Matrix(N, N, "center");
    const XX = C.dot(druid.Matrix.from(X)).to2dArray
    const YY = C.dot(druid.Matrix.from(Y)).to2dArray
    console.log(XX, YY)
    const DX = distance_matrix(XX);
    const DY = distance_matrix(YY);

    const max_DX = d3.max(DX.flat(), d => Math.abs(d));
    const max_DY = d3.max(DY.flat(), d => Math.abs(d));

    for (let i = 0; i < N; ++i) {
        for (let j = 0; j < N; ++j) {
            DX[i][j] /= max_DX;
            DY[i][j] /= max_DY
        }
    }

    const idx = d3.range(0, N)
    const a = d3.mean(DX.flat());
    const b = d3.mean(DY.flat());
    console.log(a, b)
    let DC = new Array(N).fill(0)
    DC = DC.map(() => new Array(N).fill(0));
    
    for (let i = 0; i < N; ++i) {
        let ai_ = d3.mean(DX[i])
        let bi_ = d3.mean(DY[i])
        let a_j = d3.mean(idx.map(j => DX[i][j]));
        let b_j = d3.mean(idx.map(j => DY[i][j]));
        for (let j = 0; j < N; ++j) {
            const aij = +DX[i][j];
            const bij = +DY[i][j];
            const Aij = aij - ai_ - a_j + a;
            const Bij = bij - bi_ - b_j + b;
            DC[i][j] = Aij * Bij
            //console.log(Aij, Bij)
        }
    }

    return {
        "value": Math.sqrt((druid.neumair_sum(DC.flat()) / N**2)),
        "measures": DC.map(druid.neumair_sum).map(d => d / N),
    }
}
