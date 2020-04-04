export { version } from '../package.json';
export * from "./hilbert.js";
export * from "./gosper.js";

/*export default function() {
    console.log(`version ${version}`);
    hilbert_encode([1,1], 1);
    hilbert_decode;
    gosper_encode;
    gosper_decode;
}*/

/*
export function gridify(X, plus_level = 0) {
    const N = X.length;
    const l = Math.ceil(Math.log(N) / Math.log(4)) + plus_level;
    const curveSize = 2 ** l;
    const P = new Map();
    const Y = new Array(N).fill(0);
    const x = d3.scaleLinear()
        .domain(d3.extent(X, d => d[0]))
        .range([0, curveSize])
    
    const y = d3.scaleLinear()
        .domain(d3.extent(X, d => d[1]))
        .range([0, curveSize])
    
    const XX = X.map(([x1, y1]) => [x(x1), y(y1)]);
    
    const times = [];
    XX.forEach((d, i) => {
        const start = performance.now()
        let p = hilbert_encode(d, curveSize);
        if (P.has(p)) {
            let q = P.get(p);
            while(P.has(p)) {
                p = (p + 1) % (curveSize ** 2) 
            }
        }
        P.set(p, d);
        Y[i] = hilbert_decode(p, curveSize).map(d => d / curveSize);

        const end = performance.now()
        times.push(end-start)
    })
    console.log(`mean ${d3.mean(times.slice(20))} and deviation ${d3.deviation(times.slice(20))}`)
    
    return Y;
}*/