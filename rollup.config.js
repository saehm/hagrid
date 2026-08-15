import typescript from "@rollup/plugin-typescript";
import resolve from "@rollup/plugin-node-resolve";
import terser from "@rollup/plugin-terser";
import replace from "@rollup/plugin-replace";
import dts from "rollup-plugin-dts";
import pkg from "./package.json" with { type: "json" };

const banner = `/* ${pkg.name} v${pkg.version} — ${pkg.homepage} */`;

// d3-delaunay is a runtime dependency, so it must stay out of the ESM/CJS
// bundles. The browser build is the exception: a <script> tag has no resolver,
// so that one bundles it in.
const external = [...Object.keys(pkg.dependencies ?? {})];

const plugins = (opts = {}) => [
    replace({
        preventAssignment: true,
        __VERSION__: JSON.stringify(pkg.version),
    }),
    resolve(),
    typescript({
        tsconfig: "./tsconfig.json",
        // Types are emitted once, by the dts bundle below.
        declaration: false,
        declarationMap: false,
        outDir: undefined,
        ...opts,
    }),
];

export default [
    {
        input: "src/main.ts",
        external,
        output: [
            { file: pkg.module, format: "es", banner, sourcemap: true },
            { file: pkg.main, format: "cjs", banner, sourcemap: true, exports: "named" },
        ],
        plugins: plugins(),
    },
    {
        input: "src/main.ts",
        output: {
            file: "dist/hagrid.umd.js",
            format: "umd",
            name: "hagrid",
            banner,
            sourcemap: true,
            plugins: [terser()],
        },
        plugins: plugins(),
    },
    {
        // Two identical declaration bundles. Under node16/nodenext resolution a
        // `require` of this package looks for types next to the CJS entry and in
        // CJS form, so the ESM .d.ts alone leaves those consumers untyped.
        input: "src/main.ts",
        external,
        output: [
            { file: "dist/hagrid.d.ts", format: "es", sourcemap: true },
            { file: "dist/hagrid.d.cts", format: "es", sourcemap: true },
        ],
        plugins: [dts()],
    },
];
