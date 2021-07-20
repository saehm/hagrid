import json from '@rollup/plugin-json';
import { terser } from 'rollup-plugin-terser';
import resolve from 'rollup-plugin-node-resolve';

export default {
  input: 'src/main.js',
  output: [{
    file: 'dist/hagrid.js',
    name: 'hagrid',
    format: 'umd'
  },
  {
    file: 'dist/hagrid.es.js',
    name: 'hagrid',
    format: 'es'
  },
  {
    file: 'dist/hagrid.min.js',
    format: 'umd',
    name: 'hagrid',
    plugins: [terser()]
  }],
  plugins: [ json(), resolve() ],
  external:[ "d3Delaunay" ],
};