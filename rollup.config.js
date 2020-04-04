import json from '@rollup/plugin-json';
import {terser} from 'rollup-plugin-terser';

export default {
  input: 'src/main.js',
  output: [{
    file: 'dist/hagrid.js',
    name: 'hagrid',
    format: 'umd'
  },
  {
    file: 'dist/hagrid.min.js',
    format: 'iife',
    name: 'hagrid',
    plugins: [terser()]
  }],
  plugins: [ json() ]
};