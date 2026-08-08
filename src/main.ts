declare const __VERSION__: string;

/** The library version, matching the published package. */
export const version: string = __VERSION__;

export * from "./types.js";
export * from "./hilbert.js";
export * from "./gosper.js";
export * from "./dgrid.js";
export * from "./nmap.js";
export * from "./cmds.js";
export * from "./gridfit.js";
export * from "./gridify.js";
export * as utils from "./utils.js";
