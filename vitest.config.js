import { defineConfig } from "vitest/config";
import pkg from "./package.json" with { type: "json" };

export default defineConfig({
    define: {
        // Mirrors the rollup replace plugin so src/main.ts resolves under test.
        __VERSION__: JSON.stringify(pkg.version),
    },
    test: {
        include: ["test/**/*.test.ts"],
        environment: "node",
    },
});
