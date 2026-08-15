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
        coverage: {
            provider: "v8",
            // The library only — the rollup config and this file are build
            // scaffolding, and counting them would flatter the number.
            include: ["src/**/*.ts"],
            reporter: ["text", "lcov"],
        },
    },
});
