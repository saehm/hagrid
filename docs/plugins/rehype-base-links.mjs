/**
 * Prefixes root-relative links in Markdown content with the site's `base`.
 *
 * The site deploys as a GitHub Pages *project* site, at `/hagrid/` rather than
 * at a domain root. Astro's `base` handles every URL the framework generates —
 * nav, sidebar, pagination — but not the ones written by hand inside content,
 * so `[Gilbert](/techniques/gilbert/)` would resolve to the domain root and
 * 404. The generated API reference is the larger half of the problem: it is
 * rewritten by TypeDoc on every build with around ninety absolute `/api/…`
 * cross-links, so fixing those by hand is not a thing that stays fixed.
 *
 * Doing it here means content can keep writing plain root-relative links, which
 * are also what works in the editor and on GitHub's own Markdown preview.
 *
 * Hand-rolled rather than using `unist-util-visit`: it is a six-line walk, and
 * the plugin should not add a dependency to the docs build for that.
 */
export function rehypeBaseLinks({ base = "/" } = {}) {
    const prefix = base.replace(/\/+$/, "");

    // A root deploy needs no rewriting, and the walk should not cost anything.
    if (!prefix) return () => () => {};

    const rewrite = (value) => {
        if (typeof value !== "string") return value;
        // Protocol-relative (`//host`) is absolute, not root-relative.
        if (!value.startsWith("/") || value.startsWith("//")) return value;
        // Idempotent: a link already carrying the prefix is left alone, so this
        // is safe if the plugin ever runs twice over the same tree.
        if (value === prefix || value.startsWith(`${prefix}/`)) return value;
        return prefix + value;
    };

    const ATTRS = { a: "href", img: "src", source: "src", video: "src", audio: "src" };

    return () => (tree) => {
        const walk = (node) => {
            if (node.type === "element") {
                const attr = ATTRS[node.tagName];
                if (attr && node.properties && attr in node.properties) {
                    node.properties[attr] = rewrite(node.properties[attr]);
                }
            }
            for (const child of node.children ?? []) walk(child);
        };
        walk(tree);
    };
}
