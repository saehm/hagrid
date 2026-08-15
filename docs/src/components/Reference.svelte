<script>
    /**
     * One bibliographic reference.
     *
     * Everything but `title` is optional, so this suits a journal article, a
     * talk or a bare link equally. A `doi` becomes the link; otherwise `url`
     * does, and with neither the entry is plain text.
     */
    let {
        /** "Surname, A. and Surname, B." — written however you cite. */
        authors = "",
        title,
        /** Journal, conference or venue. */
        venue = "",
        volume = "",
        issue = "",
        pages = "",
        year = "",
        doi = "",
        url = "",
        /** One line on why it is cited here. */
        note = "",
    } = $props();

    let href = $derived(doi ? `https://doi.org/${doi}` : url);
    let locator = $derived(
        [volume && issue ? `${volume}(${issue})` : volume, pages, year]
            .filter(Boolean)
            .join(", "),
    );
</script>

<div class="reference">
    <p class="cite">
        {#if authors}<span class="authors">{authors}</span>{/if}
        {#if href}
            <a class="title" {href} rel="noreferrer">{title}</a>
        {:else}
            <span class="title">{title}</span>
        {/if}
        <!-- The closing period belongs to the venue/locator clause, so an entry
             with neither — a bare repository, say — must not emit it: the title
             has already supplied one and the reader would see two. -->
        {#if venue}<span class="venue">{venue}</span>{/if}{#if locator}<span class="locator">
                {" "}{locator}</span>{/if}{#if venue || locator}.{/if}
        {#if doi}<a class="doi" href="https://doi.org/{doi}" rel="noreferrer">doi:{doi}</a>{/if}
    </p>
    {#if note}<p class="note">{note}</p>{/if}
</div>

<style>
    .reference {
        margin: 0.75rem 0;
        padding-left: 0.85rem;
        border-left: 2px solid var(--sl-color-gray-5);
    }
    .cite {
        margin: 0;
        font-size: var(--sl-text-sm);
        line-height: 1.55;
    }
    .authors::after { content: ". "; }
    .title {
        font-weight: 600;
        color: var(--sl-color-white);
    }
    a.title { text-decoration: none; }
    a.title:hover { text-decoration: underline; }
    .title::after { content: ". "; }
    .venue { font-style: italic; }
    .locator { font-variant-numeric: tabular-nums; }
    .doi {
        margin-left: 0.35rem;
        font-size: var(--sl-text-xs);
        color: var(--sl-color-gray-3);
    }
    .note {
        margin: 0.25rem 0 0;
        font-size: var(--sl-text-xs);
        color: var(--sl-color-gray-3);
    }
</style>
