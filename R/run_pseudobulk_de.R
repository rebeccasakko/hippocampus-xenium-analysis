# ============================================================
# Pseudobulk DESeq2 differential expression (secondary analysis)
#
# Mouse-level replication (n = 3 per group): all cells of one cell type,
# from one mouse, are aggregated into a single expression profile before
# testing. This is the conservative, standard-practice unit of biological
# replication -- see protocols/differential_expression_protocol.md for
# the full rationale and how this compares to the Wilcoxon/FindMarkers
# primary analysis in R/run_wilcox_de.R.
# ============================================================

# Builds the pseudobulk DESeqDataSet for one cell type from the merged,
# classified Seurat object. `cell_type_filter` should match a value in
# `xenium$cell_type` (e.g. "Astrocyte (Active)") or `xenium$atlas_identity`
# (e.g. "Dentate gyrus granule") depending on which grouping is wanted.
build_pseudobulk_dds <- function(xenium, cell_type_filter, identity_col = "atlas_identity",
                                 min_row_sum = 10) {
    Idents(xenium) <- identity_col
    subset_obj <- subset(xenium, idents = cell_type_filter)

    pseudobulk <- AggregateExpression(subset_obj, group.by = "sample_id",
                                      assays = "Xenium", slot = "counts")$Xenium

    sample_meta <- data.frame(
        sample_id = colnames(pseudobulk),
        group = ifelse(grepl("^gi-", colnames(pseudobulk)), "gi", "gii")
    )

    dds <- DESeqDataSetFromMatrix(countData = round(as.matrix(pseudobulk)),
                                   colData = sample_meta, design = ~ group)
    dds <- dds[rowSums(counts(dds)) >= min_row_sum, ]
    DESeq(dds)
}

# Full-panel results table (every gene passing the pseudobulk count
# filter), BH-corrected across the whole panel. This is what the
# genome-wide null result (zero significant genes, every cell type) was
# established from -- see the protocol doc for that diagnosis.
get_pseudobulk_results <- function(dds) {
    results(dds, contrast = c("group", "gi", "gii")) %>%
        as.data.frame() %>%
        tibble::rownames_to_column("gene") %>%
        arrange(padj)
}

# Candidate-panel-only results: re-corrects BH across just the candidate
# gene list (see R/de_candidate_panels.R) instead of the whole
# transcriptome, since testing thousands of untested genes' worth of
# multiple-comparison burden against a handful of pre-specified genes is
# both statistically wrong and unnecessarily conservative.
run_pseudobulk_panel_de <- function(res_full, candidate_genes, dds) {
    genes_present <- intersect(candidate_genes, rownames(dds))
    missing <- setdiff(candidate_genes, genes_present)
    if (length(missing) > 0) {
        message("Not present after the pseudobulk count filter (rowSums >= 10): ",
                paste(missing, collapse = ", "))
    }

    res_full %>%
        filter(gene %in% genes_present) %>%
        mutate(padj_panel = p.adjust(pvalue, method = "BH")) %>%
        select(gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj_panel) %>%
        arrange(padj_panel)
}

# ============================================================
# Usage per cell type (repeat for dg, gaba, astro, microglia):
#
#   dds_astro <- build_pseudobulk_dds(combined, "Astrocyte")
#   res_astro <- get_pseudobulk_results(dds_astro)
#   res_astro_panel <- run_pseudobulk_panel_de(res_astro, astro_candidates, dds_astro)
# ============================================================
