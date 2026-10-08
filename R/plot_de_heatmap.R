# ============================================================
# Heatmap visualization for both DE methods (R/run_pseudobulk_de.R,
# R/run_wilcox_de.R). Both methods use the SAME display convention --
# z-scored, per-mouse pseudobulk expression (not per-cell), gi/gii
# annotation bar, significance asterisks, free (unforced) column
# clustering -- so the two analyses' figures are directly comparable.
# Even the Wilcoxon-selected genes are displayed at the mouse level here:
# gene SELECTION can come from the cell-level test, but display always
# happens at the level a human can actually read and that means
# something biologically (one column per animal), never one column per
# cell. See protocols/differential_expression_protocol.md for why free
# (not forced-block) clustering is used: it's independent corroboration
# -- clustering that has no knowledge of group labels either confirms or
# fails to confirm the hypothesized gi/gii split, which forced ordering
# would hide.
# ============================================================

# Tiered significance asterisks matching the lab's convention:
# * padj<=0.05, ** padj<=0.01, *** padj<=0.001, **** padj<=0.0001
sig_label <- function(p) {
    dplyr::case_when(
        p <= 0.0001 ~ "****",
        p <= 0.001  ~ "***",
        p <= 0.01   ~ "**",
        p <= 0.05   ~ "*",
        TRUE        ~ ""
    )
}

# Draws AND saves in one call (pheatmap's `filename=` argument suppresses
# on-screen drawing, so this redraws the returned grob to a png() device
# instead -- keeps both the RStudio Plots pane and the saved file).
save_pheatmap_png <- function(x, filename, width = 6, height = 5, res = 300) {
    png(filename, width = width, height = height, units = "in", res = res)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
}

# Shared heatmap builder. `genes` should already be the final, filtered
# gene list to display (significant Wilcoxon genes, or a candidate
# panel) -- this function does not do any further gene selection.
plot_de_heatmap <- function(dds, genes, sig_lookup, cell_type_label, filename) {
    genes_present <- intersect(genes, rownames(dds))
    if (length(genes_present) == 0) {
        message("No genes present in the pseudobulk matrix for ", cell_type_label, " -- skipping.")
        return(invisible(NULL))
    }

    norm_counts <- log2(counts(dds, normalized = TRUE)[genes_present, , drop = FALSE] + 1)
    z <- t(scale(t(norm_counts)))

    rownames(z) <- ifelse(sig_lookup[rownames(z)] == "" | is.na(sig_lookup[rownames(z)]),
                          rownames(z), paste0(rownames(z), " ", sig_lookup[rownames(z)]))

    sample_order <- colnames(z)
    # Group from the DESeqDataSet's own colData (set from group_id metadata in
    # build_pseudobulk_dds(); audit finding B5), with the name-based rule as the
    # fallback. NOTE for the fallback: "^gi-" not "^gi" -- "gii" also starts with "gi".
    grp_all  <- setNames(as.character(dds$group), colnames(dds))
    group_of <- unname(grp_all[sample_order])
    if (anyNA(group_of)) group_of <- ifelse(grepl("^gi-", sample_order), "gi", "gii")
    ordered_samples <- c(sample_order[group_of == "gi"], sample_order[group_of == "gii"])
    z <- z[, ordered_samples, drop = FALSE]

    group_anno <- data.frame(Group = group_of[match(ordered_samples, sample_order)])
    rownames(group_anno) <- ordered_samples
    ann_colors <- list(Group = c(gi = "#5BC8C2", gii = "#F2A8A0"))

    p <- pheatmap::pheatmap(z,
        color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
        breaks = seq(-1.5, 1.5, length.out = 101),
        annotation_col = group_anno,
        annotation_colors = ann_colors,
        cluster_rows = length(genes_present) > 1,  # can't cluster a single row
        cluster_cols = TRUE,                       # free clustering -- see file header
        clustering_method = "complete",
        fontsize_row = 10,
        main = paste0(cell_type_label, ": gi vs. gii"))

    save_pheatmap_png(p, filename, height = max(4, 1.5 + 0.22 * length(genes_present)))
    invisible(p)
}

# ---- Wilcoxon primary analysis: genes are whatever cleared padj_cutoff ----
plot_wilcox_heatmap <- function(dds, wilcox_res, cell_type_label, filename, padj_cutoff = 0.05) {
    wilcox_res$sig_label <- sig_label(wilcox_res$p_val_adj)
    sig_genes <- wilcox_res$gene[wilcox_res$p_val_adj < padj_cutoff]
    sig_lookup <- setNames(wilcox_res$sig_label, wilcox_res$gene)
    plot_de_heatmap(dds, sig_genes, sig_lookup, paste0(cell_type_label, " (FindMarkers/Wilcoxon sig. genes)"), filename)
}

# ---- Pseudobulk secondary analysis: genes are the whole candidate panel ----
# (significant or not -- the point of this panel is to show where every
# pre-specified candidate landed, not just the ones that happened to clear
# the bar).
plot_panel_heatmap <- function(dds, panel_res, cell_type_label, filename,
                               gene_col = "gene", padj_col = "padj_panel") {
    sig_lookup <- setNames(sig_label(panel_res[[padj_col]]), panel_res[[gene_col]])
    plot_de_heatmap(dds, panel_res[[gene_col]], sig_lookup, paste0(cell_type_label, " candidate panel"), filename)
}

# ============================================================
# Usage:
#   plot_wilcox_heatmap(dds_astro, wilcox_astro, "Astrocyte", "outputs/wilcox_heatmap_astro.png")
#   plot_panel_heatmap(dds_astro, res_astro_panel, "Astrocyte", "outputs/pseudobulk_heatmap_astro.png")
# ============================================================
