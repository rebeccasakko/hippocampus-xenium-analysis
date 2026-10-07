# ============================================================
# Secondary (exploratory) DE: Seurat FindMarkers, Wilcoxon, full gene panel
#
# Cell-level replication: every individual cell is treated as an
# independent observation, tested across the FULL gene panel (no
# candidate-list restriction). This is the method described in the
# Qiu et al. (2022, Front. Neurosci.), matched to its exact thresholds (logfc.threshold
# = 0.25, min.pct = 0.1 -- Seurat's own defaults). Much higher power than
# the pseudobulk approach in R/run_pseudobulk_de.R, at the cost of
# pseudoreplication risk -- see protocols/differential_expression_
# protocol.md for the full comparison and citation (Squair et al., 2021,
# Nature Communications, "Confronting false discoveries in single-cell
# differential expression").
# ============================================================

run_wilcox_degs <- function(xenium, cell_type_filter, identity_col = "atlas_identity",
                            group_col = "group_id") {
    Idents(xenium) <- identity_col
    subset_obj <- subset(xenium, idents = cell_type_filter)
    Idents(subset_obj) <- group_col

    FindMarkers(
        subset_obj,
        ident.1         = "gi",
        ident.2         = "gii",
        test.use        = "wilcox",
        logfc.threshold = 0.25,
        min.pct         = 0.1,
        only.pos        = FALSE
    ) %>%
        tibble::rownames_to_column("gene") %>%
        arrange(p_val_adj)
}

# ============================================================
# Usage per cell type (repeat for dg, gaba, astro, microglia):
#
#   wilcox_astro <- run_wilcox_degs(combined, "Astrocyte")
#   sum(wilcox_astro$p_val_adj < 0.05)   # significant gene count
# ============================================================
