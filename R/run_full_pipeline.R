# ============================================================
# Full pipeline: raw samples -> classification -> composition ->
# visualization -> differential expression (both methods)
#
# Run with the working directory set to the repository root (e.g. open
# hippocampus-xenium-analysis.Rproj, or setwd() here first) so the
# source() calls below resolve correctly.
#
# This script is the map of how the individual R/ functions fit
# together end to end -- it is not meant to be run unattended in one
# click the first time. EDIT THE CONFIG SECTION below with your actual
# data paths before running, and check each stage's printed/plotted
# output (cell_type counts, gate-cluster review messages, QC checks)
# before moving to the next stage, same as you would running it by hand.
#
# Required packages: Seurat, dplyr, tidyr, ggplot2, scales, sp, DESeq2,
# tibble, pheatmap
# ============================================================

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(sp)
library(DESeq2)
library(tibble)
library(pheatmap)

source("R/marker_panels.R")
source("R/classify_atlas_index.R")
source("R/process_sample_template.R")
source("R/build_composition_summary.R")
source("R/plot_composition_barchart.R")
source("R/plot_umap_celltype.R")
source("R/de_candidate_panels.R")
source("R/run_pseudobulk_de.R")
source("R/run_wilcox_de.R")
source("R/plot_de_heatmap.R")

dir.create("outputs", showWarnings = FALSE)

# ============================================================
# CONFIG -- edit these paths for your machine before running.
# One row per sample; group_id must be exactly "gi" or "gii".
# data_dir = the Xenium output folder LoadXenium() expects;
# lasso_csv = the hippocampus ROI boundary CSV for that sample
# (the ones read via `read.csv(lasso_csv, skip = 2)`).
# ============================================================
sample_config <- tibble::tribble(
    ~sample_id, ~group_id, ~data_dir,                          ~lasso_csv,
    "gi_1_b",   "gi",      "/path/to/raw/gi_1_b/",             "~/Desktop/hippocampus_gi_1_b.csv",
    "gi_2_d",   "gi",      "/path/to/raw/gi_2_d/",             "~/Desktop/hippocampus_gi_2_d.csv",
    "gi_3_f",   "gi",      "/path/to/raw/gi_3_f/",             "~/Desktop/hippocampus_gi_3_f.csv",
    "gii_1_a",  "gii",     "/path/to/raw/gii_1_a/",            "~/Desktop/hippocampus_gii_1_a.csv",
    "gii_2_c",  "gii",     "/path/to/raw/gii_2_c/",            "~/Desktop/hippocampus_gii_2_c.csv",
    "gii_3_e",  "gii",     "/path/to/raw/gii_3_e/",            "~/Desktop/hippocampus_gii_3_e.csv"
)
# NOTE: data_dir (the raw Xenium output folder LoadXenium() reads) is
# still a placeholder -- those live on "Joe" or elsewhere, not Desktop.
# Fill those in before running.

# ============================================================
# Stage 1: per-sample load -> QC -> normalize -> lasso-restrict -> classify
# ============================================================
sample_objects <- lapply(seq_len(nrow(sample_config)), function(i) {
    row <- sample_config[i, ]
    message("\n==== Processing ", row$sample_id, " ====")
    process_hippocampus_sample(
        data_dir  = row$data_dir,
        lasso_csv = row$lasso_csv,
        sample_id = row$sample_id,
        group_id  = row$group_id
    )
})
names(sample_objects) <- sample_config$sample_id

# ============================================================
# Stage 2: merge into one object for cross-sample analyses
# (JoinLayers() collapses the per-sample Seurat v5 layers so downstream
# functions like AggregateExpression/FindMarkers see one unified counts
# matrix rather than one layer per sample)
# ============================================================
combined <- merge(sample_objects[[1]], y = sample_objects[-1],
                  add.cell.ids = names(sample_objects))
combined <- JoinLayers(combined)

saveRDS(combined, "outputs/hippocampus_combined.rds")  # gitignored (*.rds) -- local checkpoint only

# ============================================================
# Stage 3: composition summary + bar chart
# ============================================================
composition <- build_composition_summary(sample_objects)
group_colors <- c(gi = "#5BC8C2", gii = "#F2A8A0")
composition_plot <- plot_composition_barchart(composition$summary_df, group_colors,
                                              composition_df = composition$composition_df,
                                              show_points = TRUE)
print(composition_plot)
ggsave("outputs/composition_barchart.png", composition_plot, width = 9, height = 6, dpi = 200, bg = "white")

# ============================================================
# Stage 4: UMAP visualization + microglia state QC check
# (classify_atlas_index() does not compute a reduction -- see
# protocols/classification_protocol.md -- so PCA/UMAP are run here,
# only for display, never feeding back into the classification itself)
# ============================================================
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined)
combined <- RunPCA(combined, npcs = 30)
combined <- RunUMAP(combined, dims = 1:30)

umap_celltype <- plot_umap_celltype(combined)
print(umap_celltype)
ggsave("outputs/umap_celltype_gi_vs_gii.png", umap_celltype, width = 10, height = 6, dpi = 200, bg = "white")

umap_microglia_check <- plot_umap_microglia_check(combined)
print(umap_microglia_check)
ggsave("outputs/umap_microglia_check_gi_vs_gii.png", umap_microglia_check, width = 10, height = 6, dpi = 200, bg = "white")

# ============================================================
# Stage 5: differential expression, both methods, all four cell types
# See protocols/differential_expression_protocol.md for why both
# methods are run and how to read them together.
# ============================================================
de_cell_types <- list(
    dg        = list(label = "Dentate gyrus granule", candidates = dg_candidates),
    gaba      = list(label = "GABAergic interneuron", candidates = gaba_candidates),
    astro     = list(label = "Astrocyte",              candidates = astro_candidates),
    microglia = list(label = "Microglia",               candidates = microglia_candidates)
)

de_results <- list()

for (key in names(de_cell_types)) {
    label <- de_cell_types[[key]]$label
    message("\n==== Differential expression: ", label, " ====")

    # ---- Secondary: pseudobulk DESeq2, candidate panel ----
    dds  <- build_pseudobulk_dds(combined, label)
    res_full  <- get_pseudobulk_results(dds)
    res_panel <- run_pseudobulk_panel_de(res_full, de_cell_types[[key]]$candidates, dds)
    plot_panel_heatmap(dds, res_panel, label, sprintf("outputs/pseudobulk_heatmap_%s.png", key))

    # ---- Primary: Wilcoxon FindMarkers, full gene panel ----
    wilcox_res <- run_wilcox_degs(combined, label)
    plot_wilcox_heatmap(dds, wilcox_res, label, sprintf("outputs/wilcox_heatmap_%s.png", key))

    de_results[[key]] <- list(dds = dds, res_full = res_full, res_panel = res_panel, wilcox_res = wilcox_res)
}

# Quick summary table matching protocols/differential_expression_protocol.md
sig_counts <- sapply(de_results, function(x) sum(x$wilcox_res$p_val_adj < 0.05))
cat("\n==== Wilcoxon significant gene counts (padj < 0.05) ====\n")
print(sig_counts)
