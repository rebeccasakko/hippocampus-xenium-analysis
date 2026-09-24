# ============================================================
# Full pipeline: raw samples -> classification -> composition ->
# visualization -> differential expression (both methods)
#
# Run with the working directory set to the repository root (e.g. open
# the .Rproj, or setwd() here first) so the source() calls below
# resolve correctly.
#
# Structured as one reusable function (run_region_pipeline()) called
# once per brain region, because the underlying data doesn't change
# per region -- same 6 mice, same raw Xenium output -- only the lasso
# ROI polygon (and therefore which cells end up in the object) does.
# `sample_base` below is the shared, region-independent sample lineup;
# each region supplies its own lasso CSVs and its own DE cell-type list
# (a region's atlas_identity composition determines which cell types
# actually have enough cells to test -- e.g. optic tract, being white
# matter, won't have dentate gyrus granule or GABAergic interneuron
# cells the way hippocampus does).
#
# Classification itself (classify_atlas_index()) needs NO changes
# between regions: it's a per-cell argmax over the same 13 atlas
# panels regardless of region, and panels for cell types genuinely
# absent from a region's tissue just don't win the argmax for any
# cell there (the same mechanism that already gives the Purkinje
# sanity check ~0% in hippocampus). See protocols/classification_
# protocol.md.
#
# This script is a map of how the pieces fit together, not something
# meant to run unattended in one click. EDIT THE CONFIG SECTIONS below
# with your actual data paths before running, and check each stage's
# printed/plotted output before moving to the next, same as you would
# running it by hand.
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

# ============================================================
# run_region_pipeline(): stages 1-5 for one brain region.
#
#   sample_base  -- tibble with sample_id, group_id, data_dir (shared
#                   across all regions -- same 6 mice, same raw data)
#   lasso_map    -- named character vector, names = sample_id, values
#                   = that region's lasso ROI CSV path for that sample
#   region_label -- e.g. "Hippocampus", "Optic Tract" (goes into
#                   obj$region and plot titles)
#   de_cell_types -- named list of list(label=, candidates=) for
#                    whichever cell types this region has enough of to
#                    run differential expression on (see
#                    R/de_candidate_panels.R -- not every region will
#                    have a candidate panel defined for every cell
#                    type yet; add new ones there as needed)
#   output_dir   -- e.g. "outputs/hippocampus" -- kept separate per
#                   region so runs don't overwrite each other
# ============================================================
run_region_pipeline <- function(sample_base, lasso_map, region_label, de_cell_types, output_dir) {

    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

    sample_config <- sample_base %>%
        mutate(lasso_csv = unname(lasso_map[sample_id]))
    if (any(is.na(sample_config$lasso_csv))) {
        stop("lasso_map is missing an entry for: ",
             paste(sample_config$sample_id[is.na(sample_config$lasso_csv)], collapse = ", "))
    }

    # ---- Stage 1: per-sample load -> QC -> normalize -> lasso-restrict -> classify ----
    sample_objects <- lapply(seq_len(nrow(sample_config)), function(i) {
        row <- sample_config[i, ]
        message("\n==== ", region_label, ": processing ", row$sample_id, " ====")
        process_hippocampus_sample(
            data_dir  = row$data_dir,
            lasso_csv = row$lasso_csv,
            sample_id = row$sample_id,
            group_id  = row$group_id,
            region    = region_label
        )
    })
    names(sample_objects) <- sample_config$sample_id

    # ---- Stage 2: merge (JoinLayers() collapses Seurat v5's per-sample layers) ----
    combined <- merge(sample_objects[[1]], y = sample_objects[-1],
                      add.cell.ids = names(sample_objects))
    combined <- JoinLayers(combined)
    saveRDS(combined, file.path(output_dir, "combined.rds"))  # gitignored (*.rds) -- local checkpoint only

    # ---- Stage 3: composition summary + bar chart ----
    composition <- build_composition_summary(sample_objects)
    group_colors <- c(gi = "#5BC8C2", gii = "#F2A8A0")
    composition_plot <- plot_composition_barchart(
        composition$summary_df, group_colors,
        title = paste0(region_label, " cell-type composition: gii vs gi"),
        composition_df = composition$composition_df, show_points = TRUE)
    print(composition_plot)
    ggsave(file.path(output_dir, "composition_barchart.png"), composition_plot,
           width = 9, height = 6, dpi = 200, bg = "white")

    # ---- Stage 4: UMAP + microglia state QC check ----
    combined <- FindVariableFeatures(combined)
    combined <- ScaleData(combined)
    combined <- RunPCA(combined, npcs = 30)
    combined <- RunUMAP(combined, dims = 1:30)

    umap_celltype <- plot_umap_celltype(combined)
    print(umap_celltype)
    ggsave(file.path(output_dir, "umap_celltype_gi_vs_gii.png"), umap_celltype,
           width = 10, height = 6, dpi = 200, bg = "white")

    umap_microglia_check <- plot_umap_microglia_check(combined)
    print(umap_microglia_check)
    ggsave(file.path(output_dir, "umap_microglia_check_gi_vs_gii.png"), umap_microglia_check,
           width = 10, height = 6, dpi = 200, bg = "white")

    # ---- Stage 5: differential expression, both methods, per this region's cell types ----
    de_results <- list()
    for (key in names(de_cell_types)) {
        label <- de_cell_types[[key]]$label
        message("\n==== ", region_label, " differential expression: ", label, " ====")

        dds       <- build_pseudobulk_dds(combined, label)
        res_full  <- get_pseudobulk_results(dds)
        res_panel <- run_pseudobulk_panel_de(res_full, de_cell_types[[key]]$candidates, dds)
        plot_panel_heatmap(dds, res_panel, label, file.path(output_dir, sprintf("pseudobulk_heatmap_%s.png", key)))

        wilcox_res <- run_wilcox_degs(combined, label)
        plot_wilcox_heatmap(dds, wilcox_res, label, file.path(output_dir, sprintf("wilcox_heatmap_%s.png", key)))

        de_results[[key]] <- list(dds = dds, res_full = res_full, res_panel = res_panel, wilcox_res = wilcox_res)
    }

    sig_counts <- sapply(de_results, function(x) sum(x$wilcox_res$p_val_adj < 0.05))
    cat("\n==== ", region_label, ": Wilcoxon significant gene counts (padj < 0.05) ====\n")
    print(sig_counts)

    list(sample_objects = sample_objects, combined = combined,
         composition = composition, de_results = de_results, sig_counts = sig_counts)
}

# ============================================================
# CONFIG -- sample lineup shared by every region (same 6 mice, same
# raw Xenium output). Fill in the real data_dir paths before running.
# ============================================================
sample_base <- tibble::tribble(
    ~sample_id, ~group_id, ~data_dir,
    "gi_1_b",   "gi",      "/path/to/raw/gi_1_b/",
    "gi_2_d",   "gi",      "/path/to/raw/gi_2_d/",
    "gi_3_f",   "gi",      "/path/to/raw/gi_3_f/",
    "gii_1_a",  "gii",     "/path/to/raw/gii_1_a/",
    "gii_2_c",  "gii",     "/path/to/raw/gii_2_c/",
    "gii_3_e",  "gii",     "/path/to/raw/gii_3_e/"
)

# ============================================================
# CONFIG -- Hippocampus: lasso ROI CSVs + DE cell-type list
# ============================================================
hippocampus_lassos <- c(
    gi_1_b  = "~/Desktop/hippocampus_gi_1_b.csv",
    gi_2_d  = "~/Desktop/hippocampus_gi_2_d.csv",
    gi_3_f  = "~/Desktop/hippocampus_gi_3_f.csv",
    gii_1_a = "~/Desktop/hippocampus_gii_1_a.csv",
    gii_2_c = "~/Desktop/hippocampus_gii_2_c.csv",
    gii_3_e = "~/Desktop/hippocampus_gii_3_e.csv"
)

hippocampus_de_cell_types <- list(
    dg        = list(label = "Dentate gyrus granule", candidates = dg_candidates),
    gaba      = list(label = "GABAergic interneuron", candidates = gaba_candidates),
    astro     = list(label = "Astrocyte",              candidates = astro_candidates),
    microglia = list(label = "Microglia",               candidates = microglia_candidates)
)

# ============================================================
# CONFIG -- Optic Tract: lasso ROI CSVs + DE cell-type list
#
# Same 6 mice/data_dir as hippocampus (sample_base above) -- only the
# lasso polygon differs. Fill in the real paths once the new lassos
# are drawn. White matter is oligodendrocyte/astrocyte/microglia/
# vasculature-dominated, so dg/gaba are dropped from the DE list here
# (hippocampus-specific cell types that won't have meaningful numbers
# of cells in optic tract) -- astro and microglia reuse the SAME
# candidate panels as hippocampus (R/de_candidate_panels.R), since
# those panels are about the cell type's biology, not the region.
# Oligodendrocyte lineage has no candidate panel defined yet -- add
# one to R/de_candidate_panels.R and a new entry below once you have
# literature-backed candidates for it.
# ============================================================
optic_tract_lassos <- c(
    gi_1_b  = "/path/to/lasso/optic_tract_gi_1_b.csv",
    gi_2_d  = "/path/to/lasso/optic_tract_gi_2_d.csv",
    gi_3_f  = "/path/to/lasso/optic_tract_gi_3_f.csv",
    gii_1_a = "/path/to/lasso/optic_tract_gii_1_a.csv",
    gii_2_c = "/path/to/lasso/optic_tract_gii_2_c.csv",
    gii_3_e = "/path/to/lasso/optic_tract_gii_3_e.csv"
)

optic_tract_de_cell_types <- list(
    astro     = list(label = "Astrocyte", candidates = astro_candidates),
    microglia = list(label = "Microglia", candidates = microglia_candidates)
)

# ============================================================
# RUN -- one call per region. Comment out whichever isn't ready yet.
# ============================================================
hippocampus_results <- run_region_pipeline(
    sample_base, hippocampus_lassos, "Hippocampus",
    hippocampus_de_cell_types, "outputs/hippocampus")

optic_tract_results <- run_region_pipeline(
    sample_base, optic_tract_lassos, "Optic Tract",
    optic_tract_de_cell_types, "outputs/optic_tract")
