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
# Classification itself (classify_atlas_index()) mostly needs no changes
# between regions: it's a per-cell argmax over the same atlas panels
# regardless of region, and panels for cell types genuinely absent from a
# region's tissue mostly just don't win the argmax for any cell there (the
# same mechanism that gives the Purkinje sanity check ~0% in hippocampus).
# One confirmed exception: Hippocampal pyramidal and Dentate gyrus granule
# had to be dropped from the candidate list specifically for Optic Tract
# (optic_tract_markers in R/marker_panels.R) after a 2026-09-28 diagnostic
# showed they were winning the argmax there via a shared non-specific gene
# (Pou3f1, Prox1), not via real signal -- see protocols/classification_
# protocol.md for the full writeup and diagnostic.
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
source("R/apply_confidence_filter.R")
source("R/build_cell_count_table.R")
source("R/process_sample_template.R")
source("R/build_composition_summary.R")
source("R/plot_composition_barchart.R")
source("R/plot_umap_celltype.R")
source("R/plot_spatial_celltype.R")
source("R/de_candidate_panels.R")
source("R/run_pseudobulk_de.R")
source("R/run_wilcox_de.R")
source("R/plot_de_heatmap.R")
source("R/build_stats_table.R")
source("R/check_marker_spillover.R")
source("R/plot_lasso_overview.R")
source("R/check_driver_genes.R")

# Area (um^2) of the lasso polygon(s) in a Xenium Explorer selection CSV
# (2 comment lines, then Selection, X, Y, Class, Color). Returns NA if the
# file can't be read in that format, so it never stops a run.
lasso_polygon_area <- function(lasso_csv) {
  tryCatch({
    d <- read.csv(lasso_csv, skip = 2)
    poly_area <- function(x, y) abs(sum(x * c(y[-1], y[1]) - c(x[-1], x[1]) * y)) / 2
    sum(vapply(split(d, d$Selection), function(s) poly_area(s$X, s$Y), numeric(1)))
  }, error = function(e) NA_real_)
}

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
#   confidence_percentile -- if supplied, derive the confidence-margin
#                   threshold as this percentile of the region's pooled
#                   margins (NULL = no pooled filter)
#   confidence_categories -- which cell types the pooled filter may
#                   reclassify as Unclassified (e.g. "BAM"). NULL = all
#                   cell types (used for the global sensitivity run)
#   run_label    -- "primary" or "sens_global"; written into the
#                   cell-count table so runs can be told apart
#   prefilter_objects -- optional: the pre-filter per-sample objects from
#                   an earlier run (its $sample_objects_prefilter). If
#                   supplied, Stage 1 is skipped, so a sensitivity run
#                   doesn't have to reload the raw data
# ============================================================
run_region_pipeline <- function(sample_base, lasso_map, region_label, de_cell_types, output_dir,
                                marker_list = atlas_index_markers,
                                confidence_percentile = NULL,
                                confidence_categories = NULL,
                                run_label = "primary",
                                prefilter_objects = NULL) {
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  sample_config <- sample_base %>%
    mutate(lasso_csv = unname(lasso_map[sample_id]))
  if (any(is.na(sample_config$lasso_csv))) {
    stop("lasso_map is missing an entry for: ",
         paste(sample_config$sample_id[is.na(sample_config$lasso_csv)], collapse = ", "))
  }
  
  # Lasso overview (primary run only; sensitivity runs use the same lassos)
  if (is.null(prefilter_objects)) {
    plot_lasso_overview(sample_base, lasso_map, region_label,
                        file.path(output_dir, "lasso_overview.png"))
  }
  
  # ---- Stage 1: per-sample load -> QC -> normalize -> lasso-restrict -> classify ----
  if (is.null(prefilter_objects)) {
    sample_objects <- lapply(seq_len(nrow(sample_config)), function(i) {
      row <- sample_config[i, ]
      ckpt <- file.path(output_dir, paste0(row$sample_id, "_classified_prefilter.rds"))
      if (file.exists(ckpt)) {
        message("\n==== ", region_label, ": loading saved ", row$sample_id, " ====")
        return(readRDS(ckpt))
      }
      message("\n==== ", region_label, ": processing ", row$sample_id, " ====")
      obj <- process_hippocampus_sample(
        data_dir    = row$data_dir,
        lasso_csv   = row$lasso_csv,
        sample_id   = row$sample_id,
        group_id    = row$group_id,
        region      = region_label,
        marker_list = marker_list
      )
      saveRDS(obj, ckpt)   # gitignored (*.rds) -- lets an interrupted run resume
      gc()
      obj
    })
    names(sample_objects) <- sample_config$sample_id
    # Local checkpoint (gitignored, *.rds) so a sensitivity run or a crash
    # doesn't mean reloading the raw data
    saveRDS(sample_objects, file.path(output_dir, "sample_objects_prefilter.rds"))
  } else {
    message("\nUsing supplied pre-filter sample objects; skipping Stage 1 load.")
    sample_objects <- prefilter_objects
  }
  sample_objects_prefilter <- sample_objects
  
  # ---- Stage 1b: confidence-margin filter (optional, region-derived) ----
  # Only applied when confidence_percentile is supplied -- hippocampus
  # is called without it and stays on its existing BAM-specific fix.
  # See R/apply_confidence_filter.R for the rationale on deriving the
  # threshold from this region's own pooled margin distribution rather
  # than reusing a fixed number across regions.
  # confidence_categories limits WHICH cell types can be reclassified
  # (primary run: only types with a diagnosed near-tie problem, e.g. BAM);
  # NULL applies the filter to every cell (global sensitivity run).
  if (!is.null(confidence_percentile)) {
    sample_objects <- apply_confidence_filter_pooled(
      sample_objects, percentile = confidence_percentile,
      categories = confidence_categories)
  }
  
  # ---- Stage 1c: cells-per-sample table (standard output of every run) ----
  lasso_area_um2 <- setNames(vapply(sample_config$lasso_csv, lasso_polygon_area, numeric(1)),
                             sample_config$sample_id)
  cell_counts <- build_cell_count_table(
    sample_objects, region = region_label, run = run_label,
    lasso_area_um2 = lasso_area_um2,
    out_file = file.path(output_dir, "cell_counts_per_sample.csv"))
  print(cell_counts)
  
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
  
  list(sample_objects = sample_objects, sample_objects_prefilter = sample_objects_prefilter,
       combined = combined, composition = composition, cell_counts = cell_counts,
       de_results = de_results, sig_counts = sig_counts)
}

# ============================================================
# CONFIG -- sample lineup shared by every region (same 6 mice, same
# raw Xenium output). Fill in the real data_dir paths before running.
# ============================================================
sample_base <- tibble::tribble(
  ~sample_id, ~group_id, ~data_dir,
  "gi_1_b",   "gi",      "/Volumes/Joe/Becca and Sakthi Files/Xenium_Data/Group_i/0071470_gi_1_b/",
  "gi_2_d",   "gi",      "/Volumes/Joe/Becca and Sakthi Files/Xenium_Data/Group_i/0071470_gi_2_d/",
  "gi_3_f",   "gi",      "/Volumes/Joe/Becca and Sakthi Files/Xenium_Data/Group_i/0071470_gi_3_f/",
  "gii_1_a",  "gii",     "/Volumes/Joe/Becca and Sakthi Files/Xenium_Data/Group_ii/0071470_gii_1_a/",
  "gii_2_c",  "gii",     "/Volumes/Joe/Becca and Sakthi Files/Xenium_Data/Group_ii/0071470_gii_2_c/",
  "gii_3_e",  "gii",     "/Volumes/Joe/Becca and Sakthi Files/Xenium_Data/Group_ii/0071470_gii_3_e/"
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
# lasso polygon differs. White matter is oligodendrocyte/astrocyte/
# microglia/vasculature-dominated, so dg/gaba are dropped from the DE
# list here (hippocampus-specific cell types that won't have meaningful
# numbers of cells in optic tract) -- astro and microglia reuse the SAME
# candidate panels as hippocampus (R/de_candidate_panels.R), since
# those panels are about the cell type's biology, not the region.
# Oligodendrocyte lineage has no candidate panel defined yet -- add
# one to R/de_candidate_panels.R and a new entry below once you have
# literature-backed candidates for it.
#
# Classification itself also uses a region-specific candidate list here:
# optic_tract_markers (R/marker_panels.R) drops Hippocampal pyramidal and
# Dentate gyrus granule, confirmed 2026-09-28 to be gene-overlap false
# positives in this region rather than real signal -- see
# protocols/classification_protocol.md.
#
# LASSOS: from 2026-10-05 the optic tract uses redrawn, border-inclusive
# lassos (Desktop/B6_Xenium_Opt). The file-to-sample pairing was checked by
# overlaying each outline on its section (docs/optic_tract_new_lassos_overview.png).
# The old lassos are kept below for reference only and are not used.
# ============================================================
optic_tract_lassos <- c(
  gi_1_b  = "~/Desktop/B6_Xenium_Opt/GI-1_coordinates.csv",
  gi_2_d  = "~/Desktop/B6_Xenium_Opt/GI-2_coordinates.csv",
  gi_3_f  = "~/Desktop/B6_Xenium_Opt/GI-3_coordinates.csv",
  gii_1_a = "~/Desktop/B6_Xenium_Opt/GII-1_coordinates.csv",
  gii_2_c = "~/Desktop/B6_Xenium_Opt/GII-2_coordinates.csv",
  gii_3_e = "~/Desktop/B6_Xenium_Opt/GII-3_coordinates.csv"
)

optic_tract_lassos_old <- c(   # reference only -- NOT used
  gi_1_b  = "~/Desktop/optic_tract_gi_1_b.csv",
  gi_2_d  = "~/Desktop/optic_tract_gi_2_d.csv",
  gi_3_f  = "~/Desktop/optic_tract_gi_3_f.csv",
  gii_1_a = "~/Desktop/optic_tract_gii_1_a.csv",
  gii_2_c = "~/Desktop/optic_tract_gii_2_c.csv",
  gii_3_e = "~/Desktop/optic_tract_gii_3_e.csv"
)

optic_tract_de_cell_types <- list(
  astro     = list(label = "Astrocyte", candidates = astro_candidates),
  microglia = list(label = "Microglia", candidates = microglia_candidates)
)

# ============================================================
# The analysis runs live in R/run_analyses.R.
# This file is setup only and is safe to source():
#   source("R/run_full_pipeline.R")
# ============================================================