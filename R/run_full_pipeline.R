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
# meant to run unattended in one click. Machine-specific paths live in
# config.local.R (gitignored; copy config.local.example.R), read in the
# SETTINGS block at the bottom. Check each stage's printed/plotted output
# before moving to the next, same as you would running it by hand.
#
# Required packages: Seurat, dplyr, tidyr, ggplot2, scales, sp, DESeq2,
# tibble, pheatmap, Matrix (checked by R/_setup.R before anything loads)
# ============================================================

# Everything downstream is resolved from the working directory, so check it first.
if (!file.exists("hippocampus-xenium-analysis.Rproj"))
  stop("Set the working directory to the repository root - open ",
       "hippocampus-xenium-analysis.Rproj in RStudio, or setwd() there. ",
       "Currently: ", getwd(), call. = FALSE)
dir.repo <- normalizePath(getwd(), winslash = "/")

source("R/_setup.R")      # base R only: package check, input manifest, input check
fn.check_packages()       # stops with the install line(s) if anything is missing

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
                                prefilter_objects = NULL,
                                skip_embedding_and_de = FALSE) {
  # skip_embedding_and_de: TRUE skips stages 4 (UMAP) and 5 (DE). The confidence
  # filter only rewrites cell_type while DE selects on atlas_identity, so a
  # sensitivity run's DE is identical to its primary run's (audit finding M2);
  # the default FALSE keeps the original behaviour and outputs.

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  stage_times <- list(start = Sys.time())   # written to run_info.txt at the end (audit T1)
  
  sample_config <- sample_base %>%
    mutate(lasso_csv = unname(lasso_map[sample_id]))
  if (any(is.na(sample_config$lasso_csv))) {
    stop("lasso_map is missing an entry for: ",
         paste(sample_config$sample_id[is.na(sample_config$lasso_csv)], collapse = ", "))
  }

  # Input check before any slow work: a missing bundle or lasso file stops here with a
  # table naming it, instead of inside LoadXenium() minutes in. A sensitivity run that
  # reuses loaded objects only needs the lasso files (for the area column).
  manifest <- fn.data_manifest(setNames(sample_config$data_dir, sample_config$sample_id),
                               setNames(sample_config$lasso_csv, sample_config$sample_id),
                               region = tolower(gsub("[^A-Za-z]+", "_", region_label)))
  if (!is.null(prefilter_objects)) manifest <- manifest[!grepl("^xenium_", manifest$key), ]
  fn.check_inputs(manifest)

  # Lasso overview (primary run only; sensitivity runs use the same lassos)
  if (is.null(prefilter_objects)) {
    plot_lasso_overview(sample_base, lasso_map, region_label,
                        file.path(output_dir, "lasso_overview.png"))
  }
  
  # ---- Stage 1: per-sample load -> QC -> normalize -> lasso-restrict -> classify ----
  # What a checkpoint was made with, so a rerun into the same folder after a
  # panel or lasso change does not silently load a stale classification
  # (audit finding R7). Stored beside the .rds as <ckpt>.meta.txt.
  ckpt_meta <- function(lasso_csv) {
    c(paste0("region=", region_label),
      paste0("panels=", paste(names(marker_list), collapse = ";")),
      paste0("lasso_md5=", unname(tools::md5sum(lasso_csv))),
      "qc=nCount_Xenium>=10;nFeature_Xenium>=5")
  }
  if (is.null(prefilter_objects)) {
    sample_objects <- lapply(seq_len(nrow(sample_config)), function(i) {
      row <- sample_config[i, ]
      ckpt <- file.path(output_dir, paste0(row$sample_id, "_classified_prefilter.rds"))
      meta_file <- paste0(ckpt, ".meta.txt")
      expected  <- ckpt_meta(row$lasso_csv)
      if (file.exists(ckpt)) {
        if (!file.exists(meta_file)) {
          message("\n==== ", region_label, ": loading saved ", row$sample_id,
                  " (checkpoint made before the sidecar existed; panels and lasso not verifiable) ====")
          return(readRDS(ckpt))
        } else if (identical(readLines(meta_file, warn = FALSE), expected)) {
          message("\n==== ", region_label, ": loading saved ", row$sample_id, " ====")
          return(readRDS(ckpt))
        } else {
          warning(sprintf("Checkpoint %s was made with a different region, panel set, lasso file or QC (see %s); recomputing it",
                          basename(ckpt), basename(meta_file)), call. = FALSE)
        }
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
      writeLines(expected, meta_file)
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
  stage_times$stage1_classified <- Sys.time()

  # Which panels were scored in which mouse (the 0.5% detection floor is per
  # sample), recorded by classify_atlas_index() in obj@misc$panel_status
  # (audit finding B6). Checkpoints made before that exist carry no table.
  panel_status <- do.call(rbind, lapply(names(sample_objects), function(s) {
    ps <- sample_objects[[s]]@misc[["panel_status"]]
    if (is.null(ps)) NULL else cbind(sample = s, ps, stringsAsFactors = FALSE)
  }))
  if (is.null(panel_status)) {
    message("panel_status_per_sample.csv not written: the loaded objects predate the panel-status record.")
  } else {
    write.csv(panel_status, file.path(output_dir, "panel_status_per_sample.csv"), row.names = FALSE)
  }

  # Per-sample depth of the ROI cells that passed QC (audit finding R3; the
  # previously committed hippocampus_per_sample_depth.csv had no producing code).
  depth <- data.frame(
    sample = names(sample_objects),
    group  = vapply(sample_objects, function(o) ifelse(as.character(unique(o$group_id))[1] == "gi", "TBI", "sham"), character(1)),
    median_transcripts = vapply(sample_objects, function(o) as.numeric(median(o$nCount_Xenium)), numeric(1)),
    median_genes       = vapply(sample_objects, function(o) as.numeric(median(o$nFeature_Xenium)), numeric(1)),
    row.names = NULL)
  write.csv(depth, file.path(output_dir, "per_sample_depth.csv"), row.names = FALSE)
  
  # ---- Stage 1b: confidence-margin filter (optional, region-derived) ----
  # Only applied when confidence_percentile is supplied: the primary run of
  # each region passes nothing and stays unfiltered; the two sensitivity
  # runs pass 0.15 (see R/run_analyses.R). See R/apply_confidence_filter.R
  # for the rationale on deriving the threshold from this region's own
  # pooled margin distribution rather than reusing a fixed number across
  # regions. confidence_categories limits WHICH cell types can be
  # reclassified (sens_neartie: the categories with a diagnosed near-tie
  # problem); NULL applies the filter to every cell (sens_global).
  confidence_threshold <- NA_real_
  if (!is.null(confidence_percentile)) {
    sample_objects <- apply_confidence_filter_pooled(
      sample_objects, percentile = confidence_percentile,
      categories = confidence_categories)
    confidence_threshold <- attr(sample_objects, "threshold")
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
  if (interactive()) print(composition_plot)
  ggsave(file.path(output_dir, "composition_barchart.png"), composition_plot,
         width = 9, height = 6, dpi = 200, bg = "white")
  # The per-sample table behind the chart (audit finding R1; previously written
  # by hand). Same columns as the committed composition_per_sample.csv files.
  write.csv(composition$composition_df, file.path(output_dir, "composition_per_sample.csv"),
            row.names = FALSE)
  stage_times$stage3_composition <- Sys.time()

  de_results <- list()
  sig_counts <- integer(0)
  if (isTRUE(skip_embedding_and_de)) {
    message("\nskip_embedding_and_de = TRUE: stages 4 (UMAP) and 5 (DE) skipped for this run.")
  } else {
  # ---- Stage 4: UMAP + microglia state QC check ----
  # Randomness: no set.seed() anywhere in this repository; RunPCA() and
  # RunUMAP() use Seurat's default seed.use = 42, AddModuleScore() (stage 1)
  # its default seed = 1, so results repeat for a given Seurat/uwot version.
  # UMAP coordinates can still differ between package versions (audit S1).
  combined <- FindVariableFeatures(combined)
  combined <- ScaleData(combined)
  combined <- RunPCA(combined, npcs = 30)
  combined <- RunUMAP(combined, dims = 1:30)

  umap_celltype <- plot_umap_celltype(combined)
  if (interactive()) print(umap_celltype)
  ggsave(file.path(output_dir, "umap_celltype_gi_vs_gii.png"), umap_celltype,
         width = 10, height = 6, dpi = 200, bg = "white")

  umap_microglia_check <- plot_umap_microglia_check(combined)
  if (interactive()) print(umap_microglia_check)
  ggsave(file.path(output_dir, "umap_microglia_check_gi_vs_gii.png"), umap_microglia_check,
         width = 10, height = 6, dpi = 200, bg = "white")
  stage_times$stage4_umap <- Sys.time()

  # ---- Stage 5: differential expression, both methods, per this region's cell types ----
  # Result tables and DESeq2 diagnostics go to <output_dir>/de/ (audit findings
  # R2 and M8; the tables were previously written by hand). Diagnostics are
  # wrapped so a plotting problem never stops a run.
  de_dir <- file.path(output_dir, "de")
  dir.create(de_dir, showWarnings = FALSE, recursive = TRUE)
  for (key in names(de_cell_types)) {
    label <- de_cell_types[[key]]$label
    message("\n==== ", region_label, " differential expression: ", label, " ====")

    dds       <- build_pseudobulk_dds(combined, label)
    res_full  <- get_pseudobulk_results(dds)
    res_panel <- run_pseudobulk_panel_de(res_full, de_cell_types[[key]]$candidates, dds)
    plot_panel_heatmap(dds, res_panel, label, file.path(output_dir, sprintf("pseudobulk_heatmap_%s.png", key)))

    wilcox_res <- run_wilcox_degs(combined, label)
    plot_wilcox_heatmap(dds, wilcox_res, label, file.path(output_dir, sprintf("wilcox_heatmap_%s.png", key)))

    write.csv(res_full,   file.path(de_dir, paste0(key, "_res_full.csv")),   row.names = FALSE)
    write.csv(res_panel,  file.path(de_dir, paste0(key, "_res_panel.csv")),  row.names = FALSE)
    write.csv(wilcox_res, file.path(de_dir, paste0(key, "_wilcox_res.csv")), row.names = FALSE)
    tryCatch({
      png(file.path(de_dir, paste0(key, "_pvalue_hist.png")), width = 6, height = 4, units = "in", res = 150)
      hist(res_full$pvalue, breaks = 40, main = paste0(label, ": pseudobulk DESeq2 p-values"),
           xlab = "p-value (gi vs gii)", col = "grey80")
      dev.off()
      png(file.path(de_dir, paste0(key, "_dispersion.png")), width = 6, height = 5, units = "in", res = 150)
      DESeq2::plotDispEsts(dds, main = paste0(label, ": DESeq2 dispersion estimates"))
      dev.off()
    }, error = function(e) {
      try(dev.off(), silent = TRUE)
      message("DESeq2 diagnostic plots not written for ", label, ": ", conditionMessage(e))
    })

    de_results[[key]] <- list(dds = dds, res_full = res_full, res_panel = res_panel, wilcox_res = wilcox_res)
  }

  sig_counts <- sapply(de_results, function(x) sum(x$wilcox_res$p_val_adj < 0.05))
  cat("\n==== ", region_label, ": Wilcoxon significant gene counts (padj < 0.05) ====\n")
  print(sig_counts)
  stage_times$stage5_de <- Sys.time()
  }

  # ---- Run record: what this run was, and the environment it ran in (audit R5, B4, T1) ----
  stage_times$end <- Sys.time()
  stage_names <- names(stage_times)
  stage_lines <- vapply(seq_along(stage_times)[-1], function(k)
    sprintf("time_%s=%s (%.1f min since %s)", stage_names[k],
            format(stage_times[[k]], "%Y-%m-%d %H:%M:%S"),
            as.numeric(difftime(stage_times[[k]], stage_times[[k - 1]], units = "mins")),
            stage_names[k - 1]), character(1))
  run_info <- c(
    paste0("run_date=", format(stage_times$start, "%Y-%m-%d %H:%M:%S")),
    paste0("region=", region_label),
    paste0("run_label=", run_label),
    paste0("output_dir=", output_dir),
    paste0("git_commit=", fn.git_head(dir.repo)),
    paste0("marker_panels=", paste(names(marker_list), collapse = ";")),
    paste0("confidence_percentile=", if (is.null(confidence_percentile)) "none" else confidence_percentile),
    paste0("confidence_categories=", if (is.null(confidence_categories)) "all" else paste(confidence_categories, collapse = ";")),
    paste0("confidence_threshold=", confidence_threshold),
    paste0("skip_embedding_and_de=", skip_embedding_and_de),
    paste0("de_cell_types=", paste(vapply(de_cell_types, `[[`, character(1), "label"), collapse = ";")),
    paste0("n_cells=", paste(sprintf("%s:%d", cell_counts$sample, cell_counts$n_cells), collapse = ";")),
    paste0("stage1_source=", if (is.null(prefilter_objects)) "loaded or checkpoint" else "prefilter_objects supplied"),
    paste0("total_minutes=", round(as.numeric(difftime(stage_times$end, stage_times$start, units = "mins")), 1)),
    stage_lines,
    paste0("R=", R.version.string),
    paste0("Seurat=", as.character(utils::packageVersion("Seurat"))),
    paste0("DESeq2=", as.character(utils::packageVersion("DESeq2"))))
  writeLines(run_info, file.path(output_dir, "run_info.txt"))
  writeLines(utils::capture.output(utils::sessionInfo()), file.path(output_dir, "session_info.txt"))

  list(sample_objects = sample_objects, sample_objects_prefilter = sample_objects_prefilter,
       combined = combined, composition = composition, cell_counts = cell_counts,
       de_results = de_results, sig_counts = sig_counts,
       confidence_threshold = confidence_threshold, stage_times = stage_times)
}

# ============================================================
# SETTINGS -- repo-relative defaults, then local overrides.
#
# The six mice are the same for every region; only the lasso polygon
# differs. Machine-specific paths do not belong in this tracked file:
# copy config.local.example.R to config.local.R (gitignored) and set
# xenium_dirs / hippocampus_lassos / optic_tract_lassos there. The
# defaults below expect the inputs under data/ (gitignored); see README
# "Getting the data".
# ============================================================
sample_ids    <- c("gi_1_b", "gi_2_d", "gi_3_f", "gii_1_a", "gii_2_c", "gii_3_e")
sample_groups <- c("gi", "gi", "gi", "gii", "gii", "gii")   # gi = TBI, gii = sham (unblinded 2026-09-30)

dir.data <- file.path(dir.repo, "data")
dir.out  <- file.path(dir.repo, "outputs")
xenium_dirs <- setNames(file.path(dir.data, "xenium", paste0("0071470_", sample_ids)), sample_ids)
hippocampus_lassos <- setNames(
  file.path(dir.data, "lassos", "hippocampus", paste0("hippocampus_", sample_ids, ".csv")), sample_ids)
optic_tract_lassos <- setNames(
  file.path(dir.data, "lassos", "optic_tract",
            paste0(c("GI-1", "GI-2", "GI-3", "GII-1", "GII-2", "GII-3"), "_coordinates.csv")), sample_ids)

if (file.exists("config.local.R")) {
  source("config.local.R", local = TRUE)   # evaluated where this file is being evaluated
  message("Applied local overrides from config.local.R")
}

# Sample lineup shared by every region, built from the settings above so
# run_region_pipeline() and every function signature are unchanged.
sample_base <- tibble::tibble(
  sample_id = sample_ids,
  group_id  = sample_groups,
  data_dir  = unname(xenium_dirs[sample_ids])
)

# ============================================================
# CONFIG -- Hippocampus: DE cell-type list (lassos: hippocampus_lassos above)
# ============================================================
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
# lassos (optic_tract_lassos in the SETTINGS block above; the author's
# files are GI-1 .. GII-3_coordinates.csv). The file-to-sample pairing was
# checked by overlaying each outline on its section
# (docs/optic_tract_new_lassos_overview.png). The lassos used before
# 2026-10-05 are listed, for the record only, in config.local.example.R.
# ============================================================
optic_tract_de_cell_types <- list(
  astro     = list(label = "Astrocyte", candidates = astro_candidates),
  microglia = list(label = "Microglia", candidates = microglia_candidates)
)

# ============================================================
# The analysis runs live in R/run_analyses.R.
# This file is setup only and is safe to source():
#   source("R/run_full_pipeline.R")
# ============================================================
