# ============================================================
# run_analyses.R -- the analysis runs behind the reported results.
# Run ONE block at a time (each is heavy). In each region the PRIMARY
# run must come first: the sensitivity runs reuse its loaded samples
# ($sample_objects_prefilter), so they don't reload raw data.
#
# Design (both regions):
#   primary       -- unfiltered (no confidence filter)
#   sens_global   -- cells below the region's pooled 15th-percentile
#                    classification_margin -> Unclassified, all categories
#   sens_neartie  -- same threshold, only in categories where >15% of
#                    cells fall below it in the primary run
# Thresholds and near-tie lists were derived from pooled margins
# (all 6 samples) without group labels.
# ============================================================

# Working directory: open the repository folder in RStudio (hippocampus-xenium-analysis.Rproj
# is the root anchor); R/run_full_pipeline.R stops with a message if it is not the root.
source("R/run_full_pipeline.R")   # setup only: functions, panels, lassos

# ---- Hippocampus (12 panels; Choroid plexus excluded 2026-10-06) ----
# Choroid plexus removed after driver-gene diagnostics (see Diagnostics
# section below). The earlier 13-panel runs are kept for comparison in
# outputs/hippocampus_*_withCP/ (threshold 0.186; near-tie: BAM 73.0%,
# Ependymal/ciliated 42.3%, Microglia (Active) 20.5%; Choroid plexus 14.9%).
# 12-panel runs: pooled 15th-percentile margin = 0.1745; near-tie list
# unchanged: BAM 73.8%, Ependymal/ciliated 42.3%, Microglia (Active) 22.7%
# (next highest: Fibroblast/meningeal 12.8%).

hippocampus_results <- run_region_pipeline(
  sample_base, hippocampus_lassos, "Hippocampus",
  hippocampus_de_cell_types, file.path(dir.out, "hippocampus_unfiltered"),
  marker_list = hippocampus_markers,
  run_label = "primary")

hippocampus_sens_global <- run_region_pipeline(
  sample_base, hippocampus_lassos, "Hippocampus",
  hippocampus_de_cell_types, file.path(dir.out, "hippocampus_sens_global"),
  marker_list = hippocampus_markers,
  confidence_percentile = 0.15,
  confidence_categories = NULL,
  run_label = "sens_global",
  prefilter_objects = hippocampus_results$sample_objects_prefilter)

# Near-tie categories: more than 15% of a category's cells below the pooled
# 15th-percentile margin in the primary run. The recorded list (12-panel run,
# 2026-10-06) is used; the same rule is now computed from the primary run and
# any difference is reported (audit finding R6, decision 10).
hippocampus_near_tie_cats <- c("BAM", "Ependymal/ciliated", "Microglia (Active)")
hippocampus_near_tie_computed <- near_tie_categories(hippocampus_results$sample_objects_prefilter, 0.15, 0.15)
if (!setequal(hippocampus_near_tie_computed, hippocampus_near_tie_cats))
  warning("Hippocampus: computed near-tie list (", paste(hippocampus_near_tie_computed, collapse = ", "),
          ") differs from the recorded list (", paste(hippocampus_near_tie_cats, collapse = ", "),
          "); the recorded list is used. See issue #6, finding R6.", call. = FALSE)

hippocampus_sens_neartie <- run_region_pipeline(
  sample_base, hippocampus_lassos, "Hippocampus",
  hippocampus_de_cell_types, file.path(dir.out, "hippocampus_sens_neartie"),
  marker_list = hippocampus_markers,
  confidence_percentile = 0.15,
  confidence_categories = hippocampus_near_tie_cats,
  run_label = "sens_neartie",
  prefilter_objects = hippocampus_results$sample_objects_prefilter)

# ---- Optic tract (16 panels; redrawn lassos, 2026-10-05) ----

optic_tract_results <- run_region_pipeline(
  sample_base, optic_tract_lassos, "Optic Tract",
  optic_tract_de_cell_types, file.path(dir.out, "optic_tract_newlasso_unfiltered"),
  marker_list = optic_tract_markers,
  run_label = "primary")

optic_tract_sens_global <- run_region_pipeline(
  sample_base, optic_tract_lassos, "Optic Tract",
  optic_tract_de_cell_types, file.path(dir.out, "optic_tract_newlasso_sens_global"),
  marker_list = optic_tract_markers,
  confidence_percentile = 0.15,
  confidence_categories = NULL,
  run_label = "sens_global",
  prefilter_objects = optic_tract_results$sample_objects_prefilter)

optic_tract_near_tie_cats <- c("Choroid plexus", "Neutrophil", "GABAergic interneuron",
                               "Peripheral macrophage", "B cell", "BAM")
optic_tract_near_tie_computed <- near_tie_categories(optic_tract_results$sample_objects_prefilter, 0.15, 0.15)
if (!setequal(optic_tract_near_tie_computed, optic_tract_near_tie_cats))
  warning("Optic tract: computed near-tie list (", paste(optic_tract_near_tie_computed, collapse = ", "),
          ") differs from the recorded list (", paste(optic_tract_near_tie_cats, collapse = ", "),
          "); the recorded list is used. See issue #6, finding R6.", call. = FALSE)

optic_tract_sens_neartie <- run_region_pipeline(
  sample_base, optic_tract_lassos, "Optic Tract",
  optic_tract_de_cell_types, file.path(dir.out, "optic_tract_newlasso_sens_neartie"),
  marker_list = optic_tract_markers,
  confidence_percentile = 0.15,
  confidence_categories = optic_tract_near_tie_cats,
  run_label = "sens_neartie",
  prefilter_objects = optic_tract_results$sample_objects_prefilter)

# ---- Robustness tables: TBI (gi) minus sham (gii), percentage points ----
# "flips" = direction changes in at least one sensitivity run. A category
# absent from one of the three runs gives NA in that column; such rows are
# reported with a message rather than sorted away silently (audit finding B8).
robustness_table <- function(primary, global, neartie) {
  diff_tbl <- function(res) {
    cd <- res$composition$composition_df
    a <- tapply(cd$pct, list(cd$cell_type, cd$group_id), mean)
    round(a[, "gi"] - a[, "gii"], 2)
  }
  d1 <- diff_tbl(primary); d2 <- diff_tbl(global); d3 <- diff_tbl(neartie)
  cats <- union(union(names(d1), names(d2)), names(d3))
  out <- data.frame(cell_type = cats, unfiltered = d1[cats], global = d2[cats],
                    neartie = d3[cats], row.names = NULL)
  sign_flip <- function(a, b) !is.na(a) & !is.na(b) & sign(a) != sign(b)
  out$flips <- with(out, sign_flip(global, unfiltered) | sign_flip(neartie, unfiltered))
  out$max_shift <- with(out, pmax(abs(global - unfiltered), abs(neartie - unfiltered), na.rm = TRUE))
  incomplete <- out$cell_type[is.na(out$unfiltered) | is.na(out$global) | is.na(out$neartie)]
  if (length(incomplete) > 0)
    message("robustness_table: category absent from at least one run (NA in the table): ",
            paste(incomplete, collapse = ", "))
  out[order(-out$max_shift), ]
}

write.csv(robustness_table(hippocampus_results, hippocampus_sens_global, hippocampus_sens_neartie),
          file.path(dir.out, "hippocampus_near_tie_sensitivity.csv"), row.names = FALSE)
write.csv(robustness_table(optic_tract_results, optic_tract_sens_global, optic_tract_sens_neartie),
          file.path(dir.out, "optic_tract_near_tie_sensitivity.csv"), row.names = FALSE)

# ---- Composition stats: Welch t-test + BH, sample as the unit ----
write.csv(build_stats_table(hippocampus_results$composition$composition_df),
          file.path(dir.out, "hippocampus_unfiltered_stats_table.csv"), row.names = FALSE)

# ---- Diagnostics: single-gene driver checks, hippocampus (2026-10-06) ----
# Choroid plexus: Kl alone (66-98% of calls; Folr1/Otx2 ~1%), ~74% neuronal
#   runner-up, median margin 2.41, cells lie along pyramidal/DG layers ->
#   Kl-positive neurons, not choroid plexus. Excluded from the hippocampus panel.
# Ependymal/ciliated: mostly Foxj1-driven, but some calls line the outer
#   (ventricular) edge -> mixed real edge cells + scattered noise. Kept;
#   reported as not robust.
# Uses the 13-panel run (before Choroid plexus was excluded). That run is not
# produced by this script any more, so its checkpoint exists only on the disk
# where it was made; see issue #9 (audit finding R4) for regenerating it.
hc_ckpt <- file.path(dir.out, "hippocampus_unfiltered_withCP", "sample_objects_prefilter.rds")
if (!file.exists(hc_ckpt)) {
  message("Diagnostics skipped: 13-panel checkpoint not found at ", hc_ckpt, " (see issue #9).")
} else {
  hc_objs <- readRDS(hc_ckpt)

  write.csv(driver_check(hc_objs, "Choroid plexus", c("Folr1", "Kl", "Otx2")),
            file.path(dir.out, "hippocampus_choroid_plexus_gene_detection.csv"), row.names = FALSE)
  write.csv(runner_up_summary(hc_objs, "Choroid plexus"),
            file.path(dir.out, "hippocampus_choroid_plexus_runner_up.csv"), row.names = FALSE)
  plot_category_spatial(hc_objs, "Choroid plexus", file.path(dir.out, "hippocampus_choroid_plexus_spatial.png"))

  write.csv(driver_check(hc_objs, "Ependymal/ciliated", c("Foxj1", "Pifo", "Dnah12")),
            file.path(dir.out, "hippocampus_ependymal_gene_detection.csv"), row.names = FALSE)
  write.csv(runner_up_summary(hc_objs, "Ependymal/ciliated"),
            file.path(dir.out, "hippocampus_ependymal_runner_up.csv"), row.names = FALSE)
  plot_category_spatial(hc_objs, "Ependymal/ciliated", file.path(dir.out, "hippocampus_ependymal_spatial.png"))
}
