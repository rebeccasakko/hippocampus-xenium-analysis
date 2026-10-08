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

setwd("~/Documents/GitHub/hippocampus-xenium-analysis")
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
  hippocampus_de_cell_types, "outputs/hippocampus_unfiltered",
  marker_list = hippocampus_markers,
  run_label = "primary")

hippocampus_sens_global <- run_region_pipeline(
  sample_base, hippocampus_lassos, "Hippocampus",
  hippocampus_de_cell_types, "outputs/hippocampus_sens_global",
  marker_list = hippocampus_markers,
  confidence_percentile = 0.15,
  confidence_categories = NULL,
  run_label = "sens_global",
  prefilter_objects = hippocampus_results$sample_objects_prefilter)

hippocampus_near_tie_cats <- c("BAM", "Ependymal/ciliated", "Microglia (Active)")

hippocampus_sens_neartie <- run_region_pipeline(
  sample_base, hippocampus_lassos, "Hippocampus",
  hippocampus_de_cell_types, "outputs/hippocampus_sens_neartie",
  marker_list = hippocampus_markers,
  confidence_percentile = 0.15,
  confidence_categories = hippocampus_near_tie_cats,
  run_label = "sens_neartie",
  prefilter_objects = hippocampus_results$sample_objects_prefilter)

# ---- Optic tract (15 panels; redrawn lassos, 2026-10-05) ----
# Neutrophil excluded 2026-10-08 after driver-gene diagnostics (see
# Diagnostics section below). The earlier 16-panel runs are kept for
# comparison in outputs/optic_tract_newlasso_*_withNeut/.
# 15-panel runs: pooled 15th-percentile margin = 0.0257; near-tie (% below):
# Choroid plexus 54.5, GABAergic interneuron 27.4, Pericyte/vascular SMC 22.0,
# B cell 21.7, Peripheral macrophage 19.3 (BAM 14.5 = borderline, not included).

optic_tract_results <- run_region_pipeline(
  sample_base, optic_tract_lassos, "Optic Tract",
  optic_tract_de_cell_types, "outputs/optic_tract_newlasso_unfiltered",
  marker_list = optic_tract_markers,
  run_label = "primary")

optic_tract_sens_global <- run_region_pipeline(
  sample_base, optic_tract_lassos, "Optic Tract",
  optic_tract_de_cell_types, "outputs/optic_tract_newlasso_sens_global",
  marker_list = optic_tract_markers,
  confidence_percentile = 0.15,
  confidence_categories = NULL,
  run_label = "sens_global",
  prefilter_objects = optic_tract_results$sample_objects_prefilter)

optic_tract_near_tie_cats <- c("Choroid plexus", "GABAergic interneuron", "Pericyte/vascular SMC",
                               "B cell", "Peripheral macrophage")

optic_tract_sens_neartie <- run_region_pipeline(
  sample_base, optic_tract_lassos, "Optic Tract",
  optic_tract_de_cell_types, "outputs/optic_tract_newlasso_sens_neartie",
  marker_list = optic_tract_markers,
  confidence_percentile = 0.15,
  confidence_categories = optic_tract_near_tie_cats,
  run_label = "sens_neartie",
  prefilter_objects = optic_tract_results$sample_objects_prefilter)

# ---- Robustness tables: TBI (gi) minus sham (gii), percentage points ----
# "flips" = direction changes in at least one sensitivity run.
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
  out$flips <- with(out, sign(global) != sign(unfiltered) | sign(neartie) != sign(unfiltered))
  out$max_shift <- with(out, pmax(abs(global - unfiltered), abs(neartie - unfiltered)))
  out[order(-out$max_shift), ]
}

write.csv(robustness_table(hippocampus_results, hippocampus_sens_global, hippocampus_sens_neartie),
          "outputs/hippocampus_near_tie_sensitivity.csv", row.names = FALSE)
write.csv(robustness_table(optic_tract_results, optic_tract_sens_global, optic_tract_sens_neartie),
          "outputs/optic_tract_near_tie_sensitivity.csv", row.names = FALSE)

# ---- Composition stats: Welch t-test + BH, sample as the unit ----
write.csv(build_stats_table(hippocampus_results$composition$composition_df),
          "outputs/hippocampus_unfiltered_stats_table.csv", row.names = FALSE)

# ---- Diagnostics: single-gene driver checks, hippocampus (2026-10-06) ----
# Choroid plexus: Kl alone (66-98% of calls; Folr1/Otx2 ~1%), ~74% neuronal
#   runner-up, median margin 2.41, cells lie along pyramidal/DG layers ->
#   Kl-positive neurons, not choroid plexus. Excluded from the hippocampus panel.
# Ependymal/ciliated: mostly Foxj1-driven, but some calls line the outer
#   (ventricular) edge -> mixed real edge cells + scattered noise. Kept;
#   reported as not robust.
# Uses the 13-panel run (before Choroid plexus was excluded):
hc_objs <- readRDS("outputs/hippocampus_unfiltered_withCP/sample_objects_prefilter.rds")

write.csv(driver_check(hc_objs, "Choroid plexus", c("Folr1", "Kl", "Otx2")),
          "outputs/hippocampus_choroid_plexus_gene_detection.csv", row.names = FALSE)
write.csv(runner_up_summary(hc_objs, "Choroid plexus"),
          "outputs/hippocampus_choroid_plexus_runner_up.csv", row.names = FALSE)
plot_category_spatial(hc_objs, "Choroid plexus", "outputs/hippocampus_choroid_plexus_spatial.png")

write.csv(driver_check(hc_objs, "Ependymal/ciliated", c("Foxj1", "Pifo", "Dnah12")),
          "outputs/hippocampus_ependymal_gene_detection.csv", row.names = FALSE)
write.csv(runner_up_summary(hc_objs, "Ependymal/ciliated"),
          "outputs/hippocampus_ependymal_runner_up.csv", row.names = FALSE)
plot_category_spatial(hc_objs, "Ependymal/ciliated", "outputs/hippocampus_ependymal_spatial.png")

# ---- Diagnostics: single-gene driver checks, optic tract (2026-10-08) ----
# Neutrophil: Csf3r alone (96-100% of calls in 5/6 mice; Ly6g 0% in all
#   mice), median margin 0.21, ~half myeloid runner-up, scattered evenly;
#   Csf3r loads with microglia/macrophage genes on PC1 -> Csf3r-positive
#   myeloid cells, not neutrophils. Excluded from the optic tract panel.
# Choroid plexus: Folr1 + Otx2 (Kl low), median margin 0.36, mixed
#   runner-ups; cluster at one tip of the lasso + scattered calls -> mixed
#   real and noisy signal. Kept; reported as not robust.
# Uses the 16-panel run (before Neutrophil was excluded):
ot_objs <- lapply(c("gi_1_b", "gi_2_d", "gi_3_f", "gii_1_a", "gii_2_c", "gii_3_e"), function(s)
  readRDS(file.path("outputs/optic_tract_newlasso_unfiltered_withNeut", paste0(s, "_classified_prefilter.rds"))))
names(ot_objs) <- c("gi_1_b", "gi_2_d", "gi_3_f", "gii_1_a", "gii_2_c", "gii_3_e")
ctx <- "Oligodendrocyte (Mature)"

write.csv(driver_check(ot_objs, "Neutrophil", c("Ly6g", "Csf3r", "Cxcr2")),
          "outputs/optic_tract_neutrophil_gene_detection.csv", row.names = FALSE)
write.csv(runner_up_summary(ot_objs, "Neutrophil"),
          "outputs/optic_tract_neutrophil_runner_up.csv", row.names = FALSE)
plot_category_spatial(ot_objs, "Neutrophil", "outputs/optic_tract_neutrophil_spatial.png", context = ctx)

write.csv(driver_check(ot_objs, "Choroid plexus", c("Folr1", "Kl", "Otx2")),
          "outputs/optic_tract_choroid_plexus_gene_detection.csv", row.names = FALSE)
write.csv(runner_up_summary(ot_objs, "Choroid plexus"),
          "outputs/optic_tract_choroid_plexus_runner_up.csv", row.names = FALSE)
plot_category_spatial(ot_objs, "Choroid plexus", "outputs/optic_tract_choroid_plexus_spatial.png", context = ctx)