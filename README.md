# Xenium Spatial Transcriptomics: Cell-Type Composition and Differential Expression

Analysis pipeline for a Xenium in situ dataset (10x Xenium Prime Mouse 5K Pan
Tissue & Pathways panel) from 6 mouse brains: 3 TBI (group "gi") and 3 sham
(group "gii"). The same 6 mice are used for every brain region; only the
region-of-interest (ROI) lasso differs. Samples were processed and analyzed
under blinded labels (gi/gii) and unblinded on 2026-09-30. Abisambra Lab,
University of Florida.

Regions analyzed: Hippocampus and Optic Tract.

## Quick start

1. Open the repository folder in RStudio (`hippocampus-xenium-analysis.Rproj` is the root
   anchor; the setup file stops with a message if the working directory is not the root).
2. Install the packages if needed: Seurat, dplyr, tidyr, ggplot2, scales, sp, tibble, pheatmap,
   Matrix from CRAN and DESeq2 from Bioconductor. The setup file prints the exact install lines
   for anything missing.
3. Get the data (next section) and copy `config.local.example.R` to `config.local.R`, pointing
   `xenium_dirs`, `hippocampus_lassos` and `optic_tract_lassos` at your copies. `config.local.R`
   is gitignored.
4. `source("R/run_full_pipeline.R")`: loads packages, functions, panels and settings; runs
   nothing heavy. Then run `R/run_analyses.R` one block at a time (primary run of a region
   first; each block is heavy).

The first thing a run prints is the input check. `OK` on every row means it will not fail for a
missing file later; `MISSING` stops the run before any slow work and says where to get the input:

```
Input check
----------------------------------------------------------------------------------------------------
  input                      stage   need     status  path
  xenium_gi_1_b              0,1     required OK      /Volumes/Joe/.../0071470_gi_1_b
  lasso_hippocampus_gii_3_e  0,1,1c  required MISSING ~/Desktop/hippocampus_gii_3_e.csv
----------------------------------------------------------------------------------------------------
Missing input(s):
  lasso_hippocampus_gii_3_e  Xenium Explorer selection export (...) for gii_3_e, hippocampus. Request from the author; ...
```

## Getting the data

Nothing the pipeline reads ships with the repository. The inputs are unpublished lab data and
stay on the lab drive; request them from the Abisambra Lab.

| Input | Size | Used by | How to get it | Where the code looks by default |
|---|---|---|---|---|
| Six Xenium output bundles, one folder per mouse (`0071470_gi_1_b` ... `0071470_gii_3_e`), read by `Seurat::LoadXenium()`; the code also reads `cells.csv.gz` by name | several GB each | stages 0 and 1 of a primary run | by request from the lab | `data/xenium/0071470_<sample>/` (`xenium_dirs`) |
| Six hippocampus lasso CSVs (Xenium Explorer selection exports: two comment lines, then `Selection, X, Y, Class, Color`) | KB | stages 0, 1 and the area column | by request from the author | `data/lassos/hippocampus/hippocampus_<sample>.csv` (`hippocampus_lassos`) |
| Six optic tract lasso CSVs (redrawn 2026-10-05; `GI-1` ... `GII-3_coordinates.csv`) | KB | same | by request from the author | `data/lassos/optic_tract/<name>_coordinates.csv` (`optic_tract_lassos`) |

`data/` is gitignored. To use copies that already sit elsewhere (a shared drive, the lab drive),
set the three keys in `config.local.R` instead of copying the files; the original author's
layout is the commented example in `config.local.example.R`.

## Configuration

`R/run_full_pipeline.R` sets repo-relative defaults in its SETTINGS block, then sources
`config.local.R` if it exists (and says so), then builds the sample table. Keys:

| Key | Default | Purpose |
|---|---|---|
| `dir.data` | `file.path(dir.repo, "data")` | Folder holding the inputs that cannot ship |
| `xenium_dirs` | `file.path(dir.data, "xenium", "0071470_<sample>")`, named by sample id | The six Xenium bundles |
| `hippocampus_lassos` | `file.path(dir.data, "lassos", "hippocampus", "hippocampus_<sample>.csv")` | Hippocampus lasso files |
| `optic_tract_lassos` | `file.path(dir.data, "lassos", "optic_tract", "<GI-1..GII-3>_coordinates.csv")` | Optic tract lasso files |
| `dir.out` | `file.path(dir.repo, "outputs")` | Where runs write; the tracked deliverables stay in the repository by default |

Pipeline parameters (QC cut-offs, detection floor, confidence percentile, DE thresholds) are
function defaults and call-site values documented in `docs/review/workflow-review.md` §2.2; they
are the author's methods choices and are not configuration.

## Cell-type classification

Every cell is classified individually; no clustering is used. Each cell is
scored against a set of marker panels with Seurat::AddModuleScore(), scores
are z-scored within each panel, and the panel with the highest score
(argmax) sets the cell's identity. Astrocyte, Microglia and Oligodendrocyte
lineage cells are then split into states (Active/Homeostatic;
Mature/Precursor) using state-specific panels.

- Panel source: the identity panels were derived from the Allen Institute
  whole-mouse-brain taxonomy (Yao et al. 2023, Nature) and cross-checked gene
  by gene against the Xenium 5K panel. BAM and the state panels are
  literature-based curated deviations. Optic-tract immune panels:
  [TODO: source].
- Hippocampus uses 12 candidate panels: the 13 original panels minus Choroid
  plexus (excluded 2026-10-06 after driver-gene diagnostics). The Optic Tract uses 16: the 13
  plus six immune panels (Peripheral macrophage, T cell, Neutrophil,
  NK cell, B cell, Dendritic cell), minus four hippocampus- or
  tract-inappropriate panels (Hippocampal pyramidal, Dentate gyrus granule,
  Ependymal/ciliated, RGC). Reasons for each exclusion (single-gene driver
  false positives) are documented in R/marker_panels.R and
  protocols/classification_protocol.md.
- Detection floor: a panel is excluded from scoring in a sample if its best
  marker is detected in fewer than 0.5% of cells.
- Each cell records classification_margin (winning minus runner-up score)
  and runner_up.
- Purkinje (Pcp2) is a global negative control and should read about 0%.

Full protocol: protocols/classification_protocol.md.

## Region of interest

ROIs are lasso polygons drawn in Xenium Explorer and exported as CSV; cells
are kept by point-in-polygon. Optic Tract lassos were replaced on 2026-10-05
by new lassos drawn by a lab neuroscientist (one polygon per sample). Lasso
areas range from about 49,000 to 227,000 square micrometres, so cell counts
differ by sample; counts and cells per 1000 square micrometres are reported
with every run. Each primary run also saves lasso_overview.png, showing every
lasso drawn over its sample's cells, as a check of the file-to-sample mapping.

## Confidence filtering

Optic Tract. The primary analysis applies no confidence filter. Two
sensitivity analyses reassign cells below the pooled 15th-percentile
classification margin to Unclassified: (1) across all categories
(sens_global), and (2) only within categories where more than 15% of cells
fell below that threshold in the unfiltered run (sens_neartie: Choroid
plexus, Neutrophil, GABAergic interneuron, Peripheral macrophage, B cell,
BAM). The threshold and category list were derived from the margin
distribution pooled across all six samples; group labels were not used to
select either. A finding is considered robust to near-tie handling if the
direction of the TBI-versus-sham difference is unchanged across all three
runs; the size of any change is reported.

Hippocampus. Earlier versions applied a fixed margin > 1 rule to BAM only.
This was replaced (2026-10-06) by the same design as the Optic Tract: an
unfiltered primary analysis plus two sensitivity runs, using the region's own
pooled 15th-percentile margin (0.1745 with the 12-panel set; sens_neartie:
BAM, Ependymal/ciliated, Microglia (Active)). A fixed threshold is not
transferable across regions because the margin depends on the number of
competing panels.

## Differential expression (TBI vs sham)

Two methods per cell type, reported together: pseudobulk DESeq2 (mouse-level
replication, n = 3 per group) as the primary analysis, with a genome-wide
pass and a pre-specified candidate-panel pass; and Seurat FindMarkers
Wilcoxon on the full gene panel (cell-level replication) as a secondary,
exploratory analysis. Cell-level tests treat cells as replicates, which
inflates false positives (Squair et al., 2021); the pseudobulk analysis
avoids this at lower power. Full details:
protocols/differential_expression_protocol.md.

Composition differences are tested per category with Welch's t-test and
Benjamini-Hochberg correction (R/build_stats_table.R), with the sample as the
unit (n = 3 per group).

## Outputs

Every run writes into `outputs/<region>_<run>/` (tracked deliverables, regenerated in place so a
rerun shows in `git diff`; `.rds` checkpoints and their `.meta.txt` sidecars are gitignored):

- Per-cell (in the Seurat objects, not in files): atlas_identity, cell_type,
  classification_margin, runner_up, one idx_<Panel> score per panel.
- `cell_counts_per_sample.csv`: sample, n_cells, n_unclassified, pct_unclassified,
  pct_pass_qc20_10, group, region, run, lasso_area_um2, cells_per_1000um2. pct_pass_qc20_10 is
  the share of cells (already passing nCount >= 10 and nFeature >= 5) that would also pass a
  stricter nCount >= 20 and nFeature >= 10 cutoff.
- `composition_per_sample.csv` (counts and percentages per mouse and category),
  `composition_barchart.png` (mean and SEM bars with the six mice as points).
- `per_sample_depth.csv` (median transcripts and genes per cell, ROI cells after QC) and
  `panel_status_per_sample.csv` (which panels were scored in which mouse, the best gene and its
  detection rate, cells called; the 0.5% floor is applied per mouse).
- `umap_celltype_gi_vs_gii.png`, `umap_microglia_check_gi_vs_gii.png`.
- Per DE cell type: `pseudobulk_heatmap_<ct>.png`, `wilcox_heatmap_<ct>.png`, and under `de/`
  the tables `<ct>_res_full.csv` (genome-wide DESeq2), `<ct>_res_panel.csv` (candidate panel,
  BH within the panel), `<ct>_wilcox_res.csv` (FindMarkers, Bonferroni `p_val_adj`) and the
  diagnostics `<ct>_pvalue_hist.png`, `<ct>_dispersion.png`.
- `run_info.txt` (date, settings, confidence threshold, commit, cell counts, stage timings,
  R/Seurat/DESeq2 versions) and `session_info.txt`.
- `lasso_overview.png` in primary runs: every lasso drawn over its own section.
- At the top of `outputs/`: `<region>_near_tie_sensitivity.csv` (TBI minus sham per category
  under the three runs), `hippocampus_unfiltered_stats_table.csv` (Welch and BH), and the
  driver-gene diagnostics `hippocampus_{choroid_plexus,ependymal}_*`.
- Per-sample spatial cell-type maps are not written by the pipeline; `plot_spatial_celltype()`
  (R/plot_spatial_celltype.R) is an interactive helper for drawing one.
- Folders ending in `_withCP` hold the earlier 13-panel hippocampus runs (before Choroid plexus
  was excluded), kept for comparison; they and the driver-gene diagnostics depend on a
  checkpoint that is not regenerated by the current script (issue #9).

## Known limitations

- n = 3 per group: every mouse-level test (Welch t-test on composition,
  pseudobulk DESeq2) has very little power with three mice per group, and
  the cell-level Wilcoxon test treats correlated cells as independent;
  composition and DE results are exploratory. The Wilcoxon `p_val_adj` from
  Seurat `FindMarkers()` is Bonferroni-adjusted over all genes tested; the
  pseudobulk and composition adjustments are Benjamini-Hochberg.
- The smallest Optic Tract lasso (gii_3_e) has 156 cells; rare categories in
  that sample are unreliable.
- Astrocyte Active/Homeostatic (Gfap, Serpina3n) is compromised in white
  matter: optic nerve astrocytes have a high baseline of Gfap, and
  Serpina3n is largely oligodendrocyte-expressed in white matter. The
  buckets are best read as marker-high versus marker-low. The combined
  astrocyte total is reported alongside the split.
- Sensitivity to near-tie handling (Optic Tract; see
  outputs/optic_tract_near_tie_sensitivity.csv): the TBI-minus-sham difference
  in composition changed sign when near-tie cells were reassigned to
  Unclassified for Unclassified, Neutrophil, B cell, Choroid plexus and
  Peripheral macrophage (shifts of 2.5 to 6.7 percentage points). Group
  differences in these five categories are not robust and are not
  interpreted. All other categories kept their direction, and the largest
  (Mature oligodendrocyte, Homeostatic microglia, Homeostatic astrocyte,
  Precursor oligodendrocyte, Endothelial) were unchanged between runs.
  (Pericyte/vascular SMC changed sign in the all-category run by 0.2
  percentage points, which is negligible.)
- Sensitivity to near-tie handling (Hippocampus, 12 panels; see
  outputs/hippocampus_near_tie_sensitivity.csv): Ependymal/ciliated and
  Fibroblast/meningeal changed sign, and the Pericyte/vascular SMC
  difference fell to about zero in the all-category run. BAM kept its sign,
  but its difference fell from 1.29 to 0.08 percentage points, so the
  unfiltered BAM difference is attributed to near-tie calls. Hippocampal
  pyramidal, Dentate gyrus granule, Homeostatic astrocyte and Mature
  oligodendrocyte were stable across all three runs.
- Hippocampal Ependymal/ciliated calls are mixed: mostly Foxj1-driven, with
  some calls lining the ventricular edge of the lasso. The category is kept
  and reported as not robust.
- Marker spillover from neighbouring cells (segmentation) can contribute to
  signal for some genes.
- Sample quality: in the hippocampus, median transcripts per cell were
  comparable between groups (TBI 613-639, sham 597-726), so composition
  differences there are unlikely to reflect sample quality.

## Open items

- Optic Tract: check Neutrophil (Csf3r) and Choroid plexus for single-gene
  drivers, as was done for hippocampal Choroid plexus (R/check_driver_genes.R).
- Optic Tract: compare per-sample transcript depth between groups.
- UMAP coloured by sample, to check whether any island comes from a single
  mouse.
- gi_2_d BAM: 0 cells under the detection floor in an earlier run, 264 in
  the current run. The floor is calculated per mouse on the cells inside the
  lasso that passed QC (classification runs after the spatial restriction);
  `panel_status_per_sample.csv` now records the decision in every run. Whether
  a per-mouse floor is the right design is a question in issue #8 (M3).
- Robustness to the percentile choice (10th and 25th, alongside the 15th).
- Optic Tract composition bar charts: redraw with the corrected plotting code
  (fixed 2026-10-07) when the optic tract runs are next repeated (issue #9).
- Sources for the optic-tract immune panels.
- Exploratory astrocyte-only state analysis using a pre-specified,
  literature-derived gene set (sample as unit).

## Planned analysis: mixed-model differential expression (FLASH-MM)

**Status:** proposed, not yet implemented. For biostatistician review.

Current DE runs within every pipeline run, for astrocytes, microglia,
dentate gyrus granule cells and GABAergic interneurons (TBI vs sham):

- **Pseudobulk** (primary): counts summed per mouse, so the mouse is the
  unit of analysis (n = 3 vs 3). Statistically appropriate but low-powered.
- **Wilcoxon** (exploratory): cell-level test that treats cells as
  independent, which inflates significance because cells from the same
  mouse are correlated.

**Proposed:** a cell-level linear mixed model with mouse as a random
effect (e.g. FLASH-MM), for the same cell types and comparison, with BH
correction. This keeps cell-level resolution while accounting for
within-mouse correlation.

**Questions for review:** model specification (fixed and random effects,
covariates such as library size), normalisation, and whether to report
it alongside or instead of the Wilcoxon results.

## Repository structure

- R/run_full_pipeline.R: setup only (root-anchor check, packages, functions,
  marker panels, SETTINGS block and local overrides). Safe to source():
  source("R/run_full_pipeline.R"). Defines run_region_pipeline(), the full
  pipeline end to end (input check, classification, composition, UMAP, both
  DE analyses, run record).
- R/_setup.R: base R only; package check, input manifest and input check,
  commit reader for run_info.txt.
- config.local.example.R: every configuration key with its default and the
  original author's paths as commented examples; copy to config.local.R.
- R/run_analyses.R: the exact runs behind the reported results (primary and
  both sensitivity runs per region, near-tie list check, robustness and stats
  tables, driver-gene diagnostics). Run one block at a time; each run is
  heavy. Per-sample checkpoints (*_classified_prefilter.rds, with a .meta.txt
  sidecar recording panels, lasso and region) are gitignored.
- R/marker_panels.R: marker panel definitions and region-specific candidate sets
- R/classify_atlas_index.R: per-cell classification
- R/process_sample_template.R: per-sample load, QC, normalize, lasso, classify
- R/apply_confidence_filter.R: margin filter (fixed and pooled-percentile versions)
- R/build_cell_count_table.R: per-sample cell-count table
- R/build_composition_summary.R: composition percentages and group summaries
- R/build_stats_table.R: Welch t-test with Benjamini-Hochberg correction
- R/check_marker_spillover.R: marker spillover diagnostic
- R/check_driver_genes.R: single-gene driver diagnostics (per-gene detection,
  runner-up labels, spatial location of a category's cells)
- R/plot_lasso_overview.R: lasso outlines over each sample's cells
- R/plot_composition_barchart.R, R/plot_umap_celltype.R,
  R/plot_spatial_celltype.R, R/plot_de_heatmap.R: figures
- R/de_candidate_panels.R, R/run_pseudobulk_de.R, R/run_wilcox_de.R: DE analyses
- protocols/: written classification and differential expression protocols
- outputs/: generated figures and tables
- docs/review/: the audit documents (science review, workflow review,
  findings log for issue #6); CLAUDE.md: conventions for anyone working here

## For collaborators

- **Your paths are yours.** Tracked code holds no machine paths; `config.local.R` (gitignored)
  holds yours, and `config.local.example.R` documents every key. If you add an input, add it to
  `fn.data_manifest()` in R/_setup.R, to the config example and to "Getting the data".
- **Nothing scientific changes silently.** Marker panels, candidate gene lists, thresholds and
  tests are the author's decisions; open questions are tracked in issue #8.
- **Contributing.** Branch from an issue (`N-short-name`), commits split by type (code / docs /
  outputs / housekeeping / record) with subjects ending `(issue #N)` and bodies stating the
  problem, the change and the evidence; pull request with `Fixes #N`.
- **Review documents.** `docs/review/science-review.md` (what the pipeline does and why),
  `docs/review/workflow-review.md` (how to run and maintain it),
  `docs/review/issue6-audit-repo-log.md` (findings, decisions, change log).

## Reproducibility

- Tested environment: the author's macOS machine (R and package versions recorded in
  `session_info.txt` of each run from the audit branch onward; the versions behind the committed
  outputs predate that record). The analyst's Windows 11, R 4.5.3 machine ran the static checks
  (parse, anchor, package and input checks, table recomputations) without the data.
- No `set.seed()` is used: `AddModuleScore()` (seed 1), `RunPCA()` and `RunUMAP()` (seed 42) run
  on Seurat defaults, so a run repeats for a given Seurat/uwot version; UMAP coordinates and PNG
  bytes can differ between package versions and machines, tables should not.
- The committed tables reproduce from each other: the near-tie sensitivity table is byte-identical
  when recomputed from the three committed composition tables, and the Welch/BH table matches to
  1.6e-14. End-to-end regeneration on a second machine is recorded in
  `docs/review/workflow-review.md` §9 as it happens.
- To fingerprint the outputs before and after a run:
  `f <- list.files("outputs", recursive = TRUE, full.names = TRUE); write.csv(data.frame(file = f, md5 = unname(tools::md5sum(f))), "md5.csv", row.names = FALSE)`.

## Data

Raw Xenium output and lasso CSV files are not tracked in this repository (unpublished research
data). See "Getting the data" above for what they are, where the code looks for them and how to
point it at existing copies.