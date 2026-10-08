# Workflow review: hippocampus-xenium-analysis

## 0. How to read this document

Written last (Phase 9). Until then: this is the plain-language guide to running and maintaining
the pipeline; the science review (`docs/review/science-review.md`) says what it does and why; the
findings log (`docs/review/issue6-audit-repo-log.md`) is the audit record and holds every
decision. Line numbers refer to commit 869af28. Last updated 7 Oct 2026
(Phase 2 draft: §1 and §2 only).

## 1. The pipeline in one page

Six Xenium output bundles (one per mouse) and six region polygons go in. For each mouse the code
loads the bundle, drops near-empty cells, keeps the cells inside the polygon and gives every cell a
cell-type label by scoring it against small marker-gene sets. The six labelled samples are then
summarised (how many cells of each type per mouse, TBI versus sham), merged and drawn on a UMAP (a
two-dimensional map in which similar cells sit together), and, for a short list of cell types,
tested for genes that differ between TBI and sham in two ways: once with each mouse as one data
point (pseudobulk DESeq2) and once with each cell as one data point (Wilcoxon). Each region is run
three times: unfiltered, and twice with low-confidence labels set to "Unclassified". Everything is
driven by one function, `run_region_pipeline()`, and a run script that calls it six times.

| Stage | What goes in | What comes out | How long | Where in the code (file and main function) |
|---|---|---|---|---|
| 0 Lasso overview (primary runs only) | `cells.csv.gz` per mouse, six lasso CSVs | `lasso_overview.png` | pending | `R/plot_lasso_overview.R`, `plot_lasso_overview()` |
| 1 Load, QC, normalise, restrict, classify (per mouse) | Xenium bundle, lasso CSV, marker list | one labelled Seurat object per mouse; checkpoint `<sample>_classified_prefilter.rds` (ignored) | pending ("heavy", `R/run_analyses.R:3`) | `R/process_sample_template.R`, `process_hippocampus_sample()`; `R/classify_atlas_index.R`, `classify_atlas_index()` |
| 1b Confidence filter (sensitivity runs) | the six objects, percentile, category list | same objects with some `cell_type` set to Unclassified | pending | `R/apply_confidence_filter.R`, `apply_confidence_filter_pooled()` |
| 1c Cell-count table | the six objects, lasso areas | `cell_counts_per_sample.csv` | pending | `R/build_cell_count_table.R`, `build_cell_count_table()` |
| 2 Merge | the six objects | one merged object; `combined.rds` (ignored) | pending | `R/run_full_pipeline.R:190-193`, `merge()` + `JoinLayers()` |
| 3 Composition | the six objects | composition tables in memory; `composition_barchart.png` | pending | `R/build_composition_summary.R`, `build_composition_summary()`; `R/plot_composition_barchart.R` |
| 4 UMAP | merged object | `umap_celltype_gi_vs_gii.png`, `umap_microglia_check_gi_vs_gii.png` | pending | `R/run_full_pipeline.R:207-220`; `R/plot_umap_celltype.R` |
| 5 Differential expression (per DE cell type) | merged object, candidate panel | `pseudobulk_heatmap_<ct>.png`, `wilcox_heatmap_<ct>.png`; result tables in memory | pending | `R/run_pseudobulk_de.R`, `build_pseudobulk_dds()`; `R/run_wilcox_de.R`, `run_wilcox_degs()`; `R/plot_de_heatmap.R` |
| 6 Robustness and stats tables (run script) | the three runs' composition tables | `<region>_near_tie_sensitivity.csv`, `hippocampus_unfiltered_stats_table.csv` | seconds | `R/run_analyses.R:86-108`, `robustness_table()`; `R/build_stats_table.R`, `build_stats_table()` |
| 7 Driver-gene diagnostics (run script) | the 13-panel hippocampus checkpoint (`sample_objects_prefilter.rds`, ignored) | `hippocampus_<category>_gene_detection.csv`, `_runner_up.csv`, `_spatial.png` | pending | `R/run_analyses.R:118-130`; `R/check_driver_genes.R` |

Run times are "pending" because no run has been timed: the repository holds no log, timer or stated
duration beyond "each run is heavy" (README "Repository structure"). The baseline run on the
author's machine (Phase 4 hand-off) fills this column.

Three tracked output groups have no stage in this table because no tracked code writes them:
`composition_per_sample.csv` in each 12-panel hippocampus folder, the twelve DE tables in
`outputs/hippocampus_unfiltered/de/`, and `outputs/hippocampus_per_sample_depth.csv` (findings log,
R lens).

## 2. Setting up and running

### 2.1 Software

- **Language.** R; the version the author used is not recorded anywhere in the repository (no
  `sessionInfo()`, lock file or README line). The analyst's machine has R 4.5.3 on Windows 11.
- **Packages attached by the setup file** (`R/run_full_pipeline.R:42-50`): Seurat, dplyr, tidyr,
  ggplot2, scales, sp, DESeq2 (Bioconductor), tibble, pheatmap. Also attached when the setup file
  is sourced: `tidyverse` and Seurat again, by `R/check_marker_spillover.R:36-37`, which the setup
  file sources at line 66; `tidyverse` is not in the README's or the header's package list. Used by
  namespace without being attached: Matrix (`R/classify_atlas_index.R:12`), grid
  (`R/plot_de_heatmap.R:35`). `Seurat::LoadXenium()` needs the packages Seurat uses to read Xenium
  bundles (at least `arrow` for parquet in recent Seurat versions); which ones depends on the
  Seurat version, to be confirmed from the author's session.
- **System libraries.** None beyond what the packages themselves need; DESeq2 is installed through
  `BiocManager`.
- **What a missing package does today.** `source("R/run_full_pipeline.R")` stops at the first
  failing `library()` call with R's own "there is no package called" message; there is no package
  check that lists everything missing at once.

### 2.2 Configuration

There is no configuration file. Every machine-specific value is a literal in the setup file:

| Value | Where | What it is | Who it works for |
|---|---|---|---|
| `sample_base$data_dir` (six paths) | `R/run_full_pipeline.R:252-260` | Folder of each mouse's Xenium bundle on the lab drive, `/Volumes/Joe/Becca and Sakthi Files/Xenium_Data/Group_i/0071470_gi_1_b/` and five more | the author on macOS with the drive mounted |
| `hippocampus_lassos` (six paths) | `R/run_full_pipeline.R:265-272` | `~/Desktop/hippocampus_<sample>.csv` | the author |
| `optic_tract_lassos` (six paths) | `R/run_full_pipeline.R:306-313` | `~/Desktop/B6_Xenium_Opt/<GI-1 .. GII-3>_coordinates.csv` | the author |
| `optic_tract_lassos_old` (six paths) | `R/run_full_pipeline.R:315-322` | superseded lassos, unused | nobody (reference only) |
| working directory | `R/run_analyses.R:17` | `setwd("~/Documents/GitHub/hippocampus-xenium-analysis")` | the author |
| output folders | `R/run_analyses.R:31, 37, 48, 59, 65, 77` | `outputs/<region>_<run>` relative to the repository root | anyone, once the working directory is right |

Pipeline parameters with their defaults (all set in function signatures or call sites, none in
one place): QC cut-offs 10 transcripts / 5 genes (`R/process_sample_template.R:4`); detection
floor 0.5 % and `min_z = 0` (`R/classify_atlas_index.R:2`); confidence percentile 0.15 and the two
near-tie category lists (`R/run_analyses.R:39, 44, 50, 67, 72, 79`); pseudobulk minimum row sum 10
(`R/run_pseudobulk_de.R:17`); Wilcoxon `logfc.threshold = 0.25`, `min.pct = 0.1`
(`R/run_wilcox_de.R:27-28`); 30 principal components (`R/run_full_pipeline.R:209-210`).

The audit's fix design (findings log §4) proposes one gitignored `config.local.R` read after
repo-relative defaults, with the author's paths above kept as commented examples so their layout
keeps working.

### 2.3 Entry point and run modes

As found, the pipeline is interactive only:

1. Open RStudio in the repository (there is no tracked `.Rproj` file; the analyst's clone has an
   untracked one created by RStudio).
2. `source("R/run_full_pipeline.R")`: attaches packages, sources the 17 function files, defines
   the sample table, lasso maps and DE cell-type lists. Safe to source; runs nothing heavy
   (`R/run_full_pipeline.R:330-333`).
3. Run `R/run_analyses.R` one block at a time (its header says so, line 3). The primary run of a
   region must come first because the two sensitivity runs reuse its loaded samples
   (`prefilter_objects = <primary>$sample_objects_prefilter`, lines 42, 53, 70, 82).

`Rscript R/run_analyses.R` would not work as an entry point today: line 17 sets the working
directory to the author's home folder, and line 118 reads a checkpoint file that only exists after
a 13-panel run that the script no longer performs. Figures are also `print()`ed to the plot pane
(`R/run_full_pipeline.R:202, 213, 218`), which under `Rscript` lands in `Rplots.pdf`.

As found, `source("R/run_full_pipeline.R")` itself fails at commit 869af28: `R/build_stats_table.R`
does not parse (a misplaced parenthesis in the `pick(everything())` edit of commit 869af28, line
31), so line 65 of the setup file stops with `unexpected '}'`. Verified on the analyst's machine
with the skill's parse-only sweep (18 of 19 files parse). Findings log B lens; the fix is one
character and does not change any output (the committed stats table recomputes from the committed
composition table to 1.6e-14).

### 2.4 The three run targets

| Target | Machine and layout | What it needs today | What it will need after the audit |
|---|---|---|---|
| The analyst (holden-harris) | Windows 11, R 4.5.3, no Xenium data, Seurat and DESeq2 not installed, `gh` not installed | the six bundles and twelve lasso CSVs copied locally, six `data_dir` edits and twelve lasso-path edits in the setup file, the `setwd()` line removed, Seurat/DESeq2 installed | a `config.local.R` pointing at the data; the package check's install line |
| The author (rebeccasakko) | macOS, lab drive mounted at `/Volumes/Joe`, lassos on the Desktop, repository under `~/Documents/GitHub/` | nothing beyond the parse fix | nothing: the defaults keep matching, or a two-line `config.local.R` with the Desktop paths |
| A fresh clone | any OS, data placed where the documentation says | undocumented: no manifest says what the bundles must contain or where to put them | the README "Getting the data" table, `config.local.R`, and a missing-input table printed before any slow work |

Section 9 (verification record) holds the evidence for each target once runs have happened; until
then its cells read "pending: run on the author's machine".

### 2.5 What a first run prints

Today: nothing before the slow work. With the data missing, the first message is R's own error
from `LoadXenium()` inside the first sample's stage 1, after the packages have loaded; with a lasso
file missing, `read.csv()`'s "cannot open file" after that sample's bundle has been loaded and
normalised. The audit proposes an input-check table (findings log §4) that reports every bundle
and lasso file as OK or MISSING before stage 1.

## 3. Data inputs

To be filled in Phase 5.

## 4. Pipeline and data transformations

To be filled in Phase 5.

## 5. Statistical models and estimators

To be filled in Phase 5.

## 6. Randomness and reproducibility

To be filled in Phase 5.

## 7. Outputs

To be filled in Phase 5.

## 8. Run times and efficiency

To be filled in Phases 4 and 5.

## 9. Verification record

### 9.1 As-found baseline

**Result: pending, run on the author's machine** (hand-off block in the findings log §3, posted on
PR #7 on 2026-10-07). The analyst's machine has neither the data nor Seurat and DESeq2, so no run
can be made here.

Baseline notes, written before the run:

- The baseline is run on branch `6-audit-repo`, which differs from `main` (869af28) in exactly one
  character of code: the misplaced parenthesis in `R/build_stats_table.R:31` (finding B1, decision
  4). Without it `source("R/run_full_pipeline.R")` stops at line 65 with `unexpected '}'`, so no run
  of any kind is possible at 869af28. Evidence that the fix changes nothing: with the corrected
  function, the committed `outputs/hippocampus_unfiltered/composition_per_sample.csv` gives a stats
  table whose 15 p-values and BH-adjusted p-values match the committed
  `outputs/hippocampus_unfiltered_stats_table.csv` to 1.6e-14 and 5.1e-15 in the same row order;
  the CSV is not byte-identical only because the last printed digit of some doubles and the line
  endings differ between the analyst's Windows R 4.5.3 and the author's machine.
- Smoke checks on 2026-10-07 (analyst's machine, parse-only): root anchor absent (P3); 19 of 19
  files parse after the fix (18 of 19 before it, B1); the package check and the input check cannot
  run here (no Seurat/DESeq2, no data) and are handed to the author with the baseline.
- The hand-off asks for the hippocampus block first (three runs plus the robustness and stats
  tables), the optic tract block if time allows, and the diagnostics block only if the 13-panel
  checkpoint still exists on the author's disk (R4). Fingerprints of `outputs/` are taken before
  and after with `tools::md5sum()`; the comparison table goes here when the CSVs arrive.

### 9.3 As-left verification

Pending (Phase 8).

## 10. Appendices

To be filled in Phases 5 and 8.
