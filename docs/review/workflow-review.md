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

Nothing the pipeline reads ships with the repository, and there is no manifest, so this section is
the first written record of the inputs. Sizes and the exact folder contents are pending the author's
answer (findings log §8, question for the author); everything else below is read from the code.

### 3.1 Xenium output bundles (six, one per mouse)

What it is: the folder the 10x Xenium instrument and analysis software produce for one tissue
section. It holds the cell-by-gene count matrix, a table of cell centroids (`cells.csv.gz`) and
segmentation and transcript files. `Seurat::LoadXenium(data_dir, fov = "fov")`
(`R/process_sample_template.R:7`) reads the matrix and the centroids, and depending on the Seurat
version also the segmentation and transcript files; `R/plot_lasso_overview.R:6` reads only
`cells.csv.gz`. Who provides it: the lab; it is unpublished research data and cannot be
redistributed (README "Data"). Where it sits: `/Volumes/Joe/Becca and Sakthi Files/Xenium_Data/
Group_i/0071470_gi_1_b/` and five siblings (`R/run_full_pipeline.R:254-259`), typed into
`sample_base$data_dir`. Which stages read it: stage 0 and stage 1 (primary runs only; sensitivity
runs reuse the loaded objects). What happens without it: `LoadXenium()` fails inside stage 1 with
R's own file error, after the packages have loaded and after any earlier sample has been processed;
nothing checks the six folders up front. Vintage: unknown; the Xenium software version is not
recorded anywhere (finding R5). Size: pending (a 5K-panel bundle is typically several gigabytes
per section).

### 3.2 Lasso polygons (six per region)

What it is: the hand-drawn outline of the region of interest, exported from Xenium Explorer as a
CSV with two comment lines followed by the columns `Selection, X, Y, Class, Color`
(`R/run_full_pipeline.R:71-72`); X and Y are in micrometres in the same frame as the cell
centroids. Who provides it: the author (hippocampus) and a lab neuroscientist (optic tract,
2026-10-05). Where it sits: `~/Desktop/hippocampus_<sample>.csv` and
`~/Desktop/B6_Xenium_Opt/<GI-1 .. GII-3>_coordinates.csv` on the author's Mac
(`R/run_full_pipeline.R:265-313`). Not tracked, although each is a small text file; why is
unexplained (findings log §8). Which stages read it: stage 0 (overview), stage 1 (membership test,
`R/process_sample_template.R:11-13`) and stage 1c (area, `R/run_full_pipeline.R:73-79`). What
happens without it: `read.csv()` fails after that sample's bundle has been loaded and normalised.
Vintage: the optic tract set was replaced on 2026-10-05; the old set is still listed at
`R/run_full_pipeline.R:315-322` and unused. The membership test treats the whole file as one
polygon; the area function splits by `Selection` (finding B2); the hand-off asks the author how many
selections each file holds.

### 3.3 Marker panels and candidate gene lists

Both ship with the code: `R/marker_panels.R` (20 identity panels, three state pairs, two
region-specific subsets) and `R/de_candidate_panels.R` (four DE candidate lists). They are inputs
in the scientific sense (a change to either changes every output) but need no download. Their
provenance is in the science review §4.

### 3.4 Inputs created by earlier runs

- `outputs/<run>/<sample>_classified_prefilter.rds` and `outputs/<run>/sample_objects_prefilter.rds`:
  per-run checkpoints (ignored by `.gitignore:48`) that let an interrupted primary run resume and
  let the sensitivity runs skip stage 1. Keyed on the output folder and the sample id only (finding
  R7).
- `outputs/hippocampus_unfiltered_withCP/sample_objects_prefilter.rds`: the 13-panel checkpoint
  the diagnostics block reads (`R/run_analyses.R:118`); no current code path creates it (finding R4).

### 3.5 Summary table

| Input | Size | Used by | Source | How a user gets it |
|---|---|---|---|---|
| Six Xenium bundles | pending (GB each) | stages 0, 1 | lab instrument run `0071470` | by request from the lab; place anywhere and point `sample_base$data_dir` at them (after the audit: `config.local.R`) |
| Six hippocampus lasso CSVs | KB each | stages 0, 1, 1c | author, Xenium Explorer | by request from the author |
| Six optic tract lasso CSVs | KB each | stages 0, 1, 1c | lab neuroscientist, 2026-10-05 | by request from the author |
| Marker panels, candidate lists | in `R/` | stage 1, stage 5 | literature (science review §4) | ships with the repository |
| Per-run checkpoints (`*.rds`) | pending | stage 1 skip; diagnostics | produced by a primary run | regenerate by running the primary block |

## 4. Pipeline and data transformations

### 4.1 Stage 0: lasso overview

Reads every mouse's `cells.csv.gz` and every lasso file, draws the outline over the centroids, one
panel per mouse, y axis reversed to match Xenium Explorer, and saves `lasso_overview.png`
(`R/plot_lasso_overview.R`). Purpose: a visual check that each lasso file is paired with the right
section. Runs only when `prefilter_objects` is `NULL` (`R/run_full_pipeline.R:128-132`). Failure
behaviour: a missing file stops the run here, before any heavy work, which is the one place the
pipeline fails early.

### 4.2 Stage 1: per mouse, load to classify

In order (`R/process_sample_template.R`): `LoadXenium()`; keep cells with `nCount_Xenium >= 10`
and `nFeature_Xenium >= 5`; `NormalizeData()` (log-normalise, Seurat default scale factor 10,000);
read the lasso; `point.in.polygon()` on the centroids from `GetTissueCoordinates(image = "fov")`;
keep cells with a result above 0 (inside or on the edge); tag `sample_id`, `group_id`, `region`;
free the full object; classify. Classification (`R/classify_atlas_index.R`) runs on the
lasso-restricted, QC-passed object, so the 0.5 % detection floor and the within-sample z-scores are
computed on the cells inside the region, not on the whole section (this answers the README "Open
items" question about which cells the floor is calculated on). Output: a Seurat object with
`atlas_identity`, `cell_type`, `classification_margin`, `runner_up` and one `idx_<Panel>` score
column per panel. A checkpoint is saved per mouse (`R/run_full_pipeline.R:151`) and the list of
six is saved again at line 158.

Invariants that should hold and where they are checked: every cell has exactly one
`atlas_identity` (by construction); `cell_type` equals `atlas_identity` except for the three
refined lineages (by construction); panels excluded by the floor appear as a `warning()` only
(`R/classify_atlas_index.R:15`), nowhere in an output file (finding B6); the Purkinje negative
control percentage is printed to the console (lines 65-70), not saved.

Failure behaviour: a missing bundle or lasso stops the run inside this stage with R's own message;
a panel with no gene on the Xenium panel warns and is skipped; a sample with zero cells inside the
lasso would reach `AddModuleScore()` with an empty object and fail with a Seurat error.

### 4.3 Stage 1b: confidence filter (sensitivity runs)

Pools `classification_margin` over the six objects, takes the requested percentile as the
threshold, and sets `cell_type` (not `atlas_identity`) to "Unclassified" for cells below it, in all
categories or in the listed ones (`R/apply_confidence_filter.R`). One `message()` reports the
threshold, the number of cells pooled and the number reclassified; that line is the only record of
the threshold value and lives in the console (the README quotes 0.1745 for the 12-panel
hippocampus from such a line). Invariant: total cell counts are unchanged by the filter (verified on
the committed tables: `n_cells` identical across the three hippocampus runs).

### 4.4 Stage 1c: cell-count table

Per mouse: cells, Unclassified count and percentage, the 20/10 QC pass rate, group, region, run
label, lasso area and cells per 1,000 square micrometres (`R/build_cell_count_table.R`). Written to
`cell_counts_per_sample.csv` and printed. The group is inferred from the sample name (`grepl("gii")`),
not from `group_id`.

### 4.5 Stage 2: merge

`merge()` of the six objects with `add.cell.ids`, then `JoinLayers()` to collapse Seurat v5's
per-sample layers; saved as `combined.rds` (ignored). Nothing is recomputed; the per-sample
normalised data and metadata are carried over.

### 4.6 Stage 3: composition

Counts per `cell_type` per mouse over the union of categories seen in any mouse, zero-filled;
percentage of that mouse's cells; group mean and SEM; categories ordered by total mean percentage
(`R/build_composition_summary.R`). Bar chart with per-mouse points, filled for TBI and hollow for
sham, dodged into the bar of their group (`R/plot_composition_barchart.R`). Invariant: percentages
sum to 100 within each mouse (verified on the committed tables: exactly 100 in all six). The
composition tables are returned in memory only; the committed `composition_per_sample.csv` files
were written interactively (finding R1).

### 4.7 Stage 4: UMAP

`FindVariableFeatures()` (default 2,000 features), `ScaleData()`, `RunPCA(npcs = 30)`,
`RunUMAP(dims = 1:30)` on the merged object; two `DimPlot` figures split by group. Display only.

### 4.8 Stage 5: differential expression

For each entry of the region's DE list: subset the merged object to cells whose `atlas_identity`
equals the label; sum counts per mouse (`AggregateExpression`, `slot = "counts"`); build a
`DESeqDataSet` with design `~ group`, drop genes with fewer than 10 counts in total, run `DESeq()`;
results for `gi` vs `gii`; the candidate-panel subset re-adjusted within the panel; then
`FindMarkers()` Wilcoxon on the same cells with `group_id` as identity. Heatmaps for the panel and
for the Wilcoxon genes at `p_val_adj < 0.05`. Failure behaviour: a cell type with no cells in one
mouse gives a pseudobulk matrix with five columns and DESeq2 still runs; a panel with no gene
surviving the count filter skips the heatmap with a message. The four result tables per cell type
stay in memory (returned in `de_results`); the committed `de/*.csv` were written interactively
(finding R2). Note that `FindMarkers()` adjusts p-values by Bonferroni over all genes tested, so
`p_val_adj` is not a BH value (finding D9).

### 4.9 Stage 6: robustness and stats tables (run script)

`robustness_table()` takes the three runs' composition tables, computes mean TBI minus mean sham
percentage per category under each run, flags a sign change in either sensitivity run and the
largest absolute shift (`R/run_analyses.R:86-99`). `build_stats_table()` runs the Welch test and BH
on the primary run (`R/build_stats_table.R`). Both reproduce exactly from the committed composition
tables (findings log §3).

### 4.10 Stage 7: driver-gene diagnostics (run script)

For a named category: per mouse, the percentage of its cells with each panel gene detected, against
all other cells (`driver_check()`); the pooled runner-up distribution and median margin
(`runner_up_summary()`); a spatial plot of the category's cells against two context categories
(`plot_category_spatial()`). Reads the 13-panel checkpoint (findings R4).

## 5. Statistical models and estimators

### 5.1 Per-cell identity by marker-panel module score and argmax

Question: which of the candidate cell types does this cell most resemble? Inputs: the log-normalised
counts of the cells inside the lasso and the panel gene lists. Method: `AddModuleScore()` computes,
per cell, the mean expression of the panel genes minus the mean of control genes drawn from the
same expression bins (Seurat default 100 controls per gene, 24 bins, `seed = 1`); the score is
z-scored across the cells of that sample within each panel; argmax over panels. Assumptions, not
stated in the code: that a short panel captures the identity; that z-scoring within a sample is
appropriate even though it removes any between-sample difference in overall signal (findings M4);
that the detection floor is a property of the sample rather than of the tissue (finding M3).
Diagnostics saved: none in files; `classification_margin` and `runner_up` are per-cell columns;
the Purkinje percentage is printed. Parameters: `min_detect_rate = 0.005`, `min_z = 0`, in the
function signature.

### 5.2 State refinement

Question: within astrocytes, microglia and oligodendrocyte lineage, is this cell in the "Active"
(or "Mature") or "Homeostatic" (or "Precursor") state? Method: two module scores on the subset,
larger wins, ties to the second (`R/classify_atlas_index.R:41-47`). Assumption not stated: that
every cell of the lineage is in one of the two states (no undetermined state; finding B3). For
astrocytes and microglia the Homeostatic panel is the identity panel, so the comparison is between
a fresh score and the genes the cell was selected on (finding M10).

### 5.3 Pooled-percentile confidence threshold

Question: which calls are too close to a tie to trust? Method: 15th percentile of the pooled margin
distribution; cells below become Unclassified in the sensitivity runs. Assumption: a fixed share of
calls (15 %) is noise in every region; the percentile was chosen on the optic tract RGC diagnostic
and reused (finding M6). Diagnostic: the console message with the threshold and count.

### 5.4 Composition test: Welch t-test and Benjamini-Hochberg

Question: does the mean percentage of a category differ between three TBI and three sham mice?
Method: `t.test(pct ~ group_id)` per category, two-sided, unequal variances; `p.adjust(method =
"BH")` across the categories tested (15 in the 12-panel hippocampus). Assumptions, not stated in
the code: approximate normality of percentages within group at n = 3 (untestable); independence
between categories, which does not hold because percentages within a mouse sum to 100;
"Unclassified" treated as a category (finding M5). Result on the committed run: smallest raw p
0.125 (Astrocyte Homeostatic), all BH-adjusted p at or above 0.74. Parameters: none configurable.

### 5.5 Pseudobulk DESeq2

Question: for a cell type, which genes differ in summed expression between three TBI and three sham
mice? Method: negative-binomial generalised linear model per gene on per-mouse summed counts,
design `~ group`, DESeq2 defaults (median-of-ratios size factors, dispersion shrinkage, Wald test,
independent filtering at alpha 0.1, Cook's distance outlier flagging, which with three replicates
per group sets flagged genes' p-values to NA rather than replacing counts). Two multiplicity
treatments: BH over all genes passing `rowSums >= 10` (genome-wide), and BH over the pre-specified
candidate panel using the same per-gene p-values (panel). Assumptions: the mouse is the unit of
replication (stated and correct); the three mice per group are exchangeable; no covariates (lasso
area, cell count per mouse, section depth) enter the model, which is unexplained. Diagnostics: none
saved (the protocol cites p-value histograms seen interactively; finding M8). Parameters:
`min_row_sum = 10`, contrast `c("group", "gi", "gii")`, fold changes TBI relative to sham.

### 5.6 Wilcoxon FindMarkers

Question: for a cell type, which genes differ between TBI cells and sham cells treated as
independent observations? Method: Seurat `FindMarkers(test.use = "wilcox", logfc.threshold = 0.25,
min.pct = 0.1)`, Bonferroni adjustment over all genes tested. Assumption that does not hold:
independence of cells within a mouse (pseudoreplication; stated in the code header and the
protocol). The result feeds only the Wilcoxon heatmap and the printed count of genes below 0.05.

### 5.7 Robustness criterion

A composition difference is "robust to near-tie handling" if its sign is unchanged across the
three runs (README "Confidence filtering"). This is a direction check, not a test; a category with a
difference near zero flips sign on tiny shifts (Pericyte/vascular SMC, 0.30 to -0.01 percentage
points), and the README already reads those as negligible.

## 6. Randomness and reproducibility

Random draws in the code, each with the seed that governs it:

| Draw | Where | Seed | Scope of the seed |
|---|---|---|---|
| Control-gene sampling in `AddModuleScore()` | `R/classify_atlas_index.R:23, 44-45` | Seurat default `seed = 1`, set inside the function on every call | per call; also resets the session's random state as a side effect |
| `RunPCA()` (irlba start) | `R/run_full_pipeline.R:209` | Seurat default `seed.use = 42` | per call |
| `RunUMAP()` | `R/run_full_pipeline.R:210` | Seurat default `seed.use = 42` | per call |
| Point jitter on the composition chart | `R/plot_composition_barchart.R:24` | `position_jitterdodge(seed = 1)` | per plot |

No `set.seed()` exists anywhere in the repository (static sweep, 0 hits), so reproducibility rests
on package defaults, which are deterministic for a given package version. What changes between
runs with identical inputs and packages: nothing, in principle; the baseline comparison (§9.1) is
the test. What still moves across machines or versions: UMAP coordinates (uwot and Seurat
versions), DESeq2 results in the last digits (BLAS), and PNG bytes (graphics device, font
rendering), so figures are expected to differ byte-wise between the author's Mac and any other
machine even when the tables are identical. Same-seed and cross-seed evidence: pending (the
author's baseline run gives the same-input comparison against the committed outputs; a second run
on the same machine would give the same-seed evidence and is requested if time allows).

## 7. Outputs

Tracked deliverables are the figures (PNG) and the tables (CSV) under `outputs/`; `.gitignore:49`
ignores every CSV and `:54-55` re-admits those under `outputs/`; `.gitignore:48` ignores every
`.rds`, which keeps the checkpoints and the merged object out of the repository. Everything is
overwritten in place by a rerun into the same folder, so a regeneration shows in `git diff`.

| Output | Format and units | Written by | Consumer and sanity check |
|---|---|---|---|
| `outputs/<run>/cell_counts_per_sample.csv` | CSV, 6 rows; cells, Unclassified %, QC pass %, lasso area in square micrometres, cells per 1,000 square micrometres | stage 1c | first thing to read after a run; `n_cells` must be identical across the three runs of a region |
| `outputs/<run>/composition_per_sample.csv` | CSV, 6 x categories rows; counts and percentages | no tracked code (R1) | percentages sum to 100 per mouse |
| `outputs/<run>/composition_barchart.png` | PNG, mean and SEM bars with per-mouse points, categories ordered by abundance | stage 3 | each category shows six points, three filled and three hollow |
| `outputs/<run>/umap_celltype_gi_vs_gii.png`, `umap_microglia_check_gi_vs_gii.png` | PNG, split by group | stage 4 | the two panels should show the same islands; the microglia check highlights two states against grey |
| `outputs/<run>/pseudobulk_heatmap_<ct>.png` | PNG, candidate genes by mouse, z-scored log2 normalised counts, clipped at +-1.5, asterisks for panel padj | stage 5 | identical across the three runs of a region (M2) |
| `outputs/<run>/wilcox_heatmap_<ct>.png` | PNG, Wilcoxon genes at Bonferroni-adjusted p < 0.05, same display | stage 5 | identical across the three runs of a region |
| `outputs/<run>/lasso_overview.png` | PNG, six panels | stage 0, primary runs | each red outline sits on the hippocampus (or tract) of its own section |
| `outputs/hippocampus_unfiltered/de/<ct>_{res_full,res_panel,wilcox_res}.csv` | CSV; DESeq2 columns (`baseMean, log2FoldChange, lfcSE, stat, pvalue, padj`), the panel table with `padj_panel`, FindMarkers columns with `p_val_adj`; leading unnamed index column | no tracked code (R2) | `padj_panel` equals BH over the panel's rows (verified) |
| `outputs/<region>_near_tie_sensitivity.csv` | CSV, one row per category; TBI minus sham percentage points under each run, `flips`, `max_shift` | stage 6 | reproduces from the three composition tables (verified) |
| `outputs/hippocampus_unfiltered_stats_table.csv` | CSV, 15 rows; `p_value`, `p_adj_BH` | stage 6 | reproduces to 1.6e-14 (verified) |
| `outputs/hippocampus_{choroid_plexus,ependymal}_{gene_detection,runner_up}.csv`, `_spatial.png` | CSV and PNG | stage 7 | needs the 13-panel checkpoint (R4) |
| `outputs/hippocampus_per_sample_depth.csv` | CSV, 6 rows; median transcripts and genes | no tracked code (R3) | |
| `outputs/*_withCP/` | the 13-panel hippocampus runs, kept for comparison (README "Outputs") | superseded code state | not regenerable from `main` (R4) |

Output tree at 869af28 (tracked): 11 files at the top of `outputs/`, 13 per run folder (14 plus a
`de/` folder of 12 for `hippocampus_unfiltered`), 8 per optic tract folder (no `lasso_overview.png`,
no composition table); 12.4 MB of PNG and CSV in all. Ignored and present only on the author's
disk: per-run `*.rds` checkpoints and `combined.rds`.

## 8. Run times and efficiency

### 8.1 Timing table

Pending: the author's baseline log (hand-off in the findings log §3). No measurement exists in the
repository; README says only "each run is heavy".

| Stage | Wall time | What dominates | Cache or skip candidate |
|---|---|---|---|
| 1 load and classify (six mice) | pending | `LoadXenium()` and `AddModuleScore()` per panel | already checkpointed per mouse (`*_classified_prefilter.rds`) |
| 4 UMAP | pending | `RunUMAP()` on about 36,000 cells | skip in sensitivity runs (the embedding does not depend on the filter) |
| 5 DE (four cell types) | pending | `DESeq()` and `FindMarkers()` | skip in sensitivity runs (results identical, M2) |
| Total, one region block | pending | | |

### 8.2 Existing caches and rerun behaviour

- Per-mouse checkpoint: if `outputs/<run>/<sample>_classified_prefilter.rds` exists it is loaded
  instead of recomputed (`R/run_full_pipeline.R:138-141`). It is not invalidated by a change to the
  marker list, the lasso file or the QC cut-offs (finding R7).
- `prefilter_objects`: the sensitivity runs reuse the primary run's objects from memory, so they
  skip stage 1 but repeat stages 2 to 5.
- A rerun overwrites every output in place; nothing fails on a second run.

### 8.3 Proposed caches and toggles

Skip stages 4 and 5 when the run is a sensitivity run (an argument defaulting to the current
behaviour), saving two of three UMAP and DE passes per region with byte-identical results for the
figures that remain; measured savings pending the baseline timings.

### 8.4 Code inventory

From the parse-only sweep of 2026-10-07 (19 files, 33 function definitions, no duplicates):
`plot_spatial_celltype()` is never called (E1); `check_marker_spillover()` and its three helpers
are a diagnostic module sourced on every setup and called from nowhere in the pipeline (E4);
`apply_confidence_filter()` and `lasso_polygon_area()` are called by name through `lapply()` and
`vapply()` and show as zero-call leads only. No file is unsourced. No function carries a roxygen
header. Tracked artifacts no code reads: `docs/hippocampus_lassos_overview.png` (a byte-identical
copy of a run output, E2) and `docs/optic_tract_new_lassos_overview.png` (referenced from a comment
at `R/run_full_pipeline.R:303`).

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
