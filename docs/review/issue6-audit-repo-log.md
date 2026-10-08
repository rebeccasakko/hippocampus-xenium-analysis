# hippocampus-xenium-analysis issue #6: audit, make runnable by three people, document methods and reproducibility

## 0. Status

```
Issue: #6 (https://github.com/rebeccasakko/hippocampus-xenium-analysis/issues/6)
Branch: 6-audit-repo        Pull request: #7 (draft, https://github.com/rebeccasakko/hippocampus-xenium-analysis/pull/7)
Mode: audit (degraded: this session cannot run the code; runs are handed to the author)
Org profile: default
Reviewer: Holden Harris (holden-harris), with Claude Code        Original author: rebeccasakko
Line numbers refer to commit: 869af28 (code unchanged by 2a1df70, which added the two review drafts)
Last completed step: §5 step 5 (portability code)        Next step: pause 7, then §5 step 6; step 2 baseline result pending from the author        Last updated: 2026-10-07
```

## 1. Scope and targets

The three run targets:

1. **The analyst**, Holden Harris: Windows 11 Home, R 4.5.3 (`C:/Program Files/R/R-4.5.3`), RStudio.
   No Xenium data on the machine; Seurat and DESeq2 not installed. Can run the static checks and,
   once the data is placed, the pipeline.
2. **The original author**, rebeccasakko: macOS, repository under `~/Documents/GitHub/`, Xenium
   bundles on the lab drive at `/Volumes/Joe/Becca and Sakthi Files/Xenium_Data/`, lasso CSVs on
   the Desktop, runs interactively from RStudio one block of `R/run_analyses.R` at a time. Their
   layout and habits must keep working unchanged.
3. **A fresh clone** on any OS with only the documented data placed and a local config file added,
   run from the command line.

In scope: the whole repository (19 R files, two protocols, README, tracked outputs), all eight
lenses. Out of scope: the proposed FLASH-MM mixed-model DE (README "Planned analysis"; noted in the
science review, not reviewed). Statistical and methods choices (M lens) are examined in depth and
recorded as questions for the author in §8 and as sub-issues; none is changed in this pull request
(decision 1).

Reference outputs: the committed `outputs/` tree at 869af28 (decision 2). Baseline: run by the
author on their Mac from the hand-off block in §3 (decision 3); until the log arrives every run
measurement reads "pending".

Documents produced: `docs/review/science-review.md`, `docs/review/workflow-review.md`, this findings
log. Commit trailer: `Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>`. Commit subject
suffix `(issue #6)`.

Related: issue #1 (the analyst's orientation checklist); PR #3 (the 5 Oct 2026 read-through review,
`docs/review_2026-10-05.md`, since removed from `main` in 683e6bc). Its open items are carried into
§2 with "(carry-over from PR #3 review, item X)" in the Problem cell.

## 2. Findings register

Rows below were established during Phases 2 and 4 (static). Phase 5 completes each table. Every
row cites 869af28.

### 2.1 Portability (P)

| ID | Where | Problem | Kind | Severity | Status |
|---|---|---|---|---|---|
| P1 | `R/run_full_pipeline.R:252-260` | Six `data_dir` paths are literals on the author's lab drive (`/Volumes/Joe/...`); required by stage 0 and stage 1 | mechanical | high | fix: step 5 |
| P2 | `R/run_full_pipeline.R:265-272, 306-313` | Twelve lasso paths are literals under the author's `~/Desktop`; required by stage 0, 1 and 1c | mechanical | high | fix: step 5 |
| P3 | `R/run_analyses.R:17` | `setwd("~/Documents/GitHub/hippocampus-xenium-analysis")`: works for the author only; no root anchor anywhere, and the `.Rproj` file is not tracked | mechanical | high | fix: step 5 |
| P4 | `R/check_marker_spillover.R:36` | `library(tidyverse)` at source time, sourced by the setup file at `R/run_full_pipeline.R:66`; `tidyverse` is declared in neither the header list (`:38-39`) nor the README | mechanical | medium | fix: step 5 |
| P5 | `R/run_full_pipeline.R:42-50` | No package check: a missing package fails at the first `library()` with no list of everything missing; packages used by namespace (Matrix, grid) and Seurat's Xenium readers (e.g. `arrow`) are declared nowhere | mechanical | medium | fix: step 5 |
| P6 | `R/run_full_pipeline.R:202, 213, 218`; `R/plot_lasso_overview.R:25`; `R/check_driver_genes.R:53` | Figures `print()`ed unguarded; under `Rscript` they go to `Rplots.pdf` in the working directory | mechanical | low | fix: step 5 |

### 2.2 Reproducibility (R)

| ID | Where | Problem | Kind | Severity | Status |
|---|---|---|---|---|---|
| R1 | `outputs/hippocampus_{unfiltered,sens_global,sens_neartie}/composition_per_sample.csv`; `R/run_full_pipeline.R:196-205` | Tracked per-sample composition tables that no tracked code writes (`build_composition_summary()` returns them in memory only; added in commit fdde845 with no code change). `git log --all -S"composition_per_sample"` finds no writer | mechanical | medium | fix: step 6 |
| R2 | `outputs/hippocampus_unfiltered/de/*.csv` (12 files); `R/run_full_pipeline.R:223-238` | Tracked DE result tables with no writer in the code (added in 713b00d with no code change); their header carries a duplicated first column name (`"gene","gene"`), the mark of an interactive `write.csv()` with row names | mechanical | medium | fix: step 6 |
| R3 | `outputs/hippocampus_per_sample_depth.csv` | Tracked table (median transcripts and genes per sample) with no producing code anywhere in the repository | mechanical | low | fix: step 6 |
| R4 | `R/run_analyses.R:118` with `R/marker_panels.R:96` | The driver-gene diagnostics read `outputs/hippocampus_unfiltered_withCP/sample_objects_prefilter.rds`, a gitignored checkpoint of a 13-panel run; `hippocampus_markers` now has 12 panels and no code path produces the `_withCP` run, so the six `hippocampus_choroid_plexus_*` / `hippocampus_ependymal_*` outputs and the four `_withCP` folders cannot be regenerated from the script | mechanical | medium | issue #9 |
| R5 | repository-wide | No record of R, Seurat, DESeq2 or other package versions: no `sessionInfo()`, lock file or README line; `R/run_pseudobulk_de.R:26` and `R/plot_de_heatmap.R:57` depend on Seurat v5 turning `_` into `-` in aggregated sample names | mechanical | medium | fix: step 6 |
| R6 | `R/run_analyses.R:44, 72` | The near-tie category lists are typed by hand from the primary run's output ("categories where >15 % of cells fall below the threshold"); no code computes or checks them, so a rerun with changed inputs silently keeps the old list | mechanical | medium | fix: step 7 |
| R7 | `R/run_full_pipeline.R:137-141` | The per-sample checkpoint `<sample>_classified_prefilter.rds` is keyed on `output_dir` and sample only; a rerun into the same folder after a marker-list or lasso change loads the stale classification without warning | mechanical | medium | fix: step 7 |
| R8 | `outputs/optic_tract_newlasso_*/` | Three optic tract run folders hold no `lasso_overview.png` and their `composition_barchart.png` predates the plotting fix of 8d9b99f (README "Open items"); they were produced by an earlier state of the code than 869af28 | mechanical | low | issue #9 |

### 2.3 Bugs and fragility (B)

| ID | Where | Problem | Kind | Severity | Status |
|---|---|---|---|---|---|
| B1 | `R/build_stats_table.R:31` | `t.test(pct ~ group_id, data = pick(everything())$p.value, error = function(e) NA_real_)`: the closing parenthesis of `t.test()` is misplaced, so the file does not parse (`38:1: unexpected '}'`) and `source("R/run_full_pipeline.R")` stops at line 65 for every user. Found by the parse sweep on 2026-10-07 (18 of 19 files parse). The committed stats table predates the edit and recomputes from the committed composition table with the corrected call (max p difference 1.6e-14). Fixed before the baseline because it prevents any run at all (decision 4) | mechanical | high | fix: step 2a (decision 4) |
| B2 | `R/process_sample_template.R:11-13` vs `R/run_full_pipeline.R:76-78` | Point-in-polygon treats the whole lasso CSV as one polygon and ignores `Selection`; the area function splits by `Selection` and sums. A file holding more than one selection is joined into one wrong outline by the first and measured correctly by the second, with no check on the number of selections (carry-over from PR #3 review, item V) | mechanical (check) / scientific (if a file has several selections) | medium | fix: step 7 |
| B3 | `R/classify_atlas_index.R:46` | State refinement `ifelse(StateA1 > StateB1, ...)`: exact ties and cells with no signal from either state panel go to Homeostatic/Precursor; no "cannot tell" state (carry-over item K) | scientific | low | author (issue #8) |
| B4 | `R/run_full_pipeline.R:118` | `dir.create(output_dir, recursive = TRUE)` is correct, but every writer (`write.csv`, `ggsave`, `png`) overwrites silently, so a second run into a folder of a different run label mixes files; no run label or commit is written into the outputs | mechanical | low | fix: step 6 |
| B5 | `R/build_cell_count_table.R:17`; `R/run_pseudobulk_de.R:26`; `R/plot_de_heatmap.R:57, 61` | Group membership is inferred from the sample name by regex (`grepl("gii", sample)`, `"^gi-"`) instead of from `group_id`; correct for the six current names, silently wrong for any renamed or added sample, and the `"^gi-"` form also depends on Seurat v5 renaming `_` to `-` (R5) | mechanical | low | fix: step 7 |
| B6 | `R/classify_atlas_index.R:9, 15` | Panel exclusions (no gene on the Xenium panel; best gene below the 0.5 % floor) are `warning()`s only. In a six-mouse run R defers and truncates them, and nothing writes which panels were scored in which mouse, although composition depends on it (M3; README "Open items" asks exactly this about BAM in `gi_2_d`) | mechanical (write a per-sample panel table) | medium | fix: step 6 |
| B7 | `R/run_pseudobulk_de.R:22` | `AggregateExpression(..., slot = "counts")`: `slot` is deprecated in Seurat v5 in favour of `layer`; runs with a warning today and may stop in a later Seurat | mechanical | low | fix: step 7 |
| B8 | `R/run_analyses.R:93-97` | `robustness_table()`: a category present in one run and absent in another gives `NA` in `unfiltered`, `global` or `neartie`, so `flips` and `max_shift` are `NA` and the row sorts to the end without notice | mechanical | low | fix: step 7 |

### 2.4 Documentation (D)

| ID | Where | Problem | Kind | Severity | Status |
|---|---|---|---|---|---|
| D1 | `protocols/differential_expression_protocol.md:313-316` vs `:381-383` and `outputs/hippocampus_unfiltered/de/dg_res_full.csv` | "A genome-wide test ... returned zero significant genes in every cell type" while the results table of the same document and the committed table list Ddit4 at padj 0.004 in dentate gyrus granule cells | mechanical | medium | fix: step 8 |
| D2 | `R/de_candidate_panels.R:4-7`; `R/plot_de_heatmap.R:80, 88`; `R/run_pseudobulk_de.R:2`; `R/run_wilcox_de.R:2` | Comments call pseudobulk the "secondary" and Wilcoxon the "primary" analysis; README and the DE protocol (since 07b9305) say the reverse | mechanical | low | fix: step 8 |
| D3 | `protocols/classification_protocol.md:40-42, 104-108, 276-280` | Says the hippocampus uses "~13 panels" and "the full, unrestricted `atlas_index_markers`", that the optic tract drops two panels, and that the hippocampus "is not (yet) using" the percentile filter; the code uses 12 and 16 panels and the percentile design for both regions (`R/marker_panels.R:79-96`; `R/run_analyses.R:35-53`) | mechanical | medium | fix: step 8 |
| D4 | `R/run_full_pipeline.R:166-168, 172-173` | Comment says the hippocampus "stays on its existing BAM-specific fix" and the primary run filters "only types with a diagnosed near-tie problem, e.g. BAM"; neither is true at 869af28 | mechanical | low | fix: step 8 |
| D5 | `R/build_stats_table.R:16-17` | "16 for Hippocampus's full panel": the 12-panel run tests 15 categories (`outputs/hippocampus_unfiltered_stats_table.csv`) | mechanical | low | fix: step 8 |
| D6 | README "Outputs" ("spatial maps ... in outputs/"); `R/plot_spatial_celltype.R:1` | No per-sample spatial cell-type map exists in `outputs/`, and `plot_spatial_celltype()` is called from nowhere; the only spatial figures are the two driver-gene diagnostics | mechanical | low | fix: step 8 |
| D7 | README "Cell-type classification" ("Optic-tract immune panels: [TODO: source]") | Six panels have no recorded source (carry-over item W) | mechanical | low | author (issue #8) |
| D8 | README lines 1-10 | The injury model (CHIMERA), the time point (7 days) and the tissue section plane are not stated anywhere in the repository; they are known only from the author's notes quoted in the PR #3 review | mechanical | low | author (issue #8) |
| D9 | README "Known limitations" bullet 1; `R/run_wilcox_de.R:22-29`; DE protocol "Results summary" | "The smallest possible two-sided Wilcoxon p-value is 0.1" describes a 3-vs-3 rank test, but composition uses Welch's t-test and the cell-level Wilcoxon has thousands of observations. Separately, `FindMarkers()` adjusts by Bonferroni over all genes, so the Wilcoxon "padj < 0.05" counts (23, 9, 4, 1) are Bonferroni, not BH; no document says so | mechanical | low | fix: step 8 |

### 2.5 Methods and statistics (M)

Rows here are questions for the author (§8) unless marked otherwise; the arithmetic stays as it is.

| ID | Where | Problem | Kind | Severity | Status |
|---|---|---|---|---|---|
| M1 | `R/marker_panels.R:4` with `R/de_candidate_panels.R:18`; `R/marker_panels.R:3` with `:15` | Circularity: P2ry13, Tmem119 and Sall1 both define microglia and are tested within microglia; Aldh1l1 and Sox9 likewise for astrocytes. P2ry13 (padj 0.015, lower in TBI) is the headline microglia result (carry-over item C) | scientific | high | author (issue #8) |
| M2 | `R/apply_confidence_filter.R:2-5` with `R/run_pseudobulk_de.R:18-19`, `R/run_wilcox_de.R:18-19` | The confidence filter rewrites `cell_type` while DE selects on `atlas_identity`, so the sensitivity runs never test DE sensitivity; verified: all eight DE heatmaps byte-identical across the three hippocampus runs (carry-over item A) | scientific (design) | medium | author (issue #8) |
| M3 | `R/classify_atlas_index.R:12-20` | The 0.5 % detection floor is evaluated per sample, so a panel can be scored in one mouse and not another; README "Open items" records BAM in `gi_2_d` going from 0 to 264 cells between runs (carry-over item E) | scientific | medium | author (issue #8) |
| M4 | `R/classify_atlas_index.R:28` | Module scores are z-scored within each sample before the argmax, which removes between-mouse differences in overall score level and can attenuate a true composition shift (carry-over item D) | scientific | medium | author (issue #8) |
| M5 | `R/build_stats_table.R:26-38` | Welch t-tests on compositional percentages, one per category including "Unclassified", BH across 15; n = 3 per group. Alternatives built for proportions exist (e.g. `speckle::propeller`) (carry-over item G) | scientific | medium | author (issue #8) |
| M6 | `R/run_analyses.R:39, 67`; classification protocol "Confidence-margin filtering" | The 15th percentile was justified on the optic tract RGC margin distribution and reused for the hippocampus for comparability; README "Open items" already lists a 10th/25th robustness check | scientific | low | author (issue #8) |
| M7 | `R/run_pseudobulk_de.R:51-64` | Candidate-panel padj re-adjusts the genome-wide p-values within 3 to 7 pre-specified genes after the genome-wide pass has been read; the pre-specification is documented, but the two passes are reported side by side with no multiplicity accounting across them | scientific | low | author (issue #8) |
| M8 | `R/run_pseudobulk_de.R:29-31` | DESeq2 defaults: independent filtering, Cook's outlier handling and `rowSums >= 10` on three-vs-three pseudobulks; no dispersion plot or p-value histogram is saved although the protocol cites "flat p-value histograms" as evidence | mechanical (save diagnostics) | low | fix: step 6 |
| M9 | `R/plot_de_heatmap.R:50-51` | Heatmaps display `log2(normalized counts + 1)` z-scored by gene across six mice; with one mouse driving a gene the scale saturates at +-1.5 (`breaks`), which is fine for display but is not stated in the caption | mechanical (caption) | low | fix: step 8 |
| M10 | `R/marker_panels.R:3-4, 98-101` with `R/classify_atlas_index.R:41-58` | The Astrocyte and Microglia Homeostatic state panels are the identity panels (protocol "double duty"), so the Active/Homeostatic split compares Gfap/Serpina3n (or Trem2/Cst7/Lpl) against the genes the cells were selected on; "Active" then means "the state genes beat an already-high baseline", which biases the split towards Homeostatic by construction and makes the Active share hard to interpret across regions (README already flags the white-matter astrocyte case) | scientific | medium | author (issue #8) |
| M11 | `R/run_analyses.R:11-12, 44, 72` | Near-tie categories are those with more than 15 % of cells below the pooled 15th percentile; the pooled share is 15 % by construction, so the rule flags categories noisier than average and its outcome depends on how many categories there are. Stated as a heuristic, never tested; the README's 10th/25th percentile check would exercise it | scientific | low | author (issue #8) |

### 2.6 Stochasticity (S)

| ID | Where | Problem | Kind | Severity | Status |
|---|---|---|---|---|---|
| S1 | `R/classify_atlas_index.R:23, 44-45` (`AddModuleScore`, control-gene sampling, Seurat default `seed = 1`); `R/run_full_pipeline.R:209-210` (`RunPCA`, `RunUMAP`, Seurat default `seed.use = 42`); `R/plot_composition_barchart.R:24` (`position_jitterdodge(seed = 1)`) | No `set.seed()` anywhere (`static_sweep.R`: 0 hits); every random step relies on a package default seed, which is deterministic for a given package version but is nowhere stated. UMAP coordinates will still move with uwot/Seurat versions | mechanical (document) | low | fix: step 8 |

### 2.7 Run time (T)

| ID | Where | Problem | Kind | Severity | Status |
|---|---|---|---|---|---|
| T1 | README "Repository structure" ("each run is heavy"); `R/run_analyses.R:3` | No measured run time, timer or log anywhere; claimed, not measured. Baseline pending on the author's machine | mechanical | low | fix: step 6 (timing lines), step 9 (measured) |
| T2 | `R/run_full_pipeline.R:207-238` | The sensitivity runs recompute PCA, UMAP and both DE analyses although DE is identical to the primary run (M2) and the UMAP embedding depends only on the merged counts; a skip toggle would save two of three heavy DE passes per region | mechanical (if outputs identical) | low | fix: step 7 |

### 2.8 Efficiency and artifacts (E)

| ID | Where | Problem | Kind | Severity | Status |
|---|---|---|---|---|---|
| E1 | `R/plot_spatial_celltype.R:1` | `plot_spatial_celltype()` is defined, sourced by the setup file, and called from nowhere (sweep: 0 call sites; README lists it as a figure function) | mechanical | low | wontfix (decision 7) |
| E2 | `docs/hippocampus_lassos_overview.png` | Byte-identical (md5 `2086421d...`) to `outputs/hippocampus_unfiltered/lasso_overview.png` and to the three `_withCP` copies; five copies of one 363 KB figure tracked | mechanical | low | wontfix (decision 7) |
| E3 | `R/run_full_pipeline.R:315-322` | `optic_tract_lassos_old`: six unused paths kept "for reference only"; documented, so a candidate for the config example's comments rather than live code | mechanical | low | fix: step 5 |
| E4 | `R/check_marker_spillover.R` | Diagnostic module sourced on every setup although nothing in the pipeline calls it (sweep: `check_marker_spillover()` 0 call sites); it is the only reason `tidyverse` is attached (P4) | mechanical | low | fix: step 5 |

## 3. Evidence notes

- **B1** (2026-10-07, analyst's machine, R 4.5.3, dplyr 1.2.0): `Rscript static_sweep.R` reports
  `R/build_stats_table.R:38:1: unexpected '}'`. Recomputing Welch p-values and BH from
  `outputs/hippocampus_unfiltered/composition_per_sample.csv` with `t.test(pct ~ group_id,
  data = pick(everything()))$p.value` reproduces `outputs/hippocampus_unfiltered_stats_table.csv`
  to 1.6e-14 (p) and 5.1e-15 (padj), 15 rows.
- **M2 / T2** (2026-10-07): `md5sum` of the eight `pseudobulk_heatmap_*.png` and
  `wilcox_heatmap_*.png` in `hippocampus_unfiltered`, `hippocampus_sens_global` and
  `hippocampus_sens_neartie`: eight distinct hashes, each appearing exactly three times.
- **Robustness table** (2026-10-07): `outputs/hippocampus_near_tie_sensitivity.csv` reproduces
  exactly from the three committed `composition_per_sample.csv` files with the `robustness_table()`
  arithmetic (`R/run_analyses.R:86-99`). Percentages sum to 100 in every sample; 15 categories in
  every sample.
- **DE tables** (2026-10-07): the committed `de/` tables give 23, 9, 4 and 1 Wilcoxon genes at padj
  < 0.05 (DG, astrocyte, microglia, GABA), Ddit4 as the only genome-wide hit (DG, padj 0.0044,
  log2FC +0.46), Id2 (DG, panel padj 0.048, +0.37) and P2ry13 (microglia, panel padj 0.015,
  -0.45), matching the DE protocol results table. `padj_panel` equals `p.adjust(pvalue, "BH")` over
  the panel's own rows in all four cell types.
- **R1-R3**: `git log --all -S"composition_per_sample"` and `-S"_res_full.csv"` return nothing;
  the files were added in fdde845 and 713b00d alongside README and protocol edits only.
- **E2**: `md5sum docs/*.png outputs/*/lasso_overview.png`.
- **Baseline run**: pending, run on the author's machine; hand-off block below.

### Hand-off block for the baseline (pause 4)

Posted as a comment on PR #7 for the author on 2026-10-07. The branch they run is `6-audit-repo`
after the B1 parse fix (step 2a); no other code differs from `main`. Everything the block produces
lands in `../hippocampus-xenium-analysis-audit/<date>/`, outside the repository.

```
Run this on your Mac, then paste back the four items at the end (as a comment on PR #7).

1. Get the branch (your working tree must be clean; commit or stash anything pending):
     cd ~/Documents/GitHub/hippocampus-xenium-analysis
     git fetch origin && git checkout 6-audit-repo
   The only code difference from main is a one-character fix in R/build_stats_table.R
   that lets source("R/run_full_pipeline.R") work again.

2. In RStudio, open the project folder as usual, then in the console (base R only):
     stamp <- format(Sys.Date(), "%Y%m%d")
     A <- file.path("..", "hippocampus-xenium-analysis-audit", stamp); dir.create(A, recursive = TRUE)
     f <- list.files("outputs", recursive = TRUE, full.names = TRUE)
     write.csv(data.frame(file = f, md5 = unname(tools::md5sum(f))), file.path(A, "md5_before.csv"), row.names = FALSE)
     sink(file.path(A, "run_hippocampus.log"), split = TRUE)
     t0 <- Sys.time(); print(t0)
     source("R/run_full_pipeline.R")

3. Run the HIPPOCAMPUS block of R/run_analyses.R exactly as you normally do (lines 29-53:
   hippocampus_results, hippocampus_sens_global, hippocampus_sens_neartie), then the two
   write.csv lines for the hippocampus robustness and stats tables (lines 101-102 and 107-108).
   Do not run setwd() on line 17 (you are already in the folder). Print the time after each run:
     print(Sys.time() - t0)

4. When the hippocampus block has finished:
     print(Sys.time() - t0); sink()
     f <- list.files("outputs", recursive = TRUE, full.names = TRUE)
     write.csv(data.frame(file = f, md5 = unname(tools::md5sum(f))), file.path(A, "md5_after.csv"), row.names = FALSE)
     writeLines(capture.output(sessionInfo()), file.path(A, "sessionInfo.txt"))
     file.copy("outputs", A, recursive = TRUE)
     # one quick check for the audit (B2): how many selections does each lasso file hold?
     sapply(c(hippocampus_lassos, optic_tract_lassos), function(p) if (file.exists(p)) length(unique(read.csv(p, skip = 2)$Selection)) else NA)

5. Optional, if you have the time: the OPTIC TRACT block (lines 57-82 and 103-104) the same way,
   into a second log (run_optic_tract.log). The diagnostics block (lines 118-130) only if
   outputs/hippocampus_unfiltered_withCP/sample_objects_prefilter.rds still exists on your disk.

6. Do NOT commit the regenerated outputs/; after step 4 you can restore them with
     git checkout -- outputs
   (the copies are safe in the audit folder).

Paste back:
  a. run_hippocampus.log (whole file if under ~300 lines, else the last 60 lines plus every line
     containing "Error", "Warning" or "apply_confidence_filter_pooled"), and the printed times.
  b. sessionInfo.txt.
  c. md5_before.csv and md5_after.csv.
  d. The lasso selection counts from step 4.
```

Until these arrive, the workflow review §9.1 and every run-time cell read "pending: run on the
author's machine".

## 4. Recommended changes (fix design)

Constraints (decision 5): no new package dependency; every existing function keeps its signature,
new arguments default to the current behaviour; nothing in `R/marker_panels.R`,
`R/de_candidate_panels.R` or any statistical call changes; the author's interactive habit (open the
folder in RStudio, `source("R/run_full_pipeline.R")`, run blocks of `R/run_analyses.R`) keeps
working with one short untracked file; base R for everything that runs before packages are
attached. Evidence for each code step is static on this machine (parse, and where possible a base-R
or dplyr-only recomputation) and dynamic on the author's machine in Phase 8.

### 4.1 Configuration

One gitignored `config.local.R` in the repository root, sourced by `R/run_full_pipeline.R` after a
SETTINGS block of repo-relative defaults and before any input is read; a tracked
`config.local.example.R` documents every key with the author's real paths as commented examples.

| Key | Default | Purpose |
|---|---|---|
| `dir.data` | `file.path(dir.repo, "data")` | Folder holding the inputs that cannot ship; gitignored |
| `xenium_dirs` | `file.path(dir.data, "xenium", "0071470_<sample>")`, named by sample id | The six Xenium bundles. The author's layout: `/Volumes/Joe/Becca and Sakthi Files/Xenium_Data/Group_i/0071470_gi_1_b/` and siblings (commented example) |
| `hippocampus_lassos` | `file.path(dir.data, "lassos", "hippocampus", "hippocampus_<sample>.csv")` | Six hippocampus lasso files. Author: `~/Desktop/hippocampus_<sample>.csv` |
| `optic_tract_lassos` | `file.path(dir.data, "lassos", "optic_tract", "<GI-1 .. GII-3>_coordinates.csv")` | Six optic tract lasso files. Author: `~/Desktop/B6_Xenium_Opt/...` |
| `dir.out` | `file.path(dir.repo, "outputs")` | Where runs write; tracked deliverables stay in the repository by default |

`sample_base` is then built from `xenium_dirs`, so `run_region_pipeline()` and every function
signature are unchanged. The superseded `optic_tract_lassos_old` vector moves into the example file
as a comment (E3).

### 4.2 Setup checks

- **Root anchor** at the top of `R/run_full_pipeline.R`: stop with a readable message naming
  `hippocampus-xenium-analysis.Rproj` and printing `getwd()` when the file is absent. The `.Rproj`
  file (RStudio's default, created by the analyst's clone) becomes tracked so a fresh clone has the
  anchor. `setwd()` leaves `R/run_analyses.R` (P3).
- **`R/_setup.R`** (base R only): `fn.check_packages()` over the required list (Seurat, dplyr,
  tidyr, ggplot2, scales, sp, DESeq2, tibble, pheatmap, Matrix) with one `install.packages()` /
  `BiocManager::install()` line for what is missing (P5); `fn.data_manifest(region)` with one row
  per bundle (presence tested on `cells.csv.gz` inside the folder) and per lasso file;
  `fn.check_inputs(manifest)` printing an OK/MISSING table and stopping before stage 1 when a
  required input is absent. `run_region_pipeline()` calls the check for its region first, so the
  author's flow gets it too, and the README quick start shows the standalone call.
- `library(tidyverse)` leaves `R/check_marker_spillover.R`; the file attaches dplyr, tibble and
  ggplot2, which are already in the required list (P4; this also closes E4 without changing how
  the file is sourced, decision 9).
- `print()` of figures wrapped in `if (interactive())` (P6).

### 4.3 Code changes per finding

| ID | Change | Backward compatibility |
|---|---|---|
| R1 | `run_region_pipeline()` writes `composition_per_sample.csv` (the `composition_df`, `row.names = FALSE`) | new file; same columns as the committed one |
| R2 | writes `de/<key>_res_full.csv`, `_res_panel.csv`, `_wilcox_res.csv` with `row.names = FALSE` | new files; the committed ones carry an extra unnamed index column, so the regenerated header differs by that column (recorded for pause 8) |
| R3 | writes `per_sample_depth.csv` (median `nCount_Xenium`, median `nFeature_Xenium` per sample, with group) | new file |
| R5, B4, T1 | writes `run_info.txt` (date, region, run label, marker panel names, percentile, categories, threshold, stage timings) and `session_info.txt` (`sessionInfo()`) per run | new files |
| R6 | `near_tie_categories(primary_result, percentile, min_share)` computes the list from the primary run; `R/run_analyses.R` prints it and warns if it differs from the recorded list, which stays in use | outputs unchanged |
| R7 | checkpoint gets a sidecar `<sample>_classified_prefilter.meta.txt` (panel names, lasso md5, QC cut-offs); on load a mismatch warns and recomputes; a missing sidecar (old checkpoints) warns and loads | the author's existing checkpoints still load |
| B2 | `process_hippocampus_sample()` warns when a lasso file holds more than one `Selection` (stop would change the author's run before their counts are known; decision 10) | outputs unchanged |
| B5 | group from `group_id` metadata in `build_cell_count_table()`, `build_pseudobulk_dds()` and `plot_de_heatmap()` (via `colData(dds)$group`) | identical for the six current names |
| B6 | `classify_atlas_index()` stores a per-panel status table (`scored`, `no_gene`, `below_floor`, best gene and rate, cells called) in `obj@misc$panel_status`; the pipeline writes `panel_status_per_sample.csv` | new file; the warnings stay |
| B7 | `slot = "counts"` becomes `layer = "counts"` in `AggregateExpression()` | Seurat v5 (already required by `JoinLayers()`) |
| B8 | `robustness_table()` handles `NA` differences explicitly | identical for complete tables |
| M8 | per DE cell type, save `de/<key>_pvalue_hist.png` and `de/<key>_dispersion.png` | new files |
| T2 | `run_region_pipeline(..., skip_embedding_and_de = FALSE)`; `R/run_analyses.R` keeps the default, so outputs are unchanged; README documents the toggle | default reproduces current behaviour |
| D1-D6, D9, M9, S1 | prose: protocol corrections, code comments, README limitation sentence, seeds and heatmap scale notes | no code |

### 4.4 Housekeeping

`.gitignore` gains `data/` (the default input folder). Nothing moves or leaves the index: E1 and
E2 are wontfix (decision 7), the `_withCP` folders stay as the README describes them (R4 to a
sub-issue).

### 4.5 Alternatives considered

- Making the author's paths the defaults, so they need no config file: rejected, because tracked
  code would still carry machine paths and the fresh-clone target would need edits to tracked
  files.
- A per-run `renv.lock`: rejected as a new dependency; `session_info.txt` per run records the same
  facts.
- Stopping (not warning) on a multi-selection lasso and on a near-tie list mismatch: deferred until
  the author's baseline answers whether either case occurs (decision 10).

### 4.6 Evidence commands

Parse sweep (`static_sweep.R`); dplyr-only recomputation of the stats and robustness tables from the
committed CSVs; `git grep` convention checks; and, on the author's machine, the Phase 8 runs with
`tools::md5sum()` before and after.

## 5. Implementation plan

Pause points (decision 8): after step 5 (the portability step), then run through; pause 8 before
any regenerated outputs are committed; pause 9 before any `git rm --cached` (none planned); pause
10 before the PR is marked ready.

- [x] 0. GitHub scaffold: issue #6, branch `6-audit-repo`, `.gitignore` additions (`config.local.R`, `hoard/`, `Rplots.pdf`). (housekeeping, part of step 1's commit)
- [x] 1. Commit this findings log and the science and workflow review drafts; push; open the draft pull request. (docs; `2f414b9`)
- [x] 2a. Minimal pre-baseline fix: the misplaced parenthesis in `R/build_stats_table.R:31`, because it blocks every run (B1; decision 4). Evidence: 19 of 19 files parse; the recomputed stats table matches the committed one. (code; `d17429b`)
- [ ] 2. Baseline: hand-off block posted on PR #7 for the author; record the result in the workflow review §9.1 when it arrives. (record)
- [x] 3. Phase 5 static review: registers completed; workflow review §3-§8; science review §5-§7. (record; `9dfcd45`)
- [x] 4. Pause 5 statuses, §4 fix design, this plan, decisions 5-10, §8 entries, sub-issues #8 and #9 opened and linked under #6. (record)
- [x] 5. Portability (`26d5d98`): track `.Rproj`; root anchor, SETTINGS block and `config.local.R` sourcing in `R/run_full_pipeline.R`; `R/_setup.R`; `config.local.example.R`; `setwd()` removed; `tidyverse` removed; prints guarded; old lassos to the example; `.gitignore` `data/`. (code; P1-P6, E3, E4) **Pause 7 after this step.**
- [ ] 6. Provenance outputs: composition, DE, depth and panel-status tables; `run_info.txt` with stage timings; `session_info.txt`; DESeq2 diagnostics. (code; R1, R2, R3, R5, B4, B6, M8, T1)
- [ ] 7. Robustness: group from metadata; `layer=`; NA-safe flips; lasso selection warning; checkpoint sidecar; near-tie list computed and compared; skip toggle. (code; B2, B5, B7, B8, R6, R7, T2)
- [ ] 8. Documentation fixes in protocols, comments and README. (docs; D1-D6, D9, M9, S1)
- [ ] 9. Baseline result recorded when the author's log arrives: workflow review §9.1 result line, stage table, comparison table; T1 measured. (record)
- [ ] 10. Phase 8 verification on the author's machine and, if the data can be placed, on the analyst's: fresh clone without data; fresh clone with data and `config.local.R`; interactive run; fingerprints compared; convention checks. (record) **Pause 8 if any tracked output is to be regenerated and committed.**
- [ ] 11. README and `CLAUDE.md` sync from the review documents; PR description; pre-ready check; §11 closing summary. (docs) **Pause 10.**

## 6. Change log

| Commit | Kind | IDs closed | Evidence | Documents updated |
|---|---|---|---|---|
| `2a1df70` | docs | none | none | science review §0-5 and workflow review §1-2 drafts (committed by the analyst on `main`) |
| `2f414b9` | docs | none | none | findings log created; `.gitignore` (`config.local.R`, `hoard/`, `Rplots.pdf`); issue number filled in |
| `d17429b` | code | B1 | 19/19 files parse; stats table recomputed from the committed composition table (see §3) | findings log §2.3, §5, §7 |
| `d767858` | record | none | none | findings log §0, §3 hand-off, §5, §6, §7 (decision 4), §10; workflow review §9.1 baseline notes |
| `9dfcd45` | record | none | static evidence in §3 (md5 of DE heatmaps; robustness and stats tables recomputed; DE tables recounted) | findings log §2 (B5-B8, D9, M10, M11); workflow review §3-§8; science review §5-§7 |
| `72a2a34` | record | none | none | findings log §2 statuses, §4, §5, §7 (decisions 5-10), §8, §10 (sub-issues #8, #9) |
| `26d5d98` | code | P1-P6, E3, E4 | 21/21 files parse; anchor stop from a temp dir; package-check stop at the root; OK/MISSING table and stop on a temp tree; override test; author examples identical() to the removed literals (commit body) | findings log §5, §6; workflow review §2 (step 8) |

## 7. Decisions log

1. **2026-10-07 (Holden Harris):** Mode `audit`, degraded; the analyst is a collaborator and the author reviews the PR; the whole repository is in scope; FLASH-MM out of scope; methods and statistics are reviewed in depth and recorded as questions for the author, never changed in this PR. (scope)
2. **2026-10-07 (Holden Harris):** The committed `outputs/` at 869af28 are the reference outputs. (scope)
3. **2026-10-07 (Holden Harris):** The baseline run is made by the author on their Mac from a hand-off block; nothing is recorded as measured until their log arrives. (scope)
4. **2026-10-07 (Holden Harris, with Claude Code):** Fix B1 (one misplaced parenthesis in `R/build_stats_table.R:31`) before the baseline, because at 869af28 `source("R/run_full_pipeline.R")` fails for everyone and no run of any kind is possible; the skill allows the minimum fix in that case. Mechanical: the recomputed stats table equals the committed one. Recorded in the workflow review §9.1 baseline notes. (B1)
5. **2026-10-07 (Holden Harris), pause 5:** Every mechanical row is fixed in this PR: P1-P6, R1-R3, R5-R7, B2 (as a warning), B4-B8, D1-D6, D9, M8-M9, S1, T1-T2, E3-E4. Each must show identical outputs or add new files only; any that turns out behavioural returns for a decision. (statuses in §2)
6. **2026-10-07 (Holden Harris), pause 5:** The scientific rows M1-M7, M10, M11, B3, D7 and D8 go to the author in one combined sub-issue, #8, with the arithmetic unchanged. (§8)
7. **2026-10-07 (Holden Harris), pause 5:** R4 and R8 (reruns the audit cannot perform) go to sub-issue #9. E1 (`plot_spatial_celltype()` unused) and E2 (five identical copies of the lasso overview, one under `docs/`) are wontfix: the function stays as an interactive helper and is documented as such; the `docs/` copy stays because code comments and the README refer to it. Nothing is moved or deleted.
8. **2026-10-07 (Holden Harris), pause 6:** Phase 7 pauses once, after step 5 (portability), then runs through to pause 8 (regenerated outputs), pause 9 (none planned) and pause 10 (ready for review).
9. **2026-10-07 (Holden Harris, with Claude Code):** E4 is closed by removing `library(tidyverse)` from `R/check_marker_spillover.R` (the file attaches dplyr, tibble and ggplot2 instead) rather than by sourcing it on demand, so the author's interactive calls keep working unchanged. (P4, E4)
10. **2026-10-07 (Holden Harris, with Claude Code):** B2 (a lasso file with more than one `Selection`) and R6 (a computed near-tie list that differs from the recorded one) warn rather than stop until the author's baseline shows whether either case occurs in the real files; converting either to a stop is a follow-up decision. (B2, R6)

## 8. Needs examination (for the original author)

All items are collected in sub-issue #8 (decision 6). The arithmetic stays as it is in this PR so
the outputs remain comparable; each entry gives the question, the evidence, the options and what
was deliberately not changed.

- **M1, circular gene selection (high).** P2ry13, Tmem119 and Sall1 both define microglia
  (`R/marker_panels.R:4`) and are tested within microglia (`R/de_candidate_panels.R:18`); Aldh1l1
  and Sox9 play the same double role for astrocytes. Evidence: the headline microglia result is
  P2ry13 (panel padj 0.015, log2FC -0.45), and the same gene helped decide which cells are
  microglia. Options: re-classify microglia with a panel that excludes P2ry13 and re-test; report
  the result as conditional on the selection; keep. Not changed: the panels and the candidate lists.
- **M10, state split against identity genes (medium).** The Astrocyte and Microglia Homeostatic
  state panels are the identity panels (`R/marker_panels.R:3-4, 98-101`), so "Active" means the
  state genes beat an already-high baseline. Options: disjoint state panels; document the split as
  marker-high versus marker-low in every region; keep.
- **M3, per-mouse detection floor (medium).** `R/classify_atlas_index.R:12-20` evaluates the 0.5 %
  floor per sample on the ROI cells; README "Open items" records BAM in `gi_2_d` at 0 then 264
  cells. The new `panel_status_per_sample.csv` (step 6) records what happened in each mouse.
  Options: pooled floor; fixed panel set per region; keep with the record.
- **M4, within-sample z-scoring (medium).** `scale(raw)` at `R/classify_atlas_index.R:28` removes
  between-mouse differences in overall score level. Option: a comparison run with z-scores pooled
  across the six mice.
- **M5, composition test (medium).** Welch per category on percentages that sum to 100 within a
  mouse, "Unclassified" as a category, BH across 15 (`R/build_stats_table.R`). Options: a
  proportions-aware test as a comparison; keep with the caveat stated in the README.
- **M2, sensitivity runs do not reach DE (medium).** `apply_confidence_filter()` rewrites
  `cell_type`; `build_pseudobulk_dds()` and `run_wilcox_degs()` select on `atlas_identity`; the
  eight DE heatmaps are byte-identical across the three hippocampus runs. Decision: intended, or
  should DE use the filtered labels in the sensitivity runs?
- **M6 and M11, sensitivity heuristics (low).** 15th percentile chosen on the optic tract RGC
  margins and reused; near-tie categories defined as more than 15 % below the threshold. The
  README's 10th/25th percentile check would exercise both.
- **M7, two-pass multiplicity (low).** Genome-wide BH and within-panel BH on the same p-values
  reported side by side (`R/run_pseudobulk_de.R:39-64`).
- **B3, no undetermined state (low).** `ifelse(StateA1 > StateB1, ...)` sends ties and no-signal
  cells to Homeostatic/Precursor (`R/classify_atlas_index.R:46`).
- **B2, lasso files with several selections.** Membership treats the file as one polygon
  (`R/process_sample_template.R:11-13`); the area sums per selection (`R/run_full_pipeline.R:
  73-79`). Step 7 adds a warning; the baseline hand-off asks for the count per file; a stop or
  per-selection membership is the author's call once the counts are known (decision 10).
- **D7 and D8, facts only the author has.** Sources of the six optic tract immune panels; the
  injury model, time point and section plane for the README.

## 9. Acceptance criteria

- [ ] Baseline run recorded before the first code change (workflow review §9 as-found).
- [ ] A fresh clone with only the local config added completes from the command line (workflow review §9 as-left).
- [ ] Run from the wrong folder stops with the anchor message; a missing input stops before any slow work with a table naming the file and where to get it (workflow review §9 as-left).
- [ ] `git grep -nE "[A-Za-z]:/|/Users/|/home/|/Volumes/|OneDrive|AppData" -- '*.R'` matches only commented examples; `git grep -n "windows(" -- '*.R'` matches nothing.
- [ ] Outputs regenerate identically given the same inputs, or every difference has a cause (workflow review §9 as-left table).
- [ ] `git ls-files --cached --ignored --exclude-standard` is empty.
- [ ] The science review, the workflow review and the README describe the stages, inputs, outputs and tested environment.
- [ ] The original author runs the branch with their local config and confirms the outputs.

## 10. GitHub record

Issue #6; branch `6-audit-repo` created from the issue's Development link (`gh issue develop`);
draft PR #7 against `main` with `Fixes #6` and the §9 checklist, opened 2026-10-07 after commit
`2f414b9`. Merge policy: merge commit, delete branch. Sub-issues, linked under #6 through the
sub-issues API on 2026-10-07: #8 (methods questions, §8; label `question`) and #9 (reruns R4 and
R8; label `enhancement`). The baseline hand-off is PR #7 comment `issuecomment-6050572603`.

## 11. Closing summary

Written at Phase 9.
