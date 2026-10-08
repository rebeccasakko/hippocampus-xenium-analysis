# hippocampus-xenium-analysis issue #6: audit, make runnable by three people, document methods and reproducibility

## 0. Status

```
Issue: #6 (https://github.com/rebeccasakko/hippocampus-xenium-analysis/issues/6)
Branch: 6-audit-repo        Pull request: pending (draft, opened in Phase 3 step 5)
Mode: audit (degraded: this session cannot run the code; runs are handed to the author)
Org profile: default
Reviewer: Holden Harris (holden-harris), with Claude Code        Original author: rebeccasakko
Line numbers refer to commit: 869af28 (code unchanged by 2a1df70, which added the two review drafts)
Last completed step: §5 step 1 (docs commit)        Next step: §5 step 2 (baseline hand-off)        Last updated: 2026-10-07
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
| P1 | `R/run_full_pipeline.R:252-260` | Six `data_dir` paths are literals on the author's lab drive (`/Volumes/Joe/...`); required by stage 0 and stage 1 | mechanical | high | open |
| P2 | `R/run_full_pipeline.R:265-272, 306-313` | Twelve lasso paths are literals under the author's `~/Desktop`; required by stage 0, 1 and 1c | mechanical | high | open |
| P3 | `R/run_analyses.R:17` | `setwd("~/Documents/GitHub/hippocampus-xenium-analysis")`: works for the author only; no root anchor anywhere, and the `.Rproj` file is not tracked | mechanical | high | open |
| P4 | `R/check_marker_spillover.R:36` | `library(tidyverse)` at source time, sourced by the setup file at `R/run_full_pipeline.R:66`; `tidyverse` is declared in neither the header list (`:38-39`) nor the README | mechanical | medium | open |
| P5 | `R/run_full_pipeline.R:42-50` | No package check: a missing package fails at the first `library()` with no list of everything missing; packages used by namespace (Matrix, grid) and Seurat's Xenium readers (e.g. `arrow`) are declared nowhere | mechanical | medium | open |
| P6 | `R/run_full_pipeline.R:202, 213, 218`; `R/plot_lasso_overview.R:25`; `R/check_driver_genes.R:53` | Figures `print()`ed unguarded; under `Rscript` they go to `Rplots.pdf` in the working directory | mechanical | low | open |

### 2.2 Reproducibility (R)

| ID | Where | Problem | Kind | Severity | Status |
|---|---|---|---|---|---|
| R1 | `outputs/hippocampus_{unfiltered,sens_global,sens_neartie}/composition_per_sample.csv`; `R/run_full_pipeline.R:196-205` | Tracked per-sample composition tables that no tracked code writes (`build_composition_summary()` returns them in memory only; added in commit fdde845 with no code change). `git log --all -S"composition_per_sample"` finds no writer | mechanical | medium | open |
| R2 | `outputs/hippocampus_unfiltered/de/*.csv` (12 files); `R/run_full_pipeline.R:223-238` | Tracked DE result tables with no writer in the code (added in 713b00d with no code change); their header carries a duplicated first column name (`"gene","gene"`), the mark of an interactive `write.csv()` with row names | mechanical | medium | open |
| R3 | `outputs/hippocampus_per_sample_depth.csv` | Tracked table (median transcripts and genes per sample) with no producing code anywhere in the repository | mechanical | low | open |
| R4 | `R/run_analyses.R:118` with `R/marker_panels.R:96` | The driver-gene diagnostics read `outputs/hippocampus_unfiltered_withCP/sample_objects_prefilter.rds`, a gitignored checkpoint of a 13-panel run; `hippocampus_markers` now has 12 panels and no code path produces the `_withCP` run, so the six `hippocampus_choroid_plexus_*` / `hippocampus_ependymal_*` outputs and the four `_withCP` folders cannot be regenerated from the script | mechanical | medium | open |
| R5 | repository-wide | No record of R, Seurat, DESeq2 or other package versions: no `sessionInfo()`, lock file or README line; `R/run_pseudobulk_de.R:26` and `R/plot_de_heatmap.R:57` depend on Seurat v5 turning `_` into `-` in aggregated sample names | mechanical | medium | open |
| R6 | `R/run_analyses.R:44, 72` | The near-tie category lists are typed by hand from the primary run's output ("categories where >15 % of cells fall below the threshold"); no code computes or checks them, so a rerun with changed inputs silently keeps the old list | mechanical | medium | open |
| R7 | `R/run_full_pipeline.R:137-141` | The per-sample checkpoint `<sample>_classified_prefilter.rds` is keyed on `output_dir` and sample only; a rerun into the same folder after a marker-list or lasso change loads the stale classification without warning | mechanical | medium | open |
| R8 | `outputs/optic_tract_newlasso_*/` | Three optic tract run folders hold no `lasso_overview.png` and their `composition_barchart.png` predates the plotting fix of 8d9b99f (README "Open items"); they were produced by an earlier state of the code than 869af28 | mechanical | low | open |

### 2.3 Bugs and fragility (B)

| ID | Where | Problem | Kind | Severity | Status |
|---|---|---|---|---|---|
| B1 | `R/build_stats_table.R:31` | `t.test(pct ~ group_id, data = pick(everything())$p.value, error = function(e) NA_real_)`: the closing parenthesis of `t.test()` is misplaced, so the file does not parse (`38:1: unexpected '}'`) and `source("R/run_full_pipeline.R")` stops at line 65 for every user. Found by the parse sweep on 2026-10-07 (18 of 19 files parse). The committed stats table predates the edit and recomputes from the committed composition table with the corrected call (max p difference 1.6e-14) | mechanical | high | open |
| B2 | `R/process_sample_template.R:11-13` vs `R/run_full_pipeline.R:76-78` | Point-in-polygon treats the whole lasso CSV as one polygon and ignores `Selection`; the area function splits by `Selection` and sums. A file holding more than one selection is joined into one wrong outline by the first and measured correctly by the second, with no check on the number of selections (carry-over from PR #3 review, item V) | mechanical (check) / scientific (if a file has several selections) | medium | open |
| B3 | `R/classify_atlas_index.R:46` | State refinement `ifelse(StateA1 > StateB1, ...)`: exact ties and cells with no signal from either state panel go to Homeostatic/Precursor; no "cannot tell" state (carry-over item K) | scientific | low | open |
| B4 | `R/run_full_pipeline.R:118` | `dir.create(output_dir, recursive = TRUE)` is correct, but every writer (`write.csv`, `ggsave`, `png`) overwrites silently, so a second run into a folder of a different run label mixes files; no run label or commit is written into the outputs | mechanical | low | open |

### 2.4 Documentation (D)

| ID | Where | Problem | Kind | Severity | Status |
|---|---|---|---|---|---|
| D1 | `protocols/differential_expression_protocol.md:313-316` vs `:381-383` and `outputs/hippocampus_unfiltered/de/dg_res_full.csv` | "A genome-wide test ... returned zero significant genes in every cell type" while the results table of the same document and the committed table list Ddit4 at padj 0.004 in dentate gyrus granule cells | mechanical | medium | open |
| D2 | `R/de_candidate_panels.R:4-7`; `R/plot_de_heatmap.R:80, 88`; `R/run_pseudobulk_de.R:2`; `R/run_wilcox_de.R:2` | Comments call pseudobulk the "secondary" and Wilcoxon the "primary" analysis; README and the DE protocol (since 07b9305) say the reverse | mechanical | low | open |
| D3 | `protocols/classification_protocol.md:40-42, 104-108, 276-280` | Says the hippocampus uses "~13 panels" and "the full, unrestricted `atlas_index_markers`", that the optic tract drops two panels, and that the hippocampus "is not (yet) using" the percentile filter; the code uses 12 and 16 panels and the percentile design for both regions (`R/marker_panels.R:79-96`; `R/run_analyses.R:35-53`) | mechanical | medium | open |
| D4 | `R/run_full_pipeline.R:166-168, 172-173` | Comment says the hippocampus "stays on its existing BAM-specific fix" and the primary run filters "only types with a diagnosed near-tie problem, e.g. BAM"; neither is true at 869af28 | mechanical | low | open |
| D5 | `R/build_stats_table.R:16-17` | "16 for Hippocampus's full panel": the 12-panel run tests 15 categories (`outputs/hippocampus_unfiltered_stats_table.csv`) | mechanical | low | open |
| D6 | README "Outputs" ("spatial maps ... in outputs/"); `R/plot_spatial_celltype.R:1` | No per-sample spatial cell-type map exists in `outputs/`, and `plot_spatial_celltype()` is called from nowhere; the only spatial figures are the two driver-gene diagnostics | mechanical | low | open |
| D7 | README "Cell-type classification" ("Optic-tract immune panels: [TODO: source]") | Six panels have no recorded source (carry-over item W) | mechanical | low | open |
| D8 | README lines 1-10 | The injury model (CHIMERA), the time point (7 days) and the tissue section plane are not stated anywhere in the repository; they are known only from the author's notes quoted in the PR #3 review | mechanical | low | open |

### 2.5 Methods and statistics (M)

Rows here are questions for the author (§8) unless marked otherwise; the arithmetic stays as it is.

| ID | Where | Problem | Kind | Severity | Status |
|---|---|---|---|---|---|
| M1 | `R/marker_panels.R:4` with `R/de_candidate_panels.R:18`; `R/marker_panels.R:3` with `:15` | Circularity: P2ry13, Tmem119 and Sall1 both define microglia and are tested within microglia; Aldh1l1 and Sox9 likewise for astrocytes. P2ry13 (padj 0.015, lower in TBI) is the headline microglia result (carry-over item C) | scientific | high | open |
| M2 | `R/apply_confidence_filter.R:2-5` with `R/run_pseudobulk_de.R:18-19`, `R/run_wilcox_de.R:18-19` | The confidence filter rewrites `cell_type` while DE selects on `atlas_identity`, so the sensitivity runs never test DE sensitivity; verified: all eight DE heatmaps byte-identical across the three hippocampus runs (carry-over item A) | scientific (design) | medium | open |
| M3 | `R/classify_atlas_index.R:12-20` | The 0.5 % detection floor is evaluated per sample, so a panel can be scored in one mouse and not another; README "Open items" records BAM in `gi_2_d` going from 0 to 264 cells between runs (carry-over item E) | scientific | medium | open |
| M4 | `R/classify_atlas_index.R:28` | Module scores are z-scored within each sample before the argmax, which removes between-mouse differences in overall score level and can attenuate a true composition shift (carry-over item D) | scientific | medium | open |
| M5 | `R/build_stats_table.R:26-38` | Welch t-tests on compositional percentages, one per category including "Unclassified", BH across 15; n = 3 per group. Alternatives built for proportions exist (e.g. `speckle::propeller`) (carry-over item G) | scientific | medium | open |
| M6 | `R/run_analyses.R:39, 67`; classification protocol "Confidence-margin filtering" | The 15th percentile was justified on the optic tract RGC margin distribution and reused for the hippocampus for comparability; README "Open items" already lists a 10th/25th robustness check | scientific | low | open |
| M7 | `R/run_pseudobulk_de.R:51-64` | Candidate-panel padj re-adjusts the genome-wide p-values within 3 to 7 pre-specified genes after the genome-wide pass has been read; the pre-specification is documented, but the two passes are reported side by side with no multiplicity accounting across them | scientific | low | open |
| M8 | `R/run_pseudobulk_de.R:29-31` | DESeq2 defaults: independent filtering, Cook's outlier handling and `rowSums >= 10` on three-vs-three pseudobulks; no dispersion plot or p-value histogram is saved although the protocol cites "flat p-value histograms" as evidence | mechanical (save diagnostics) | low | open |
| M9 | `R/plot_de_heatmap.R:50-51` | Heatmaps display `log2(normalized counts + 1)` z-scored by gene across six mice; with one mouse driving a gene the scale saturates at +-1.5 (`breaks`), which is fine for display but is not stated in the caption | mechanical (caption) | low | open |

### 2.6 Stochasticity (S)

| ID | Where | Problem | Kind | Severity | Status |
|---|---|---|---|---|---|
| S1 | `R/classify_atlas_index.R:23, 44-45` (`AddModuleScore`, control-gene sampling, Seurat default `seed = 1`); `R/run_full_pipeline.R:209-210` (`RunPCA`, `RunUMAP`, Seurat default `seed.use = 42`); `R/plot_composition_barchart.R:24` (`position_jitterdodge(seed = 1)`) | No `set.seed()` anywhere (`static_sweep.R`: 0 hits); every random step relies on a package default seed, which is deterministic for a given package version but is nowhere stated. UMAP coordinates will still move with uwot/Seurat versions | mechanical (document) | low | open |

### 2.7 Run time (T)

| ID | Where | Problem | Kind | Severity | Status |
|---|---|---|---|---|---|
| T1 | README "Repository structure" ("each run is heavy"); `R/run_analyses.R:3` | No measured run time, timer or log anywhere; claimed, not measured. Baseline pending on the author's machine | mechanical | low | open |
| T2 | `R/run_full_pipeline.R:207-238` | The sensitivity runs recompute PCA, UMAP and both DE analyses although DE is identical to the primary run (M2) and the UMAP embedding depends only on the merged counts; a skip toggle would save two of three heavy DE passes per region | mechanical (if outputs identical) | low | open |

### 2.8 Efficiency and artifacts (E)

| ID | Where | Problem | Kind | Severity | Status |
|---|---|---|---|---|---|
| E1 | `R/plot_spatial_celltype.R:1` | `plot_spatial_celltype()` is defined, sourced by the setup file, and called from nowhere (sweep: 0 call sites; README lists it as a figure function) | mechanical | low | open |
| E2 | `docs/hippocampus_lassos_overview.png` | Byte-identical (md5 `2086421d...`) to `outputs/hippocampus_unfiltered/lasso_overview.png` and to the three `_withCP` copies; five copies of one 363 KB figure tracked | mechanical | low | open |
| E3 | `R/run_full_pipeline.R:315-322` | `optic_tract_lassos_old`: six unused paths kept "for reference only"; documented, so a candidate for the config example's comments rather than live code | mechanical | low | open |
| E4 | `R/check_marker_spillover.R` | Diagnostic module sourced on every setup although nothing in the pipeline calls it (sweep: `check_marker_spillover()` 0 call sites); it is the only reason `tidyverse` is attached (P4) | mechanical | low | open |

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

To be pasted to the author after the parse fix question is settled (B1 blocks `source()`; see §5
step 2 and decision 4).

## 4. Recommended changes (fix design)

To be written in Phase 6 after pause 5.

## 5. Implementation plan

- [x] 0. GitHub scaffold: issue #6, branch `6-audit-repo`, `.gitignore` additions (`config.local.R`, `hoard/`, `Rplots.pdf`). (housekeeping, part of step 1's commit)
- [x] 1. Commit this findings log and the science and workflow review drafts; push; open the draft pull request. (docs)
- [ ] 2. Baseline: hand-off block to the author; record the result in the workflow review §9 as-found when it arrives. (record)
- [ ] 3. Phase 5 review: complete the eight registers; fill the workflow review §3-§8; correct the science review §5-§7. (record)
- [ ] 4. Pause 5 and Phase 6: statuses per row; §4 fix design; this plan rewritten with commit types and IDs. (record)
- [ ] 5+. Implementation steps, from Phase 6.

## 6. Change log

| Commit | Kind | IDs closed | Evidence | Documents updated |
|---|---|---|---|---|
| pending | docs | none | none | science review §0-5, workflow review §1-2, findings log created |

## 7. Decisions log

1. **2026-10-07 (Holden Harris):** Mode `audit`, degraded; the analyst is a collaborator and the author reviews the PR; the whole repository is in scope; FLASH-MM out of scope; methods and statistics are reviewed in depth and recorded as questions for the author, never changed in this PR. (scope)
2. **2026-10-07 (Holden Harris):** The committed `outputs/` at 869af28 are the reference outputs. (scope)
3. **2026-10-07 (Holden Harris):** The baseline run is made by the author on their Mac from a hand-off block; nothing is recorded as measured until their log arrives. (scope)

## 8. Needs examination (for the original author)

Populated at Phase 5 from the M rows and B2/B3; each item gets a sub-issue.

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

Issue #6; branch `6-audit-repo` created from the issue's Development link; draft PR against `main`
with `Fixes #6` and the §9 checklist (number recorded when opened). Merge policy: merge commit, delete
branch. Sub-issues for §8 items opened in Phase 6.

## 11. Closing summary

Written at Phase 9.
