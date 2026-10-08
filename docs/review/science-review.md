# Science review: hippocampus-xenium-analysis

## 0. Header

```
Repository: rebeccasakko/hippocampus-xenium-analysis (main at 869af28).
Original author: rebeccasakko (research assistant, Abisambra Lab, University of Florida).
Reviewers: Holden Harris (holden-harris), with Claude Code. Date: 7 Oct 2026.
Status: DRAFT (Phase 2 orientation), for the author's review. Not yet agreed with the author.
How to cite: "Science review of hippocampus-xenium-analysis, draft of 7 Oct 2026, docs/review/science-review.md".
Version history: 2026-10-07 first draft (issue #6).
```

Every statement below cites where it comes from: a file and line against commit 869af28, a README
section, a protocol section, a commit, or the earlier review document of 5 Oct 2026 (merged as
PR #3, removed from `main` in commit 683e6bc, still readable on branch
`origin/claude/project-thread-c42vt7` as `docs/review_2026-10-05.md`; it in turn quotes the
author's project notes, which are not in the repository). Statements the code does not explain are
marked "unexplained".

## 1. Summary

The pipeline turns 10x Xenium in situ transcriptomics data from six mouse brains (three traumatic
brain injury, three sham) into two products per brain region: a cell-type composition comparison
between groups, and a differential-expression (DE) comparison between groups within selected cell
types. Every cell inside a hand-drawn region-of-interest polygon is classified individually by
scoring it against small marker-gene panels and taking the best-scoring panel; no clustering is
used. Two regions have been run, the hippocampus (12 candidate panels) and the optic tract (16).
For each region three runs are reported: an unfiltered primary run and two sensitivity runs that
reassign low-confidence calls to "Unclassified". DE is tested with pseudobulk DESeq2 at the mouse
level (primary, n = 3 per group) and with a cell-level Wilcoxon test (secondary, exploratory). A
reader must know two things before using the outputs: with n = 3 per group every comparison is
exploratory, and the raw data and region polygons are not in the repository, so nothing here can be
rerun from a clone without them. (README lines 1-10, 12-45, 84-103; `R/run_full_pipeline.R:111-246`.)

## 2. Objectives

**Questions answered.** (a) Does the proportion of each cell type inside the region differ between
TBI and sham mice? (b) Within astrocytes, microglia, dentate gyrus granule cells and GABAergic
interneurons (hippocampus) or astrocytes and microglia (optic tract), which genes differ between
TBI and sham? (README "Cell-type classification", "Differential expression";
`R/run_full_pipeline.R:274-279, 324-327`.)

**Products.** Per region and run: a per-sample cell-count table, a per-sample composition table, a
composition bar chart, two UMAP figures, and per cell type a pseudobulk candidate-panel heatmap and
a Wilcoxon heatmap; per region: a near-tie sensitivity table; for the hippocampus: a Welch t-test
table, four sets of DE result tables, driver-gene diagnostic tables and figures, and a per-sample
transcript-depth table (README "Outputs"; `outputs/` tree; workflow review §7).

**Consumers.** The lab (a presentation on 2026-09-30 is cited at `R/marker_panels.R:61` and in
`protocols/classification_protocol.md:199`), a planned methods paper, and a planned hand-off to a
biostatistician (earlier review §1, from the author's notes; README "Planned analysis").

**What "correct" means.** For this audit: the committed files in `outputs/` at 869af28 are the
reference, and the tracked code should regenerate them from the same inputs, identically or with
every difference explained (decision at pause 1, findings log §1).

**Non-goals.** The proposed FLASH-MM mixed-model DE (README "Planned analysis") is noted and not
reviewed. The B cell / Cd22 line of work described in the earlier review is not in this repository.

## 3. Introduction and context

Six mice received either a CHIMERA closed-head impact (TBI) or sham handling and were sampled seven
days later (earlier review §1, from the author's notes; the repository itself does not state the
injury model or the time point: README lines 1-10 say only "TBI" and "sham"). Sections were
profiled with the Xenium Prime 5K Mouse Pan Tissue and Pathways panel (README lines 2-3). The
analysis ran under blinded labels `gi` and `gii` and was unblinded on 2026-09-30, with `gi` = TBI
and `gii` = sham (README lines 5-6).

The code grew out of an interactive analysis and was consolidated into one reusable function,
`run_region_pipeline()`, called once per region, with the setup (`R/run_full_pipeline.R`) separated
from the actual runs (`R/run_analyses.R`) on 2026-10-06 (commit a612776). Dated decisions are
written into the code comments and the two protocol documents; the main ones, in order:

| Date | Decision | Where recorded |
|---|---|---|
| 2026-09-24 | Six immune panels and an RGC sanity-check panel added for the optic tract | `R/marker_panels.R:17-26`; classification protocol table |
| 2026-09-28 | Hippocampal pyramidal, Dentate gyrus granule, RGC and Ependymal/ciliated dropped from the optic tract candidate list (single-gene drivers; bimodal RGC margins) | `R/marker_panels.R:29-83`; classification protocol "Region-specific candidate panel set" |
| 2026-09-28 | Confidence threshold changed from a fixed margin > 1 (hippocampus BAM only) to a percentile of the region's pooled margin distribution | `R/apply_confidence_filter.R:8-30`; classification protocol "Confidence-margin filtering" |
| 2026-09-28 | Composition rows zero-filled so an absent category counts as 0 % rather than being dropped from the mean | `R/build_composition_summary.R:4-19` |
| 2026-10-05 | Optic tract lassos redrawn by a lab neuroscientist; old lassos kept for reference only | `R/run_full_pipeline.R:301-322`; README "Region of interest" |
| 2026-10-06 | Hippocampus moved to the same unfiltered + two-sensitivity-run design as the optic tract | README "Confidence filtering"; `R/run_analyses.R:7-14` |
| 2026-10-06 | Choroid plexus excluded from the hippocampus candidate list after driver-gene diagnostics (Kl alone drove 66-98 % of calls) | `R/marker_panels.R:91-96`; `outputs/hippocampus_choroid_plexus_*` |
| 2026-10-07 | Pseudobulk DESeq2 made the primary DE analysis, Wilcoxon secondary | commit 07b9305; DE protocol "Why two methods" |

An earlier read-through review (5 Oct 2026, PR #3) listed 24 lettered issues; the ones since
addressed in the code are the hippocampus-specific marker list, the percentile filter for the
hippocampus, saved result tables, TBI/Sham plot labels and the `cur_data()` deprecation. The ones
still open are carried into this audit's findings log §2 and §8 with their origin noted.

## 4. Data sources

| Source | What it is | Coverage and vintage | Terms | Where cited |
|---|---|---|---|---|
| Xenium output bundles, one per mouse | 10x Xenium Analyzer output (cell-by-gene matrix, `cells.csv.gz` centroids, segmentation files) read by `Seurat::LoadXenium()` | 6 mice, one section each; run id prefix `0071470_`; panel: Xenium Prime 5K Mouse Pan Tissue and Pathways (about 5,000 genes) | Unpublished research data; lab external drive; not redistributable (README "Data") | `R/run_full_pipeline.R:252-260`; `R/process_sample_template.R:7`; `R/plot_lasso_overview.R:6` |
| Lasso polygons | Xenium Explorer selection exports (two comment lines, then `Selection, X, Y, Class, Color`), one file per mouse per region | Hippocampus: six files; optic tract: six files redrawn 2026-10-05 (`GI-1` to `GII-3`), older set retained unused | Not tracked; stored on the author's Desktop. Why they are not tracked is unexplained (they are small text files) | `R/run_full_pipeline.R:265-322`; README "Region of interest" |
| Identity marker panels | 20 panels of 1 to 8 genes; 12 used for the hippocampus, 16 for the optic tract | Derived from the Allen Institute whole-mouse-brain taxonomy (Yao et al. 2023, Nature) and checked against the 5K panel; BAM and the state panels are literature-based; the six immune panels' sources are a README TODO | Public literature | `R/marker_panels.R`; README "Cell-type classification"; classification protocol "Marker panels" |
| DE candidate panels | 5 to 7 genes per cell type, pre-specified | Astrocyte and microglia panels predate the file; DG from Arneson et al. 2018 and Truettner et al. 1999; GABA from eNeuro 2020, Cerebral Cortex 2022 and Front. Cell. Neurosci. 2025 | Public literature | `R/de_candidate_panels.R:14-40` |
| Methods references | Qiu et al. 2022 (Wilcoxon thresholds); Squair et al. 2021 (pseudoreplication) | | | `R/run_wilcox_de.R:7`; DE protocol |

Sizes, on-disk layout and how to obtain each input belong in the workflow review §3 (pending the
author's answers; see findings log §8).

## 5. Scientific methods

Stages are in run order as `run_region_pipeline()` executes them (`R/run_full_pipeline.R:111-246`).
Parameters are listed with the justification found in the repository or marked "unexplained".

### 5.1 Per-sample load and quality control

Each mouse's bundle is loaded with `LoadXenium(data_dir, fov = "fov")` and cells with fewer than
10 transcripts or fewer than 5 genes are removed (`R/process_sample_template.R:7-8`). Justification:
at very low counts normalisation amplifies noise and such cells are prone to spillover
(classification protocol step 1). A stricter 20/10 cut was evaluated and found to remove 0.1 to
0.2 % of cells; its pass rate is reported in every cell-count table as `pct_pass_qc20_10`
(`R/build_cell_count_table.R:11-16`; README "Outputs"). Counts are then log-normalised with
Seurat's default `NormalizeData()` (`R/process_sample_template.R:9`).

### 5.2 Region restriction

Cells whose centroid falls inside the lasso polygon are kept, using `sp::point.in.polygon()` on the
whole CSV as one polygon (`R/process_sample_template.R:11-16`). The lasso area used for cell
density is computed separately per `Selection` group and summed (`R/run_full_pipeline.R:73-79`),
so a file with more than one selection is treated as one polygon by the membership test and as
several by the area calculation. Methods question (findings log §8).

### 5.3 Per-cell identity classification

For each candidate panel, genes present in the object are kept; a panel is skipped in that sample
if its best single gene is detected in fewer than 0.5 % of cells (`min_detect_rate`,
`R/classify_atlas_index.R:12-20`); justification: a module score on so sparse a signal is noise
dominated (protocol step 4). For the remaining panels `Seurat::AddModuleScore()` is run on the
detected genes, the resulting scores are z-scored across cells within each panel and within each
sample (`scale(raw)`, line 28), and the panel with the highest z-score becomes `atlas_identity`
(line 40; `min_z = 0`, so a cell whose best z-score is negative is "Unclassified"). The gap between
first and second z-score is `classification_margin` and the second panel is `runner_up` (lines
33-39). Classification runs on the lasso-restricted, QC-passed cells (`R/process_sample_template.R:
16-23`), so the detection floor and the z-scores are computed on the cells inside the region, not
on the whole section. Choices that are methods questions rather than facts: the 0.5 % floor is per sample, so a
panel can be on in one mouse and off in another (README "Open items" records a BAM case in
`gi_2_d`); z-scoring within each sample removes between-sample differences in overall score level
before comparison; `min_z = 0` is unexplained.

### 5.4 State refinement

Astrocyte, Microglia and Oligodendrocyte lineage cells are rescored against two state panels each
(Active vs Homeostatic; Mature vs Precursor; `R/marker_panels.R:98-103`) and labelled by the larger
module score, ties going to the second label (`R/classify_atlas_index.R:41-63`). The Astrocyte and
Microglia identity panels are identical to their Homeostatic state panels by design (protocol
"double duty"), so for those two lineages the state score competes against the genes the cells were
selected on (methods question M10). No "cannot tell" state exists.

### 5.5 Confidence filtering (sensitivity runs only)

The primary run applies no filter. The pooled 15th percentile of `classification_margin` across all
six samples sets a threshold; cells below it are relabelled "Unclassified" in `cell_type` only,
either for every category (`sens_global`) or for categories where more than 15 % of cells fell
below the threshold in the primary run (`sens_neartie`; `R/apply_confidence_filter.R`;
`R/run_analyses.R:35-53`). The 15th percentile was chosen on the optic tract as the smallest
percentile capturing the RGC near-tie population (classification protocol "Confidence-margin
filtering"); its reuse for the hippocampus is a design choice for comparability (README
"Confidence filtering"). The near-tie category list is typed into the run script by hand from the
primary run's output (`R/run_analyses.R:44, 72`); no code derives it. Because the filter rewrites
`cell_type` and DE selects cells by `atlas_identity`, the sensitivity runs change composition but
not DE (verified: the eight DE heatmaps are byte-identical across the three hippocampus runs).

### 5.6 Composition summary and group test

Per sample, counts and percentages for every `cell_type` seen in any sample, with absent categories
zero-filled (`R/build_composition_summary.R:19-33`); per group, mean and SEM over n = 3
(lines 35-38). Welch's t-test on the percentage per category, Benjamini-Hochberg across the
categories tested (15 for the 12-panel hippocampus, including "Unclassified";
`R/build_stats_table.R`). Justification for Welch over Student: equal variances cannot be checked at
n = 3 (file header). Percentages are compositional (they sum to 100 within a sample), which the
per-category test ignores; "Unclassified" is tested as if it were a cell type. Methods question.
The near-tie robustness table reports the TBI-minus-sham difference in mean percentage under each
run and flags sign changes (`R/run_analyses.R:86-104`).

### 5.7 Embedding for display

`FindVariableFeatures`, `ScaleData`, `RunPCA(npcs = 30)`, `RunUMAP(dims = 1:30)` on the merged
object (`R/run_full_pipeline.R:207-210`); used only to draw cell types and the microglia state
split, never to assign labels (protocol step 6). `npcs = 30` is unexplained.

### 5.8 Pseudobulk DESeq2 (primary DE)

For each DE cell type, cells with that `atlas_identity` are aggregated to one count profile per
mouse (`AggregateExpression`), genes with total count below 10 are dropped, and DESeq2 is fit with
design `~ group`, contrast `gi` vs `gii` (`R/run_pseudobulk_de.R:16-44`). Two result sets: the
genome-wide table (BH across all genes passing the filter) and a candidate-panel table that keeps
the genome-wide p-values and re-adjusts BH within the pre-specified panel (lines 51-64; verified on
the committed tables: `padj_panel` equals BH over the panel's own p-values in all four cell types).
The protocol describes the genome-wide pass as a null result everywhere, but the committed table
lists Ddit4 in dentate gyrus granule cells at padj 0.004 (DE protocol "Method 1" vs "Results
summary"); see findings log D lens. Several candidate genes are also identity-panel genes for the
same cell type (P2ry13, Tmem119, Sall1 for microglia; Aldh1l1, Sox9 for astrocytes); selecting
cells partly on a gene and then testing it is a circularity question for the author (findings log
§8).

### 5.9 Wilcoxon FindMarkers (secondary DE)

Seurat `FindMarkers` with `test.use = "wilcox"`, `logfc.threshold = 0.25`, `min.pct = 0.1`, cells
as replicates, full gene panel (`R/run_wilcox_de.R:16-33`), matching Qiu et al. 2022. The
pseudoreplication caveat is stated in the file header and the protocol. `FindMarkers()` adjusts
p-values by Bonferroni over all genes tested, so the "padj < 0.05" gene counts reported for this
method are Bonferroni-adjusted, which no document states (finding D9).

### 5.10 Heatmaps

Both DE methods are displayed as per-mouse z-scored log2 normalised pseudobulk counts with free
column clustering (`R/plot_de_heatmap.R:43-78`); the protocol explains why clustering is left free.

### 5.11 Diagnostics

Driver-gene check (per-gene detection among a category's calls versus all other cells, runner-up
distribution, spatial location; `R/check_driver_genes.R`) and a marker-spillover check
(`R/check_marker_spillover.R`; run once interactively for Gfap in microglia, header lines 29-33).
Neither is part of `run_region_pipeline()`.

## 6. Results

Everything in this section is read from the committed outputs at 869af28 (the reference, decision
2). Verification against a fresh run is pending the author's baseline (workflow review §9.1).

### 6.1 Cells captured

| Region | Mouse | Group | Cells | Unclassified (%) | Lasso area (mm²) | Cells per 1,000 µm² |
|---|---|---|---|---|---|---|
| Hippocampus | gi_1_b | TBI | 5,425 | 3.8 | 2.13 | 2.55 |
| | gi_2_d | TBI | 6,023 | 2.7 | 2.24 | 2.69 |
| | gi_3_f | TBI | 5,644 | 4.8 | 2.18 | 2.59 |
| | gii_1_a | sham | 5,827 | 5.3 | 2.01 | 2.89 |
| | gii_2_c | sham | 6,283 | 4.5 | 2.42 | 2.60 |
| | gii_3_e | sham | 7,270 | 1.7 | 3.22 | 2.26 |
| Optic tract | gi_1_b | TBI | 609 | 6.4 | 0.13 | 4.57 |
| | gi_2_d | TBI | 763 | 10.1 | 0.23 | 3.36 |
| | gi_3_f | TBI | 866 | 14.8 | 0.21 | 4.20 |
| | gii_1_a | sham | 494 | 3.0 | 0.14 | 3.46 |
| | gii_2_c | sham | 310 | 11.3 | 0.08 | 4.10 |
| | gii_3_e | sham | 156 | 2.6 | 0.05 | 3.20 |

Source: `cell_counts_per_sample.csv` of the two unfiltered runs. All hippocampus samples pass the
stricter 20/10 QC at 99.8 to 99.9 %. The three sham optic tract lassos are smaller than the three
TBI ones (0.05 to 0.14 mm² against 0.13 to 0.23 mm²), so sham optic tract percentages rest on 156
to 494 cells (README "Known limitations").

**Figure: lasso overview** (`outputs/hippocampus_unfiltered/lasso_overview.png`). What to look for:
one red outline per section, sitting on the hippocampal formation; the sections differ in
orientation and in extent (gi_3_f carries a second piece of tissue above the main section, gii_1_a
a strip at the top edge), and none of the outlines touches those extra pieces.

### 6.2 Composition, hippocampus (12 panels, unfiltered run)

| Cell type | TBI mean % | Sham mean % | TBI minus sham (pp) | Welch p | BH p | Robust to near-tie handling |
|---|---|---|---|---|---|---|
| Dentate gyrus granule | 24.9 | 25.7 | -0.79 | 0.546 | 0.744 | yes |
| Hippocampal pyramidal | 18.1 | 19.2 | -1.17 | 0.243 | 0.744 | yes |
| GABAergic interneuron | 8.5 | 8.8 | -0.29 | 0.157 | 0.744 | yes |
| Oligodendrocyte (Mature) | 6.4 | 5.7 | +0.72 | 0.356 | 0.744 | yes |
| Pericyte/vascular SMC | 5.2 | 4.9 | +0.30 | 0.536 | 0.744 | sign change in `sens_global` (to -0.01) |
| Astrocyte (Active) | 4.7 | 4.8 | -0.06 | 0.864 | 0.925 | yes |
| Astrocyte (Homeostatic) | 4.9 | 4.0 | +0.84 | 0.125 | 0.744 | yes |
| Fibroblast/meningeal | 4.1 | 4.5 | -0.35 | 0.691 | 0.797 | sign change in `sens_global` |
| Endothelial | 4.4 | 3.9 | +0.46 | 0.252 | 0.744 | yes |
| Microglia (Homeostatic) | 4.1 | 4.0 | +0.14 | 0.645 | 0.797 | yes |
| BAM | 4.6 | 3.4 | +1.29 | 0.446 | 0.744 | sign kept, difference falls to +0.08 |
| Unclassified | 3.8 | 3.9 | -0.07 | 0.958 | 0.958 | |
| Oligodendrocyte (Precursor) | 3.1 | 3.3 | -0.23 | 0.457 | 0.744 | yes |
| Ependymal/ciliated | 2.2 | 3.2 | -1.01 | 0.529 | 0.744 | sign change in both sensitivity runs |
| Microglia (Active) | 0.9 | 0.7 | +0.22 | 0.243 | 0.744 | yes |

Sources: `composition_per_sample.csv` (means recomputed), `hippocampus_unfiltered_stats_table.csv`,
`hippocampus_near_tie_sensitivity.csv`. No category differs between groups after correction; the
smallest raw p is 0.125. The two largest raw differences (BAM, Ependymal/ciliated) are each driven
by one sham mouse (gii_2_c: BAM 0.8 %, Ependymal 5.7 %) and are the two that the sensitivity runs
remove or reverse.

**Figure: composition bar chart** (`outputs/hippocampus_unfiltered/composition_barchart.png`).
What to look for: bars are mean and SEM over three mice, filled points are the three TBI mice and
hollow points the three sham mice; in every category the six points overlap except BAM and
Ependymal/ciliated, where one hollow square (Sham 2) sits far from the other five.

**Figure: UMAP by cell type** (`outputs/hippocampus_unfiltered/umap_celltype_gi_vs_gii.png`). What
to look for: the labels were assigned without this embedding, so coherent islands of one colour
(dentate gyrus granule, hippocampal pyramidal, GABAergic, mature oligodendrocytes, microglia) mean
the per-cell calls agree with the expression structure; the central cloud mixes astrocytes,
vascular types, fibroblasts and Unclassified, which is where the near-tie calls live. The two
panels (sham left, TBI right) show the same islands.

### 6.3 Composition, optic tract (16 panels, redrawn lassos)

No per-sample composition table or stats table is committed for the optic tract; the README
reports the near-tie sensitivity only. From `optic_tract_near_tie_sensitivity.csv`: five categories
change sign between the unfiltered and the sensitivity runs (Unclassified, Neutrophil, B cell,
Choroid plexus, Peripheral macrophage; shifts of 2.5 to 6.7 percentage points) and are reported as
not robust; the largest differences (Mature oligodendrocyte -2.15, Homeostatic microglia -2.38,
Homeostatic astrocyte -1.23, Precursor oligodendrocyte +1.13, Endothelial -1.04) keep their sign.
The optic tract bar charts predate the plotting fix of 8d9b99f (finding R8).

### 6.4 Differential expression, hippocampus (unfiltered run; identical in the sensitivity runs)

| Cell type | Genes tested (pseudobulk) | Pseudobulk genome-wide, BH < 0.05 | Candidate panel, BH within panel < 0.05 | Wilcoxon, Bonferroni < 0.05 |
|---|---|---|---|---|
| Dentate gyrus granule | 3,939 | Ddit4 (padj 0.0044, log2FC +0.46) | Id2 (0.048, +0.37); Bdnf, Ntf3 null; Npy, Penk, Ptprn below the count filter | 23 |
| GABAergic interneuron | 3,689 | none | none (Pvalb +0.60, p 0.28) | 1 |
| Astrocyte | 3,621 | none | none (Sox9 p 0.028, padj 0.19; Mapt -0.07, p 0.61) | 9 |
| Microglia | 3,407 | none | P2ry13 (0.015, -0.45); Tmem119 (0.083, -0.30) | 4 |

Source: `outputs/hippocampus_unfiltered/de/`. Fold changes are TBI relative to sham. Genes
significant under both the mouse-level and the cell-level test: Ddit4 and Id2 (dentate gyrus,
higher in TBI) and P2ry13 (microglia, lower in TBI), as the DE protocol states. The original
question of the project, whether astrocyte Mapt changes, is answered null at the mouse level.

**Figure: microglia candidate-panel heatmap**
(`outputs/hippocampus_unfiltered/pseudobulk_heatmap_microglia.png`). What to look for: one column
per mouse, rows z-scored across the six; the asterisk on P2ry13 is the panel-level result; the
unforced column clustering does not separate the groups (gi-2-d pairs with gii-2-c, gi-1-b with the
two other shams), which is the honest reading of a one-gene result at n = 3.

### 6.5 Diagnostics

Choroid plexus (13-panel run): Kl detected in 66 to 98 % of the 245 to 428 cells called per mouse,
Folr1 and Otx2 in 0.6 to 2 %; runner-up Hippocampal pyramidal 46 % and Dentate gyrus granule 24 %;
the panel was removed (`hippocampus_choroid_plexus_*`). Ependymal/ciliated: Foxj1 23 to 96 %, Pifo
and Dnah12 0.3 to 19 %; runner-up spread across BAM, oligodendrocyte and neuronal panels; kept and
reported as not robust (`hippocampus_ependymal_*`). Per-sample depth: median transcripts per cell
613 to 639 (TBI) and 597 to 726 (sham), so depth does not separate the groups
(`hippocampus_per_sample_depth.csv`).

### 6.6 QC checks the code performs

Purkinje negative control: printed to the console per sample, not saved; "about 0 %" per the README.
The 20/10 QC pass rate is saved per sample (99.8 to 99.9 %). Panel exclusions by the detection floor
are printed as warnings only (finding B6). Verification against a fresh run: pending.

## 7. Discussion

### 7.1 Intended use

The outputs support an exploratory comparison of six mice. The repository's own documents say so,
and the results bear it out: no composition difference survives correction, and the mouse-level DE
yields one genome-wide gene and two panel genes. The pipeline's strengths are that every cell is
labelled by a stated rule, every decision is dated, and two DE methods with different replication
units are reported together.

### 7.2 Caveats, ranked by how much they could move a result

1. **Sample size.** Three mice per group bound every test; the smallest Welch p observed is 0.125
   and the composition tests have essentially no power for differences under a few percentage
   points (§5.6; finding M5).
2. **Circular gene selection in the headline DE result.** P2ry13, Tmem119 and Sall1 define the
   microglia population and are then tested within it; the same holds for Aldh1l1 and Sox9 in
   astrocytes (§5.8; M1). The P2ry13 result may be real, but it is not independent of how microglia
   were chosen.
3. **Per-sample panel switching.** The 0.5 % floor turns a panel on or off per mouse, and nothing
   records which panels were scored where (§5.3; M3, B6). A category present in two mice and absent
   in the third enters the mean as 0 % for the third.
4. **Within-sample z-scoring.** Comparing each mouse's cells to that mouse's own score
   distribution before the argmax can attenuate real between-group composition shifts (§5.3; M4).
5. **The sensitivity design does not reach DE.** The three runs differ only in `cell_type`, so the
   DE tables are the same in all three (§5.5; M2); the robustness claim applies to composition only.
6. **State labels for astrocytes and microglia** are a contrast against the identity genes, not a
   free classification (§5.4; M10); the README already treats the white-matter astrocyte split as
   marker-high versus marker-low.
7. **Compositional statistics.** Percentages that sum to 100 are tested one at a time, with
   "Unclassified" as a category (§5.6; M5).
8. **Threshold and category rules** for the sensitivity runs are heuristics chosen on one region
   (§5.5; M6, M11).
9. **Lasso files with more than one selection** would be handled inconsistently (§5.2; B2); the
   author's baseline hand-off checks the count per file.
10. **Reproducibility of the record.** Three output groups, including the DE tables that the
    results above are read from, were written by code that is not in the repository (findings
    R1-R3), and the 13-panel diagnostics cannot be regenerated from the script (R4). No package
    versions are recorded (R5).

### 7.3 Methods questions for the author

Mirrored in the findings log §8 with options and sub-issues; the arithmetic stays as it is until the
author decides: M1 (circularity: re-classify microglia without P2ry13 and re-test, or report the
result as conditional), M2 (whether DE should use the filtered labels), M3 (a pooled or fixed panel
set per region instead of a per-mouse floor), M4 (z-score pooled across the six mice as a
comparison), M5 (a proportions-aware composition test), M10 (state panels disjoint from identity
panels, or a documented "marker-high" reading), B2 (one-selection check or per-selection
membership), B3 (an undetermined state).

### 7.4 Suggested follow-ups

Save the per-sample panel-inclusion table and the DESeq2 diagnostics the protocol cites (B6, M8);
make the composition, DE and depth tables outputs of the pipeline (R1-R3); record the environment
(R5); run the 10th and 25th percentile sensitivity the README lists (M6, M11). Each is tracked in
the findings log §5 or a sub-issue.

## 8. References

- Yao Z. et al. (2023). A high-resolution transcriptomic and spatial atlas of cell types in the whole mouse brain. Nature 624, 317-332.
- Squair J. W. et al. (2021). Confronting false discoveries in single-cell differential expression. Nature Communications 12, 5692.
- Qiu et al. (2022). Front. Neurosci. (full citation to be confirmed with the author; cited at `R/run_wilcox_de.R:7`).
- Arneson D. et al. (2018). Single cell molecular alterations reveal target cells and pathways of concussive brain injury. Nature Communications 9, 3894.
- Truettner J. et al. (1999). J. Neurotrauma (BDNF/NT-3 mRNA after fluid percussion injury).
- Love M. I., Huber W., Anders S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology 15, 550.
- Hao Y. et al. (2024). Dictionary learning for integrative, multimodal and scalable single-cell analysis. Nature Biotechnology 42, 293-304 (Seurat v5).
- Earlier review: `docs/review_2026-10-05.md` on branch `claude/project-thread-c42vt7` (PR #3).

## 9. Appendices

9.1 Glossary and 9.2 figure list: to be added in Phase 5.
