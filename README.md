# Xenium Spatial Transcriptomics: Cell-Type Composition and Differential Expression

Analysis pipeline for a Xenium in situ dataset (10x Xenium Prime Mouse 5K Pan
Tissue & Pathways panel) from 6 mouse brains: 3 TBI (group "gi") and 3 sham
(group "gii"). The same 6 mice are used for every brain region; only the
region-of-interest (ROI) lasso differs. Samples were processed and analyzed
under blinded labels (gi/gii) and unblinded on 2026-09-30. Abisambra Lab,
University of Florida.

Regions analyzed: Hippocampus and Optic Tract.

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
  plexus (excluded 2026-10-06 after driver-gene diagnostics). The Optic Tract uses 15: the 13
  plus six immune panels (Peripheral macrophage, T cell, Neutrophil,
  NK cell, B cell, Dendritic cell) and an RGC panel, minus five panels
  (Hippocampal pyramidal, Dentate gyrus granule, Ependymal/ciliated, RGC,
  Neutrophil). Reasons for each exclusion (single-gene driver
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
fell below that threshold in the unfiltered run (sens_neartie: Choroid plexus, GABAergic interneuron, Pericyte/vascular SMC,
B cell, Peripheral macrophage; threshold 0.0257 with the 15-panel set). The threshold and category list were derived from the margin
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

- Per-cell: atlas_identity, cell_type, classification_margin, runner_up.
- Per-sample cell-count table (cell_counts_per_sample.csv, every run):
  sample, n_cells, n_unclassified, pct_unclassified, pct_pass_qc20_10,
  group, region, run, lasso_area_um2, cells_per_1000um2. pct_pass_qc20_10 is
  the share of cells (already passing nCount >= 10 and nFeature >= 5) that
  would also pass a stricter nCount >= 20 and nFeature >= 10 cutoff.
- Composition tables and bar charts, UMAP, spatial maps and DE heatmaps in
  outputs/.
- Folders ending in _withCP hold the earlier 13-panel hippocampus runs
  (before Choroid plexus was excluded), kept for comparison.
- Folders ending in _withNeut hold the earlier 16-panel optic tract runs
  (before Neutrophil was excluded), kept for comparison.

## Known limitations

- n = 3 per group: with 3 versus 3, the smallest possible two-sided Wilcoxon
  p-value is 0.1; composition and DE results are exploratory.
- The smallest Optic Tract lasso (gii_3_e) has 156 cells; rare categories in
  that sample are unreliable.
- Astrocyte Active/Homeostatic (Gfap, Serpina3n) is compromised in white
  matter: optic nerve astrocytes have a high baseline of Gfap, and
  Serpina3n is largely oligodendrocyte-expressed in white matter. The
  buckets are best read as marker-high versus marker-low. The combined
  astrocyte total is reported alongside the split.
- Sensitivity to near-tie handling (Optic Tract, 15 panels; see
  outputs/optic_tract_near_tie_sensitivity.csv): the TBI-minus-sham
  difference changed sign for B cell and Unclassified. Peripheral
  macrophage kept its direction, but its difference grew from 0.60 to 3.05
  percentage points, so its size is not robust; Choroid plexus kept its
  direction but halved (-1.19 to -0.53). Mature oligodendrocyte,
  Homeostatic microglia, Homeostatic astrocyte, Precursor oligodendrocyte
  and Endothelial were essentially unchanged across all three runs. Before
  Neutrophil was excluded (16 panels; outputs/
  optic_tract_near_tie_sensitivity_withNeut.csv), five categories changed sign.
- Optic tract Neutrophil calls were Csf3r-driven (Ly6g 0% in all mice) and
  the panel was excluded; optic tract Choroid plexus calls are mixed (a
  cluster at one end of the lasso plus scattered low-confidence calls) and
  are kept but reported as not robust.
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

- Optic Tract: compare per-sample transcript depth between groups.
- UMAP coloured by sample, to check whether any island comes from a single
  mouse.
- gi_2_d BAM: 0 cells under the detection floor in an earlier run, 264 in
  the current run; confirm which cells the 0.5% floor is calculated on.
- Robustness to the percentile choice (10th and 25th, alongside the 15th).
- Optic Tract composition bar charts: redraw with the corrected plotting code
  (fixed 2026-10-07) when the optic tract runs are next repeated.
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

- R/run_full_pipeline.R: setup only (packages, functions, marker panels,
  sample paths and lassos). Safe to source():
  source("R/run_full_pipeline.R"). Defines run_region_pipeline(), the full
  pipeline end to end (classification, composition, UMAP, both DE analyses).
- R/run_analyses.R: the exact runs behind the reported results (primary and
  both sensitivity runs per region, robustness and stats tables, driver-gene
  diagnostics). Run one block at a time; each run is heavy. Per-sample
  checkpoints (*_classified_prefilter.rds) are gitignored.
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

## Data

Raw Xenium output and lasso CSV files are not tracked in this repository
(unpublished research data; see .gitignore). They are stored on the lab's
external drive.