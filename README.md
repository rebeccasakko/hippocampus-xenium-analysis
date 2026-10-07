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
- Hippocampus uses 13 candidate panels. The Optic Tract uses 16: the 13
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
with every run.

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

Hippocampus. A fixed margin > 1 rule was applied to BAM only (BAM cells
below the threshold were reassigned to Unclassified). [TODO: this rule is not
yet in the repository; add and verify before sharing.] The threshold is not
transferable across regions because the margin depends on the number of
competing panels.

## Differential expression (TBI vs sham)

Two methods per cell type, reported together: pseudobulk DESeq2 (mouse-level
replication, n = 3 per group, candidate-panel restricted) as the conservative
secondary analysis, and Seurat FindMarkers Wilcoxon on the full gene panel
(cell-level replication) as the primary, more sensitive analysis. Cell-level
tests treat cells as replicates, which carries a pseudoreplication risk; the
pseudobulk analysis addresses it at lower power. Full details:
protocols/differential_expression_protocol.md.

Composition differences are tested per category with Welch's t-test and
Benjamini-Hochberg correction (R/build_stats_table.R), with the sample as the
unit (n = 3 per group).

## Outputs

- Per-cell: atlas_identity, cell_type, classification_margin, runner_up.
- Per-sample cell-count table (cell_counts_per_sample.csv, every run):
  sample, n_cells, n_unclassified, pct_unclassified, group, region, run,
  lasso_area_um2, cells_per_1000um2.
- Composition tables and bar charts, UMAP, spatial maps and DE heatmaps in
  outputs/.

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
- Marker spillover from neighbouring cells (segmentation) can contribute to
  signal for some genes.

## Open items

- Add the hippocampus BAM margin rule to the repository and verify that it
  reproduces the reported BAM percentages.
- Check Neutrophil (Csf3r) and Choroid plexus for single-gene drivers, as was
  done for the panels excluded from the optic tract.
- Re-derive the hippocampus confidence threshold on the same percentile basis
  as the optic tract.
- Exploratory astrocyte-only state analysis using a pre-specified,
  literature-derived gene set (sample as unit).

## Repository structure

- R/run_full_pipeline.R: run_region_pipeline(), the full pipeline end to end
  (classification, composition, UMAP, both DE analyses), called once per
  region. Edit the CONFIG/RUN sections (data paths) before running; run
  section by section the first time. Per-sample checkpoints
  (*_classified_prefilter.rds) are gitignored.
- R/marker_panels.R: marker panel definitions and region-specific candidate sets
- R/classify_atlas_index.R: per-cell classification
- R/process_sample_template.R: per-sample load, QC, normalize, lasso, classify
- R/apply_confidence_filter.R: margin filter (fixed and pooled-percentile versions)
- R/build_cell_count_table.R: per-sample cell-count table
- R/build_composition_summary.R: composition percentages and group summaries
- R/build_stats_table.R: Welch t-test with Benjamini-Hochberg correction
- R/check_marker_spillover.R: marker spillover diagnostic
- R/plot_composition_barchart.R, R/plot_umap_celltype.R,
  R/plot_spatial_celltype.R, R/plot_de_heatmap.R: figures
- R/de_candidate_panels.R, R/run_pseudobulk_de.R, R/run_wilcox_de.R: DE analyses
- protocols/: written classification and differential expression protocols
- outputs/: generated figures and tables

## Data

Raw Xenium output and lasso CSV files are not tracked in this repository
(unpublished research data; see .gitignore). They are stored on the lab's
external drive.
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
  