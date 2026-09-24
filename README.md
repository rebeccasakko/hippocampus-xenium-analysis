# Hippocampus Xenium Analysis

Spatial transcriptomics analysis pipeline for hippocampal cell-type composition in a Xenium in situ dataset (mouse brain, blinded groups "gi" and "gii", n = 3 samples per group). Abisambra Lab, University of Florida.

## Methodology (v3: atlas-informed marker-index classification)

Every cell is individually scored against ~13 atlas-informed marker panels via `Seurat::AddModuleScore()`; the panel with the highest z-scored module score (argmax) determines top-level identity. Astrocyte, Microglia, and Oligodendrocyte-lineage cells are further sub-split into Active/Homeostatic or Mature/Precursor states using dedicated state-specific marker panels. No clustering (FindClusters/Louvain) is used.

A minimum marker-detection-rate floor (0.5% of cells) excludes any panel whose best single marker gene falls below reliable detection in a given sample, preventing module scores built on sparse/noisy signal from spuriously winning the argmax.

See `protocols/classification_protocol.md` for the full written protocol (per-step rationale, marker panel tables, and QC checks).

Note: a separate, earlier sub-analysis of the optic tract region (astrocyte/microglia/oligodendrocyte only) uses actual Seurat clustering (`FindClusters`/Louvain) as part of a different, cluster-then-gate-then-substate workflow. That is not part of this pipeline -- see the protocol doc for how the two are distinguished.

## Repository structure

- `R/marker_panels.R` — marker gene panel definitions
- `R/classify_atlas_index.R` — core per-cell classification function
- `R/process_sample_template.R` — `process_hippocampus_sample()`: per-sample load → QC → normalize → lasso-restrict → classify pipeline
- `R/build_composition_summary.R` — `build_composition_summary()`: aggregates per-sample counts into composition percentages and group-level summaries
- `R/plot_composition_barchart.R` — `plot_composition_barchart()`: the gii-vs-gi composition bar chart
- `R/plot_umap_celltype.R` — `plot_umap_celltype()` / `plot_umap_microglia_check()`: UMAP visualization of `cell_type` (not part of classification itself -- PCA/UMAP are computed separately, only for display) and the post-classification microglia state QC check
- `protocols/classification_protocol.md` — full written classification protocol
- `outputs/` — generated figures (e.g. UMAP PNGs); not raw data, safe to track

## Data

Raw Xenium output and lasso ROI boundary files are not tracked in this repository (see `.gitignore`) — too large, and this is unpublished research data. Data lives on the "Joe" external drive and the analyst's local Desktop.
