# Hippocampus Xenium Analysis

Spatial transcriptomics analysis pipeline for hippocampal cell-type composition and differential expression in a Xenium in situ dataset (mouse brain, blinded groups "gi" and "gii", n = 3 samples per group). Abisambra Lab, University of Florida.

## Cell-type classification (v3: atlas-informed marker-index)

Every cell is individually scored against ~13 atlas-informed marker panels via `Seurat::AddModuleScore()`; the panel with the highest z-scored module score (argmax) determines top-level identity. Astrocyte, Microglia, and Oligodendrocyte-lineage cells are further sub-split into Active/Homeostatic or Mature/Precursor states using dedicated state-specific marker panels. No clustering (FindClusters/Louvain) is used.

A minimum marker-detection-rate floor (0.5% of cells) excludes any panel whose best single marker gene falls below reliable detection in a given sample, preventing module scores built on sparse/noisy signal from spuriously winning the argmax.

Full protocol: `protocols/classification_protocol.md`.

Note: a separate, earlier sub-analysis of the optic tract region (astrocyte/microglia/oligodendrocyte only) uses actual Seurat clustering (`FindClusters`/Louvain) as part of a different, cluster-then-gate-then-substate workflow. That is not part of this pipeline -- see the protocol doc for how the two are distinguished.

## Differential expression (gi vs. gii)

Two methods are run per cell type and reported together: a pseudobulk DESeq2 analysis (mouse-level replication, n=3/group, candidate-panel restricted) as the secondary/conservative analysis, and a Seurat `FindMarkers` Wilcoxon test on the full gene panel (cell-level replication) as the primary/sensitive analysis, matching the method from the paper Sakthi shared. Full protocol, results summary, and method comparison (benefits/drawbacks, pseudoreplication caveat): `protocols/differential_expression_protocol.md`.

## Repository structure

- `R/marker_panels.R` — cell-type classification marker gene panel definitions
- `R/classify_atlas_index.R` — core per-cell classification function
- `R/process_sample_template.R` — `process_hippocampus_sample()`: per-sample load → QC → normalize → lasso-restrict → classify pipeline
- `R/build_composition_summary.R` — `build_composition_summary()`: aggregates per-sample counts into composition percentages and group-level summaries
- `R/plot_composition_barchart.R` — `plot_composition_barchart()`: the gii-vs-gi composition bar chart
- `R/plot_umap_celltype.R` — `plot_umap_celltype()` / `plot_umap_microglia_check()`: UMAP visualization of `cell_type` (not part of classification itself -- PCA/UMAP are computed separately, only for display) and the post-classification microglia state QC check
- `R/de_candidate_panels.R` — literature-backed candidate gene panels per cell type, used by the pseudobulk DE analysis
- `R/run_pseudobulk_de.R` — pseudobulk DESeq2 differential expression (secondary analysis): `build_pseudobulk_dds()`, `get_pseudobulk_results()`, `run_pseudobulk_panel_de()`
- `R/run_wilcox_de.R` — `run_wilcox_degs()`: Seurat FindMarkers Wilcoxon differential expression on the full gene panel (primary analysis)
- `R/plot_de_heatmap.R` — shared heatmap visualization for both DE methods: `plot_wilcox_heatmap()`, `plot_panel_heatmap()`
- `protocols/classification_protocol.md` — full written classification protocol
- `protocols/differential_expression_protocol.md` — full written DE protocol: both methods, results summary, benefits/drawbacks
- `outputs/` — generated figures (e.g. UMAP and heatmap PNGs); not raw data, safe to track

## Data

Raw Xenium output and lasso ROI boundary files are not tracked in this repository (see `.gitignore`) — too large, and this is unpublished research data. Data lives on the "Joe" external drive and the analyst's local Desktop.
