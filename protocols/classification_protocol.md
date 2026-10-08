# Cell-Type Classification Protocol (v3: Atlas-Informed Marker-Index)

Abisambra Lab, University of Florida. Hippocampal Xenium in situ dataset,
blinded groups "gi" and "gii", n = 3 mice per group. This document is the
written protocol behind the pipeline implemented in `R/`; the code is the
source of truth, this is the explanation of what it does and why.

## Summary

Every cell is classified individually. There is no clustering step
(`FindNeighbors`/`FindClusters`/Louvain) anywhere in this pipeline --
identity is decided per cell, not per cluster. This is a deliberate
methodological choice, not an oversight: see "Why no clustering?" below.

(Note: a separate, earlier sub-analysis restricted to the optic tract
region -- covering only astrocytes, microglia, and oligodendrocytes --
*does* use real Seurat clustering (`FindClusters`, `seurat_clusters`) as
part of a cluster-then-gate-then-substate workflow. That is a distinct
methodology for a distinct analysis and is not part of this pipeline;
don't conflate the two.)

## Pipeline steps

1. **Load + QC** (`process_hippocampus_sample()`, `R/process_sample_template.R`).
   `LoadXenium()`, then filter to `nCount_Xenium >= 10 & nFeature_Xenium >= 5`.
   Cells below this floor are excluded before anything else: at very low
   counts, normalization amplifies noise from a tiny denominator, and these
   cells are disproportionately vulnerable to transcript spillover from
   neighboring cells during segmentation.

2. **Normalize.** Standard `NormalizeData()`.

3. **Spatial restriction.** The sample is intersected with a lasso ROI
   polygon (`point.in.polygon()` against the boundary CSV) to keep only
   cells inside the region of interest (e.g. hippocampus).

4. **Per-cell atlas classification** (`classify_atlas_index()`,
   `R/classify_atlas_index.R`, marker panels in `R/marker_panels.R`).
   For each atlas cell-type panel in the candidate list passed as
   `marker_list` (`hippocampus_markers`, 12 panels, since Choroid plexus was
   excluded on 2026-10-06; `optic_tract_markers`, 16 panels, which drops
   four panels and adds six immune panels -- see "Region-specific candidate
   panel set" below and `R/marker_panels.R`), every cell gets a
   `Seurat::AddModuleScore()` value for that panel. Scores are z-scored
   across the cells of that sample *within each panel*, and each cell's
   top-level identity (`atlas_identity`) is the panel with the highest
   z-scored score (argmax). The margin between the best and second-best
   panel is recorded (`classification_margin`), and the second-best panel
   is recorded too (`runner_up`) -- both useful for spotting ambiguous calls.

   A panel is excluded from scoring in a given sample if its single
   best-detected marker gene is expressed in fewer than 0.5% of cells
   (`min_detect_rate`) -- a module score built on that sparse a signal is
   noise-dominated and can win the argmax spuriously. The floor and the
   z-scores are computed on the cells inside the lasso that passed QC
   (classification runs after the spatial restriction), so the decision is
   per sample and per region; which panels were scored in which sample is
   written to `panel_status_per_sample.csv` in every run folder.

5. **State refinement for Astrocyte / Microglia / Oligodendrocyte lineage.**
   Cells classified into these three top-level identities are re-scored
   against a second, state-specific pair of panels (Active vs. Homeostatic
   for astrocyte/microglia; Mature vs. Precursor for oligodendrocyte), and
   `cell_type` records the winning state as a suffix, e.g.
   `"Microglia (Active)"`. Every other cell type's `cell_type` equals its
   `atlas_identity` unchanged.

6. **Visualization (not part of classification itself).**
   `classify_atlas_index()` does not compute PCA or UMAP -- identity is
   decided from raw module scores, not from an embedding. `RunPCA()` /
   `RunUMAP()` are run separately, only for visualizing the already-decided
   `cell_type` labels (`R/plot_umap_celltype.R`), never to inform the
   labels themselves.

## Marker panels

| Top-level identity | Marker panel (argmax gate) | Region |
|---|---|---|
| Astrocyte | Aldh1l1, Sox9 | Hippocampus + Optic Tract |
| Microglia | Tmem119, Sall1, P2ry13 | Hippocampus + Optic Tract |
| Oligodendrocyte lineage | Sox10, Olig1, Olig2 | Hippocampus + Optic Tract |
| Endothelial | Pecam1, Id1 | Hippocampus + Optic Tract |
| Pericyte/vascular SMC | Kcnj8, Myh11, Pdgfrb | Hippocampus + Optic Tract |
| Fibroblast/meningeal | Col1a1, Col1a2, Lum, Fn1 | Hippocampus + Optic Tract |
| Ependymal/ciliated | Foxj1, Pifo, Dnah12 | Hippocampus only -- excluded from Optic Tract candidate list, 2026-09-28 -- see below |
| Choroid plexus | Folr1, Kl, Otx2 | Hippocampus + Optic Tract |
| BAM | Mrc1, F13a1, Cd163 | Hippocampus + Optic Tract |
| Hippocampal pyramidal | Slc17a7, Fibcd1, Spink8, Hes1, Rab38, Pou3f1, Hcn1, Map3k15 | Hippocampus only -- excluded from Optic Tract candidate list, 2026-09-28 |
| Dentate gyrus granule | Slc17a7, Prox1 | Hippocampus only -- excluded from Optic Tract candidate list, 2026-09-28 |
| GABAergic interneuron | Gad1, Gad2, Slc32a1 | Hippocampus + Optic Tract |
| Purkinje (sanity check) | Pcp2 -- negative control; hippocampus should classify ~0% of cells here | Hippocampus + Optic Tract |
| Peripheral macrophage | Ccr2, Cd14, Ms4a7 -- distinguishes infiltrating from resident (BAM) macrophage | Added for Optic Tract, 2026-09-24 |
| T cell | Cd3e, Cd3d, Cd3g -- pan-T-cell (TCR complex) | Added for Optic Tract, 2026-09-24 |
| Neutrophil | Ly6g, Csf3r, Cxcr2 | Added for Optic Tract, 2026-09-24 |
| NK cell | Ncr1, Klrb1c, Gzmb, Klrk1, Prf1 | Added for Optic Tract, 2026-09-24 |
| B cell | Cd19, Cd79a, Cd79b, Ms4a1, Cd22, Pax5 | Added for Optic Tract, 2026-09-24 |
| Dendritic cell | Itgax, Cd83, Batf3, Xcr1, Clec9a | Added for Optic Tract, 2026-09-24 |
| RGC (sanity check) | Pou4f1, Sncg -- negative control; optic tract is axonal, no RGC soma present, should read ~0% | Added for Optic Tract, 2026-09-24; excluded from the Optic Tract candidate list, 2026-09-28 -- see below |

Every panel listed here is checked against the Xenium Prime Mouse 5K
pan-tissue panel before being added (the same citation-and-verification
discipline for every row). Most panels apply to every region unchanged --
"Region" mostly records when/why a panel was added, not a restriction. The
two exceptions are Hippocampal pyramidal and Dentate gyrus granule, which
are excluded from the Optic Tract candidate list specifically (see below).

### Region-specific candidate panel set (Optic Tract)

`optic_tract_markers` (`R/marker_panels.R`) is `atlas_index_markers` with
Hippocampal pyramidal, Dentate gyrus granule, RGC (sanity check) and
Ependymal/ciliated removed (the latter two for the reasons in the two
subsections that follow), and is what `run_region_pipeline()` passes as
`marker_list` when classifying Optic Tract samples. The hippocampus uses
`hippocampus_markers`: `atlas_index_markers` without the six immune panels
and RGC added for the optic tract, and, since 2026-10-06, without Choroid
plexus (12 panels; see `R/marker_panels.R:85-96` and
`outputs/hippocampus_choroid_plexus_*` for the driver-gene evidence).

This was not the first response to the problem. Both panels were initially
left in as candidates everywhere, on the assumption that the argmax
mechanism itself would suppress them in a region where they don't belong --
the same mechanism that gives the Purkinje and RGC sanity-check panels
~0% in the wrong region. Instead, the first Optic Tract composition chart
(2026-09-24) showed both reading at non-trivial, non-incidental
percentages (Dentate gyrus granule ~8.9-9.7%, Hippocampal pyramidal
~5.7-6.7%) -- too high to be background noise given neither cell type has
any anatomical basis for being in a white matter tract.

A two-part diagnostic (2026-09-28) distinguished two possible causes --
near-tie noise (the same failure mode already found and being addressed
for the RGC/BAM panels) versus genuine marker cross-reactivity (the same
failure mode already found for RGC's Sncg gene) -- and found the second:

- **classification_margin ruled out noise.** Cells classified as either
  panel had a median margin of 1.9-2.4, and 0% of them would flip to
  Unclassified under any candidate confidence-filter threshold up to 1.0.
  These were among the most confident calls the classifier made, not
  borderline ones -- so the pending classification_margin recalibration
  (needed for RGC/BAM) does not touch this problem at all.
- **Per-gene detection rates identified the driver gene in each panel.**
  For Hippocampal pyramidal, Fibcd1/Spink8/Rab38 were at 0% detection in
  every single Optic Tract sample, and Hes1/Hcn1/Map3k15 were low/spotty --
  but Pou3f1 was detected in 50-100% of these cells in every sample.
  Pou3f1 (Oct-6/SCIP) has a well-established role in oligodendrocyte/
  Schwann cell myelination and is biologically expected to be active in a
  heavily myelinated white matter tract; it is not hippocampus-specific.
  For Dentate gyrus granule, Slc17a7 was near-absent but Prox1 was at
  100% detection in every sample -- since it's a 2-gene panel, Prox1 alone
  was generating every one of these calls. Prox1 is not DG-restricted
  brain-wide (it also marks e.g. Prox1+ interneuron populations found
  scattered through adult white matter).
- **runner_up confirmed the biological direction of the error.**
  Oligodendrocyte lineage was the dominant runner-up for both panels
  (29% of Hippocampal pyramidal calls, 45% of Dentate gyrus granule
  calls) -- consistent with real oligodendrocyte-lineage cells losing the
  argmax to these false-positive calls, which also explains why dropping
  the two panels was expected to (and did) raise Oligodendrocyte lineage's
  share of the Optic Tract composition.

GABAergic interneuron was evaluated against the same question and kept as
an Optic Tract candidate: Gad1/Gad2/Slc32a1 are not hippocampus-restricted
markers to begin with (inhibitory interneurons are present throughout the
CNS, including along white matter tracts), and it did not show the
low-diversity single-gene-driver pattern found above.

### Ependymal/ciliated -- excluded from Optic Tract, 2026-09-28

Same failure mode as Hippocampal pyramidal/Dentate gyrus granule, caught
the same way. Foxj1 is detected in 65-100% of the cells this panel calls
in every Optic Tract sample, while Pifo and Dnah12 -- the genes that would
actually confirm ciliated-ependymal identity -- sit at 0-2.9% in every
sample. Real ependymal tissue expresses all three cilia-program genes
together; Foxj1 alone firing is more consistent with its broader tie to
motile ciliogenesis generally than to ependymal identity specifically.
This also argues against a competing hypothesis raised during review (that
the Optic Tract lasso was simply catching adjacent periventricular
tissue): true lasso contamination would be expected to carry all three
genes together as real tissue, not one gene in isolation. Choroid plexus
was evaluated against the same two-gene-coherence question and left
unresolved -- Otx2 and Folr1 both contribute somewhat (not one gene doing
all the work), so it doesn't show the same clean single-driver signature;
revisit with more time.

### RGC (sanity check) -- excluded from Optic Tract, 2026-09-28

Unlike Hippocampal pyramidal and Dentate gyrus granule, RGC's problem
turned out not to be a single clean explanation. Its classification_margin
distribution among Optic Tract cells labeled RGC is cleanly bimodal:
roughly half sit at near-zero margin (real near-tie noise against panels
like BAM/Dendritic cell -- the same failure mode already fixed for BAM in
hippocampus), and the other half sit at a very large margin (4.9-6.9) --
confidently, unambiguously called RGC, not a tie. There is a clean gap
between the two clusters, so no single confidence-margin threshold can
resolve both: any threshold that clears the noise half either does
nothing more for the confident half (if too low) or starts costing real
cells in unrelated panels like Astrocyte/Microglia for no additional RGC
benefit (if raised further).

The confident half is suspected, but not yet confirmed, to be driven by
Sncg alone -- Sncg has documented cross-reactivity with vascular/renal/
smooth-muscle cell types in the reference atlas, and Optic Tract has real
vasculature (Endothelial, Pericyte/vascular SMC, Fibroblast/meningeal)
that could plausibly supply confident false positives the same way Pou3f1
and Prox1 did for the two panels excluded above. This has not been
verified with a per-gene detection-rate split by margin group the way the
Hippocampal pyramidal/Dentate gyrus granule cause was.

Given time constraints ahead of a 2026-09-30 presentation, RGC was
excluded from the Optic Tract candidate list rather than left half-fixed
-- this is a deliberate, time-boxed decision, not a resolution. It gives
up the sanity check's original purpose for Optic Tract (catching a lasso
that strayed toward the retina/optic nerve head). **Revisit before final
submission**: confirm the Sncg driver-gene hypothesis, decide on a
corrected RGC panel (e.g. requiring both Pou4f1 and Sncg rather than
either alone, if Pou4f1 independently clears the detection floor), and
reinstate it as a sanity check once fixed.

| State refinement | Active / Mature panel | Homeostatic / Precursor panel |
|---|---|---|
| Astrocyte | Gfap, Serpina3n | Aldh1l1, Sox9 |
| Microglia | Trem2, Cst7, Lpl | Tmem119, Sall1, P2ry13 |
| Oligodendrocyte lineage | Mog, Mag, Cnp, Opalin | Pdgfra, Cspg4, Gpr17, Vcan |

Note that the Astrocyte and Microglia gate panels are identical to their
Homeostatic-state panels ("double duty") -- no marker survived verification
as truly state-independent for these two families, so the same genes serve
both roles. Oligodendrocyte lineage is the exception: Sox10/Olig1/Olig2 are
pan-lineage and state-independent, so the gate panel is kept separate from
both the Mature and Precursor state panels.

## Why no clustering?

This is a deliberate departure from the more common cluster-then-annotate
workflow (Louvain cluster, then label each cluster by its dominant marker
expression). Per-cell argmax classification was chosen instead so that
identity is decided independently for every cell, rather than inherited
from cluster membership -- this avoids diluting or misassigning cells at
cluster boundaries, and avoids the resolution parameter as a hidden
analyst-facing lever on the final cell-type counts. The tradeoff is that
per-cell scores are noisier than cluster-averaged scores, which is why
step 4's detection-rate floor and the classification-margin diagnostic
exist -- they're the safeguards that compensate for not having the
noise-averaging effect of clustering.

## Outputs

- `atlas_identity` -- top-level cell-type call (12 categories in the hippocampus, 16 in the optic tract, plus Unclassified)
- `cell_type` -- `atlas_identity`, with Astrocyte/Microglia/Oligodendrocyte
  lineage further split into their state-refined subtypes
- `classification_margin` -- gap between the winning and runner-up panel's
  z-scored module score (larger = more confident call)
- `runner_up` -- the second-place panel, for reviewing ambiguous cells

## Confidence-margin filtering (region-derived threshold)

`apply_confidence_filter_pooled()` (`R/apply_confidence_filter.R`) derives
a `classification_margin` floor as a percentile of the POOLED margin
distribution across all samples in a region, rather than reusing a fixed
raw number across regions. This matters because `classification_margin`
is mechanically panel-count-dependent: more competing candidate panels
shrinks typical winning margins regardless of classification quality, so
a threshold calibrated in one region's panel-count context doesn't
transfer to a region with a different panel count.

This was discovered 2026-09-28 when Hippocampus's validated BAM threshold
(`margin > 1`, from the detection-floor noise fix documented above) was
first reused as-is for Optic Tract and badly over-corrected: Unclassified
jumped to ~50%, and even confidently-called categories like Astrocyte and
Microglia lost ~37% of their cells at that threshold -- well inside their
normal, non-ambiguous margin range (1st-quartile margins of 0.67-0.83),
not a genuine confidence floor for Optic Tract's panel count. A percentile
of the region's own observed distribution was adopted instead so the
DERIVATION RULE stays uniform across regions even though the resulting
raw threshold differs by region -- analogous to how a multiple-testing
correction threshold is expected to shift with the number of tests
performed, not stay fixed.

Optic Tract uses the 15th percentile of its own pooled margin distribution
(`run_region_pipeline(..., confidence_percentile = 0.15)`) -- chosen as
the smallest percentile that already captures the full near-tie noise
population found in the RGC (sanity check) diagnostic before RGC was
excluded (RGC's margin distribution was cleanly bimodal, with the
noise-driven half fully separated from the confident half by a wide gap;
15% was the smallest percentile threshold that captured all of the
noise-driven half with minimal cost elsewhere). Since 2026-10-06 the
hippocampus uses the same design: an unfiltered primary run plus two
sensitivity runs at its own pooled 15th percentile (0.1745 with the
12-panel set; the earlier BAM-only `margin > 1` rule is retired). In both
regions the filter is applied only in the sensitivity runs, never in the
primary run; the "near-tie" category list for `sens_neartie` is the set of
categories in which more than 15% of cells fall below the threshold in the
primary run, computed by `near_tie_categories()` and recorded in
`R/run_analyses.R`. The filter rewrites `cell_type` only; the differential
expression step selects cells by `atlas_identity`, so the sensitivity runs
change composition but not DE (audit issue #6, finding M2; whether that is
the intended design is a question for the author in issue #8). Note that
the 15% figure was chosen on the optic tract and reused for the hippocampus
for comparability; a 10th/25th percentile check is listed under README
"Open items".

## QC checks

- **Purkinje sanity check**: the hippocampus should contain ~0% Purkinje
  cells; a nonzero percentage flags a contaminated ROI or a marker-panel
  problem.
- **Microglia UMAP highlight check** (`plot_umap_microglia_check()`,
  `R/plot_umap_celltype.R`): after classification, Microglia (Active) and
  Microglia (Homeostatic) cells are highlighted against every other cell
  type (dimmed grey) on the UMAP, split by group, as a visual check that
  the state split looks spatially coherent rather than scattered/noisy
  before it's trusted in any downstream analysis (e.g. differential
  expression panels).
