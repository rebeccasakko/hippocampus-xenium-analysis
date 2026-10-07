# Differential Expression Protocol: gi vs. gii

Abisambra Lab, University of Florida. Applies to the classified hippocampal
Xenium dataset (`atlas_identity` / `cell_type`, see
`protocols/classification_protocol.md`). Two independent methods are run
per cell type and reported together, by design -- neither is treated as a
replacement for the other.

## Why two methods

The Seurat `FindMarkers` Wilcoxon method follows Qiu et al. (2022, *Front. Neurosci.*).
Pseudobulk DESeq2 is the **primary** analysis and the Wilcoxon method the
**secondary** (exploratory) analysis. They differ on exactly one axis that matters a
great deal at n = 3 mice per group: the unit of biological replication.

## Method 1 (primary): Pseudobulk DESeq2 -- `R/run_pseudobulk_de.R`

All cells of one cell type, from one mouse, are aggregated into a single
expression profile before testing (`Seurat::AggregateExpression()`). The
real sample size the statistics see is n = 3 vs. n = 3 -- the standard,
defensible unit of replication. A genome-wide test (every gene passing
`rowSums(counts) >= 10`) returned **zero significant genes in every cell
type** -- confirmed to be a true null (flat p-value histograms, not a
correction artifact) and confirmed to hold even when isolating the Active
astrocyte subpopulation specifically (ruling out dilution-by-pooling as an
explanation).

Because genome-wide testing is underpowered at this sample size, a second
pass restricts to a small, pre-specified, literature-justified candidate
gene list per cell type (`R/de_candidate_panels.R`) and re-corrects BH
across just that panel (`run_pseudobulk_panel_de()`). This is hypothesis-
driven, not data-dredged -- every gene in every panel has a citation.

**Strength:** conservative, defensible unit of replication; the standard
a reviewer or statistician expects for a "does the treatment have an
effect" claim.
**Weakness:** very low power at n = 3/group; a real effect can easily go
undetected (type II error), which is likely part of why the genome-wide
pass came back null everywhere.

## Method 2 (secondary, exploratory): Seurat FindMarkers, Wilcoxon, full gene panel -- `R/run_wilcox_de.R`

Every individual cell is treated as an independent observation -- no
aggregation, no candidate-list restriction, run across every gene passing
Seurat's own default thresholds (`logfc.threshold = 0.25`, `min.pct =
0.1`, matching the paper's method exactly). Effective sample size balloons
into the tens of thousands of cells, even though it is still only 3
animals per group biologically.

**Strength:** far more sensitive; surfaces genes no one would have thought
to pre-specify (e.g. Ddit4 recurring across several cell types); doesn't
require deciding on a candidate list in advance.
**Weakness:** pseudoreplication -- treating thousands of cells nested
within just 3 animals as independent samples is a well-documented source
of false positives in single-cell differential expression (Squair et al.,
2021, *Nature Communications*, "Confronting false discoveries in
single-cell differential expression"). A handful of cells from one animal
can drive a "significant" result that does not generalize across the
group. More hits does not mean more true positives.

## Visualization -- `R/plot_de_heatmap.R`

Both methods share one display convention so their figures are directly
comparable: z-scored (`log2(normalized counts + 1)`, then per-gene
z-score), one column per MOUSE (never per cell, even for
Wilcoxon-selected genes -- gene selection can come from the cell-level
test, but display always happens at a level a person can read and that
means something biologically), gi/gii annotation bar, free (unforced)
column clustering, and asterisk significance labels (* padj<=0.05, **
<=0.01, *** <=0.001, **** <=0.0001).

Free clustering is deliberate: it is independent statistical
corroboration. A cluster split that has no knowledge of group labels can
either confirm or fail to confirm the hypothesized gi/gii split. Forcing
columns into gi/gii blocks would hide that information and risks reading
as confirmation bias. A heatmap built on very few genes (e.g. a single
significant gene) has correspondingly little clustering evidence to offer
-- six values sorted along a line will always produce *some* split,
whether or not it means anything -- so read the clustering pattern's
weight in proportion to how many genes are actually driving it.

## Results summary (padj < 0.05)

Hippocampus, 12-panel primary run (2026-10-07). Full tables in
`outputs/hippocampus_unfiltered/de/`. Fold changes are TBI relative to
sham (positive = higher in TBI), confirmed against group means.

| Cell type | Wilcoxon (full panel) sig. genes | Pseudobulk genome-wide | Pseudobulk candidate panel |
|---|---|---|---|
| Dentate gyrus granule | 23 | **Ddit4** (padj = 0.004, higher in TBI) | **Id2** (padj = 0.048, higher in TBI) |
| Astrocyte | 9 | Null | Null (Sox9 raw p = 0.028, padj = 0.19) |
| Microglia | 4 | Null | **P2ry13** (padj = 0.015, lower in TBI) |
| GABAergic interneuron | 1 | Null | Null |

**Findings that survive both methods:** Ddit4 (dentate gyrus) is the only
gene significant in the genome-wide mouse-level test, and is also
significant under Wilcoxon. P2ry13 (microglia) is significant in the
candidate panel and under Wilcoxon, lower in TBI, consistent with reduced
expression of a homeostatic microglial marker after injury. Id2 (dentate
gyrus) is significant in the candidate panel (padj = 0.048, just below the
threshold) and under Wilcoxon, higher in TBI. Every other Wilcoxon hit has
no independent corroboration from the mouse-level test, consistent with
the pseudoreplication caveat above.

## Bottom line

The two methods are not in conflict so much as answering different
questions at different levels of rigor: "does this hold up mouse-to-mouse"
(pseudobulk) versus "is there any detectable cell-level shift, however it
arose" (Wilcoxon). Report both, and lead with the genes that survive both:
Ddit4, P2ry13 and Id2.