# Xenium Pipeline: Project Context for Claude Code

*Compiled 2026-10-05 from the "Xenium Pipeline" Cowork session and the running project notes. Use this as starting context for coding work in Claude Code. See the last section for what this document does not capture.*

---

## 1. Project at a glance

- **Who / where:** Bec, research assistant in Dr. Abisambra's neuroscience lab, University of Florida. Day-to-day direction comes from Sakthivel Ravi ("Sakthi"). Dr. Abisambra is the PI.
- **Data:** 10x Genomics Xenium spatial transcriptomics, mouse brain, **Xenium Prime 5K Mouse Pan Tissue & Pathways panel (~5000 genes)**.
- **Design:** 6 samples, 3 TBI and 3 sham.
  - Analysis ran **blinded** as `gi` / `gii`.
  - **Unblinded 2026-09-30: `gi` = TBI, `gii` = sham.** Group identity can now be referenced directly going forward.
- **Sample IDs:** `gi_1_b`, `gi_2_d`, `gi_3_f`, `gii_1_a`, `gii_2_c`, `gii_3_e`. Raw folders are prefixed `0071470_`.
- **Tools:** R, Seurat, tidyverse, RStudio on a Mac. Xenium Explorer is used to view raw data and draw ROI lassos. Code lives in a private GitHub repo, `rebeccasakko/hippocampus-xenium-analysis`, managed through GitHub Desktop.
- **Raw data location:** lab external drive "Joe", at `/Volumes/Joe/Becca and Sakthi Files/Xenium_Data/Group_i` (or `Group_ii`)`/<sample folder>`. All 6 folders sit flat, with no nested duplicates.
- **Goal / deadline:** Submit for publication by end of 2026. Methodological decisions have to survive peer review. Eventual handoff of the repo to a biostatistician.
- **Original question:** identify astrocytes and determine whether they express Tau (*Mapt*). It has since grown into region-wide cell-type composition and DE analysis for **Hippocampus** and **Optic Tract (OT)**.

## 2. Repo layout (as of last known state)

```
R/marker_panels.R               # atlas_index_markers + state panels + optic_tract_markers
R/classify_atlas_index.R        # classify_atlas_index(), min_detect_rate = 0.005
R/process_sample_template.R     # process_hippocampus_sample(); later refactored to run_region_pipeline()
R/build_composition_summary.R   # build_composition_summary() (zero-fill bug fixed, see section 5)
R/plot_composition_barchart.R   # plot_composition_barchart() (labels rotated 90 degrees)
R/build_stats_table.R           # build_stats_table(): Welch t-test + BH correction
run_pseudobulk_de.R, run_wilcox_de.R, plot_de_heatmap.R   # generic across regions
run_full_pipeline.R             # includes optic_tract_de_cell_types list
README.md                       # documents the v3 methodology
```

- `apply_bam_confidence_filter` and `audit_detection` existed in Bec's R session but were not yet added to the repo as of the last record. Check whether they were added since.
- `apply_confidence_filter_pooled()` was added for OT and is wired in via `run_region_pipeline(confidence_percentile = ...)`.
- Repo not yet shared. Plan is to add Sakthi and Dr. Abisambra as collaborators.

## 3. Methodology (v3, the "atlas-informed marker index")

**Approach:** Full reference-atlas label transfer was rejected (too large, Python-only for her setup). Instead, every cell is classified by module scores (`AddModuleScore`) against a marker "index" built from the Allen Institute whole-mouse-brain taxonomy (Yao et al. 2023, `cl.df_CCN202307220.xlsx`) and checked gene by gene against her actual Xenium panel. The class with the highest score wins (argmax).

**QC:** `nCount_Xenium >= 10 & nFeature_Xenium >= 5` is the baseline. A strict regime of `>= 20 & >= 10` retains 99.8-99.9% of cells and gives essentially identical composition. QC stringency was ruled out as a driver of findings. The one exception is a BAM detection-floor gap in `gi_2_d`, which strict QC resolved (256 BAM cells recovered vs 0).

**Established identity panels (Hippocampus, in production):**

| Identity | Genes |
|---|---|
| Astrocyte | Aldh1l1, Sox9 |
| Microglia | Tmem119, Sall1, P2ry13 |
| Oligodendrocyte lineage | Sox10, Olig1, Olig2 |
| Endothelial | Pecam1, Id1 |
| Pericyte / vascular SMC | Kcnj8, Myh11, Pdgfrb |
| Fibroblast / meningeal | Col1a1, Col1a2, Lum, Fn1 |
| Ependymal / ciliated | Foxj1, Pifo, Dnah12 |
| Choroid plexus | Folr1, Kl, Otx2 |
| BAM | **Mrc1, F13a1, Cd163** (fixed from Mrc1+Maf, since Maf is non-specific) |
| Hippocampal pyramidal | Slc17a7, Fibcd1, Spink8, Hes1, Rab38, Pou3f1, Hcn1, Map3k15 |
| Dentate gyrus granule | Slc17a7, Prox1 |
| GABAergic interneuron | Gad1, Gad2, Slc32a1 |
| Purkinje (sanity check, should be ~0%) | Pcp2 |

**State refinement (only within Astrocyte, Microglia and Oligodendrocyte lineage):**
- Astrocyte Active = Gfap, Serpina3n. Homeostatic = Aldh1l1, Sox9.
- Microglia Active (DAM) = **Trem2, Cst7, Lpl**. The generic myeloid genes Cd68, Itgam and Aif1 were dropped. Homeostatic = Tmem119, Sall1, P2ry13.
- Oligodendrocyte Mature = Mog, Mag, Cnp, Opalin. Precursor = Pdgfra, Cspg4, Gpr17, Vcan.
- **Excluded from the astrocyte panel on purpose:** Gja1, Aqp4, Slc1a3 (reactivity/TBI-associated, not state-independent). Sakthi also asked for Gja1 to be dropped.

**Optic Tract additions (confirmed present in the 5K panel):**
- Peripheral macrophage: Ccr2, Cd14, Ms4a7
- T cell: Cd3e, Cd3d, Cd3g
- Neutrophil: Ly6g, Csf3r, Cxcr2
- NK cell: Ncr1, Klrb1c, Gzmb, Klrk1, Prf1
- B cell: Cd19, Cd79a, Cd79b, Ms4a1, Cd22, Pax5
- Dendritic cell: Itgax, Cd83, Batf3, Xcr1, Clec9a

**Detection-rate floor:** In `classify_atlas_index()` (`min_detect_rate = 0.005`), a panel's best single marker is the gene with the highest fraction of cells with a nonzero count. If that rate is below 0.5%, the whole panel is excluded from scoring for that sample. Only one gene has to clear the bar. A panel whose markers have zero detected transcripts must be excluded before scoring, otherwise `AddModuleScore` gives a degenerate score that can spuriously win the argmax.

**Confidence filter:** `classification_margin` is the winning score minus the runner-up score. A near-zero margin means a near-tie.
- **Hippocampus:** BAM kept only if margin > 1, otherwise reclassified to Unclassified. This fixed BAM overcalling (3-6% of cells vs <1-2% expected).
- **Optic Tract:** the fixed margin > 1 cutoff does not transfer, because margins shrink mechanically with panel count (13 vs ~20 panels). Applying it pushed Unclassified to ~50%. OT uses a **region-derived percentile threshold** (15th percentile of that region's pooled margin distribution) via `apply_confidence_filter_pooled()`, set with `confidence_percentile = 0.15`. Hippocampus is still on the fixed rule, with a to-do to re-derive it percentile-style before submission.

## 4. Region-specific decisions

**Hippocampus:** the full 13-panel set above plus state refinement. Lasso CSVs live on the Desktop as `hippocampus_<sample>.csv`. The gi_3_f file has a typo in its name (`hipocampus_gi_3_f.csv`).

**Optic Tract:**
- Drop **Hippocampal pyramidal** and **Dentate gyrus granule** as candidates. They are single-gene false positives: Pou3f1 (myelination-associated) and Prox1 (not DG-restricted) drive them alone, while the truly region-specific genes never fire. Margins were high, so this is not near-tie noise.
- Drop **Ependymal/ciliated**. Foxj1 drives it alone (65-100% detection), while Pifo and Dnah12 sit at 0-3%.
- Drop **RGC (sanity check)**. Its margin distribution was bimodal (near-tie noise plus a confidently called half, suspected Sncg cross-reactivity), so no single threshold could fix it. This was a deliberate, time-boxed tradeoff. It loses the check for a lasso straying toward the retina/optic nerve head, so **revisit before final submission**.
- Keep **GABAergic interneuron** (no problem seen).
- **Choroid plexus** is murky (Otx2 and Folr1 contribute, Kl lags) and unresolved.
- **White-matter astrocyte panel:** Gfap and Serpina3n are both compromised in white matter. Vimentin is not a safe substitute. Serpina3n is mostly oligodendrocyte-expressed there, so signal may be segmentation spillover. Candidates C3, Lcn2 and Cxcl10 are unvalidated for white matter. The panel redesign is parked, and the limitation is documented. Region-specific marker panels are now an explicit goal.
- Astrocyte and Microglia calls in OT are not meaningfully affected by margin thresholds in the useful range (~1-6% flipped at the 15th-25th percentile).

## 5. Bugs found and fixed (do not reintroduce)

1. **Zero-fill bug in `build_composition_summary()`.** `table(obj$cell_type)` returns only levels present in a sample, so absent categories got no row. This shrank n_samples and biased group means upward, while the chart label still said "n = 3". It is fixed by taking the union of labels across all samples and tabulating each sample against the full factor levels. Verified it never affected Hippocampus (`filter(n_samples < 3)` returned 0 rows), so the reported Hippocampus p-values stand.
2. **Degenerate module score** when a marker panel has zero detected counts (see section 3).
3. **BAM silently dropped** in `gi_1_b` and `gi_2_d` because the best marker fell under the 0.5% floor. Fixed by the better panel (Mrc1/F13a1/Cd163) plus strict QC.
4. **Unrotated % labels** collided on the OT chart (34 labels). They are now rotated 90 degrees with extra bottom margin, in the shared function.
5. **Stale labels in lasso CSVs** (`gi_2_b`, a leftover "Optic Tract gi_3_f", and two "Unclassified"). All were confirmed harmless.
6. **Memory:** `LoadXenium` can hit R's 24GB vector limit (`R_MAX_VSIZE`). Workaround: fully quit and reopen RStudio. Intermediate objects were once lost this way, so **save `.rds` to disk after every step**. Per-sample classified objects are in `~/Desktop/R Files/`.

## 6. Results so far

**Hippocampus (gi/TBI vs gii/sham, n = 3 per group):**
- Welch t-test with BH correction across 16 categories. **Nothing is significant** (all p_adj_BH >= 0.614). The closest raw p-values were Hippocampal pyramidal (0.066) and Astrocyte Homeostatic (0.078).
- Active microglia p = 0.51. The p-value history is 0.65 (original), then 0.089 (an interim result with a known-flawed BAM panel, not adopted), then **0.51 (both corrected panels, adopted)**. The defensible panel was chosen over the more favorable p-value.
- Unclassified is ~3.8-6.1%. The Purkinje check reads ~0%. BAM is <1% post-fix.
- Files: `hippocampus_group_comparison_barchart_v6.png` (latest seen) and `hippocampus_group_comparison_stats_table_v2.csv` on the Desktop. The chart has % labels and per-sample points, as Sakthi requested.

**Optic Tract (old lassos):**
- Mature oligodendrocytes are dominant (~17-18% in the first pass), as expected for white matter.
- **B cell: 0% in all 3 gii vs 2.78-7.74% in gi.**
- **Captured cells per sample: gi 620 / 692 / 612; gii 425 / 602 / 168.** `gii_3_e` is a clear low-count outlier. Sample size alone does not fully explain the B cell result for gii_1_a and gii_2_c.
- OT lassos were harder to draw because the boundary is less anatomically defined.

**Differential expression (Wilcoxon `FindMarkers`, plus pseudobulk as a secondary method):**
- Hippocampus and OT heatmaps exist for Astrocyte and Microglia. OT Mature Oligodendrocyte shows real signal.
- OT GABAergic interneuron: no significant genes (underpowered).
- OT Oligodendrocyte Precursor: one hit, Gfap (log2FC 1.43, p_adj 0.004), probably segmentation spillover. Astrocyte Gfap detection is near ceiling in both groups (98.8% vs 92.8%) with no fold change. Excluded from the deck. Gfap is not in the astrocyte classification panel, so there is no circularity.

## 7. Current thread: B cells and Cd22 (latest work)

- Research focus has shifted to **B cells / Cd22**, building on the OT B cell finding and open **Cd22/microglia spillover** questions.
- **New refined OT lassos** from Sakthi (meant to include reactive/infiltrating cells at the border) are in the Desktop folder `B6_Xenium_Opt`. Files are `GI-1/2/3_coordinates.csv` and `GII-1/2/3_coordinates.csv`, in Xenium Explorer export format (2 comment lines, then `Selection,X,Y,Class,Color`). **The GI-N/GII-N to sample ID mapping is not yet confirmed with Sakthi.**
- The OT pipeline is **not to be rebuilt from the old lassos**. Wait for the new ones.
- **Immediate planned step:** get raw and normalized Cd22 transcript counts per sample (gi vs gii) via point-in-polygon subsetting plus QC, **without** full atlas-index classification yet.

**Framing agreed in the latest discussion:**
- Xenium here is a **screening tool** for a candidate, spatially grounded signal. It is not meant to carry the B cell result to publication-grade statistical weight on its own at n = 3.
- **Two-paper plan:**
  1. A **methods paper** that is squarely Bec's: pipeline, marker-panel validation, QC and confidence-margin work, and the lasso-size and spillover diagnostics. The method should stand on its own without the B cell result.
  2. The **B cell / immune finding** is handed off to someone with immunology expertise. It should be documented with its limitations: n = 3, the gi vs gii lasso-size disparity, the open spillover question, and what would be needed next (more samples, orthogonal validation such as flow or IHC on an independent cohort).
- An offered next step: draft a short handoff summary for the B cell signal, in Bec's own voice, separate from the methods paper.

## 8. Open items

- Rebuild the OT pipeline on the **new lassos**, after confirming the GI/GII sample mapping.
- Simple Cd22 raw and normalized count per sample (first), then the full classification.
- Refine the OT lasso with Sakthi to include border reactive/infiltrating cells.
- Build **region-specific marker panels** (e.g. different Astrocyte markers for OT vs Hippocampus).
- Second DE pipeline from the paper Sakthi flagged, and report overlap between pseudobulk and Wilcoxon. Two candidate routes, undecided:
  - External validation against a mouse TBI Visium dataset (TBI-seq; Zenodo 10.5281/zenodo.15313202, GEO GSE282909). It is a different platform, so not like-for-like.
  - Internal upgrade with FLASH-MM, a mixed-effects DE method that treats mouse as a random effect.
- Re-derive the Hippocampus confidence threshold percentile-style.
- Revisit the OT RGC sanity check, Choroid plexus, the white-matter astrocyte panel, and the oligodendrocyte secondary-panel literature pass.
- Add `apply_bam_confidence_filter` and `audit_detection` to the repo if wanted. Keep the repo current for the biostatistician handoff.
- Draft the B cell handoff summary and the methods-paper write-up.

## 9. Working conventions for Claude Code

- **Paste R code directly in chat**, not as a downloadable file. Bec prefers this.
- Bec has limited R background and likes the *why* explained alongside the code (scripts vs functions, saving and reloading, etc.).
- **Transparency matters.** Methods should be visible enough that another researcher could catch a mistake. State assumptions and limitations plainly. Prefer the methodologically defensible choice over the more favorable p-value.
- Never guess or decode group identity from the data. Treat both groups symmetrically in any rule or threshold.
- Save intermediate `.rds` files after each step. Restart RStudio between heavy sessions.
- Present both Hippocampus and OT numbers apples-to-apples (same stats function, same thresholds logic).

## 10. What this document does not capture

The detailed turn-by-turn record of the most recent session (the Cd22 / B cell / spillover diagnostics and the two-paper discussion) could not be recovered after context compaction. Section 7 reflects the conclusions and framing only. Any specific Cd22 counts, spillover test results, code snippets, or figure versions from that stretch should be re-checked in the R session, the repo's commit history, and Desktop files before being relied on. The latest Hippocampus chart version (`v6`) and the OT chart versions (`v2`-`v4`) are the ones visible on the Desktop, and newer ones may exist.
