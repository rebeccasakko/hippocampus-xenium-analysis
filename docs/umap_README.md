# UMAP plots: what they are, what they show, and how they are made

<!--
=====================================================================
UPDATE HERE FIRST
When the UMAP figures are re-made, change the lines in this box and the
"Figure register" table below. The images further down the page all read
their file paths from the [fig-...] lines here, so they update on their own.
Then follow "How to update this document" at the bottom of the page.
=====================================================================
-->

[fig-celltype]: ../outputs/optic_tract/umap_celltype_gi_vs_gii.png
[fig-microglia]: ../outputs/optic_tract/umap_microglia_check_gi_vs_gii.png

| Document status | |
|---|---|
| Last updated | 2026-10-06 |
| Figures describe | Optic Tract, **old lassos** (provisional, see issue 1) |
| Figures made | 2026-10-02 (commit `d751537`) |
| Code described | `R/plot_umap_celltype.R`, Stage 4 of `R/run_full_pipeline.R` (as of commit `66bbd8d`) |
| Hippocampus UMAPs | Not in the repository yet |

### Figure register

One row per UMAP image in the repository. When a figure is replaced, edit its row; when a new one is added, add a row (and a `[fig-...]` line in the box above if it should appear on this page).

| ID | File | Region | Lassos | Run | Labels shown | Made | Status |
|---|---|---|---|---|---|---|---|
| fig-celltype | `outputs/optic_tract/umap_celltype_gi_vs_gii.png` | Optic Tract | Old | Old-lasso run (15% filter on all cell types, *inferred*) | `cell_type` | 2026-10-02 | Provisional, to be replaced |
| fig-microglia | `outputs/optic_tract/umap_microglia_check_gi_vs_gii.png` | Optic Tract | Old | Same run as above | Microglia states vs everything else | 2026-10-02 | Provisional, to be replaced |

---

## Contents

1. [Why we make these plots](#1-why-we-make-these-plots)
2. [What the plots show](#2-what-the-plots-show)
3. [Where the data comes from and how the code works](#3-where-the-data-comes-from-and-how-the-code-works)
4. [Issues to review](#4-issues-to-review)
5. [How to update this document](#5-how-to-update-this-document)
6. [Change log](#6-change-log)

---

## 1. Why we make these plots

### The study

We want to know how gene activity in the mouse brain changes 7 days after a traumatic brain injury (TBI). Six mice were used: three received a closed-head injury with the CHIMERA model (group **gi** = TBI) and three were handled the same way without injury (group **gii** = sham). The groups were analysed under blinded labels and unblinded on 30 September 2026.

Each brain section was read with **10x Xenium** using the **Xenium Prime 5K mouse panel**. For every cell, Xenium counts how many copies of each of about 5,000 genes it carries, *and* records where the cell sits in the tissue. We draw an outline (a "lasso") around a brain region and keep only the cells inside it. The regions so far are the **hippocampus** and the **optic tract** (a bundle of nerve fibres carrying visual signals, which CHIMERA is known to injure).

### How cells get their labels

Before any plot is drawn, each cell is given a cell-type label (astrocyte, microglia, oligodendrocyte and so on). This is done one cell at a time by scoring the cell against short lists of "name-tag" genes (marker genes) and picking the best match. No clustering is used. Full details are in [`protocols/classification_protocol.md`](../protocols/classification_protocol.md).

### What a UMAP is, in plain words

Each cell has about 5,000 numbers (one per gene). That is impossible to look at directly. A **UMAP** squashes those 5,000 numbers down to just two, so every cell can be drawn as one dot on a flat page. It tries to keep cells with **similar overall gene activity close together**. So:

- **Dots close together** = cells whose gene activity looks alike.
- **Separate islands** = groups of cells that look different from the rest.
- **The axes (`umap_1`, `umap_2`) have no units or meaning.** Only "near" versus "far" matters, and even then only locally: the distance between two far-apart islands is not meaningful.
- A UMAP is a **picture, not a measurement**. It is not a statistical test.

### Why we make them here

In this project the UMAP is **a visual check, not part of the analysis**. The cell labels are decided first, from marker genes. The UMAP is then computed separately from *all the genes* and coloured by those labels. If the labels are sensible, cells with the same label should mostly sit together, because cells of one type should look alike across thousands of genes, not just their few marker genes. ([`protocols/classification_protocol.md`](../protocols/classification_protocol.md), step 6: the UMAP is "never" used "to inform the labels themselves".)

There are two plots:

1. **Cell-type UMAP:** every cell coloured by its label. Purpose: a broad sanity check of the labels.
2. **Microglia check UMAP:** only microglia coloured (active vs homeostatic, i.e. "reacting" vs "resting"), everything else grey. Purpose: check that the active/resting split looks coherent before the microglia gene-change results are trusted. Microglia are the brain's immune cells and gave the clearest TBI signal in the optic tract.

Both plots are split into two side-by-side panels, **gii (sham) on the left and gi (TBI) on the right**. Both panels share one UMAP layout (it was computed on all six mice together), so the same spot on the left and right panels means the same kind of cell.

---

## 2. What the plots show

> **Version note.** The observations below describe the figures listed in the register at the top (Optic Tract, old lassos, made 2026-10-02). They were read by eye from the images. When the figures are replaced, rewrite this section.

### Figure 1: every cell, coloured by cell type

![Figure 1. Optic Tract UMAP coloured by cell type, sham (gii) left, TBI (gi) right][fig-celltype]

**What you are looking at:** each dot is one cell from inside the optic tract outline. Left: the three sham mice pooled. Right: the three TBI mice pooled. Colour = the cell's label (`cell_type`), including the active/resting splits for astrocytes and microglia and the mature/precursor split for oligodendrocytes.

**What it shows:**

- **One large island on the right** holds most cells. It is a **mix** of mature oligodendrocytes (purple), microglia (blues), astrocytes and Unclassified cells (reds and pinks). These cell types do not form clearly separate islands of their own.
- **A small island at the top** is mostly precursor oligodendrocytes (lilac).
- **Smaller islands on the left** are mostly blood-vessel and lining cells (endothelial, fibroblast/meningeal, pericyte; greens) and a few others.
- **An island at the bottom (around `umap_2` = -10)** is large in TBI and almost empty in sham. It contains a mix of labels, including microglia, immune cells and Unclassified.
- The TBI panel looks **denser** because the TBI outlines captured more cells (see issue 5).

**What it means (and does not mean):**

- Some labels group well (precursor oligodendrocytes, the vessel cells), which supports those labels.
- The mixing in the large island means that, *when all ~5,000 genes are considered*, many astrocytes, microglia and oligodendrocytes do not look clearly different from each other here. This could be because the labels are noisy, because optic tract cells carry few RNA molecules each, or because of RNA from neighbouring cells (see issues 2 and 3). It is **not yet explained**.
- The TBI-heavy bottom island is the most striking feature. It **could** be injury-related cells, but it **could also** come mostly from one mouse or one tissue section. The current plot cannot tell these apart, because it pools mice within each group (see issue 4). It should not be read as a TBI finding yet.

### Figure 2: microglia check

![Figure 2. Optic Tract UMAP highlighting active (black) and homeostatic (green) microglia, sham left, TBI right][fig-microglia]

**What you are looking at:** the same layout as Figure 1. Only microglia are coloured: **black = active** (higher Trem2, Cst7, Lpl), **green = homeostatic / resting**. All other cells are light grey.

**What it shows:**

- **Sham:** active microglia sit mostly at the lower tip of the large island (around `umap_1` = 3, `umap_2` = -2), with resting microglia spread above them. That is a fairly tidy split.
- **TBI:** many active microglia sit in the TBI-heavy bottom island, mixed with resting ones. Fewer sit at the tip of the large island.
- In both groups, resting microglia are **spread thinly through the whole large island** rather than forming one tight group.

**What it means:** the protocol asks whether the active/resting split looks "coherent rather than scattered". The honest reading is **partly**. In sham the active cells cluster; in TBI they mostly sit in the bottom island, whose origin is not yet known. The thin spread of resting microglia across the large island fits the mixing seen in Figure 1. This is a judgement by eye, so it should be recorded as such and not used as evidence on its own.

---

## 3. Where the data comes from and how the code works

### Data provenance

| What | Where | In the repo? |
|---|---|---|
| Raw Xenium output (one folder per mouse) | Lab drive "Joe": `/Volumes/Joe/Becca and Sakthi Files/Xenium_Data/Group_i/` and `Group_ii/`, folders `0071470_<sample>` | No (unpublished data) |
| Optic tract lassos, **old** (used for the current figures) | Bec's Desktop: `optic_tract_<sample>.csv` | No (`*.csv` is ignored by git) |
| Optic tract lassos, **new** (from Sakthi, 2026-10-05) | Bec's Desktop: `B6_Xenium_Opt/GI-1..3_coordinates.csv`, `GII-1..3_coordinates.csv` | No |
| Hippocampus lassos | Bec's Desktop: `hippocampus_<sample>.csv` | No |
| Saved R objects (`*.rds` checkpoints) | Written next to the figures in the output folder | No (`*.rds` is ignored by git) |
| UMAP images | `outputs/<region folder>/umap_*.png` | Yes |

Samples: `gi_1_b`, `gi_2_d`, `gi_3_f` (TBI) and `gii_1_a`, `gii_2_c`, `gii_3_e` (sham). The sample list is in `R/run_full_pipeline.R` (the `sample_base` table).

### Where the code lives

| File | What it does for the UMAPs |
|---|---|
| [`R/run_full_pipeline.R`](../R/run_full_pipeline.R) | The whole pipeline, as one function `run_region_pipeline()`, run once per region. **Stage 4 (lines 198-212)** computes the UMAP and saves both images. |
| [`R/plot_umap_celltype.R`](../R/plot_umap_celltype.R) | Two drawing functions: `plot_umap_celltype()` (Figure 1) and `plot_umap_microglia_check()` (Figure 2). They only draw; they do not compute anything. |
| [`R/process_sample_template.R`](../R/process_sample_template.R) | Loads, cleans, normalises, trims to the lasso and labels each mouse (the inputs to the UMAP). |
| [`R/classify_atlas_index.R`](../R/classify_atlas_index.R) | Gives every cell its label (`atlas_identity` and `cell_type`). |
| [`R/apply_confidence_filter.R`](../R/apply_confidence_filter.R) | Optionally relabels low-confidence cells as "Unclassified" (changes `cell_type` only). |

### Step by step: from raw files to the picture

Steps 1 to 6 happen for every analysis; steps 7 to 11 are only for the UMAP.

| # | Step | What it does, in plain words | Code | Settings |
|---|---|---|---|---|
| 1 | Load | Reads one mouse's Xenium folder: a table of gene counts per cell, plus cell positions. | `process_sample_template.R:7` (`LoadXenium`) | |
| 2 | Clean | Drops cells with very little RNA, which are too sparse to trust. | `process_sample_template.R:8` | Keep cells with at least 10 RNA molecules and at least 5 different genes. |
| 3 | Normalise | Makes cells comparable. A big cell catches more RNA than a small one, so each cell's counts are divided by its total, multiplied by 10,000, then log-transformed (which shrinks very large numbers so a few very active genes don't dominate). | `process_sample_template.R:9` (`NormalizeData`) | Seurat default ("LogNormalize", scale 10,000). |
| 4 | Keep the region | Keeps only cells whose centre falls inside the lasso outline. | `process_sample_template.R:11-16` | One lasso file per mouse per region. |
| 5 | Label | Gives each cell a cell-type label from marker genes. | `classify_atlas_index.R` | See the classification protocol. |
| 6 | Optional filter | Relabels cells whose winning label only narrowly beat the runner-up as "Unclassified". | `run_full_pipeline.R:157-170` | Optic tract: cut-off = 15th percentile of all cells' margins. Which cell types it applies to depends on the run (see issue 7). |
| 7 | Combine | Joins all six mice into one object so they share one UMAP layout. | `run_full_pipeline.R:181-185` (`merge`, `JoinLayers`) | |
| 8 | Pick informative genes | Chooses the genes that vary most between cells, so the layout is driven by real differences rather than flat background. | `run_full_pipeline.R:199` (`FindVariableFeatures`) | Seurat default: 2,000 of the ~5,000 genes. |
| 9 | Put genes on the same scale | Rescales each chosen gene so its average is 0 and its spread is 1, so highly expressed genes don't outweigh quieter ones. | `run_full_pipeline.R:200` (`ScaleData`) | Seurat defaults. |
| 10 | Summarise (PCA) | Boils 2,000 genes down to 30 summary "directions" that capture most of the variation. This removes noise and speeds up the next step. | `run_full_pipeline.R:201` (`RunPCA`) | 30 components. |
| 11 | UMAP | Turns the 30 summary numbers per cell into 2 numbers, placing similar cells close together. | `run_full_pipeline.R:202` (`RunUMAP`) | Uses all 30 components. Other settings are Seurat defaults (30 neighbours, `min.dist` 0.3, random seed 42). |
| 12 | Draw and save | Draws the dots, coloured by label and split into sham / TBI panels, and saves the PNGs. | `plot_umap_celltype.R`, `run_full_pipeline.R:204-212` | 10 x 6 inches, 200 dpi, white background. Files: `umap_celltype_gi_vs_gii.png`, `umap_microglia_check_gi_vs_gii.png`, in the run's output folder. |

**Important:** the UMAP layout (steps 7 to 11) is computed from gene activity only. The labels and the confidence filter only change the **colours**, not where the dots are. So two runs on the same cells with different filters give the same layout with different colours.

### Where new figures will appear

`run_full_pipeline.R` now writes optic tract output to new folders, so it will **not** overwrite the figures described above:

- `outputs/optic_tract_newlasso/` (primary run, new lassos)
- `outputs/optic_tract_newlasso_sens_global/` (sensitivity run, filter on all cell types)
- `outputs/hippocampus/` (if the hippocampus run is switched back on)

When those images are added to the repository, update the register and the `[fig-...]` lines at the top.

---

## 4. Issues to review

Each item says whether it was **checked** (seen in the code or figure), or **inferred** (a reasoned concern that needs testing). Letters in brackets refer to the issues in [`docs/review_2026-10-05.md`](review_2026-10-05.md). Nothing here has been changed in the code.

### About what the current figures can tell us

1. **The figures are provisional.** *Checked.* They come from the old optic tract lassos (made 2026-10-02). The project has switched to Sakthi's new lassos, and the old-lasso results are not to be used. The figures need re-making, and this document's section 2 rewritten, once the new run is done. There are no hippocampus UMAPs in the repository.
2. **Cell types mix in the large island.** *Checked (figure); cause inferred.* Astrocytes, microglia and oligodendrocytes overlap rather than separating. Possible causes to check: few RNA molecules per optic tract cell; RNA from neighbouring cells landing in the wrong cell (spillover, see review issue F, where astrocyte genes showed up in every optic tract cell type); or noisy labels. *Suggested check:* colour the UMAP by RNA count per cell (`nCount_Xenium`) and by the confidence margin, to see whether the mixed area is the low-RNA or low-confidence cells.
3. **Few marker genes decide the labels, but 2,000 genes decide the layout.** *Checked.* This is by design and is why the plot is a useful check. It also means disagreement between the two is expected to some degree and is not automatically a labelling error.
4. **The TBI-heavy bottom island needs a per-mouse view.** *Inferred.* Each mouse is a separate tissue section, and no correction for section-to-section differences is applied before the UMAP. A large island seen mainly in one group could come from a single mouse or section rather than from the injury. One TBI mouse (`gi_1_b`) already drives much of the signal elsewhere (review L). *Suggested check:* a version of Figure 1 coloured or split by `sample_id`. This is a small new plot, so it needs Bec's go-ahead.
5. **Unequal cell numbers.** *Checked (notes and figure).* The old lassos captured 168 to 692 cells per mouse, with fewer in sham (review U). Dense and sparse areas therefore can't be compared by eye between panels. The new lassos also differ in size (about 49,000 to 227,000 µm²).

### About the methods

6. **No correction for differences between tissue sections.** *Checked (no such step in Stage 4).* This is a choice, not a mistake, for a visual check, but it should be stated. If the per-mouse view (issue 4) shows islands made of single mice, a correction method could be considered for the picture only. *Decision for Bec.*
7. **Which filter the figures use is not consistent across documents.** *Checked.* The top-level `README.md` says the primary optic tract run applies no confidence filter. `run_full_pipeline.R` runs the primary with the filter on BAM cells only, and the sensitivity run with it on all cell types. The current (old-lasso) figures were made with a 15% filter (inferred from about 18-23% Unclassified). Because the filter changes colours only (see the note in section 3), this does not move any dots, but it changes how many grey or pink "Unclassified" dots appear. The register should name the run each figure came from.
8. **The UMAP shows the filtered labels; the gene-change tests use the unfiltered ones.** *Checked* (review A). The UMAP colours come from `cell_type` (after the filter). Gene-change tests select cells by `atlas_identity` (before the filter). So the UMAP matches the composition bar chart but not exactly the cells in the gene-change heatmaps.
9. **Active vs resting astrocytes are not reliable in the optic tract.** *Notes* (review J). Gfap and Serpina3n are not good state markers in white matter. Read the two astrocyte colours as "marker-high" and "marker-low".
10. **The microglia check is a judgement by eye.** *Checked.* There is no number attached to "coherent". If it is going into a paper, consider a simple number alongside, such as the share of active microglia per mouse.

### About reproducibility and presentation

11. **The UMAP coordinates are not saved.** *Checked.* `combined.rds` is saved at line 185, *before* the UMAP is computed (lines 199-202), so the layout is recomputed each time and can't be reloaded to re-draw a figure. *Suggested fix:* save the object again after line 202. This is a code change, so it needs Bec's go-ahead.
12. **Software versions are not recorded.** *Checked* (review O). UMAP layouts can shift between Seurat or `uwot` versions even with the same seed. The seed (42) is Seurat's default, not set explicitly in the code. Saving `sessionInfo()` with each run would fix this.
13. **Colours are hard to tell apart and may change between versions.** *Checked.* Seventeen labels share an automatic rainbow palette, so neighbours look alike (for example Astrocyte (Active) and Unclassified are both red; the microglia states and Neutrophil are all blue). The colours are assigned by alphabetical order, so adding or removing a cell type shifts every colour. *Suggested fix:* a fixed, named colour for each cell type, used in every figure.
14. **Labels on the plots.** *Checked.* Panels say `gi` / `gii` rather than TBI / Sham (the groups were unblinded on 30 September). Both figures carry the title "cell_type", which is a leftover on Figure 2.
15. **Small note in the code comment.** *Checked.* The usage example at the bottom of `R/plot_umap_celltype.R` saves to `outputs/` rather than a region folder. It is only a comment, so it changes nothing, but it could mislead someone running it by hand.

---

## 5. How to update this document

When the UMAP figures are re-made:

1. **Add the new images** to the repository (for example `outputs/optic_tract_newlasso/umap_celltype_gi_vs_gii.png`).
2. **Edit the `[fig-...]` lines** at the very top so they point to the new files. The images in section 2 update automatically.
3. **Edit the "Document status" table and the figure register:** region, lassos, run, date made, commit, status.
4. **Rewrite section 2** ("What it shows" and "What it means") from the new images. Keep the version note.
5. **Check section 3** still matches the code. If Stage 4 of `run_full_pipeline.R` or `plot_umap_celltype.R` changed, update the step table (line numbers and settings).
6. **Review section 4.** Strike through issues that are fixed (`~~like this~~`) with a note of the date and commit, rather than deleting them, so the history stays visible. Add new ones at the end of the right group.
7. **Add a row to the change log** below.

To add a second region (for example hippocampus UMAPs): add rows to the register, add `[fig-hippo-celltype]`-style lines at the top, and copy the Figure 1 / Figure 2 blocks in section 2 under a new heading for that region.

---

## 6. Change log

| Date | Change | Figures | How to undo |
|---|---|---|---|
| 2026-10-06 | First version. Describes the optic tract old-lasso UMAPs and the code as of commit `66bbd8d`. | fig-celltype, fig-microglia (2026-10-02) | Delete `docs/umap_README.md`. |
