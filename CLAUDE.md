# CLAUDE.md

Guidance for anyone (person or assistant) working in this repository. Written during the code
audit of October 2026 (issue #6, pull request #7); the review documents under `docs/review/` hold
the detail.

## What this is

An R pipeline that turns 10x Xenium in situ transcriptomics from six mouse brains (three TBI,
three sham; same six sections for every region) into, per brain region, a cell-type composition
comparison and a differential-expression comparison between groups. Every cell inside a
hand-drawn lasso polygon is classified individually by scoring it against small marker-gene
panels and taking the best panel; no clustering. Two regions so far: hippocampus (12 panels) and
optic tract (16). Each region is run three times (unfiltered, and two sensitivity runs that set
low-confidence calls to Unclassified). DE is pseudobulk DESeq2 at the mouse level (primary) and
cell-level Wilcoxon (secondary). With n = 3 per group everything is exploratory. The science is
in `docs/review/science-review.md`; how to run and maintain it is in
`docs/review/workflow-review.md`.

## Running it

- Open the repository folder in RStudio (`hippocampus-xenium-analysis.Rproj` is the root anchor;
  `R/run_full_pipeline.R` stops with a message if the working directory is not the root).
- Copy `config.local.example.R` to `config.local.R` (gitignored) and point `xenium_dirs`,
  `hippocampus_lassos` and `optic_tract_lassos` at your copies of the inputs. The defaults expect
  them under `data/` (gitignored). The raw data is unpublished lab data and is not in the
  repository.
- `source("R/run_full_pipeline.R")` loads packages, functions, panels and settings; it runs
  nothing heavy. A missing package stops it with the install line; a missing input stops the
  first run with an OK/MISSING table before any slow work.
- Run `R/run_analyses.R` one block at a time. The primary run of a region must come first; the
  two sensitivity runs reuse its loaded objects. Each block is heavy (hours; see the workflow
  review §8 once the author's timings are in).
- Every run folder gets `run_info.txt` (settings, threshold, commit, cell counts, stage
  timings, package versions) and `session_info.txt`.

## Conventions that matter

- **Machine-specific paths never go in tracked files.** Defaults in `R/run_full_pipeline.R`
  are repo-relative; overrides go in `config.local.R`, documented key by key in
  `config.local.example.R`. Three people must be able to run the same tracked code: the original
  author (whose paths are the commented examples), a collaborator, and a fresh clone. If you add
  an input, add it to `fn.data_manifest()` in `R/_setup.R`, to `config.local.example.R`, and to
  README "Getting the data".
- **One file per stage under `R/`.** `R/_setup.R` is base R only and is sourced before anything
  else. `R/run_full_pipeline.R` is setup and the `run_region_pipeline()` function;
  `R/run_analyses.R` is the runs behind the reported results.
- **Nothing interactive-only runs unguarded.** Figures are saved to files; `print()` of a figure
  sits behind `if (interactive())`.
- **Nothing scientific changes silently.** Marker panels, candidate gene lists, thresholds and
  tests are the author's decisions; questions about them go to issue #8 and the findings log §8,
  and the arithmetic stays as it is until the author decides.
- **Outputs are tracked deliverables, regenerated in place.** `outputs/<run>/` holds PNG figures
  and CSV tables; `.rds` checkpoints and their `.meta.txt` sidecars stay out of git. A checkpoint
  is reused only if its sidecar matches the current panels, lasso file and region.
- **Commits are split by type** (code / docs / outputs / housekeeping / record), subject ends
  with `(issue #N)`, body says the problem, the change, the finding IDs and the evidence.

## Architecture notes that are not obvious from one file

- Classification happens per mouse on the cells inside the lasso that passed QC, so the 0.5 %
  detection floor and the within-panel z-scores depend on the region and the mouse. Which panels
  were scored in which mouse is in `panel_status_per_sample.csv`.
- The confidence filter rewrites `cell_type` only; DE selects cells by `atlas_identity`. So the
  sensitivity runs change composition and not DE; their DE heatmaps are identical to the primary
  run's. `skip_embedding_and_de = TRUE` skips that recomputation (default FALSE).
- The near-tie category list for `sens_neartie` is computed by `near_tie_categories()` from the
  primary run and compared with the recorded list in `R/run_analyses.R`; a difference warns and
  the recorded list is used.
- The Astrocyte and Microglia identity panels are also their Homeostatic state panels, and
  several DE candidate genes are identity genes for the same cell type (issue #8, M1 and M10).
- `AggregateExpression()` renames `_` to `-` in sample names (Seurat v5); group membership is
  taken from `group_id` metadata, with the name-based rule only as a fallback.
- No `set.seed()` anywhere: `AddModuleScore()` (seed 1), `RunPCA()` and `RunUMAP()` (seed 42)
  use Seurat defaults. UMAP coordinates can still move between package versions.

## Inputs, outputs, and what is gitignored

- Inputs (not tracked): six Xenium bundles (one folder per mouse, several GB each) and six
  lasso CSVs per region (Xenium Explorer exports: two comment lines, then
  `Selection, X, Y, Class, Color`). Place them under `data/` or point `config.local.R` at them.
- Tracked outputs: `outputs/<region>_<run>/` with `cell_counts_per_sample.csv`,
  `composition_per_sample.csv`, `per_sample_depth.csv`, `panel_status_per_sample.csv`,
  `composition_barchart.png`, two UMAP PNGs, per cell type `pseudobulk_heatmap_<ct>.png` and
  `wilcox_heatmap_<ct>.png`, `de/` tables and diagnostics, `run_info.txt`, `session_info.txt`,
  and `lasso_overview.png` in primary runs; at the top of `outputs/`, the near-tie sensitivity
  tables, the hippocampus stats table and the driver-gene diagnostics. Folders ending `_withCP`
  are the earlier 13-panel hippocampus runs (issue #9).
- Gitignored: `data/`, `config.local.R`, `*.rds`, `*.rds.meta.txt`, `hoard/`, `Rplots.pdf`,
  `.Rproj.user/`.
- Review documents: `docs/review/science-review.md`, `docs/review/workflow-review.md`,
  `docs/review/issue6-audit-repo-log.md` (the audit's state and decisions).
