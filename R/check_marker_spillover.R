# ============================================================
# Marker-spillover check: is a gene's signal in a target cell
# population real expression, or segmentation contamination
# ("spillover") bled in from a neighboring cell type that
# strongly expresses that marker?
#
# Written up (2026-10-02) from an interactive investigation of Gfap
# signal in the microglia panel ("real signal or spillover?"), so the
# same four-part diagnostic can be re-run on a different gene/cell-type
# pair later rather than redone by hand.
#
# Logic: if marker-positive cells in the target population are
# actually segmentation-blended with the contaminating cell type,
# they should look MORE like that contaminating type on all four
# checks below than marker-negative cells do. If all four come back
# flat between marker+ and marker- cells, that argues against
# spillover and for the signal being real.
#
#   1. Rate at which the contaminating cell type is the runner-up
#      classification call (closer/more ambiguous calls = more
#      likely boundary/blended cells)
#   2. The contaminating cell type's identity (module) score
#      (idx_<ContaminatingIdentity> in the Seurat object's metadata)
#   3. Classification confidence margin (winning score minus
#      runner-up score; smaller = less confident/more ambiguous call)
#   4. Segmentation size (nCount_Xenium; larger/noisier segmented
#      cells are more prone to picking up a neighbor's transcripts)
#
# Validated against the Gfap/Microglia/Astrocyte case: all four
# checks came back flat (p = 0.86 for nCount; none of the other
# three showed a spillover-consistent shift either), supporting
# the Gfap signal in Microglia-classified cells as real rather than
# an artifact of astrocyte segmentation bleed.
# ============================================================

library(dplyr)    # dplyr, tibble and ggplot2 are in the required package list (R/_setup.R);
library(tibble)   # tidyverse itself is not a dependency of this repository
library(ggplot2)
library(Seurat)

# Build the per-cell diagnostic table for one marker-spillover check.
# `seurat_obj` should already be subset to the target cell population
# (e.g. microglia_cells, not the full merged object). `contaminating_identity`
# must match the suffix used in that object's idx_<Identity> metadata
# columns (e.g. "Astrocyte" for idx_Astrocyte) and match a value that
# appears in the runner_up column.
build_spillover_check_data <- function(seurat_obj, gene, contaminating_identity,
                                        group_var = "group_id") {
    idx_col <- paste0("idx_", contaminating_identity)

    stopifnot(
        gene %in% rownames(seurat_obj),
        idx_col %in% colnames(seurat_obj@meta.data),
        "classification_margin" %in% colnames(seurat_obj@meta.data),
        "nCount_Xenium" %in% colnames(seurat_obj@meta.data),
        "runner_up" %in% colnames(seurat_obj@meta.data)
    )

    df <- FetchData(seurat_obj, vars = c(gene, "nCount_Xenium", group_var,
                                          idx_col, "classification_margin")) %>%
        as_tibble() %>%
        dplyr::rename(gene_expr = all_of(gene), idx_contam = all_of(idx_col)) %>%
        dplyr::mutate(
            marker_pos = gene_expr > 0,
            # runner_up is stored as a list-column in meta.data (one named
            # element per cell) -- flatten to plain character before use,
            # same fix needed interactively for the Gfap/Microglia check.
            runner_up = vapply(seurat_obj@meta.data$runner_up, function(x) {
                if (is.null(x) || length(x) == 0) NA_character_ else as.character(x)[1]
            }, character(1)) %>% unname()
        )

    df
}

# Run the four spillover-diagnostic statistical tests on the output of
# build_spillover_check_data() and summarize into one table.
run_spillover_checks <- function(data, contaminating_identity) {
    runner_up_test <- chisq.test(table(data$marker_pos, data$runner_up == contaminating_identity))
    idx_test       <- wilcox.test(idx_contam ~ marker_pos, data = data)
    margin_test    <- wilcox.test(classification_margin ~ marker_pos, data = data)
    ncount_test    <- wilcox.test(nCount_Xenium ~ marker_pos, data = data)

    means <- data %>%
        dplyr::group_by(marker_pos) %>%
        dplyr::summarise(
            pct_runner_up_contam = mean(runner_up == contaminating_identity, na.rm = TRUE) * 100,
            mean_idx_contam      = mean(idx_contam, na.rm = TRUE),
            mean_margin          = mean(classification_margin, na.rm = TRUE),
            mean_nCount          = mean(nCount_Xenium, na.rm = TRUE),
            .groups = "drop"
        )

    tibble(
        Check = c(paste0(contaminating_identity, " as runner-up call (%)"),
                  paste0(contaminating_identity, " identity score"),
                  "Classification confidence (margin)",
                  "Segmentation size (nCount)"),
        `Marker-negative` = c(means$pct_runner_up_contam[means$marker_pos == FALSE],
                               means$mean_idx_contam[means$marker_pos == FALSE],
                               means$mean_margin[means$marker_pos == FALSE],
                               means$mean_nCount[means$marker_pos == FALSE]),
        `Marker-positive` = c(means$pct_runner_up_contam[means$marker_pos == TRUE],
                               means$mean_idx_contam[means$marker_pos == TRUE],
                               means$mean_margin[means$marker_pos == TRUE],
                               means$mean_nCount[means$marker_pos == TRUE]),
        `p-value` = c(signif(runner_up_test$p.value, 3),
                      signif(idx_test$p.value, 3),
                      signif(margin_test$p.value, 3),
                      signif(ncount_test$p.value, 3))
    ) %>%
        dplyr::mutate(
            dplyr::across(c(`Marker-negative`, `Marker-positive`), ~ round(.x, 3)),
            Verdict = ifelse(`p-value` < 0.05,
                              "Differs - possible spillover signature",
                              "Flat - no spillover signature")
        )
}

# Build the 4-panel diagnostic summary plot (theme_classic, matching the
# rest of this repo's figures).
plot_spillover_checks <- function(data, contaminating_identity, gene) {
    pos_label <- paste0(gene, "+")
    neg_label <- paste0(gene, "-")

    long_df <- dplyr::bind_rows(
        data %>% dplyr::transmute(marker_pos, metric = paste(contaminating_identity, "identity score"),
                            value = idx_contam),
        data %>% dplyr::transmute(marker_pos, metric = "Classification margin",
                            value = classification_margin),
        data %>% dplyr::transmute(marker_pos, metric = "Segmentation size (nCount)",
                            value = nCount_Xenium)
    ) %>%
        dplyr::group_by(metric, marker_pos) %>%
        dplyr::summarise(mean_val = mean(value, na.rm = TRUE),
                  se_val = sd(value, na.rm = TRUE) / sqrt(n()),
                  .groups = "drop") %>%
        dplyr::bind_rows(
            data %>%
                dplyr::group_by(marker_pos) %>%
                dplyr::summarise(mean_val = mean(runner_up == contaminating_identity, na.rm = TRUE) * 100,
                          .groups = "drop") %>%
                dplyr::mutate(metric = paste(contaminating_identity, "runner-up (%)"), se_val = NA_real_)
        ) %>%
        dplyr::mutate(group_label = ifelse(marker_pos, pos_label, neg_label))

    ggplot(long_df, aes(x = group_label, y = mean_val, fill = group_label)) +
        geom_col(width = 0.6) +
        geom_errorbar(aes(ymin = mean_val - se_val, ymax = mean_val + se_val), width = 0.15, na.rm = TRUE) +
        facet_wrap(~ metric, scales = "free_y", nrow = 1) +
        scale_fill_manual(values = setNames(c("grey70", "#2C5F2D"), c(neg_label, pos_label))) +
        theme_classic() +
        theme(legend.position = "none",
              axis.title.x = element_blank(),
              strip.text = element_text(face = "bold")) +
        labs(title = paste0(pos_label, " vs ", neg_label, " segmentation-artifact checks"),
             y = NULL)
}

# Run the full marker-spillover check end to end: build the per-cell
# table, run the four tests, and build the summary plot.
check_marker_spillover <- function(seurat_obj, gene, contaminating_identity,
                                    group_var = "group_id") {
    data <- build_spillover_check_data(seurat_obj, gene, contaminating_identity, group_var)

    list(
        data = data,
        summary = run_spillover_checks(data, contaminating_identity),
        plot = plot_spillover_checks(data, contaminating_identity, gene)
    )
}

# ============================================================
# Usage:
#   result <- check_marker_spillover(microglia_cells, gene = "Gfap",
#                                     contaminating_identity = "Astrocyte")
#   result$summary   # table with p-values and verdict per check
#   result$plot      # 4-panel diagnostic figure
#   result$data      # per-cell table, for any further ad hoc digging
# ============================================================
