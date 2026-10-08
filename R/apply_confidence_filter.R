apply_confidence_filter <- function(xenium, min_margin = 1, categories = NULL) {
  low <- !is.na(xenium$classification_margin) & xenium$classification_margin < min_margin
  if (!is.null(categories)) low <- low & xenium$cell_type %in% categories
  xenium$cell_type <- ifelse(low, "Unclassified", xenium$cell_type)
  xenium
}

# apply_confidence_filter_pooled(): derives min_margin as a percentile of
# the POOLED classification_margin distribution across all samples in
# sample_list, then applies apply_confidence_filter() to each sample with
# that shared threshold, rather than reusing a raw fixed number.
#
# Why: classification_margin is mechanically panel-count-dependent -- more
# competing candidate panels shrinks typical winning margins regardless of
# classification quality, so a raw threshold calibrated in one region's
# panel-count context (e.g. hippocampus's margin > 1, validated for its
# BAM panel) does not transfer to a region with a different panel count.
# Confirmed 2026-09-28 on Optic Tract (16 candidate panels after excluding
# Hippocampal pyramidal/Dentate gyrus granule/RGC/Ependymal-ciliated):
# margin > 1 would have wiped ~37% of Astrocyte/Microglia calls, well
# above their 1st-quartile margins (0.67-0.83) -- nowhere near a genuine
# confidence floor for this panel-count context. A percentile of the
# region's own observed margin distribution keeps the DERIVATION rule
# uniform across regions even though the resulting threshold differs
# by region. See protocols/classification_protocol.md.
#
# categories: if given, only cells currently labelled with one of these
# types are eligible to be reclassified as Unclassified (e.g. "BAM").
# If NULL, the filter applies to all cells (used for the global
# sensitivity run).
apply_confidence_filter_pooled <- function(sample_list, percentile = 0.15, categories = NULL) {
  all_margins <- unlist(lapply(sample_list, function(obj) obj$classification_margin))
  threshold <- as.numeric(quantile(all_margins, percentile, na.rm = TRUE))
  out <- lapply(sample_list, apply_confidence_filter,
                min_margin = threshold, categories = categories)
  n_changed <- sum(unlist(Map(function(a, b) sum(a$cell_type != b$cell_type), sample_list, out)))
  message(sprintf(
    "apply_confidence_filter_pooled: %.0fth percentile of pooled margins = %.4f (n = %d cells, %d samples); %d cells reclassified; scope: %s",
    percentile * 100, threshold, length(all_margins), length(sample_list), n_changed,
    if (is.null(categories)) "all categories" else paste(categories, collapse = ", ")))
  attr(out, "threshold") <- threshold   # recorded in run_info.txt by run_region_pipeline() (audit R5)
  out
}

# near_tie_categories(): the categories eligible for the sens_neartie run,
# computed from the PRIMARY run's objects by the rule stated in README and
# R/run_analyses.R: categories (cell_type) in which more than `min_share` of
# cells fall below the pooled `percentile` margin threshold. "Unclassified"
# is never returned (reassigning it is a no-op). Before the audit (issue #6,
# finding R6) the list was typed into R/run_analyses.R by hand from the
# printed margins; R/run_analyses.R now computes it, prints the shares, and
# warns if it differs from the recorded list (which stays in use, decision 10).
near_tie_categories <- function(sample_list, percentile = 0.15, min_share = 0.15) {
  m  <- unlist(lapply(sample_list, function(o) o$classification_margin))
  ct <- unlist(lapply(sample_list, function(o) as.character(o$cell_type)))
  thr <- as.numeric(quantile(m, percentile, na.rm = TRUE))
  share <- sort(tapply(m < thr, ct, mean, na.rm = TRUE), decreasing = TRUE)
  message(sprintf("near_tie_categories: threshold %.4f (%.0fth percentile); share of cells below it: %s",
                  thr, percentile * 100,
                  paste(sprintf("%s %.1f%%", names(share), 100 * share), collapse = ", ")))
  setdiff(names(share)[share > min_share], "Unclassified")
}
