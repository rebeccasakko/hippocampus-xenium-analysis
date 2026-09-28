apply_confidence_filter <- function(xenium, min_margin = 1) {
    xenium$cell_type <- ifelse(
        xenium$classification_margin < min_margin,
        "Unclassified",
        xenium$cell_type
    )
    xenium
}

# apply_confidence_filter_pooled(): derives min_margin as a percentile of
# the POOLED classification_margin distribution across all samples in
# sample_list, then applies apply_confidence_filter() to each sample with
# that shared threshold, rather than reusing a fixed raw number.
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
# uniform across regions even though the resulting raw threshold differs
# by region. See protocols/classification_protocol.md.
apply_confidence_filter_pooled <- function(sample_list, percentile = 0.15) {
    all_margins <- unlist(lapply(sample_list, function(obj) obj$classification_margin))
    threshold <- as.numeric(quantile(all_margins, percentile, na.rm = TRUE))
    message(sprintf(
        "apply_confidence_filter_pooled: %.0fth percentile of pooled classification_margin = %.4f (n = %d cells across %d samples)",
        percentile * 100, threshold, length(all_margins), length(sample_list)))
    lapply(sample_list, apply_confidence_filter, min_margin = threshold)
}
