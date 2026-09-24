apply_confidence_filter <- function(xenium, min_margin = 1) {
    xenium$cell_type <- ifelse(
        xenium$classification_margin < min_margin,
        "Unclassified",
        xenium$cell_type
    )
    xenium
}
