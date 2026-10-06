# ============================================================
# Single-gene driver diagnostics for a classified category.
# A real cell type should express its panel genes together; if one
# gene is high while the others are ~0, the call is likely driven by
# that gene alone (marker cross-reactivity), even when margins are high.
# ============================================================

# % of cells called `category` (and of all other cells) with each gene detected
driver_check <- function(objs, category, genes) {
  do.call(rbind, lapply(names(objs), function(s) {
    o <- objs[[s]]; hit <- o$cell_type == category
    cnt <- Seurat::GetAssayData(o, layer = "counts")[genes, , drop = FALSE] > 0
    data.frame(sample = s, n_called = sum(hit),
               t(round(100 * rowMeans(cnt[, hit, drop = FALSE]), 1)),
               other = t(round(100 * rowMeans(cnt[, !hit, drop = FALSE]), 1)))
  }))
}

# What the `category` cells would otherwise have been called (pooled over samples)
runner_up_summary <- function(objs, category) {
  rb <- do.call(rbind, lapply(objs, function(o)
    o@meta.data[o$cell_type == category, c("runner_up", "classification_margin")]))
  message(category, ": median classification_margin = ",
          round(median(rb$classification_margin), 2), " (n = ", nrow(rb), ")")
  out <- as.data.frame(sort(round(100 * prop.table(table(rb$runner_up)), 1),
                            decreasing = TRUE))
  names(out) <- c("runner_up", "pct")
  out
}

# Where the `category` cells sit (red), against context cell types (blue)
plot_category_spatial <- function(objs, category, out_file,
                                  context = c("Hippocampal pyramidal", "Dentate gyrus granule")) {
  xy <- do.call(rbind, lapply(names(objs), function(s) {
    o <- objs[[s]]
    d <- Seurat::GetTissueCoordinates(o)
    d$cell_type <- o$cell_type[match(d$cell, colnames(o))]
    d$sample <- s
    d
  }))
  g <- ggplot2::ggplot(xy, ggplot2::aes(x, y)) +
    ggplot2::geom_point(data = subset(xy, !cell_type %in% c(context, category)),
                        colour = "grey85", size = 0.2) +
    ggplot2::geom_point(data = subset(xy, cell_type %in% context),
                        colour = "steelblue", size = 0.3) +
    ggplot2::geom_point(data = subset(xy, cell_type == category),
                        colour = "red", size = 0.6) +
    ggplot2::facet_wrap(~ sample, scales = "free") +
    ggplot2::scale_y_reverse() + ggplot2::theme_classic() +
    ggplot2::labs(title = paste0("'", category, "' calls (red) vs ",
                                 paste(context, collapse = "/"), " (blue)"))
  ggplot2::ggsave(out_file, g, width = 12, height = 7, dpi = 200)
  print(g)
  invisible(g)
}