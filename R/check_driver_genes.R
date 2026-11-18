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
                                  context = c("Hippocampal pyramidal", "Dentate gyrus granule"),
                                  scalebar_um = 500, red_size = 0.6) {
  xy <- do.call(rbind, lapply(names(objs), function(s) {
    o <- objs[[s]]
    d <- Seurat::GetTissueCoordinates(o)
    d$cell_type <- o$cell_type[match(d$cell, colnames(o))]
    d$sample <- s
    d
  }))
  nm <- c(gi_1_b = "TBI 1", gi_2_d = "TBI 2", gi_3_f = "TBI 3",
          gii_1_a = "Sham 1", gii_2_c = "Sham 2", gii_3_e = "Sham 3")
  xy$sample <- factor(xy$sample, levels = names(nm))
  
  # Same window size in every panel, centred on each section,
  # so 1 µm is the same length in all panels
  rng <- do.call(rbind, lapply(split(xy, xy$sample, drop = TRUE), function(d)
    data.frame(sample = d$sample[1],
               cx = mean(range(d$x)), cy = mean(range(d$y)),
               wx = diff(range(d$x)), wy = diff(range(d$y)))))
  half_x <- max(rng$wx) / 2 * 1.05
  half_y <- max(rng$wy) / 2 * 1.05
  band   <- 0.15 * 2 * half_y   # empty strip below the tissue for the scale bar
  
  corners <- data.frame(sample = rep(rng$sample, each = 2),
                        x = c(rbind(rng$cx - half_x, rng$cx + half_x)),
                        y = c(rbind(rng$cy - half_y, rng$cy + half_y + band)))
  
  # Scale bar in the first panel only, in the strip below the tissue
  # (y is reversed, so larger y is lower on the page)
  x0 <- rng$cx[1] - half_x * 0.95
  sb <- data.frame(sample = rng$sample[1],
                   x = x0, xend = x0 + scalebar_um,
                   y = rng$cy[1] + half_y + band * 0.7,
                   label_y = rng$cy[1] + half_y + band * 0.3,
                   label = paste(scalebar_um, "µm"))
  
  g <- ggplot2::ggplot(xy, ggplot2::aes(x, y)) +
    ggplot2::geom_blank(data = corners) +
    ggplot2::geom_point(data = subset(xy, !cell_type %in% c(context, category)),
                        colour = "grey85", size = 0.2) +
    ggplot2::geom_point(data = subset(xy, cell_type %in% context),
                        colour = "steelblue", size = 0.3) +
    ggplot2::geom_point(data = subset(xy, cell_type == category),
                        colour = "red", size = red_size) +
    ggplot2::geom_segment(data = sb,
                          ggplot2::aes(x = x, xend = xend, y = y, yend = y),
                          linewidth = 0.8) +
    ggplot2::geom_text(data = sb,
                       ggplot2::aes(x = (x + xend) / 2, y = label_y, label = label),
                       size = 3.5) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::facet_wrap(~ sample, nrow = 2, scales = "free",
                        labeller = ggplot2::as_labeller(nm)) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_reverse(expand = c(0, 0)) +
    ggplot2::theme_void() +
    ggplot2::theme(aspect.ratio = (2 * half_y + band) / (2 * half_x),
                   strip.text = ggplot2::element_text(size = 12,
                                                      margin = ggplot2::margin(b = 4)))
  ggplot2::ggsave(out_file, g, width = 12, height = 7, dpi = 300, bg = "white")
  print(g)
  invisible(g)
}
