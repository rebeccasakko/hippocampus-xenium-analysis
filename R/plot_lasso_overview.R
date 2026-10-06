# Lasso overview: each ROI outline (red) drawn over its sample's cells (grey).
# Visual check that each lasso CSV is matched to the right sample.
# Reads only cells.csv.gz (cell centroids), so it is fast and light.
plot_lasso_overview <- function(sample_base, lassos, region_label, out_file) {
  cells <- do.call(rbind, lapply(seq_len(nrow(sample_base)), function(k) {
    d <- read.csv(file.path(sample_base$data_dir[k], "cells.csv.gz"))
    data.frame(sample_id = sample_base$sample_id[k],
               x = d$x_centroid, y = d$y_centroid)
  }))
  polys <- do.call(rbind, lapply(names(lassos), function(s) {
    p <- read.csv(lassos[[s]], skip = 2)   # Xenium Explorer export: 2 comment lines
    data.frame(sample_id = s, grp = paste(s, p$Selection), x = p$X, y = p$Y)
  }))
  g <- ggplot2::ggplot(cells, ggplot2::aes(x, y)) +
    ggplot2::geom_point(size = 0.05, colour = "grey70") +
    ggplot2::geom_polygon(data = polys, ggplot2::aes(group = grp),
                          fill = NA, colour = "red", linewidth = 0.6) +
    ggplot2::facet_wrap(~ sample_id, ncol = 3) +
    ggplot2::scale_y_reverse() +          # match Xenium Explorer orientation
    ggplot2::coord_equal() +
    ggplot2::labs(title = paste(region_label, "lassos over each sample"),
                  x = "x (um)", y = "y (um)") +
    ggplot2::theme_classic()
  ggplot2::ggsave(out_file, g, width = 12, height = 8, dpi = 200)
  print(g)   # also show in the Plots pane
}