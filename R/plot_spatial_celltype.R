plot_spatial_celltype <- function(xenium, title = NULL, pt_size = 1.2) {
    if (is.null(title)) {
        title <- paste0(unique(xenium$region), " \u2014 ", unique(xenium$sample_id),
                         " \u2014 v3 atlas-informed cell types")
    }

    coords <- GetTissueCoordinates(xenium, image = "fov")
    coords$cell_type <- xenium$cell_type[match(coords$cell, colnames(xenium))]

    n_na <- sum(is.na(coords$cell_type))
    if (n_na > 0) {
        warning(sprintf(
            "%d of %d cells have no matched cell_type (cell-ID mismatch between GetTissueCoordinates() and colnames(xenium)) -- these will not be plotted",
            n_na, nrow(coords)))
    }
    message(sprintf("plot_spatial_celltype: %d cells with coordinates, %d matched to a cell_type",
                     nrow(coords), nrow(coords) - n_na))

    ggplot(coords, aes(x = x, y = y, color = cell_type)) +
        geom_point(size = pt_size, stroke = 0, alpha = 0.95) +
        coord_fixed() +
        scale_color_hue(l = 75, c = 130) +
        theme_void() +
        theme(
            plot.background  = element_rect(fill = "black", color = NA),
            panel.background = element_rect(fill = "black", color = NA),
            legend.background = element_rect(fill = "black", color = NA),
            legend.key        = element_rect(fill = "black", color = NA),
            legend.text  = element_text(color = "white"),
            legend.title = element_text(color = "white"),
            plot.title   = element_text(color = "white", face = "bold", size = 14)
        ) +
        guides(color = guide_legend(title = "cell_type", override.aes = list(size = 4))) +
        ggtitle(title)
}
