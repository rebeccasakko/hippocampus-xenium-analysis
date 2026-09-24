plot_spatial_celltype <- function(xenium, title = NULL, pt_size = 0.3) {
    if (is.null(title)) {
        title <- paste0(unique(xenium$region), " \u2014 ", unique(xenium$sample_id),
                         " \u2014 v3 atlas-informed cell types")
    }
    ImageDimPlot(xenium, fov = "fov", group.by = "cell_type",
                 size = pt_size, dark.background = TRUE, axes = FALSE) +
        ggtitle(title)
}
