process_hippocampus_sample <-
function(data_dir, lasso_csv, sample_id, group_id,
                                       region = "Hippocampus",
                                       min_count = 10, min_feature = 5,
                                       marker_list = atlas_index_markers) {
    
    xenium <- LoadXenium(data_dir, fov = "fov")
    xenium <- subset(xenium, subset = nCount_Xenium >= min_count & nFeature_Xenium >= min_feature)
    xenium <- NormalizeData(xenium)
    
    poly <- read.csv(lasso_csv, skip = 2)
    # The membership test below treats the whole file as one polygon, while the
    # area in the cell-count table is summed per Selection. A file holding more
    # than one selection would be joined into one wrong outline here, so say so
    # (audit finding B2; a stop is the author's decision once the counts are known).
    if ("Selection" %in% names(poly) && length(unique(poly$Selection)) > 1)
        warning(sprintf("%s: lasso file %s holds %d selections; they are joined into one polygon for the membership test",
                        sample_id, basename(lasso_csv), length(unique(poly$Selection))), call. = FALSE)
    all_coords <- GetTissueCoordinates(xenium, image = "fov")
    inside <- point.in.polygon(all_coords$x, all_coords$y, poly$X, poly$Y)
    cells_inside <- all_coords$cell[inside > 0]
    
    obj <- subset(xenium, cells = cells_inside)
    obj$sample_id <- sample_id
    obj$group_id  <- group_id
    obj$region    <- region
    rm(xenium); gc()
    
    classify_atlas_index(obj, marker_list = marker_list)
}
