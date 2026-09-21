process_hippocampus_sample <-
function(data_dir, lasso_csv, sample_id, group_id,
                                       region = "Hippocampus",
                                       min_count = 10, min_feature = 5) {
    
    xenium <- LoadXenium(data_dir, fov = "fov")
    xenium <- subset(xenium, subset = nCount_Xenium >= min_count & nFeature_Xenium >= min_feature)
    xenium <- NormalizeData(xenium)
    
    poly <- read.csv(lasso_csv, skip = 2)
    all_coords <- GetTissueCoordinates(xenium, image = "fov")
    inside <- point.in.polygon(all_coords$x, all_coords$y, poly$X, poly$Y)
    cells_inside <- all_coords$cell[inside > 0]
    
    obj <- subset(xenium, cells = cells_inside)
    obj$sample_id <- sample_id
    obj$group_id  <- group_id
    obj$region    <- region
    rm(xenium); gc()
    
    classify_atlas_index(obj)
}
