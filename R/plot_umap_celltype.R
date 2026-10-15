plot_umap_celltype <- function(xenium, group_levels = c("gii", "gi"), pt_size = 0.1) {
  
  group_labels <- c(gii = "Sham", gi = "TBI")
  xenium$group_id <- factor(xenium$group_id, levels = group_levels,
                            labels = group_labels[group_levels])
  
  DimPlot(xenium, reduction = "umap", group.by = "cell_type",
          split.by = "group_id", pt.size = pt_size) +
    ggtitle("cell_type")
}

# Highlights Microglia (Active) and Microglia (Homeostatic) against every
# other cell type dimmed to grey -- the QC check run after the full
# cell_type UMAP to confirm the microglia state split looks spatially
# sane before trusting it downstream (e.g. in the pseudobulk DE panels).
plot_umap_microglia_check <-
function(xenium, group_levels = c("gii", "gi"), pt_size = 0.3,
          active_color = "black", homeostatic_color = "#8ed4bb", other_color = "grey85") {

    xenium$group_id <- factor(xenium$group_id, levels = group_levels)
    xenium$microglia_highlight <- ifelse(grepl("^Microglia", xenium$cell_type),
                                          xenium$cell_type, "Other")

    DimPlot(xenium, reduction = "umap", group.by = "microglia_highlight",
            split.by = "group_id", pt.size = pt_size,
            cols = c("Microglia (Active)" = active_color,
                     "Microglia (Homeostatic)" = homeostatic_color,
                     "Other" = other_color)) +
        ggtitle("cell_type")
}

# ============================================================
# Usage (combined = merged, classify_atlas_index()'d Seurat object for
# all 6 samples). classify_atlas_index() does not compute a reduction,
# so PCA/UMAP must be run once beforehand if `combined` doesn't have one:
#
#   if (!"umap" %in% names(combined@reductions)) {
#     combined <- FindVariableFeatures(combined)
#     combined <- ScaleData(combined)
#     combined <- RunPCA(combined, npcs = 30)
#     combined <- RunUMAP(combined, dims = 1:30)
#   }
#
#   umap_celltype <- plot_umap_celltype(combined)
#   ggsave("outputs/umap_celltype_gi_vs_gii.png", umap_celltype,
#          width = 10, height = 6, dpi = 200, bg = "white")
#
#   umap_microglia_check <- plot_umap_microglia_check(combined)
#   ggsave("outputs/umap_microglia_check_gi_vs_gii.png", umap_microglia_check,
#          width = 10, height = 6, dpi = 200, bg = "white")
# ============================================================
