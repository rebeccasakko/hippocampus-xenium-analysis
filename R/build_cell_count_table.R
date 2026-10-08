build_cell_count_table <- function(sample_list, region = NA_character_, run = "primary",
                                   lasso_area_um2 = NULL, out_file = NULL) {
  tab <- data.frame(
    sample = names(sample_list),
    n_cells = vapply(sample_list, function(o) as.integer(ncol(o)), integer(1)),
    n_unclassified = vapply(sample_list,
                            function(o) as.integer(sum(o$cell_type == "Unclassified", na.rm = TRUE)),
                            integer(1)),
    row.names = NULL)
  tab$pct_unclassified <- round(100 * tab$n_unclassified / tab$n_cells, 1)
  # QC check: % of cells (already passing nCount >= 10 & nFeature >= 5) that
  # would also pass the stricter nCount >= 20 & nFeature >= 10 cutoff.
  # If this drops well below ~98%, a full 20/10 rerun is worth doing.
  tab$pct_pass_qc20_10 <- vapply(sample_list, function(o)
    round(100 * mean(o$nCount_Xenium >= 20 & o$nFeature_Xenium >= 10), 1),
    numeric(1))
  # Group from the objects' own group_id metadata rather than from the sample
  # name (audit finding B5); the name-based rule stays as the fallback for
  # objects without that column.
  grp <- vapply(sample_list, function(o)
    if ("group_id" %in% colnames(o@meta.data)) as.character(unique(o$group_id))[1] else NA_character_,
    character(1))
  tab$group  <- ifelse(is.na(grp), ifelse(grepl("gii", tab$sample), "gii", "gi"), unname(grp))
  tab$region <- region
  tab$run    <- run
  
  if (!is.null(lasso_area_um2)) {
    tab$lasso_area_um2    <- as.numeric(lasso_area_um2[tab$sample])
    tab$cells_per_1000um2 <- round(1000 * tab$n_cells / tab$lasso_area_um2, 2)
  }
  
  if (!is.null(out_file)) write.csv(tab, out_file, row.names = FALSE)
  tab
}