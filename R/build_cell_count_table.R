build_cell_count_table <- function(sample_list, region = NA_character_, run = "primary",
                                   lasso_area_um2 = NULL, out_file = NULL) {
  tab <- data.frame(
    sample = names(sample_list),
    n_cells = vapply(sample_list, ncol, integer(1)),
    n_unclassified = vapply(sample_list,
                            function(o) sum(o$cell_type == "Unclassified", na.rm = TRUE),
                            integer(1)),
    row.names = NULL)
  tab$pct_unclassified <- round(100 * tab$n_unclassified / tab$n_cells, 1)
  tab$group  <- ifelse(grepl("gii", tab$sample), "gii", "gi")
  tab$region <- region
  tab$run    <- run
  
  if (!is.null(lasso_area_um2)) {
    tab$lasso_area_um2    <- as.numeric(lasso_area_um2[tab$sample])
    tab$cells_per_1000um2 <- round(1000 * tab$n_cells / tab$lasso_area_um2, 2)
  }
  
  if (!is.null(out_file)) write.csv(tab, out_file, row.names = FALSE)
}
