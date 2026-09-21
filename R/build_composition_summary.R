build_composition_summary <-
function(sample_list, group_levels = c("gii", "gi")) {
    
    composition_df <- do.call(rbind, lapply(sample_list, function(obj) {
        tab <- table(obj$cell_type)
        data.frame(
            sample_id = unique(obj$sample_id),
            group_id  = unique(obj$group_id),
            cell_type = names(tab),
            n_cells   = as.integer(tab)
        )
    })) %>%
        group_by(sample_id) %>%
        mutate(total_cells = sum(n_cells), pct = 100 * n_cells / total_cells) %>%
        ungroup() %>%
        mutate(replicate_num = as.integer(sub("^gi+_([0-9]+)_.*$", "\\1", sample_id)))
    
    summary_df <- composition_df %>%
        group_by(group_id, cell_type) %>%
        summarise(mean_pct = mean(pct), sem_pct = sd(pct) / sqrt(n()), n_samples = n(), .groups = "drop") %>%
        mutate(group_id = factor(group_id, levels = group_levels))
    
    cell_type_order <- summary_df %>%
        group_by(cell_type) %>%
        summarise(total_pct = sum(mean_pct), .groups = "drop") %>%
        arrange(desc(total_pct)) %>%
        pull(cell_type)
    
    summary_df <- summary_df %>% mutate(cell_type = factor(cell_type, levels = cell_type_order))
    
    list(composition_df = composition_df, summary_df = summary_df)
}
