build_composition_summary <-
function(sample_list, group_levels = c("gii", "gi"), cell_type_order = NULL) {

    # Every cell_type seen in ANY sample must get an explicit row (even if
    # n_cells = 0) in EVERY sample. table(obj$cell_type) on its own only
    # returns levels actually present in that one sample -- a cell_type
    # that's absent (whether because its panel was excluded from scoring
    # entirely at the detection floor, or it just happened to get zero
    # cells) silently gets NO row at all, rather than a 0% row. That
    # silently drops that sample out of the category's mean/SEM below,
    # shrinking n_samples for that one category without changing the
    # plot's "n = 3 per group" label, and biasing mean_pct upward by
    # averaging only over the samples where the category happened to
    # appear. Confirmed 2026-09-28 via Optic Tract: BAM was silently
    # averaged over n=2 gi samples instead of 3 (10.5% shown vs the true
    # 7.0%), and B cell -- present in all 3 gi samples, absent in all 3
    # gii samples -- had NO gii row at all, rendering as a bar with only
    # one group's percentage instead of two.
    all_cell_types <- unique(unlist(lapply(sample_list, function(obj) unique(obj$cell_type))))

    composition_df <- do.call(rbind, lapply(sample_list, function(obj) {
        tab <- table(factor(obj$cell_type, levels = all_cell_types))
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
    
    # cell_type_order: pass this in to reuse another chart's category
    # ordering (e.g. levels(other_summary$summary_df$cell_type)) instead of
    # this call's own abundance ranking. Needed whenever two charts built
    # from different sample sets (a QC-sensitivity comparison, a re-run
    # after a panel fix, etc.) are meant to sit side by side -- each call
    # ranks categories by ITS OWN total_pct by default, so two charts of
    # mostly-the-same categories can still come out in a slightly different
    # order and be hard to visually compare bar-by-bar. Confirmed 2026-09-29:
    # standard-QC and strict-QC Hippocampus charts, otherwise near-identical,
    # had Unclassified/BAM/Endothelial/Microglia (Homeostatic) shuffled
    # relative to each other purely from each chart computing its own order.
    if (is.null(cell_type_order)) {
        cell_type_order <- summary_df %>%
            group_by(cell_type) %>%
            summarise(total_pct = sum(mean_pct), .groups = "drop") %>%
            arrange(desc(total_pct)) %>%
            pull(cell_type)
    }
    
    summary_df <- summary_df %>% mutate(cell_type = factor(cell_type, levels = cell_type_order))
    
    list(composition_df = composition_df, summary_df = summary_df)
}
