plot_composition_barchart <-
function(summary_df, group_colors,
                                      title = "Hippocampus cell-type composition: gii vs gi",
                                      subtitle = NULL,
                                      show_points = FALSE,
                                      composition_df = NULL) {
    
    p <- ggplot(summary_df, aes(x = cell_type, y = mean_pct, fill = group_id)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.85), width = 0.75,
                 color = "grey25", linewidth = 0.3) +
        geom_errorbar(aes(ymin = mean_pct - sem_pct, ymax = mean_pct + sem_pct),
                      position = position_dodge(width = 0.85), width = 0.25, linewidth = 0.4)
    
    if (show_points) {
        if (is.null(composition_df)) stop("composition_df is required when show_points = TRUE")
        p <- p + geom_jitter(data = composition_df,
                             aes(x = cell_type, y = pct, shape = factor(replicate_num)),
                             position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.85),
                             size = 1.8, stroke = 0.35, color = "black", alpha = 0.85, inherit.aes = FALSE) +
            scale_shape_manual(values = c("1" = 1, "2" = 0, "3" = 2), name = "Sample")
    }
    
    p +
        geom_text(aes(label = sprintf("%.1f%%", mean_pct), y = -1),
                  position = position_dodge(width = 0.85),
                  size = 2.1, fontface = "bold", color = "black", vjust = 1) +
        scale_fill_manual(values = group_colors, name = "Group") +
        scale_y_continuous(expand = expansion(mult = c(0.12, 0.05))) +
        coord_cartesian(clip = "off") +
        labs(x = NULL, y = "% of cells (mean \u00b1 SEM, n = 3 per group)",
             title = title, subtitle = subtitle) +
        theme_classic(base_size = 12) +
        theme(
            axis.text.x = element_text(angle = 40, hjust = 1, size = 9),
            axis.line = element_line(color = "grey40"),
            panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
            plot.title = element_text(face = "bold", size = 13),
            plot.subtitle = element_text(color = "grey40", size = 9.5),
            legend.title = element_text(face = "bold"),
            legend.key = element_blank(),
            plot.margin = margin(t = 10, r = 10, b = 40, l = 25)
        )
}
