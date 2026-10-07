plot_composition_barchart <-
function(summary_df, group_colors,
                                      title = "Hippocampus cell-type composition: TBI vs sham",
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
      # Match the dots' group order to the bars' order, so each dot lands in its own bar
      composition_df$group_id <- factor(composition_df$group_id,
                                        levels = if (is.factor(summary_df$group_id)) levels(summary_df$group_id) else c("gii", "gi"))
      # fill = group_id tells jitterdodge which bar each point belongs to.
      # Shape identifies the individual mouse: filled = TBI, hollow = sham.
      p <- p + geom_point(data = composition_df,
                          aes(x = cell_type, y = pct, fill = group_id, shape = sample_id),
                          position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.85,
                                                          seed = 1),
                          size = 1.8, stroke = 0.35, color = "black", alpha = 0.85,
                          inherit.aes = FALSE) +
        scale_shape_manual(values = c(16, 15, 17, 1, 0, 2), name = "Mouse",
                           labels = function(x) paste(ifelse(grepl("^gii", x), "Sham", "TBI"),
                                                      sub("^gii?_(\\d+)_.*$", "\\1", x)))
    }
    
    # % value labels: rotated 90 degrees and right-justified against the y = -1
    # baseline so each label only needs the height of a text string, not its
    # full printed width, to clear its neighbor. With up to ~20 cell-type
    # categories x 2 groups = 40 dodged bars on one chart (Optic Tract), the
    # unrotated version (size = 2.1, vjust = 1) collided into an illegible
    # smear at the bottom of the chart -- flagged 2026-09-28. Rotating fixes
    # this regardless of how many categories a given region's panel has, so
    # this is a rendering-only fix, not a per-region special case.
    p +
        geom_text(aes(label = sprintf("%.1f%%", mean_pct), y = -1),
                  position = position_dodge(width = 0.85),
                  size = 2.4, fontface = "bold", color = "black",
                  angle = 90, hjust = 1, vjust = 0.5) +
      scale_fill_manual(values = group_colors, name = "Group",
                        breaks = c("gii", "gi"), labels = c("Sham", "TBI")) +
      guides(fill = guide_legend(override.aes = list(shape = NA))) +
        scale_y_continuous(expand = expansion(mult = c(0.18, 0.05))) +
        coord_cartesian(clip = "off") +
        labs(x = NULL, y = "% of cells (mean ± SEM, n = 3 per group)",
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
            plot.margin = margin(t = 10, r = 10, b = 55, l = 25)
        )
}
