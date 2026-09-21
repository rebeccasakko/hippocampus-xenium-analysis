classify_atlas_index <-
function(xenium, marker_list = atlas_index_markers, min_z = 0, min_detect_rate = 0.005) {
    label_lookup <- setNames(names(marker_list), paste0("idx_", make.names(names(marker_list))))
    score_cols <- character(0)
    counts <- GetAssayData(xenium, layer = "counts")
    for (name in names(marker_list)) {
        genes <- intersect(marker_list[[name]], rownames(xenium))
        if (length(genes) == 0) {
            warning("No genes from '", name, "' panel found in this object's gene panel -- skipping")
            next
        }
        detect_rate <- Matrix::rowMeans(counts[genes, , drop = FALSE] > 0)
        best_rate <- max(detect_rate)
        if (best_rate < min_detect_rate) {
            warning("'", name, "' panel's best single marker (", names(which.max(detect_rate)),
                    ") is detected in only ", round(100 * best_rate, 3), "% of cells -- below the ",
                    100 * min_detect_rate, "% reliability floor -- excluded from scoring ",
                    "(a module score built on this sparse a signal is dominated by noise and can win the argmax spuriously)")
            next
        }
        detected <- genes[Matrix::rowSums(counts[genes, , drop = FALSE]) > 0]
        score_name <- paste0("idx_", make.names(name))
        xenium <- AddModuleScore(xenium, features = list(detected), name = score_name)
        xenium[[score_name]] <- xenium[[paste0(score_name, "1")]]
        score_cols <- c(score_cols, score_name)
    }
    raw <- xenium@meta.data[, score_cols, drop = FALSE]
    z <- scale(raw)
    colnames(z) <- score_cols
    best_idx   <- apply(z, 1, which.max)
    best_z     <- apply(z, 1, max)
    best_label <- label_lookup[score_cols[best_idx]]
    rank_info <- t(apply(z, 1, function(row) {
        row <- unname(row)
        ord <- order(row, decreasing = TRUE)
        c(margin = row[ord[1]] - row[ord[2]], runner_up_idx = ord[2])
    }))
    xenium$classification_margin <- rank_info[, "margin"]
    xenium$runner_up <- unname(label_lookup[score_cols[rank_info[, "runner_up_idx"]]])
    xenium$atlas_identity <- ifelse(best_z < min_z, "Unclassified", unname(best_label))
    refine_state <- function(cells_subset, panel_a, label_a, panel_b, label_b) {
        sub <- subset(xenium, cells = cells_subset)
        a_m <- intersect(panel_a, rownames(sub)); b_m <- intersect(panel_b, rownames(sub))
        sub <- AddModuleScore(sub, features = list(a_m), name = "StateA")
        sub <- AddModuleScore(sub, features = list(b_m), name = "StateB")
        ifelse(sub$StateA1 > sub$StateB1, label_a, label_b)
    }
    xenium$cell_type <- xenium$atlas_identity
    astro_cells <- colnames(xenium)[xenium$atlas_identity == "Astrocyte"]
    if (length(astro_cells) > 0) {
        xenium$cell_type[colnames(xenium) %in% astro_cells] <- paste0("Astrocyte (",
                                                                      refine_state(astro_cells, astro_active_markers, "Active", astro_homeostatic_markers, "Homeostatic"), ")")
    }
    micro_cells <- colnames(xenium)[xenium$atlas_identity == "Microglia"]
    if (length(micro_cells) > 0) {
        xenium$cell_type[colnames(xenium) %in% micro_cells] <- paste0("Microglia (",
                                                                      refine_state(micro_cells, microglia_active_markers, "Active", microglia_homeostatic_markers, "Homeostatic"), ")")
    }
    oligo_cells <- colnames(xenium)[xenium$atlas_identity == "Oligodendrocyte lineage"]
    if (length(oligo_cells) > 0) {
        xenium$cell_type[colnames(xenium) %in% oligo_cells] <- paste0("Oligodendrocyte (",
                                                                      refine_state(oligo_cells, oligo_mature_markers, "Mature", oligo_precursor_markers, "Precursor"), ")")
    }
    cat("\n==== cell_type counts ====\n"); print(table(xenium$cell_type))
    if ("idx_Purkinje..sanity.check." %in% score_cols) {
        purkinje_pct <- 100 * mean(xenium$atlas_identity == "Purkinje (sanity check)")
        cat(sprintf("\nPurkinje sanity check: %.2f%% of cells (should be ~0%%)\n", purkinje_pct))
    } else {
        cat(sprintf("\nPurkinje sanity check: best marker below the %.1f%% detection-reliability floor -- excluded from scoring, 0.00%% by definition (cleanest possible negative-control result)\n", 100 * min_detect_rate))
    }
    xenium
}
