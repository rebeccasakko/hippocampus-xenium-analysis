# ============================================================
# Group comparison statistics: Welch's t-test + BH correction
# per cell-type category (gi vs. gii)
#
# First run interactively for the Hippocampus, then written up as a
# function (2026-09-28) so the same method can be applied to Optic Tract and
# both regions stay on record as using identical statistical
# treatment, not just identical composition-summary code.
#
# Welch's t-test (R's t.test() default, var.equal = FALSE) rather
# than Student's: with n = 3 per group there's no reliable way to
# verify an equal-variance assumption, so Welch's is the safer
# default regardless of what the data happen to look like.
#
# Benjamini-Hochberg correction is applied across however many
# cell-type categories are actually tested for a given region (15 for
# the 12-panel hippocampus run: 12 identities, three of them split into
# two states, plus Unclassified; a different number for the optic tract,
# whose candidate list differs -- see R/marker_panels.R) rather than a
# fixed number, so the correction always reflects the real multiple-
# comparisons burden for that region's own candidate list. The
# percentages are compositional (they sum to 100 within a mouse) and
# "Unclassified" is tested like any category; both are known
# limitations of this per-category test (audit issue #6, finding M5).
# ============================================================

# `composition_df` is the per-sample, per-category data.frame returned
# as build_composition_summary(...)$composition_df (columns: sample_id,
# group_id, cell_type, n_cells, total_cells, pct, replicate_num).
build_stats_table <- function(composition_df) {
    composition_df %>%
        group_by(cell_type) %>%
        summarise(
            p_value = tryCatch(
                t.test(pct ~ group_id, data = pick(everything()))$p.value,
                error = function(e) NA_real_
            ),
            .groups = "drop"
        ) %>%
        mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
        arrange(p_value)
}

# ============================================================
# Usage:
#   hippocampus_stats <- build_stats_table(hippocampus_summary$composition_df)
#   optic_tract_stats  <- build_stats_table(optic_tract_summary$composition_df)
# ============================================================
