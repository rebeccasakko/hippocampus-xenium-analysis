# Local path overrides for R/run_full_pipeline.R.
#
# Copy this file to config.local.R in the repository root and edit it. config.local.R is
# gitignored, so machine-specific paths never reach the repository. This is the only file
# you should need to touch to run the pipeline on another machine.
#
# R/run_full_pipeline.R sources it AFTER setting its repo-relative defaults and BEFORE reading
# any input, so uncomment only the lines you want to change. `dir.repo` (the repository root)
# and `sample_ids` (the six sample ids, in order) are already defined at that point.
#
# Inputs that cannot ship with the repository, and how to get them, are listed in README.md
# under "Getting the data". One block per key: what it is, the default, and the assignment
# commented out.


# Folder holding the inputs that cannot ship (the six Xenium bundles and the lasso CSVs).
# The repo default is a gitignored data/ folder inside the repository.
# Default: file.path(dir.repo, "data")
# dir.data <- "D:/xenium/hippocampus-project/data"


# The six Xenium output bundles, one folder per mouse, named by sample id. Each folder is
# what the Xenium instrument wrote for that section (cells.csv.gz, cell_feature_matrix, ...).
# Default: file.path(dir.data, "xenium", paste0("0071470_", sample_ids)), named by sample_ids
# xenium_dirs <- c(                                                    # the original author's layout
#   gi_1_b  = "/Volumes/Joe/Becca and Sakthi Files/Xenium_Data/Group_i/0071470_gi_1_b/",
#   gi_2_d  = "/Volumes/Joe/Becca and Sakthi Files/Xenium_Data/Group_i/0071470_gi_2_d/",
#   gi_3_f  = "/Volumes/Joe/Becca and Sakthi Files/Xenium_Data/Group_i/0071470_gi_3_f/",
#   gii_1_a = "/Volumes/Joe/Becca and Sakthi Files/Xenium_Data/Group_ii/0071470_gii_1_a/",
#   gii_2_c = "/Volumes/Joe/Becca and Sakthi Files/Xenium_Data/Group_ii/0071470_gii_2_c/",
#   gii_3_e = "/Volumes/Joe/Becca and Sakthi Files/Xenium_Data/Group_ii/0071470_gii_3_e/")


# Hippocampus lasso polygons: Xenium Explorer selection exports, one CSV per mouse.
# Default: file.path(dir.data, "lassos", "hippocampus", paste0("hippocampus_", sample_ids, ".csv"))
# hippocampus_lassos <- c(                                             # the original author's layout
#   gi_1_b  = "~/Desktop/hippocampus_gi_1_b.csv",
#   gi_2_d  = "~/Desktop/hippocampus_gi_2_d.csv",
#   gi_3_f  = "~/Desktop/hippocampus_gi_3_f.csv",
#   gii_1_a = "~/Desktop/hippocampus_gii_1_a.csv",
#   gii_2_c = "~/Desktop/hippocampus_gii_2_c.csv",
#   gii_3_e = "~/Desktop/hippocampus_gii_3_e.csv")


# Optic tract lasso polygons: the border-inclusive set redrawn on 2026-10-05 (one polygon per
# sample; file-to-sample pairing checked in docs/optic_tract_new_lassos_overview.png).
# Default: file.path(dir.data, "lassos", "optic_tract", c("GI-1", "GI-2", "GI-3", "GII-1", "GII-2", "GII-3") + "_coordinates.csv")
# optic_tract_lassos <- c(                                             # the original author's layout
#   gi_1_b  = "~/Desktop/B6_Xenium_Opt/GI-1_coordinates.csv",
#   gi_2_d  = "~/Desktop/B6_Xenium_Opt/GI-2_coordinates.csv",
#   gi_3_f  = "~/Desktop/B6_Xenium_Opt/GI-3_coordinates.csv",
#   gii_1_a = "~/Desktop/B6_Xenium_Opt/GII-1_coordinates.csv",
#   gii_2_c = "~/Desktop/B6_Xenium_Opt/GII-2_coordinates.csv",
#   gii_3_e = "~/Desktop/B6_Xenium_Opt/GII-3_coordinates.csv")
#
# The optic tract lassos used before 2026-10-05, kept here for the record only; no run uses them:
#   gi_1_b  = "~/Desktop/optic_tract_gi_1_b.csv",   gi_2_d  = "~/Desktop/optic_tract_gi_2_d.csv",
#   gi_3_f  = "~/Desktop/optic_tract_gi_3_f.csv",   gii_1_a = "~/Desktop/optic_tract_gii_1_a.csv",
#   gii_2_c = "~/Desktop/optic_tract_gii_2_c.csv",  gii_3_e = "~/Desktop/optic_tract_gii_3_e.csv"


# Where runs write their outputs. The tracked deliverables stay inside the repository by
# default so a rerun shows in `git diff`; redirect the whole tree, or nothing.
# Default: file.path(dir.repo, "outputs")
# dir.out <- "D:/xenium/hippocampus-project/outputs"


# Keys not listed here keep their defaults. When a key is added, add its block here, its
# default to the SETTINGS block of R/run_full_pipeline.R, and its row to README > Configuration,
# all in the same commit, so the three cannot drift.
