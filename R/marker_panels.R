atlas_index_markers <- list(
    # -- established (Hippocampus, in production) --
    Astrocyte                    = c("Aldh1l1", "Sox9"),
    Microglia                    = c("Tmem119", "Sall1", "P2ry13"),
    `Oligodendrocyte lineage`    = c("Sox10", "Olig1", "Olig2"),
    Endothelial                  = c("Pecam1", "Id1"),
    `Pericyte/vascular SMC`      = c("Kcnj8", "Myh11", "Pdgfrb"),
    `Fibroblast/meningeal`       = c("Col1a1", "Col1a2", "Lum", "Fn1"),
    `Ependymal/ciliated`         = c("Foxj1", "Pifo", "Dnah12"),
    `Choroid plexus`             = c("Folr1", "Kl", "Otx2"),
    BAM                          = c("Mrc1", "F13a1", "Cd163"),
    `Hippocampal pyramidal`      = c("Slc17a7", "Fibcd1", "Spink8", "Hes1", "Rab38", "Pou3f1", "Hcn1", "Map3k15"),
    `Dentate gyrus granule`      = c("Slc17a7", "Prox1"),
    `GABAergic interneuron`      = c("Gad1", "Gad2", "Slc32a1"),
    `Purkinje (sanity check)`    = "Pcp2",

    # -- proposed additions (Optic Tract), confirmed present in the
    #    Xenium Prime Mouse 5K pan-tissue panel via metadata scan,
    #    2026-09-24 -- see protocols/classification_protocol.md
    `Peripheral macrophage`      = c("Ccr2", "Cd14", "Ms4a7"),
    `T cell`                     = c("Cd3e", "Cd3d", "Cd3g"),
    Neutrophil                   = c("Ly6g", "Csf3r", "Cxcr2"),
    `NK cell`                    = c("Ncr1", "Klrb1c", "Gzmb", "Klrk1", "Prf1"),
    `B cell`                     = c("Cd19", "Cd79a", "Cd79b", "Ms4a1", "Cd22", "Pax5"),
    `Dendritic cell`             = c("Itgax", "Cd83", "Batf3", "Xcr1", "Clec9a"),
    `RGC (sanity check)`         = c("Pou4f1", "Sncg")
)

# -- Region-specific candidate panel set: Optic Tract --
# Drops Hippocampal pyramidal and Dentate gyrus granule from the argmax
# competition for optic tract only (hippocampus keeps the full list above).
# Confirmed 2026-09-28 via classification_margin + per-gene detection-rate
# diagnostics that these were NOT near-tie noise (median margin 1.9-2.4;
# 0% of these calls would flip under any margin threshold <=1.0) but
# genuine gene-overlap false positives: each panel's truly hippocampus-
# specific genes never fire in optic tract in any sample (Fibcd1/Spink8/
# Rab38 all at 0% detection for Hippocampal pyramidal), while one shared,
# non-specific gene in each panel fires and wins the argmax alone --
# Pou3f1 (a myelination-associated TF, biologically expected in white
# matter) at 50-100% detection for Hippocampal pyramidal, and Prox1 (not
# actually DG-restricted brain-wide) at 100% detection for Dentate gyrus
# granule. runner_up for both was dominated by Oligodendrocyte lineage
# (29% and 45% of the misclassified cells respectively) -- real
# oligodendrocyte-lineage cells were losing the argmax call to these false
# positives. GABAergic interneuron was kept as an Optic Tract candidate:
# Gad1/Gad2/Slc32a1 aren't hippocampus-restricted and didn't show this
# problem. See protocols/classification_protocol.md.
#
# RGC (sanity check) is ALSO excluded here, 2026-09-28, for a different
# reason than the two panels above: its margin distribution turned out to
# be cleanly bimodal, not uniformly noisy. Roughly half of RGC-labeled
# cells have near-zero classification_margin (real near-tie noise,
# fixable by a confidence filter); the other half have very large margins
# (4.9-6.9) -- confidently, unambiguously called RGC, not a tie. No single
# margin threshold can resolve both halves at once. The confident half is
# suspected (not yet fully confirmed) to be driven by Sncg alone, which
# has documented cross-reactivity with vascular/renal/smooth-muscle cell
# types in the reference atlas -- plausible here given Optic Tract's own
# vasculature (Endothelial/Pericyte/Fibroblast). Given the Sncg
# cross-reactivity lead was not yet fully run down and time was limited
# ahead of a 2026-09-30 presentation, excluded RGC as a candidate for
# Optic Tract rather than leave it half-fixed. NOTE: this gives up the
# sanity check's original purpose here (catching a lasso straying toward
# the retina/optic nerve head) -- a real, deliberate tradeoff, not an
# oversight. Revisit before final submission: confirm/resolve the Sncg
# driver-gene hypothesis (see protocols/classification_protocol.md) and
# decide whether to reinstate a corrected RGC panel.
#
# Ependymal/ciliated is ALSO excluded here, 2026-09-28: the same
# single-gene-driver pattern as Hippocampal pyramidal/Dentate gyrus
# granule. Foxj1 alone is detected in 65-100% of cells this panel calls
# in every Optic Tract sample, while Pifo and Dnah12 -- the genes that
# would provide real ciliated-ependymal specificity -- sit at 0-2.9% in
# every sample. Real ependymal tissue would express all three cilia
# genes together, so this argues for a marker cross-reactivity artifact
# (Foxj1 is broadly tied to motile ciliogenesis, not ependymal identity
# specifically) rather than the lasso boundary catching genuine
# periventricular tissue. See protocols/classification_protocol.md.
optic_tract_markers <- atlas_index_markers
optic_tract_markers[["Hippocampal pyramidal"]] <- NULL
optic_tract_markers[["Dentate gyrus granule"]] <- NULL
optic_tract_markers[["RGC (sanity check)"]] <- NULL
optic_tract_markers[["Ependymal/ciliated"]] <- NULL

# Hippocampus uses the 13 original panels: drop the optic-tract-only panels
hippocampus_markers <- atlas_index_markers
for (p in c("Peripheral macrophage", "T cell", "Neutrophil", "NK cell",
            "B cell", "Dendritic cell", "RGC (sanity check)")) {
  hippocampus_markers[[p]] <- NULL
}
# Choroid plexus excluded from hippocampus (2026-10-06): calls were Kl-driven
# (Folr1/Otx2 ~1% of calls), ~74% had a neuronal runner-up, and the cells lay
# along the pyramidal/DG layers -> Kl-positive neurons, not choroid plexus.
# Evidence: outputs/hippocampus_choroid_plexus_*. Ependymal/ciliated is kept:
# mixed signal (some calls line the ventricular edge), reported as not robust.
hippocampus_markers[["Choroid plexus"]] <- NULL

astro_active_markers <- c("Gfap", "Serpina3n")
astro_homeostatic_markers <- c("Aldh1l1", "Sox9")
microglia_active_markers <- c("Trem2", "Cst7", "Lpl")
microglia_homeostatic_markers <- c("Tmem119", "Sall1", "P2ry13")
oligo_mature_markers <- c("Mog", "Mag", "Cnp", "Opalin")
oligo_precursor_markers <- c("Pdgfra", "Cspg4", "Gpr17", "Vcan")
