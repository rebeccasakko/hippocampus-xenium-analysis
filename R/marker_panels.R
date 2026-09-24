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

astro_active_markers <- c("Gfap", "Serpina3n")
astro_homeostatic_markers <- c("Aldh1l1", "Sox9")
microglia_active_markers <- c("Trem2", "Cst7", "Lpl")
microglia_homeostatic_markers <- c("Tmem119", "Sall1", "P2ry13")
oligo_mature_markers <- c("Mog", "Mag", "Cnp", "Opalin")
oligo_precursor_markers <- c("Pdgfra", "Cspg4", "Gpr17", "Vcan")
