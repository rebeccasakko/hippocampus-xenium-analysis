# ============================================================
# Candidate gene panels -- differential expression (gi vs. gii)
#
# Pre-specified, literature-justified gene lists used ONLY for the
# pseudobulk DESeq2 secondary analysis (see protocols/differential_
# expression_protocol.md). NOT used by the Wilcoxon/FindMarkers primary
# analysis, which tests the full gene panel with no restriction.
#
# dg and gaba panels were pulled from the literature specifically for
# this project (see citations below); astro and microglia panels predate
# this file and were confirmed against the same standard.
# ============================================================

# --- Astrocyte: reactive-astrocyte markers ---
astro_candidates <- c("Gfap", "Aqp4", "Slc1a3", "Gja1", "Mapt", "Aldh1l1", "Sox9")

# --- Microglia: homeostatic vs. disease-associated-microglia (DAM) signature ---
microglia_candidates <- c("Trem2", "Cst7", "Lpl", "Tmem119", "Sall1", "P2ry13")

# --- Dentate gyrus granule: TBI-specific single-cell hippocampal findings ---
# Source: Arneson et al., 2018, Nature Communications, "Single cell
# molecular alterations reveal target cells and pathways of concussive
# brain injury" -- Id2, Npy, Penk, Ptprn flagged as altered in DG granule
# cells post-injury; Bdnf/Ntf3 corroborated independently by Truettner et
# al., 1999, J Neurotrauma (BDNF/NT-3 mRNA after fluid percussion injury).
# Npy, Penk, Ptprn did not survive the pseudobulk rowSums >= 10 filter in
# this dataset -- kept in the candidate list below for provenance/
# transparency, filtered out automatically by run_pseudobulk_panel_de().
dg_candidates <- c("Id2", "Npy", "Bdnf", "Ntf3", "Penk", "Ptprn")

# --- GABAergic interneuron: vulnerable-subtype + inhibitory-function markers ---
# Pvalb, Sst: parvalbumin- and somatostatin-positive interneurons are the
# two most consistently reported vulnerable subtypes in hippocampal TBI
# (e.g. eNeuro 2020 ENEURO.0195-19.2020; Cerebral Cortex 2022 32(5):1093).
# Gad1/Gad2: GABA-synthesizing enzymes, standard index of inhibitory tone.
# Slc12a5 (KCC2): chloride transporter, documented dysregulation after
# TBI underlying post-traumatic inhibitory failure (Frontiers Cell
# Neurosci 2025, 10.3389/fncel.2025.1572213).
gaba_candidates <- c("Pvalb", "Sst", "Gad1", "Gad2", "Slc12a5")
