# _setup.R -------------------------------------------------------------------------------
# Session setup for R/run_full_pipeline.R: package check, input manifest, input check.
# Base R only, so it runs before any package is attached. Sourced once by
# R/run_full_pipeline.R, right after the root-anchor check and before any library() call.
# Added by the audit (issue #6); see docs/review/workflow-review.md section 2.

# Packages -------------------------------------------------------------------------------
# Everything the pipeline attaches or calls by namespace. DESeq2 comes from Bioconductor.
PACKAGES.REQUIRED <- c("Seurat", "dplyr", "tidyr", "ggplot2", "scales", "sp", "DESeq2",
                       "tibble", "pheatmap", "Matrix")
PACKAGES.BIOC     <- c("DESeq2")

#' Stop early, with the install lines, if any required package is missing.
#' Runs before any library() call so a missing package fails in the first second rather than
#' partway through a sample.
#' @param pkgs Character vector of required packages.
#' @param bioc Which of them are installed through BiocManager.
#' @return TRUE, invisibly, when everything is present.
fn.check_packages <- function(pkgs = PACKAGES.REQUIRED, bioc = PACKAGES.BIOC) {
  have <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
  missing <- pkgs[!have]
  if (length(missing) > 0) {
    cran <- setdiff(missing, bioc)
    bc   <- intersect(missing, bioc)
    msg  <- paste0("Missing R package(s): ", paste(missing, collapse = ", "), "\nInstall with:")
    if (length(cran) > 0)
      msg <- paste0(msg, "\n  install.packages(c(", paste0('"', cran, '"', collapse = ", "), "))")
    if (length(bc) > 0)
      msg <- paste0(msg, "\n  if (!requireNamespace(\"BiocManager\", quietly = TRUE)) install.packages(\"BiocManager\")",
                    "\n  BiocManager::install(c(", paste0('"', bc, '"', collapse = ", "), "))")
    stop(msg, call. = FALSE)
  }
  invisible(TRUE)
}

# Input manifest -------------------------------------------------------------------------
#' One row per input a region run reads: the six Xenium bundles and that region's six lasso
#' files. The driver, its error messages and README "Getting the data" all read from this table.
#' Presence of a bundle is tested on cells.csv.gz inside it (the one file the code reads by
#' name; Seurat::LoadXenium() reads the rest of the bundle).
#' @param xenium_dirs Named character vector, names = sample ids, values = bundle folders.
#' @param lassos Named character vector, names = sample ids, values = lasso CSV paths.
#' @param region Label used in the key column (e.g. "hippocampus").
#' @return data.frame with columns key, stage, required, how, path, test, source.
fn.data_manifest <- function(xenium_dirs, lassos, region = "region") {
  ids <- names(xenium_dirs)
  bundles <- data.frame(
    key      = paste0("xenium_", ids),
    stage    = "0,1",
    required = TRUE,
    how      = "manual",
    path     = unname(xenium_dirs),
    test     = file.path(unname(xenium_dirs), "cells.csv.gz"),
    source   = paste0("Xenium output bundle for ", ids, " (10x run 0071470; several GB). ",
                      "Unpublished lab data: request from the Abisambra Lab and point xenium_dirs at it in config.local.R."),
    stringsAsFactors = FALSE)
  lids <- names(lassos)
  lasso_rows <- data.frame(
    key      = paste0("lasso_", region, "_", lids),
    stage    = "0,1,1c",
    required = TRUE,
    how      = "manual",
    path     = unname(lassos),
    test     = unname(lassos),
    source   = paste0("Xenium Explorer selection export (CSV, two comment lines then Selection,X,Y,Class,Color) ",
                      "for ", lids, ", ", region, ". Request from the author; point ", region, "_lassos at it in config.local.R."),
    stringsAsFactors = FALSE)
  rbind(bundles, lasso_rows)
}

#' Check every input before any slow work and stop with instructions if a required one is
#' missing. Prints a status table so a run log records what was used.
#' @param manifest The data.frame from fn.data_manifest().
#' @param stop.on.missing Stop (TRUE) or only warn (FALSE) when a required input is absent.
#' @return The manifest with a logical `present` column, invisibly.
fn.check_inputs <- function(manifest, stop.on.missing = TRUE) {
  m <- manifest
  m$present <- !is.na(m$test) & nzchar(m$test) & file.exists(m$test)

  cat("\nInput check\n", strrep("-", 100), "\n", sep = "")
  cat(sprintf("  %-26s %-7s %-8s %-7s %s\n", "input", "stage", "need", "status", "path"))
  for (i in seq_len(nrow(m)))
    cat(sprintf("  %-26s %-7s %-8s %-7s %s\n", m$key[i], m$stage[i],
                if (m$required[i]) "required" else "optional",
                if (m$present[i]) "OK" else "MISSING", m$path[i]))
  cat(strrep("-", 100), "\n", sep = "")

  miss <- m[!m$present, , drop = FALSE]
  if (nrow(miss) > 0) {
    cat("\nMissing input(s):\n")
    for (i in seq_len(nrow(miss))) cat(sprintf("  %-26s %s\n", miss$key[i], miss$source[i]))
    cat("  See README.md > Getting the data, and config.local.example.R for the path overrides.\n")
  }

  short <- m$key[!m$present & m$required]
  if (length(short) > 0) {
    msg <- paste0("Missing required input(s): ", paste(short, collapse = ", "), ".")
    if (isTRUE(stop.on.missing)) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  }
  invisible(m)
}
