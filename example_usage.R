#!/usr/bin/env Rscript
# =============================================================================
# example_usage.R — How to run SigClust on your own data
# =============================================================================
# This file shows 3 different ways to use SigClust.
# Copy the section that matches your situation.
# =============================================================================

# Load the SigClust tool (adjust path to where you saved it)
source("sigclust_annotate.R")


# =============================================================================
# EXAMPLE 1: Basic usage with a Seurat object + scores CSV
# =============================================================================
# You have:
#   - A Seurat object with UMAP + clusters (from FindClusters)
#   - A CSV of UCell scores (rows = cell barcodes, columns = signature names)

obj    <- readRDS("path/to/my_clustered_object.rds")
scores <- read.csv("path/to/my_ucell_scores.csv", row.names = 1)

results <- sigclust_annotate(
  seurat_obj = obj,
  scores     = scores,
  output_dir = "sigclust_results/",
  label      = "my_sample"
)

# Access results programmatically:
results$annotations       # One row per cluster with its label
results$enrichment_df     # Full Fisher test results
results$summary           # Summary statistics


# =============================================================================
# EXAMPLE 2: Scores already in Seurat metadata (UCell columns)
# =============================================================================
# If you ran UCell::AddModuleScore_UCell() on your object, scores are already
# stored as metadata columns ending in "_UCell". Extract them:

obj <- readRDS("path/to/my_scored_object.rds")

ucell_cols <- grep("_UCell$", colnames(obj@meta.data), value = TRUE)
scores <- obj@meta.data[, ucell_cols]
colnames(scores) <- gsub("_UCell$", "", colnames(scores))  # Remove _UCell suffix

results <- sigclust_annotate(obj, scores, output_dir = "output/", label = "my_tissue")


# =============================================================================
# EXAMPLE 3: Custom parameters (more or less strict)
# =============================================================================
# Stricter: only strong enrichments (OR > 4, FDR < 0.001)
results_strict <- sigclust_annotate(
  obj, scores, output_dir = "strict_results/", label = "strict",
  params = list(min_or = 4, fdr_threshold = 0.001, log2_or_display = 2)
)

# More inclusive: top 25% cells, OR > 1.5
results_broad <- sigclust_annotate(
  obj, scores, output_dir = "broad_results/", label = "broad",
  params = list(top_pct = 0.25, min_or = 1.5, log2_or_display = 0.58)
)


# =============================================================================
# EXAMPLE 4: Multiple samples in a loop
# =============================================================================
strata <- c("FN_Diff", "FN_Progenitor", "FN_Stem", "FP_Diff")

for (s in strata) {
  obj    <- readRDS(sprintf("rds/P2_%s_harmony.rds", s))
  scores <- read.csv(sprintf("tables/ucell_scores_%s.csv", s), row.names = 1)
  sigclust_annotate(obj, scores, output_dir = "results/", label = s)
}
