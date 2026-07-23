#!/usr/bin/env Rscript
# =============================================================================
#
#   SigClust: Signature-based Cluster Annotation for Single-Cell Data
#
# =============================================================================
#
# WHAT THIS TOOL DOES:
#   Given a clustered Seurat object and pre-computed UCell signature scores,
#   SigClust identifies which gene-set signatures are statistically enriched
#   in each cluster, and produces publication-quality figures showing:
#   (A) A heatmap of enrichment strength (log2 Odds Ratio) for every
#       cluster x signature combination
#   (B) A UMAP where each cluster is colored by its dominant signature
#
# ─────────────────────────────────────────────────────────────────────────────
# REQUIRED INPUTS:
#   1. Seurat object (.rds) with UMAP embedding + cluster identities
#   2. UCell scores — either:
#      (a) A CSV file: rows = cell barcodes, columns = signature names, values = scores
#      (b) Already present in Seurat metadata (columns ending in _UCell)
#
# ─────────────────────────────────────────────────────────────────────────────
# HOW TO RUN (3 options):
#
#   OPTION A — Interactive R session (simplest):
#     source("sigclust_annotate.R")
#     results <- sigclust_annotate(
#       seurat_obj = readRDS("my_clustered_object.rds"),
#       scores     = read.csv("my_ucell_scores.csv", row.names = 1),
#       output_dir = "my_output/",
#       label      = "FN_Diff"
#     )
#
#   OPTION B — Command line:
#     Rscript sigclust_annotate.R --seurat obj.rds --scores scores.csv \
#       --output results/ --label my_sample
#
#   OPTION C — Use scores from Seurat metadata directly:
#     source("sigclust_annotate.R")
#     obj <- readRDS("my_object.rds")
#     # Extract UCell columns from metadata as the scores table:
#     ucell_cols <- grep("_UCell$", colnames(obj@meta.data), value = TRUE)
#     scores <- obj@meta.data[, ucell_cols]
#     colnames(scores) <- gsub("_UCell$", "", colnames(scores))
#     results <- sigclust_annotate(obj, scores, output_dir = "output/")
#
# ─────────────────────────────────────────────────────────────────────────────
# OUTPUTS:
#   <output_dir>/figures/
#     <label>_or_heatmap.png         — Blue-white-red enrichment heatmap
#     <label>_functional_umap.png    — UMAP colored by cluster identity
#     <label>_combined_panel.png     — Panel figure: (A) heatmap + (B) UMAP
#   <output_dir>/tables/
#     <label>_enrichment_full.csv    — All Fisher test results (every test)
#     <label>_cluster_annotations.csv — Summary: one row per cluster
#     <label>_summary_stats.csv      — Run parameters and counts
#
# ─────────────────────────────────────────────────────────────────────────────
# STATISTICAL METHOD:
#   For each signature, cells are ranked by UCell score. The top N% (default 10%)
#   are marked as "high" for that signature. A one-sided Fisher exact test asks:
#   "Are high-scoring cells over-represented in this cluster?"
#
#   The Odds Ratio (OR) quantifies how much more likely a top-scoring cell is
#   to be found in a given cluster vs. elsewhere:
#     OR = 1:   No enrichment (random distribution)
#     OR = 2:   2x more likely (log2 OR = 1) — our default significance threshold
#     OR = 4:   4x more likely (log2 OR = 2)
#     OR = 100: Nearly all top cells are in one cluster
#
#   Multiple testing correction: Benjamini-Hochberg FDR across all tests.
#
# ─────────────────────────────────────────────────────────────────────────────
# DEPENDENCIES:
#   install.packages(c("Seurat", "ggplot2", "pheatmap", "patchwork",
#                      "png", "dplyr"))
#   # Optional for CLI mode: install.packages("optparse")
#
# CITATION: If you use SigClust, please cite:
#   Khan Lab, NCI. SigClust: Signature-based Cluster Annotation (2026).
#   GitHub: [repository URL]
#
# =============================================================================


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 1: LOAD PACKAGES
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(pheatmap)
  library(patchwork)
  library(png)
  library(dplyr)
})


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 2: DEFAULT PARAMETERS (user-overridable)
# ═══════════════════════════════════════════════════════════════════════════════

DEFAULT_PARAMS <- list(

  # ─── Statistical thresholds ───────────────────────────────────────────────
  min_max_score = 0.2,
  # ^ Skip signatures whose max UCell score across ALL cells is below this.
  #   If no cell scores above 0.2 on a signature, that signature is considered
  #   "not expressed" in this dataset and is excluded from testing.
  #   Increase to 0.3 for stringent filtering; decrease to 0.1 to be inclusive.

  top_pct = 0.10,
  # ^ What fraction of cells counts as "high" for each signature.
  #   Default: top 10%. This means for each signature, the 10% of cells with
  #   the highest UCell scores are labeled "top" and tested for cluster enrichment.
  #   Decrease to 0.05 for very specific hits; increase to 0.25 for broader patterns.

  min_or = 2.0,
  # ^ Minimum Odds Ratio to label a cluster as enriched.
  #   OR=2 means top cells are 2x more likely in that cluster than expected.
  #   OR=4 (stricter) or OR=1.5 (more inclusive) are reasonable alternatives.

  fdr_threshold = 0.01,
  # ^ False Discovery Rate cutoff after Benjamini-Hochberg correction.
  #   All Fisher tests are corrected together. Only results below this FDR
  #   are considered significant.

  log2_or_display = 1.0,
  # ^ log2(OR) threshold for marking "*" in the heatmap.
  #   log2(OR) > 1 means OR > 2. Adjust to 2 (OR>4) for stricter display.

  # ─── Visualization parameters ────────────────────────────────────────────
  point_size   = 0.9,     # UMAP dot size
  point_alpha  = 0.9,     # UMAP dot opacity (0=invisible, 1=solid)
  label_size   = 3.0,     # Cluster label text size on UMAP
  heatmap_fontsize_row = 10,
  heatmap_fontsize_col = 11,
  heatmap_cellwidth    = 35,
  heatmap_cellheight   = 16,
  dpi = 300,

  # ─── Color palettes ──────────────────────────────────────────────────────
  # Heatmap: diverging blue-white-red
  heatmap_colors = c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
                     "white",
                     "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"),

  # UMAP: up to 40 distinct vivid colors
  sig_palette = c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628",
    "#F781BF", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
    "#FFD92F", "#E5C494", "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
    "#66A61E", "#E6AB02", "#A6761D", "#666666", "#8DD3C7", "#BEBADA",
    "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#BC80BD",
    "#CCEBC5", "#FFED6F", "#E31A1C", "#1F78B4", "#33A02C", "#FF7F00",
    "#6A3D9A", "#B15928", "#999999", "#000000")
)


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 3: MAIN FUNCTION
# ═══════════════════════════════════════════════════════════════════════════════

sigclust_annotate <- function(seurat_obj, scores, output_dir = "sigclust_output",
                              label = "sample", params = list()) {

  # Merge user params with defaults
  p <- DEFAULT_PARAMS
  for (nm in names(params)) p[[nm]] <- params[[nm]]

  # ─── 3.1 Validate inputs ─────────────────────────────────────────────────
  cat("=====================================================================\n")
  cat(sprintf("  SigClust: %s\n", label))
  cat("=====================================================================\n\n")

  stopifnot("umap" %in% Reductions(seurat_obj))
  stopifnot(length(levels(Idents(seurat_obj))) > 0)

  common_cells <- intersect(colnames(seurat_obj), rownames(scores))
  if (length(common_cells) == 0)
    stop("No matching cell barcodes between object and scores.")

  cat(sprintf("  Cells matched:  %d / %d (%.1f%%)\n",
              length(common_cells), ncol(seurat_obj),
              length(common_cells) / ncol(seurat_obj) * 100))
  cat(sprintf("  Signatures:     %d\n", ncol(scores)))
  cat(sprintf("  Clusters:       %d\n\n", length(levels(Idents(seurat_obj)))))

  scores <- scores[common_cells, , drop = FALSE]
  clusters <- Idents(seurat_obj)[common_cells]

  # ─── 3.2 Filter signatures ──────────────────────────────────────────────
  sig_max <- apply(scores, 2, max, na.rm = TRUE)
  retained <- names(sig_max[sig_max >= p$min_max_score])
  dropped  <- names(sig_max[sig_max < p$min_max_score])
  cat(sprintf("  Retained sigs:  %d (max >= %.2f)\n", length(retained), p$min_max_score))
  if (length(dropped) > 0)
    cat(sprintf("  Dropped sigs:   %d [%s]\n", length(dropped),
                paste(head(dropped, 5), collapse = ", ")))
  if (length(retained) == 0) stop("No signatures pass threshold.")

  # ─── 3.3 Fisher exact tests ─────────────────────────────────────────────
  cat("\n  Running Fisher exact tests...\n")
  results_list <- list()
  for (sig in retained) {
    scores_vec <- scores[, sig]
    threshold <- quantile(scores_vec, 1 - p$top_pct, na.rm = TRUE)
    is_top <- scores_vec >= threshold
    n_top <- sum(is_top)

    for (cl in levels(clusters)) {
      in_cl <- clusters == cl
      a <- sum(is_top & in_cl)
      b <- sum(is_top & !in_cl)
      cc <- sum(!is_top & in_cl)
      d <- sum(!is_top & !in_cl)
      ft <- fisher.test(matrix(c(a, b, cc, d), nrow = 2), alternative = "greater")

      results_list[[length(results_list) + 1]] <- data.frame(
        signature = sig, cluster = cl,
        n_top_in_cluster = a, n_cluster = a + cc, n_top_total = n_top,
        pct_top_in_cluster = round(a / max(a + cc, 1) * 100, 1),
        odds_ratio = round(ft$estimate, 2), p_value = ft$p.value,
        stringsAsFactors = FALSE)
    }
  }

  enrich_df <- do.call(rbind, results_list)
  enrich_df$fdr <- p.adjust(enrich_df$p_value, method = "BH")
  enrich_df$odds_ratio <- pmin(pmax(enrich_df$odds_ratio, 0.001), 1000)
  enrich_df$log2_or <- log2(enrich_df$odds_ratio)
  enrich_df$is_significant <- enrich_df$fdr < p$fdr_threshold &
                              enrich_df$log2_or > p$log2_or_display

  cat(sprintf("  Tests run:      %d\n", nrow(enrich_df)))
  cat(sprintf("  Significant:    %d\n\n", sum(enrich_df$is_significant)))

  # ─── 3.4 Assign cluster annotations ─────────────────────────────────────
  annotations <- enrich_df |>
    filter(is_significant) |>
    group_by(cluster) |>
    slice_max(log2_or, n = 1, with_ties = FALSE) |>
    ungroup() |>
    select(cluster, signature, odds_ratio, log2_or, fdr,
           n_top_in_cluster, pct_top_in_cluster)

  all_cl <- levels(clusters)
  unassigned <- setdiff(all_cl, annotations$cluster)
  if (length(unassigned) > 0) {
    annotations <- rbind(annotations, data.frame(
      cluster = unassigned, signature = "Unassigned",
      odds_ratio = NA, log2_or = NA, fdr = NA,
      n_top_in_cluster = NA, pct_top_in_cluster = NA))
  }
  annotations <- annotations[order(as.integer(annotations$cluster)), ]
  n_assigned <- sum(annotations$signature != "Unassigned")
  cat(sprintf("  Clusters assigned: %d / %d\n\n", n_assigned, length(all_cl)))

  # ─── 3.5 Create output directories ──────────────────────────────────────
  fig_dir <- file.path(output_dir, "figures")
  tab_dir <- file.path(output_dir, "tables")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

  # ─── 3.6 HEATMAP: blue-white-red, * for significant ─────────────────────
  cat("  Generating heatmap...\n")

  all_sigs <- sort(unique(enrich_df$signature))
  or_mat <- matrix(0, nrow = length(all_sigs), ncol = length(all_cl),
                   dimnames = list(all_sigs, paste0("C", all_cl)))
  sig_mat <- matrix("", nrow = length(all_sigs), ncol = length(all_cl),
                    dimnames = list(all_sigs, paste0("C", all_cl)))
  for (i in seq_len(nrow(enrich_df))) {
    rn <- enrich_df$signature[i]; cn <- paste0("C", enrich_df$cluster[i])
    if (rn %in% rownames(or_mat) && cn %in% colnames(or_mat)) {
      or_mat[rn, cn] <- enrich_df$log2_or[i]
      if (enrich_df$is_significant[i]) sig_mat[rn, cn] <- "*"
    }
  }
  or_mat[!is.finite(or_mat)] <- 0

  row_keep <- rowSums(abs(or_mat) > 0.1) > 0
  or_mat_f <- or_mat[row_keep, , drop = FALSE]
  sig_mat_f <- sig_mat[row_keep, , drop = FALSE]

  max_abs <- min(ceiling(max(abs(or_mat_f), na.rm = TRUE)), 10)
  if (!is.finite(max_abs) || max_abs < 1) max_abs <- 5
  breaks <- seq(-max_abs, max_abs, length.out = 101)
  colors <- colorRampPalette(p$heatmap_colors)(100)

  heatmap_path <- file.path(fig_dir, sprintf("%s_or_heatmap.png", label))
  png(heatmap_path, width = max(10, ncol(or_mat_f) * 1.0 + 5),
      height = max(8, nrow(or_mat_f) * 0.4 + 3), units = "in", res = p$dpi)
  pheatmap(or_mat_f, color = colors, breaks = breaks,
           display_numbers = sig_mat_f, number_color = "black", fontsize_number = 14,
           cluster_rows = TRUE, cluster_cols = TRUE, fontsize = 11,
           fontsize_row = p$heatmap_fontsize_row, fontsize_col = p$heatmap_fontsize_col,
           cellwidth = p$heatmap_cellwidth, cellheight = p$heatmap_cellheight,
           main = sprintf("%s: log2(OR) [* = FDR<%.3f & log2OR>%.1f]",
                          label, p$fdr_threshold, p$log2_or_display),
           angle_col = 45, border_color = "grey80")
  dev.off()
  cat(sprintf("    -> %s\n", basename(heatmap_path)))

  # ─── 3.7 FUNCTIONAL UMAP ────────────────────────────────────────────────
  cat("  Generating UMAP...\n")

  cluster_sig_map <- setNames(annotations$signature, annotations$cluster)
  umap_coords <- Embeddings(seurat_obj, "umap")[common_cells, ]
  umap_df <- data.frame(UMAP1 = umap_coords[,1], UMAP2 = umap_coords[,2],
                        cluster = as.character(clusters), stringsAsFactors = FALSE)
  umap_df$top_sig <- cluster_sig_map[umap_df$cluster]
  umap_df$top_sig[is.na(umap_df$top_sig)] <- "Unassigned"

  centroids <- umap_df |>
    group_by(cluster) |>
    summarize(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2), .groups = "drop") |>
    mutate(label = paste0("C", cluster, ": ", cluster_sig_map[cluster]))

  uniq_sigs <- sort(unique(umap_df$top_sig[umap_df$top_sig != "Unassigned"]))
  sig_colors <- setNames(p$sig_palette[seq_along(uniq_sigs)], uniq_sigs)
  sig_colors["Unassigned"] <- "grey80"
  all_levels <- c(uniq_sigs, "Unassigned")
  umap_df$top_sig <- factor(umap_df$top_sig, levels = all_levels)
  umap_df <- umap_df[order(umap_df$top_sig == "Unassigned", decreasing = TRUE), ]

  p_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = top_sig)) +
    geom_point(size = p$point_size, alpha = p$point_alpha, stroke = 0) +
    scale_color_manual(values = sig_colors, name = "Identity") +
    geom_text(data = centroids, aes(x = UMAP1, y = UMAP2, label = label),
              color = "white", size = p$label_size + 0.2, fontface = "bold",
              inherit.aes = FALSE, nudge_y = 0.3) +
    geom_text(data = centroids, aes(x = UMAP1, y = UMAP2, label = label),
              color = "black", size = p$label_size, fontface = "bold",
              inherit.aes = FALSE) +
    theme_minimal(base_size = 13) +
    theme(panel.grid = element_blank(), axis.title = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank(),
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA),
          legend.key.size = unit(0.5, "cm")) +
    guides(color = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 1)) +
    labs(title = sprintf("%s: Cluster Annotation", label),
         subtitle = sprintf("%d/%d clusters annotated (OR>%.0f, FDR<%.3f)",
                            n_assigned, length(all_cl), 2^p$log2_or_display, p$fdr_threshold))

  umap_path <- file.path(fig_dir, sprintf("%s_functional_umap.png", label))
  ggsave(umap_path, p_umap, width = 11, height = 9, dpi = p$dpi, bg = "white")
  cat(sprintf("    -> %s\n", basename(umap_path)))

  # ─── 3.8 COMBINED PANEL ─────────────────────────────────────────────────
  cat("  Generating combined panel...\n")
  heatmap_raster <- readPNG(heatmap_path)
  p_heat <- ggplot() +
    annotation_raster(heatmap_raster, xmin=0, xmax=1, ymin=0, ymax=1) +
    theme_void() + labs(title = "(A) Enrichment Heatmap") +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
  p_umap_b <- p_umap + labs(title = "(B) Functional UMAP")
  combined <- p_heat + p_umap_b + plot_layout(widths = c(1.2, 1)) +
    plot_annotation(title = sprintf("SigClust: %s", label),
                    theme = theme(plot.title = element_text(size = 16, face = "bold")))
  panel_path <- file.path(fig_dir, sprintf("%s_combined_panel.png", label))
  ggsave(panel_path, combined, width = 22, height = 10, dpi = p$dpi, bg = "white")
  cat(sprintf("    -> %s\n", basename(panel_path)))

  # ─── 3.9 Save tables ────────────────────────────────────────────────────
  cat("  Saving tables...\n")
  write.csv(enrich_df, file.path(tab_dir, sprintf("%s_enrichment_full.csv", label)),
            row.names = FALSE)
  write.csv(annotations, file.path(tab_dir, sprintf("%s_cluster_annotations.csv", label)),
            row.names = FALSE)
  summary_df <- data.frame(label = label, n_cells = length(common_cells),
                           n_clusters = length(all_cl), n_sigs_tested = length(retained),
                           n_sigs_dropped = length(dropped), n_assigned = n_assigned,
                           top_pct = p$top_pct, min_or = p$min_or,
                           fdr_threshold = p$fdr_threshold)
  write.csv(summary_df, file.path(tab_dir, sprintf("%s_summary_stats.csv", label)),
            row.names = FALSE)

  cat("\n=====================================================================\n")
  cat(sprintf("  DONE: %s | %d/%d clusters annotated\n", label, n_assigned, length(all_cl)))
  cat(sprintf("  Output: %s\n", output_dir))
  cat("=====================================================================\n")

  invisible(list(enrichment_df = enrich_df, annotations = annotations,
                 summary = summary_df, params = p))
}


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 4: COMMAND-LINE INTERFACE
# ═══════════════════════════════════════════════════════════════════════════════

if (!interactive() && length(commandArgs(trailingOnly = TRUE)) > 0) {
  if (!requireNamespace("optparse", quietly = TRUE))
    stop("Install optparse for CLI: install.packages('optparse')")
  library(optparse)
  opts <- list(
    make_option(c("-s", "--seurat"), help = "Seurat RDS path"),
    make_option(c("-c", "--scores"), help = "UCell scores CSV path"),
    make_option(c("-o", "--output"), default = "sigclust_output", help = "Output dir"),
    make_option(c("-l", "--label"), default = "sample", help = "Label prefix"),
    make_option("--top_pct", type = "double", default = 0.10),
    make_option("--min_or", type = "double", default = 2.0),
    make_option("--fdr", type = "double", default = 0.01))
  opt <- parse_args(OptionParser(option_list = opts))
  stopifnot(!is.null(opt$seurat), !is.null(opt$scores))
  obj <- readRDS(opt$seurat)
  sc <- read.csv(opt$scores, row.names = 1)
  sigclust_annotate(obj, sc, opt$output, opt$label,
                    list(top_pct = opt$top_pct, min_or = opt$min_or,
                         fdr_threshold = opt$fdr))
}
