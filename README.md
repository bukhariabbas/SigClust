# SigClust: Signature-based Cluster Annotation

**A standalone R tool for annotating single-cell clusters using pre-computed gene-set signature scores.**

## Overview

SigClust takes two inputs:
1. A **clustered Seurat object** (with UMAP + cluster identities)
2. **UCell signature scores** (one column per signature, one row per cell)

It then statistically identifies which signatures are enriched in which clusters, assigns functional labels, and produces publication-quality figures.

## Quick Start

```r
# Load the tool
source("sigclust_annotate.R")

# Load your data
obj    <- readRDS("my_clustered_seurat.rds")
scores <- read.csv("my_ucell_scores.csv", row.names = 1)

# Run SigClust
results <- sigclust_annotate(
  seurat_obj = obj,
  scores     = scores,
  output_dir = "my_results/",
  label      = "FN_Diff"
)
```

## How It Works

### Step 1: Filter signatures
Signatures whose maximum UCell score across all cells is below 0.2 are excluded (they're not expressed in this dataset and would produce meaningless tests).

### Step 2: Identify "high" cells
For each retained signature, cells are ranked by UCell score. The top 10% are labeled "high" for that signature.

### Step 3: Test enrichment (Fisher's exact test)
For each cluster × signature pair, a 2×2 contingency table is built:

|                | In cluster | Not in cluster |
|----------------|-----------|----------------|
| **Top 10%**    | a         | b              |
| **Bottom 90%** | c         | d              |

A one-sided Fisher's exact test asks: "Are top-scoring cells over-represented in this cluster?" The Odds Ratio (OR) quantifies the strength.

### Step 4: Correct for multiple testing
All p-values are FDR-corrected (Benjamini-Hochberg) across all tests simultaneously.

### Step 5: Assign annotations
Each cluster is labeled with the signature that has the highest OR (and passes significance thresholds). Clusters with no significant enrichment are labeled "Unassigned".

## Understanding the Odds Ratio

| log2(OR) | OR    | Meaning |
|----------|-------|---------|
| 0        | 1     | No enrichment — top cells distributed randomly |
| 1        | 2     | 2× enriched — **our default significance threshold** |
| 2        | 4     | 4× enriched |
| 3        | 8     | 8× enriched |
| 5+       | 32+   | Very strong — nearly all top cells in one cluster |
| 7+       | 128+  | Extreme — the cluster IS the signature |
| Negative | <1    | Depleted — top cells AVOID this cluster |

**Why are some ORs extremely high (100–460)?**

This happens when a cluster is essentially a pure population of one cell type. Example: If cluster 5 has 305 cells, and 282 of them are in the top-10% for "Pericyte_vSMC", then 92.5% of the cluster consists of top pericyte scorers. The Fisher test produces OR=460 because the overlap is near-complete. **This is correct behavior, not a bug** — it means Seurat's clustering and UCell's scoring are detecting the same biology independently.

## Outputs

### Figures
- **`_or_heatmap.png`** — Full log2(OR) heatmap: blue (depleted) → white (neutral) → red (enriched). Asterisks (`*`) mark statistically significant enrichments.
- **`_functional_umap.png`** — UMAP colored by each cluster's top-enriched signature.
- **`_combined_panel.png`** — Publication panel: (A) heatmap + (B) UMAP.

### Tables
- **`_enrichment_full.csv`** — Every test result (rows = cluster × signature combinations)
- **`_cluster_annotations.csv`** — One row per cluster: its identity, OR, FDR, % top cells
- **`_summary_stats.csv`** — Run parameters and summary counts

## Parameters

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `min_max_score` | 0.2 | Ignore signatures with max score below this |
| `top_pct` | 0.10 | Fraction of cells considered "high" per signature |
| `min_or` | 2.0 | Minimum OR for annotation |
| `fdr_threshold` | 0.01 | FDR cutoff for significance |
| `log2_or_display` | 1.0 | Heatmap "*" threshold (log2 OR > this) |

Override any parameter:
```r
results <- sigclust_annotate(obj, scores, params = list(top_pct = 0.25, fdr_threshold = 0.05))
```

## Command-Line Usage

```bash
Rscript sigclust_annotate.R \
  --seurat my_object.rds \
  --scores my_scores.csv \
  --output results/ \
  --label my_sample \
  --top_pct 0.10 \
  --min_or 2 \
  --fdr 0.01
```

## Dependencies

```r
install.packages(c("Seurat", "ggplot2", "pheatmap", "patchwork", "png", "dplyr"))
# Optional (for command-line mode only):
install.packages("optparse")
```

## Example with Seurat metadata scores

If UCell scores are already in your Seurat object's metadata:

```r
source("sigclust_annotate.R")
obj <- readRDS("my_object.rds")

# Extract UCell columns as the scores matrix
ucell_cols <- grep("_UCell$", colnames(obj@meta.data), value = TRUE)
scores <- obj@meta.data[, ucell_cols]
colnames(scores) <- gsub("_UCell$", "", colnames(scores))

results <- sigclust_annotate(obj, scores, output_dir = "output/", label = "my_sample")
```

## License

MIT
