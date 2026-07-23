# SigClust: Statistical Methods

## Overview

**SigClust** is a method for functionally annotating clusters in single-cell RNA-seq data using pre-computed gene-expression signature scores. It answers the question: *"What biological identity does each cluster represent?"*

The approach connects two independently-derived pieces of information:
1. **Cluster assignments** — groups of transcriptionally similar cells identified by graph-based clustering (e.g., Seurat's `FindClusters`)
2. **Signature scores** — per-cell quantification of gene-set activity (e.g., UCell, AUCell, or AddModuleScore)

For each cluster, SigClust asks: *"Are the highest-scoring cells for any given signature statistically over-represented in this cluster?"* If yes — and the enrichment is strong — that signature becomes the cluster's functional label.

This is a **post-hoc annotation method**, not a clustering method. It does not modify cluster boundaries or re-run dimensionality reduction. It simply tells you what each existing cluster *is*, using the gene programs you care about.

---

## Why This Approach?

Traditional cluster annotation methods rely on:
- **Marker genes** — identify differentially expressed genes per cluster, then manually look them up. Labor-intensive, subjective, non-reproducible.
- **Reference-based transfer** — project cells onto a labeled atlas (e.g., SingleR, Azimuth). Requires a high-quality reference for your tissue type, which may not exist.
- **Deconvolution** — estimate cell-type fractions (CIBERSORTx, xCell). Designed for bulk data, not for annotating individual clusters.

SigClust takes a different path: if you already have **curated gene signatures** that define the cell types or states you care about (from the literature, from CellMarker 2.0, from your own experiments), you can score them on every cell and let statistics identify which clusters match which signatures. This is:
- **Reproducible** — same signatures + same clusters = same annotations every time
- **Quantitative** — enrichment is measured by odds ratio with confidence intervals, not by eye
- **Multi-label capable** — a cluster can be significantly enriched for multiple signatures simultaneously (e.g., a cluster that is both "Macrophage" and "Hypoxia-high")
- **Agnostic** — works with any set of gene signatures (immune, muscle, metabolic, custom)

---

## Statistical Framework

### Input Requirements

| Input | Description | Format |
|-------|-------------|--------|
| **Seurat object** | Clustered cells with UMAP coordinates | `.rds` file with `seurat_clusters` metadata column |
| **Signature scores** | Per-cell scores for N gene signatures | Data frame: rows = cells, columns = signature names, values = scores (0–1 for UCell) |

The signature scores can come from any method that produces a per-cell continuous value:
- **UCell** (recommended) — rank-based, robust to library size differences
- **AUCell** — area under the recovery curve
- **Seurat::AddModuleScore** — z-score relative to control gene sets
- **GSVA/ssGSEA** — gene set variation analysis (less common for single cells)

### Step 1: Signature Relevance Filtering

**Problem:** Not every signature is informative for every dataset. A "B-cell" signature will score near-zero in a dataset of sorted tumor cells. Including such irrelevant signatures inflates the multiple-testing burden without adding information.

**Solution:** For each signature, compute its **maximum score** across all cells. If `max(score) < MIN_MAX_SCORE` (default: 0.2), the signature is dropped from analysis.

**Rationale for 0.2 threshold:** On UCell's 0–1 scale, a score of 0.2 means approximately 20% of the signature's genes are co-expressed above background in at least one cell. Below this, the signature is effectively absent from the dataset. The threshold is conservative — it only removes signatures with zero biological representation, not weak-but-real signals.

**This filter is applied per-dataset (or per-stratum if analyzing subsets separately).** A signature might be relevant in one subset but irrelevant in another.

### Step 2: Define "High-Scoring" Cells

For each retained signature, cells are classified as **"high"** if their score falls in the **top percentile** (default: top 10%, i.e., 90th percentile threshold).

**Why percentile-based rather than absolute threshold?**
- Different signatures have different score distributions (a T-cell signature with 200 genes will have different absolute scores than a 15-gene hypoxia signature)
- Percentile-based thresholding makes the definition comparable across signatures
- 10% is the standard in the literature (Tirosh et al. 2016, Science; Kinker et al. 2020, Nature Genetics; Gavish et al. 2023, Nature)

**Why 10% specifically?**
- Too stringent (e.g., top 1%) captures only extreme outliers, missing moderate-but-real expression
- Too permissive (e.g., top 50%) includes cells with background-level expression, diluting signal
- 10% balances sensitivity (enough cells to have statistical power per cluster) with specificity (genuinely high-scoring cells only)

The parameter `TOP_PCT` is user-configurable. For datasets with very sharp cell-type boundaries, top 5% may be more appropriate. For datasets with continuous gradients (e.g., developmental trajectories), top 25% may capture more biology.

### Step 3: Enrichment Testing (Fisher's Exact Test)

For every combination of **(signature S, cluster C)**, construct a 2×2 contingency table:

|                        | Cells in cluster C | Cells NOT in cluster C | Row total |
|------------------------|:------------------:|:----------------------:|:---------:|
| **Top cells for S**    | a                  | b                      | a + b     |
| **Non-top cells for S**| c                  | d                      | c + d     |
| **Column total**       | a + c              | b + d                  | N         |

Where:
- **a** = number of top-10% cells for signature S that are in cluster C
- **b** = number of top-10% cells for signature S that are NOT in cluster C
- **c** = number of non-top cells in cluster C
- **d** = number of non-top cells not in cluster C
- **N** = total cells in the dataset

**Test:** One-sided Fisher's exact test (`alternative = "greater"`), testing whether top-scoring cells are **over-represented** in the cluster.

**Why Fisher's exact (not chi-squared)?**
- Fisher's is exact (no large-sample approximation needed)
- Handles small cell counts correctly (some clusters have <100 cells)
- One-sided because we only care about enrichment (over-representation), not depletion

### Step 4: Odds Ratio Interpretation

The **Odds Ratio (OR)** from the Fisher test quantifies enrichment strength:

$$OR = \frac{a \cdot d}{b \cdot c}$$

**Interpretation scale:**

| OR | log₂(OR) | Meaning |
|----|----------|---------|
| 1 | 0 | No enrichment — top cells are distributed randomly across clusters |
| 2 | 1 | 2× more likely to find a top-scoring cell in this cluster than expected |
| 4 | 2 | 4× enrichment — moderate, biologically meaningful |
| 8 | 3 | Strong enrichment — the cluster is clearly dominated by this signature |
| 16+ | 4+ | Very strong — near-complete overlap between cluster and signature |
| 100+ | 6.6+ | Extreme — the cluster essentially IS the signature (see explanation below) |

**Why log₂(OR) for visualization?** Raw OR values span 0.01 to 500+ — an unmanageable range for heatmaps. Log₂ transformation compresses this to approximately [-3, +9], making patterns visible. It also makes the scale symmetric: log₂(OR) = +2 means "4× enriched" and log₂(OR) = -2 means "4× depleted."

### Step 5: Multiple-Testing Correction

With N_signatures × N_clusters tests (typically 30–40 × 10–20 = 300–800 tests per dataset), multiple-testing correction is essential.

**Method:** Benjamini-Hochberg FDR correction, applied across all tests within a dataset.

**Significance threshold:** FDR < 0.01 (default). This is deliberately strict because cluster annotations propagate to all downstream analyses — a false annotation has cascading effects. The threshold is user-configurable via `params$min_fdr`.

A cluster-signature pair is considered **significant** if:
1. FDR < 0.01, AND
2. OR ≥ 2 (at minimum 2× enrichment — avoids statistically significant but biologically trivial associations that can occur with very large cell counts)

### Step 6: Cluster Assignment

Each cluster receives a **dominant signature label** — the signature with the highest OR among all its significant associations.

**Important:** A cluster can have MULTIPLE significant signatures (visible in the heatmap), but only the top one becomes the UMAP label. The full multi-signature profile is preserved in the output CSV.

**Unassigned clusters:** If no signature passes both the OR ≥ 2 and FDR < 0.01 thresholds for a given cluster, it is labeled "Unassigned" and colored grey on the UMAP. This typically indicates:
- A transitional/intermediate population between two cell states
- A novel cell type not covered by the provided signatures
- A cluster driven by technical variation (e.g., a sample-specific batch effect)

---

## Why Are Some Odds Ratios Extremely High (100–500)?

This is **not a bug** — it is expected behavior that reflects how single-cell clustering works.

### Explanation

Graph-based clustering (Louvain/Leiden on a shared-nearest-neighbor graph) groups cells by **transcriptomic similarity**. When a dataset contains a discrete cell type (e.g., macrophages), those cells will form a tight, well-separated cluster because they share hundreds of co-expressed genes.

UCell (or any signature-scoring method) detects the **same set of co-expressed genes** that defined the cluster in the first place. So if cluster 5 was formed because 300 cells all express CD163/MRC1/STAB1/CSF1R/CD68, and your "Macrophage" signature contains those exact genes, then:
- Nearly 100% of cluster 5 cells will be in the top 10% for "Macrophage"
- Nearly 0% of non-cluster-5 cells will be in the top 10% for "Macrophage"

This produces a near-diagonal 2×2 table:

|                     | In cluster 5 (300 cells) | Not in cluster 5 (4,700 cells) |
|---------------------|:---:|:---:|
| **Top 10% Macrophage** | 282 | 218 |
| **Other 90%**          | 18  | 4,482 |

$$OR = \frac{282 \times 4482}{218 \times 18} = 322$$

### What different OR ranges mean biologically

| OR range | log₂(OR) | Biological interpretation |
|----------|----------|--------------------------|
| **100–500** | 6.6–9.0 | The cluster IS the cell type. Near-perfect 1:1 mapping. This means your clustering resolution is appropriate and your signature correctly identifies this population. |
| **10–100** | 3.3–6.6 | Strong enrichment. The cluster is dominated by this cell type but contains some additional cells (mixed population or over-clustered). |
| **2–10** | 1.0–3.3 | Moderate enrichment. Common for: (a) cell STATES (hypoxia, proliferation) that span multiple clusters, (b) mixed-identity clusters, or (c) signatures with partial gene overlap to other cell types. |
| **1** | 0 | No association. The signature has nothing to do with this cluster. |
| **<1** | <0 | Depletion. The cluster actively LACKS this signature's genes. Informative for negative identity (e.g., "this cluster is NOT immune"). |

### Mathematical intuition

OR becomes extreme when two conditions hold simultaneously:
1. **High specificity:** Most top-10% cells for a signature concentrate in ONE cluster (cell **a** dominates row 1)
2. **High purity:** That cluster contains very few cells that are NOT top-scoring (cell **c** is small)

When both specificity and purity exceed 90%, the cross-products in the OR formula produce values > 100. This is the Fisher test working correctly on a clean biological separation — not an artifact or an error.

---

## Outputs

### Figures

| Output | Description |
|--------|-------------|
| **Enrichment heatmap** | Signatures (rows) × clusters (columns), colored by log₂(OR). Blue = depletion, white = no association, red = enrichment. Cells with significant enrichment (FDR < 0.01, OR > threshold) marked with `*`. Both axes hierarchically clustered to reveal co-enrichment patterns. |
| **Functional UMAP** | Standard UMAP with each cell colored by its cluster's dominant signature. Unassigned clusters shown in grey. Cluster centroids labeled with signature name. |
| **Combined panel** | Side-by-side: (A) enrichment heatmap + (B) functional UMAP, for publication-ready display. |

### Tables

| Output | Description |
|--------|-------------|
| **Full enrichment table** (CSV) | Complete Fisher's exact results for every signature × cluster combination: n_top_in_cluster, n_cluster, n_top_total, pct_top_in_cluster, odds_ratio, p_value, fdr |
| **Cluster annotations** (CSV) | One row per cluster: cluster_id, dominant_signature, odds_ratio, fdr, n_significant_signatures |

---

## Parameters

| Parameter | Default | Description | How to adjust |
|-----------|---------|-------------|---------------|
| `TOP_PCT` | 0.10 | Percentile threshold for "high-scoring" cells | Decrease (0.05) for sharper cell types; increase (0.25) for continuous gradients |
| `MIN_MAX_SCORE` | 0.20 | Minimum max-score for a signature to be retained | Lower (0.1) to keep more signatures; raise (0.3) to be stricter |
| `MIN_OR` | 2.0 | Minimum OR for significant assignment | Raise (4.0) for stricter annotations; lower (1.5) for exploratory |
| `MIN_FDR` | 0.01 | FDR significance threshold | Standard; rarely needs changing |
| `LOG2_OR_THRESHOLD` | 1.0 | Minimum log₂(OR) shown in heatmap highlighting | Raise (2.0) to show only strong enrichments |

---

## Assumptions and Limitations

### Assumptions
1. **Signatures are biologically valid** — SigClust trusts that the gene lists you provide correctly represent the cell types/states they claim to. Garbage signatures produce garbage annotations.
2. **Clusters represent biological populations** — if clusters are driven by batch effects or doublets, annotations will be meaningless. QC and batch correction should be done before SigClust.
3. **UCell scores are comparable across signatures** — the percentile-based thresholding handles scale differences, but extremely short signatures (<5 genes) may produce noisy scores.

### Limitations
1. **Cannot discover novel cell types** — SigClust can only annotate clusters using the signatures you provide. An unlabeled cluster might be novel biology, or it might be a known cell type for which you didn't include a signature.
2. **Resolution-dependent** — annotations depend on the clustering resolution chosen upstream. Over-clustered data may split one cell type across multiple clusters; under-clustered data may merge distinct types.
3. **One dominant label per cluster** — the UMAP visualization shows only the top signature. The full multi-signature profile is in the CSV output. Always check the heatmap for multi-signature clusters.
4. **Spot-level data (Visium)** — each Visium spot contains multiple cells. SigClust can still be applied, but OR values will be lower (spots are inherently mixed) and interpretation shifts from "cell type" to "dominant program in this region."

---

## References

- **UCell:** Andreatta M & Carmona SJ (2021). UCell: Robust and scalable single-cell gene signature scoring. *Computational and Structural Biotechnology Journal* 19:3796-3798.
- **Fisher's exact test:** Fisher RA (1922). On the interpretation of χ² from contingency tables. *Journal of the Royal Statistical Society* 85:87-94.
- **Benjamini-Hochberg FDR:** Benjamini Y & Hochberg Y (1995). Controlling the false discovery rate. *Journal of the Royal Statistical Society B* 57:289-300.
- **Top-10% convention:** Tirosh I et al. (2016). Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq. *Science* 352:189-196.
- **CellMarker 2.0:** Hu C et al. (2023). CellMarker 2.0: an updated database of manually curated cell markers. *Nucleic Acids Research* 51:D1348-D1353.
- **Seurat clustering:** Hao Y et al. (2024). Dictionary learning for integrative, multimodal and scalable single-cell analysis. *Nature Biotechnology* 42:293-304.

---

## Citation

If you use SigClust in your research, please cite:

```
SigClust: Signature-based cluster annotation via enrichment testing.
https://github.com/bukhariabbas/SigClust
```
