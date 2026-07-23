# Methods: Multi-Signature Cluster Enrichment Analysis

## Overview

This analysis answers the question: **"What is each cluster on the UMAP?"** — not by
running new tools, but by leveraging the 40 gene-expression signatures already scored
via UCell to functionally annotate every Seurat cluster in every stratum.

In plain language: we already know how much each cell "looks like" a macrophage, a
myoblast, a T-cell, a pericyte, etc. (from UCell scoring). We also know which cluster
each cell belongs to (from Seurat's FindClusters). This analysis connects the two:
for each cluster, it asks "are the cells in this cluster unusually enriched for any
particular signature?" — and if so, that signature becomes the cluster's identity label.

---

## Step-by-Step Methodology

### Step 1: Input Data

For each of the 8 strata (FN/FP × Diff/Progenitor/Stem/NonMalignant):

- **Harmony-integrated Seurat object** (`P2_<stratum>_harmony.rds`)
  containing UMAP coordinates and Seurat cluster assignments (from FindClusters,
  resolution 0.5 for malignant, 0.8 for non-malignant)
- **Pre-computed UCell signature scores** (`batch*_<stratum>_scores.csv`)
  from the 64 SLURM scoring jobs — 40 signatures scored on every cell

### Step 2: Signature Filtering

Not every signature is relevant for every stratum. For example, the `CD8_Cytotoxic`
signature scores near-zero on all cells in a malignant-only stratum (because there are
no T-cells there). Keeping such signatures would add noise and inflate multiple-testing
burden.

**Filter rule:** For each signature, compute its maximum UCell score across all cells
in that stratum. If `max < 0.2`, the signature is dropped (meaning no cell in this
stratum expresses that program at even a moderate level). This is conservative — 0.2
on UCell's 0–1 scale typically means "at least some genes in the signature are
detectably co-expressed."

**Typical retention:** 33–40 out of 40 signatures pass per stratum. NonMalignant
strata retain all 40 (they contain all cell types). FP_Stem retains only 33 (immune
signatures drop out — consistent with FP-RMS being "immune cold").

### Step 3: Identify "High-Scoring" Cells

For each retained signature, we define "high-scoring" cells as the **top 10%** — cells
whose UCell score for that signature is in the 90th percentile or above within that
stratum. This is:
- Robust to signature-to-signature scale differences (each gets its own threshold)
- Biologically meaningful: top 10% captures the cells most strongly expressing the
  program, regardless of baseline expression level
- Not arbitrary: 10th percentile is a standard cutoff in scRNA-seq enrichment analyses
  (used in Tirosh 2016, Kinker 2020, Gavish 2023)

### Step 4: Fisher's Exact Test (Cluster Enrichment)

For each combination of (signature × cluster), we build a 2×2 contingency table:

|                    | In this cluster | Not in this cluster |
|--------------------|:-:|:-:|
| **Top 10% cells**  | a | b |
| **Other 90% cells**| c | d |

Then compute a **one-sided Fisher's exact test** (alternative = "greater"), asking:
"Are top-scoring cells for this signature over-represented in this cluster compared
to what you'd expect by chance?"

The **Odds Ratio (OR)** quantifies how much more likely a top-10% cell is to be in
this cluster versus elsewhere:
- OR = 1: no enrichment (top cells are evenly distributed)
- OR = 5: 5× more likely to find a top-scoring cell in this cluster than expected
- OR = 100+: extreme enrichment (almost all top cells concentrate in one cluster)

### Step 5: Multiple-Testing Correction

With ~40 signatures × ~15 clusters per stratum = ~600 tests, we apply
Benjamini-Hochberg FDR correction. A cluster-signature pair is considered significant
if `FDR < 0.01` (strict threshold to avoid false annotations).

### Step 6: Cluster Assignment

Each cluster is assigned the signature with the **highest OR** among all its
significant associations (OR ≥ 2, FDR < 0.01). The OR value is shown on the UMAP
label and in the heatmap.

---

## Why Are Some Odds Ratios Extremely High (100–460)?

This is not a bug — it reflects the biology of how scRNA-seq clusters form.

### Example: FN_Diff, Cluster 5, Pericyte_vSMC signature
- **282 of 305** cluster-5 cells (92.5%) are in the top-10% for Pericyte_vSMC
- Only **87 of 3,374** cells outside cluster 5 are top-10% Pericyte_vSMC
- This produces OR = 460

### Why this happens:

**Seurat clusters ARE cell-type-defined groups.** When FindClusters identifies
cluster 5 based on PCA/Harmony space, it's finding a group of cells that share
similar overall transcriptomes. If those cells happen to be pericytes, then nearly
ALL pericyte-signature genes are co-expressed in that cluster — and since UCell
also detects pericyte genes, the overlap between "cluster 5" and "pericyte-high
cells" is near-complete.

In other words: the cluster was already defined by the same biology the signature
measures. They're not independent — the high OR reflects a true, near-1:1 mapping
between the Seurat cluster and the signature.

**High ORs (50–500) mean:** One cluster almost perfectly captures one cell type.
This is a GOOD result — it means your clustering resolution is appropriate (not
over-splitting or under-splitting), and the signature correctly identifies that
population.

**Moderate ORs (2–20) mean:** The signature is enriched in a cluster but doesn't
fully dominate it. This happens when:
- A cluster contains a mix of cell types (e.g., cluster has both M1 and M2 macrophages)
- A signature marks a cell STATE rather than a cell TYPE (e.g., Hypoxia, Proliferation)
  that can be active in multiple clusters

**OR = 1 (not significant):** The signature has nothing to do with that cluster.

### Mathematical intuition

OR explodes when two things are true simultaneously:
1. Most top-10% cells for a signature fall in ONE cluster (high `a`)
2. That cluster contains very few NON-top cells for that signature (low `c`)

When 92% of cluster-5 is top-10% Pericyte AND 97% of non-cluster-5 is NOT top-10%
Pericyte, the 2×2 table becomes extremely unbalanced → OR > 100. This is the
Fisher test working correctly, not an artifact.

---

## Outputs

### Per stratum (8 strata × 3 outputs = 24 files):

| Output | Location | Description |
|--------|----------|-------------|
| Multi-sig UMAP | `figures/qc/P2_multisig_umap_<stratum>.png` | UMAP colored by each cluster's dominant signature. Cluster centroids labeled with `<sig_name> (OR=X.X)` |
| Enrichment heatmap | `figures/qc/P2_multisig_heatmap_<stratum>.png` | Clusters (rows) × signatures (cols), colored by log2(OR). `*` marks significant cells (FDR<0.01, OR≥2) |
| Full results table | `tables/qc/P2_multisig_enrichment_<stratum>.csv` | Complete Fisher's exact results: stratum, signature, cluster, n_top_in_cluster, n_cluster, n_top_total, pct_top_in_cluster, odds_ratio, p_value, fdr |

### Summary statistics

| Column | Meaning |
|--------|---------|
| `n_top_in_cluster` | Number of top-10% cells for this signature that fall in this cluster |
| `n_cluster` | Total cells in this cluster |
| `n_top_total` | Total top-10% cells for this signature (across all clusters) |
| `pct_top_in_cluster` | % of cluster cells that are top-10% for this signature |
| `odds_ratio` | Fisher's exact OR (one-sided, "greater") |
| `p_value` | Raw Fisher p-value |
| `fdr` | BH-adjusted p-value across all tests within a stratum |

---

## Parameters

| Parameter | Value | Justification |
|-----------|-------|---------------|
| `MIN_MAX_SCORE` | 0.2 | Drop signatures where no cell scores above 0.2 (too low to be meaningful) |
| `TOP_PCT` | 0.10 | Top 10% defines "high-scoring" cells — standard in literature |
| `MIN_OR` | 2.0 | Minimum enrichment for assignment (2× expected = meaningful) |
| `MIN_PVAL` | 0.01 | FDR threshold for significance (conservative) |
| Clustering resolution | 0.5 (malignant), 0.8 (nonmalignant) | Matches V5 convention |

---

## Relationship to Other Analyses

This analysis provides **cluster identity** that is referenced by:
- **NMF metaprograms (P3/P4):** Which clusters do MP-high cells concentrate in?
- **Survival (P6):** The protective myeloid/APC signal maps to specific clusters
- **CopyKAT validation (CK01–CK03):** Which clusters in the malignant strata
  are potential misclassification candidates (immune signatures scoring high in
  "malignant" clusters)?
- **Visium (VS05):** Spot classification mirrors this approach but at spatial resolution

---

## Script Location

```
v2/adhoc/P2_multisig_cluster_enrichment.R
```

Dependencies: requires pre-computed score CSVs from the 64 SLURM batch scoring jobs
(`v2/adhoc/umap_signature_scoring/score_and_plot_batch.R`).

---

## Interpretation Guide for Biologists

When you look at the enrichment heatmap:
- **Red columns** = signatures that concentrate strongly in one or a few clusters
  (e.g., Pan_Immune, Pericyte, Endothelial). These are "cell-type signatures" —
  they mark discrete populations.
- **Uniform/blue columns** = signatures expressed at similar levels across many
  clusters (e.g., Hypoxia, Proliferative_G2M). These are "cell-state signatures" —
  they mark transient activities, not permanent identities.
- **A cluster with OR>50 for a single signature** = that cluster IS that cell type.
  Near-perfect overlap.
- **A cluster with OR 2–10 for multiple signatures** = a mixed or transitional
  population (e.g., a cluster that's partially myoblast, partially mesenchymal).

The dominant-signature UMAP gives you a bird's-eye view: each cluster is colored
by its strongest association. The heatmap gives you the full picture: a cluster
might be "dominantly" immune but also enriched for hypoxia and proliferation
simultaneously.
