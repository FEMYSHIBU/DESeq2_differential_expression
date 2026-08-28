# DESeq2 : Theoretical Foundations

This document explains the statistical and biological reasoning behind the DESeq2 differential expression pipeline — why each step exists, not just what it does. For the actual analysis, see ['scripts/deseq2_analysis.R']

---

## 1. Why raw counts can't be compared directly

### Library size and composition bias
Two samples sequenced to different depths (library sizes) will naturally show different raw read counts for the same gene, even if the true expression is identical — a sample sequenced twice as deeply will simply have roughly twice the counts. Composition bias adds a second complication: if a handful of genes are extremely highly expressed in one sample, they consume a disproportionate share of the sequencing reads, artificially deflating the apparent counts of every other gene in that sample.

DESeq2 corrects for both using the **median-of-ratios method**: for each gene, compute the ratio of its count in a given sample to a reference value (the geometric mean of that gene's counts across all samples). Take the median of these ratios across all genes for each sample — that median becomes the sample's **size factor**. Dividing raw counts by their sample's size factor puts all samples on a comparable scale, robust to a small number of highly expressed or differentially expressed genes skewing the correction.

### Gene length bias — not applicable here
Gene length normalization (as in FPKM/TPM) matters when comparing *different genes within the same sample*, since longer transcripts naturally accumulate more reads. DESeq2 instead compares the *same gene* across different samples/conditions — the gene's length is constant across that comparison, so it cancels out in the ratio and requires no separate correction.

---

## 2. Why negative binomial, not Poisson

RNA-seq read counts are discrete, non-negative whole numbers (a gene either has 0, 1, 2... reads mapped to it — never a fraction, never negative), which rules out continuous distributions like the normal distribution as a direct model.

The Poisson distribution is the natural first choice for count data, but it makes a strict assumption: **variance equals the mean**. Real RNA-seq data violates this — variance across biological replicates is consistently *larger* than the mean, a phenomenon called **overdispersion**. This extra variability comes from genuine biological variation between replicates (different individuals, cells, or conditions), not just sequencing noise.

The **negative binomial distribution** extends the Poisson model with an extra parameter — the **dispersion parameter** — that explicitly captures this excess variance, making it the appropriate model for RNA-seq counts.

---

## 3. Filtering low-count genes

Before modeling, genes with zero (or near-zero) counts across all samples are removed. These genes carry no meaningful signal, contribute disproportionate noise to dispersion estimation, and inflate the multiple-testing burden unnecessarily (see FDR correction, below) — removing them improves both statistical power and computational efficiency.

---

## 4. Dispersion estimation and shrinkage

DESeq2 estimates a **dispersion parameter for each gene individually**, describing how much that gene's counts vary across replicates beyond what the mean alone would predict. 

With typically small sample sizes (few biological replicates per condition), per-gene dispersion estimates are inherently noisy and unreliable — genes with low mean counts in particular show artificially inflated or unstable variance estimates purely by chance. To correct this, DESeq2 applies **dispersion shrinkage**: individual gene-wise estimates are pulled ("shrunk") toward a fitted trend line that captures the typical dispersion-versus-mean relationship across all genes. Genes with similar expression levels share information with each other, giving more stable estimates. This reduces false positives caused by unreliable, small-sample variance.

---

## 5. Hypothesis testing
 
With a fitted negative binomial model and shrunk dispersion estimates in hand, DESeq2 applies a **Wald test** to each gene: it estimates the standard error of the log2 fold change and computes a p-value testing whether that fold change is significantly different from zero.
 
Testing thousands of genes simultaneously inflates the chance of false positives from multiple testing alone. DESeq2 corrects for this using the **Benjamini-Hochberg (BH) procedure**, producing an **adjusted p-value (padj / FDR)** for each gene. The FDR accounts for the total number of genes tested, so genes surviving this correction (typically padj < 0.05) are the ones considered truly significant.

---

## 6. Fold change filtering
 
Statistical significance (a low padj) doesn't guarantee *biological* relevance — a gene can be statistically significant with a trivially small expression change. A common practical filter is:
 
- **padj < 0.05** — statistically significant after FDR correction
- **|log2FoldChange| > 1** — at least a 2-fold change in expression

This log2FC cutoff is applied as a **post-hoc filter** on top of the Wald test results, not as part of the test itself.

---

## Summary of the pipeline logic

| Step | Problem it solves |
|---|---|
| Median-of-ratios normalization | Library size / composition bias between samples |
| Negative binomial model | Overdispersion in biological replicate counts |
| Low-count filtering | Removes uninformative, noise-prone genes |
| Dispersion shrinkage | Stabilizes noisy per-gene variance estimates (few replicates) |
| Wald test | Tests significance of fold change per gene |
| BH / FDR correction | Controls false positives across thousands of simultaneous tests |
| log2FC filtering / shrinkage | Separates statistical significance from biological relevance |

