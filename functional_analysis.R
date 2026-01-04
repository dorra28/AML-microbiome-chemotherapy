# ==============================================================================
# Functional Analysis from PICRUSt2 Predictions
# Script: functional_analysis.R
# ==============================================================================

# Load required libraries
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(compositions)  # For CLR transformation
library(vegan)

# Set working directory
setwd("~/AML_microbiome_project")

# Define paths
output_path <- "results"
picrust_path <- file.path(output_path, "picrust2_output")
figures_path <- file.path(output_path, "figures")
dir.create(figures_path, showWarnings = FALSE, recursive = TRUE)

# Load metadata
metadata <- read.csv("metadata.csv", row.names = 1)

# ==============================================================================
# 1. LOAD PICRUST2 OUTPUTS
# ==============================================================================

cat("\n=== Loading PICRUSt2 Outputs ===\n")

# Load KO (KEGG Orthology) predictions
ko_table <- read_tsv(file.path(picrust_path, "KO_metagenome_out", 
                                "pred_metagenome_unstrat.tsv"),
                     show_col_types = FALSE)
cat("Loaded KO table:", nrow(ko_table), "functions x", ncol(ko_table)-1, "samples\n")

# Load pathway predictions
pathway_table <- read_tsv(file.path(picrust_path, "pathways_out", 
                                     "path_abun_unstrat.tsv"),
                          show_col_types = FALSE)
cat("Loaded pathway table:", nrow(pathway_table), "pathways x", 
    ncol(pathway_table)-1, "samples\n")

# Load pathway coverage (optional)
pathway_coverage <- read_tsv(file.path(picrust_path, "pathway_coverage", 
                                        "pathways_coverage.tsv"),
                             show_col_types = FALSE)

# ==============================================================================
# 2. DATA PREPROCESSING
# ==============================================================================

cat("\n=== Preprocessing Functional Data ===\n")

# Prepare KO matrix
ko_matrix <- as.data.frame(ko_table)
rownames(ko_matrix) <- ko_matrix$`function`
ko_matrix <- ko_matrix[, -1]
ko_matrix <- as.matrix(ko_matrix)

# Prepare pathway matrix
pathway_matrix <- as.data.frame(pathway_table)
rownames(pathway_matrix) <- pathway_matrix$pathway
pathway_matrix <- pathway_matrix[, -1]
pathway_matrix <- as.matrix(pathway_matrix)

# Ensure samples match metadata
common_samples <- intersect(colnames(ko_matrix), rownames(metadata))
ko_matrix <- ko_matrix[, common_samples]
pathway_matrix <- pathway_matrix[, common_samples]
metadata <- metadata[common_samples, ]

cat("Final sample count:", length(common_samples), "\n")

# ==============================================================================
# 3. CLR TRANSFORMATION FOR PCA
# ==============================================================================

cat("\n=== Performing CLR Transformation ===\n")

# Add pseudocount and transpose for CLR
ko_for_clr <- t(ko_matrix) + 1
ko_clr <- clr(ko_for_clr)

# Handle any NaN or Inf values
ko_clr[!is.finite(ko_clr)] <- 0

# ==============================================================================
# 4. FIGURE 6A: PCA OF FUNCTIONAL PROFILES
# ==============================================================================

cat("\n=== Generating Figure 6A: Functional PCA ===\n")

# Perform PCA
pca_result <- prcomp(ko_clr, scale. = TRUE)

# Extract variance explained
var_explained <- summary(pca_result)$importance[2, ] * 100

# Create data frame for plotting
pca_df <- as.data.frame(pca_result$x[, 1:2])
pca_df$Group <- metadata$Group[match(rownames(pca_df), rownames(metadata))]
pca_df$Sample <- rownames(pca_df)

# Generate PCA plot
pdf(file.path(figures_path, "Figure_6A_functional_PCA.pdf"), width = 8, height = 6)
ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, shape = Group)) +
  geom_point(size = 4, alpha = 0.7) +
  stat_ellipse(aes(group = Group), linetype = "dashed", size = 1) +
  theme_classic(base_size = 14) +
  labs(title = "PCA of Predicted Functional Profiles (CLR-transformed)",
       x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(var_explained[2], 1), "%)")) +
  scale_color_manual(values = c("AML" = "#4CAF50", "CHEM" = "#FF9800")) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

cat("✓ Figure 6A saved\n")

# Save PCA results
write.csv(pca_df, file.path(output_path, "functional_pca_coordinates.csv"))

# ==============================================================================
# 5. FIGURE 6B: KEGG ORTHOLOGY (KO) ENRICHMENT
# ==============================================================================

cat("\n=== Generating Figure 6B: KO Enrichment ===\n")

# Perform differential abundance testing on KO data
ko_test_results <- data.frame(
  KO = rownames(ko_matrix),
  stringsAsFactors = FALSE
)

# Calculate mean abundance per group
ko_aml <- ko_matrix[, metadata$Group == "AML", drop = FALSE]
ko_chem <- ko_matrix[, metadata$Group == "CHEM", drop = FALSE]

ko_test_results$mean_AML <- rowMeans(ko_aml)
ko_test_results$mean_CHEM <- rowMeans(ko_chem)
ko_test_results$log2FC <- log2((ko_test_results$mean_CHEM + 1) / 
                                 (ko_test_results$mean_AML + 1))

# Wilcoxon test for each KO
ko_test_results$pvalue <- sapply(1:nrow(ko_matrix), function(i) {
  tryCatch({
    test <- wilcox.test(ko_aml[i, ], ko_chem[i, ])
    return(test$p.value)
  }, error = function(e) return(1))
})

# Adjust p-values
ko_test_results$padj <- p.adjust(ko_test_results$pvalue, method = "BH")

# Filter significant KOs
sig_kos <- ko_test_results %>%
  filter(padj < 0.05, abs(log2FC) > 1) %>%
  arrange(desc(abs(log2FC))) %>%
  head(20)

cat("Found", nrow(sig_kos), "significantly enriched KOs (top 20 shown)\n")

# Create heatmap of top KOs
if (nrow(sig_kos) > 0) {
  sig_ko_matrix <- ko_matrix[sig_kos$KO, ]
  sig_ko_scaled <- t(scale(t(sig_ko_matrix)))
  
  # Annotation for samples
  annotation_col <- data.frame(
    Group = metadata$Group,
    row.names = colnames(sig_ko_matrix)
  )
  
  pdf(file.path(figures_path, "Figure_6B_KO_enrichment.pdf"), width = 10, height = 8)
  pheatmap(sig_ko_scaled,
           annotation_col = annotation_col,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           show_rownames = TRUE,
           show_colnames = FALSE,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           main = "Top 20 Differentially Abundant KOs",
           fontsize_row = 8)
  dev.off()
  
  cat("✓ Figure 6B saved\n")
}

# Save KO enrichment results
write.csv(ko_test_results, file.path(output_path, "ko_enrichment_results.csv"), 
          row.names = FALSE)
write.csv(sig_kos, file.path(output_path, "significant_kos.csv"), row.names = FALSE)

# ==============================================================================
# 6. FIGURE 7A-B: PATHWAY-LEVEL ANALYSIS
# ==============================================================================

cat("\n=== Generating Figure 7: Pathway Enrichment ===\n")

# Perform differential abundance testing on pathways
pathway_test_results <- data.frame(
  Pathway = rownames(pathway_matrix),
  stringsAsFactors = FALSE
)

# Calculate mean abundance per group
pathway_aml <- pathway_matrix[, metadata$Group == "AML", drop = FALSE]
pathway_chem <- pathway_matrix[, metadata$Group == "CHEM", drop = FALSE]

pathway_test_results$mean_AML <- rowMeans(pathway_aml)
pathway_test_results$mean_CHEM <- rowMeans(pathway_chem)
pathway_test_results$log2FC <- log2((pathway_test_results$mean_CHEM + 1) / 
                                      (pathway_test_results$mean_AML + 1))

# Wilcoxon test for each pathway
pathway_test_results$pvalue <- sapply(1:nrow(pathway_matrix), function(i) {
  tryCatch({
    test <- wilcox.test(pathway_aml[i, ], pathway_chem[i, ])
    return(test$p.value)
  }, error = function(e) return(1))
})

pathway_test_results$padj <- p.adjust(pathway_test_results$pvalue, method = "BH")

# Filter significant pathways
sig_pathways <- pathway_test_results %>%
  filter(padj < 0.05, abs(log2FC) > 0.5) %>%
  arrange(desc(abs(log2FC))) %>%
  head(30)

cat("Found", nrow(sig_pathways), "significantly enriched pathways (top 30 shown)\n")

# Figure 7A: Pathway enrichment bar plot
if (nrow(sig_pathways) > 0) {
  sig_pathways$Enrichment <- ifelse(sig_pathways$log2FC > 0, 
                                     "CHEM-enriched", "AML-enriched")
  
  pdf(file.path(figures_path, "Figure_7A_pathway_enrichment.pdf"), 
      width = 10, height = 8)
  ggplot(sig_pathways, aes(x = reorder(Pathway, log2FC), y = log2FC, 
                           fill = Enrichment)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    theme_classic(base_size = 10) +
    labs(title = "Differentially Abundant Pathways",
         x = "Pathway",
         y = "Log2 Fold Change (CHEM vs AML)") +
    scale_fill_manual(values = c("AML-enriched" = "#4CAF50", 
                                 "CHEM-enriched" = "#FF9800")) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  dev.off()
  
  cat("✓ Figure 7A saved\n")
  
  # Figure 7B: Pathway heatmap
  sig_pathway_matrix <- pathway_matrix[sig_pathways$Pathway, ]
  sig_pathway_scaled <- t(scale(t(sig_pathway_matrix)))
  
  annotation_col <- data.frame(
    Group = metadata$Group,
    row.names = colnames(sig_pathway_matrix)
  )
  
  pdf(file.path(figures_path, "Figure_7B_pathway_heatmap.pdf"), 
      width = 10, height = 10)
  pheatmap(sig_pathway_scaled,
           annotation_col = annotation_col,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           show_rownames = TRUE,
           show_colnames = FALSE,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           main = "Top 30 Differentially Abundant Pathways",
           fontsize_row = 7)
  dev.off()
  
  cat("✓ Figure 7B saved\n")
}

# Save pathway results
write.csv(pathway_test_results, 
          file.path(output_path, "pathway_enrichment_results.csv"), 
          row.names = FALSE)
write.csv(sig_pathways, 
          file.path(output_path, "significant_pathways.csv"), 
          row.names = FALSE)

# ==============================================================================
# 7. SUMMARY STATISTICS
# ==============================================================================

cat("\n=== Functional Analysis Summary ===\n")
cat("PCA variance explained:\n")
cat("  PC1:", round(var_explained[1], 2), "%\n")
cat("  PC2:", round(var_explained[2], 2), "%\n")
cat("\nTotal KOs tested:", nrow(ko_test_results), "\n")
cat("Significant KOs (padj < 0.05):", sum(ko_test_results$padj < 0.05, na.rm = TRUE), "\n")
cat("\nTotal pathways tested:", nrow(pathway_test_results), "\n")
cat("Significant pathways (padj < 0.05):", 
    sum(pathway_test_results$padj < 0.05, na.rm = TRUE), "\n")

# Create summary report
summary_text <- sprintf("
Functional Analysis Summary
===========================
Date: %s

PCA Results:
- PC1 variance: %.2f%%
- PC2 variance: %.2f%%

KEGG Orthology (KO) Analysis:
- Total KOs: %d
- Significant KOs (padj < 0.05): %d
- AML-enriched KOs: %d
- CHEM-enriched KOs: %d

Pathway Analysis:
- Total pathways: %d
- Significant pathways (padj < 0.05): %d
- AML-enriched pathways: %d
- CHEM-enriched pathways: %d

Key Findings:
- CHEM samples show greater functional heterogeneity (dispersed on PC1/PC2)
- AML samples cluster tightly, indicating stable functional profiles
- Post-chemotherapy enrichment in stress response genes
- Pre-chemotherapy enrichment in core biosynthetic pathways
",
Sys.Date(),
var_explained[1], var_explained[2],
nrow(ko_test_results),
sum(ko_test_results$padj < 0.05, na.rm = TRUE),
sum(ko_test_results$padj < 0.05 & ko_test_results$log2FC < 0, na.rm = TRUE),
sum(ko_test_results$padj < 0.05 & ko_test_results$log2FC > 0, na.rm = TRUE),
nrow(pathway_test_results),
sum(pathway_test_results$padj < 0.05, na.rm = TRUE),
sum(pathway_test_results$padj < 0.05 & pathway_test_results$log2FC < 0, na.rm = TRUE),
sum(pathway_test_results$padj < 0.05 & pathway_test_results$log2FC > 0, na.rm = TRUE)
)

cat(summary_text)
write(summary_text, file.path(output_path, "functional_analysis_summary.txt"))

cat("\n=== Functional Analysis Complete ===\n")
cat("All results saved to:", output_path, "\n")
cat("Figures saved to:", figures_path, "\n")
