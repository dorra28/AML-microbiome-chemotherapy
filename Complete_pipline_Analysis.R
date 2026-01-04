# ==============================================================================
# AML Gut Microbiome Analysis: Pre- vs Post-Chemotherapy
# Data source: NCBI SRA (SRP141394) - Rashidi et al. (2022)
# ==============================================================================

# Load required libraries
library(dada2)        # v1.18
library(phyloseq)     # v1.34
library(ggplot2)      # v3.3.5 / v3.4.2
library(vegan)        # v2.5-7
library(DESeq2)       # v1.38.0
library(pheatmap)     # v1.0.12
library(readr)        # v2.1.4
library(dplyr)
library(tidyr)

# Set working directory
setwd("~/AML_microbiome_project")

# ==============================================================================
# 1. DATA RETRIEVAL AND PREPROCESSING
# ==============================================================================

# Define paths
raw_data_path <- "raw_fastq"
filtered_path <- "filtered_fastq"
output_path <- "results"

# Create directories if they don't exist
dir.create(filtered_path, showWarnings = FALSE)
dir.create(output_path, showWarnings = FALSE)

# List fastq files (assuming paired-end sequencing)
fnFs <- sort(list.files(raw_data_path, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(raw_data_path, pattern = "_R2_001.fastq.gz", full.names = TRUE))

# Extract sample names
sample_names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Define filtered file names
filtFs <- file.path(filtered_path, paste0(sample_names, "_F_filt.fastq.gz"))
filtRs <- file.path(filtered_path, paste0(sample_names, "_R_filt.fastq.gz"))

# ==============================================================================
# 2. DADA2 PIPELINE
# ==============================================================================

# Quality filtering and trimming
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen = c(240, 200),
                     maxN = 0,
                     maxEE = c(2, 2),
                     truncQ = 2,
                     rm.phix = TRUE,
                     compress = TRUE,
                     multithread = TRUE)

# Save read retention statistics
write.csv(out, file.path(output_path, "reads_filtering_stats.csv"))

# Learn error rates
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

# Plot error rates (optional)
pdf(file.path(output_path, "error_rates.pdf"))
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)
dev.off()

# Dereplicate reads
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
names(derepFs) <- sample_names
names(derepRs) <- sample_names

# Sample inference
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

# Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Remove chimeras
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", 
                                     multithread = TRUE, verbose = TRUE)

# Track reads through pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, 
               sapply(dadaFs, getN), 
               sapply(dadaRs, getN), 
               sapply(mergers, getN), 
               rowSums(seqtab_nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", 
                     "merged", "nonchim")
rownames(track) <- sample_names
write.csv(track, file.path(output_path, "track_reads.csv"))

# Create Supplementary Fig. S3 - Read retention plot
track_df <- as.data.frame(track)
track_df$sample <- rownames(track_df)
track_long <- tidyr::pivot_longer(track_df, 
                                   cols = -sample, 
                                   names_to = "step", 
                                   values_to = "reads")

pdf(file.path(output_path, "Supp_Fig_S3_read_retention.pdf"), width = 12, height = 6)
ggplot(track_long, aes(x = sample, y = reads, fill = step)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Read Retention Across DADA2 Pipeline",
       x = "Sample", y = "Number of Reads") +
  scale_fill_brewer(palette = "Set2")
dev.off()

# ==============================================================================
# 3. TAXONOMIC ASSIGNMENT
# ==============================================================================

# Assign taxonomy using SILVA v138
taxa <- assignTaxonomy(seqtab_nochim, 
                       "silva_nr99_v138.1_train_set.fa.gz", 
                       multithread = TRUE)

# Add species-level annotation (optional)
taxa <- addSpecies(taxa, "silva_species_assignment_v138.1.fa.gz")

# ==============================================================================
# 4. CREATE PHYLOSEQ OBJECT
# ==============================================================================

# Load or create metadata
# Metadata should include: Sample ID, Group (AML or CHEM), Day, etc.
metadata <- read.csv("metadata.csv", row.names = 1)

# Ensure sample names match
metadata <- metadata[rownames(seqtab_nochim), ]

# Create phyloseq object
ps <- phyloseq(otu_table(seqtab_nochim, taxa_are_rows = FALSE),
               sample_data(metadata),
               tax_table(taxa))

# Save phyloseq object
saveRDS(ps, file.path(output_path, "phyloseq_object.rds"))

# ==============================================================================
# 5. ALPHA DIVERSITY ANALYSIS
# ==============================================================================

# Calculate alpha diversity metrics
alpha_div <- estimate_richness(ps, measures = c("Shannon", "Simpson"))
alpha_div$Group <- sample_data(ps)$Group
alpha_div$Sample <- rownames(alpha_div)

# Save alpha diversity results
write.csv(alpha_div, file.path(output_path, "alpha_diversity.csv"))

# Statistical testing (Wilcoxon rank-sum test)
shannon_test <- wilcox.test(Shannon ~ Group, data = alpha_div)
simpson_test <- wilcox.test(Simpson ~ Group, data = alpha_div)

# Print summary statistics
cat("\n=== Alpha Diversity Summary ===\n")
cat("\nShannon Index:\n")
print(aggregate(Shannon ~ Group, data = alpha_div, 
                FUN = function(x) c(mean = mean(x), sd = sd(x), median = median(x))))
cat("\nSimpson Index:\n")
print(aggregate(Simpson ~ Group, data = alpha_div, 
                FUN = function(x) c(mean = mean(x), sd = sd(x), median = median(x))))

cat("\n=== Statistical Tests ===\n")
cat("Shannon p-value:", shannon_test$p.value, "\n")
cat("Simpson p-value:", simpson_test$p.value, "\n")

# Figure 1: Alpha diversity boxplots
pdf(file.path(output_path, "Figure_1_alpha_diversity.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))

# Shannon plot
p1 <- ggplot(alpha_div, aes(x = Group, y = Shannon, fill = Group)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_classic() +
  labs(title = "Shannon Index", 
       x = "", 
       y = "Shannon Diversity Index") +
  annotate("text", x = 1.5, y = max(alpha_div$Shannon), 
           label = paste0("P = ", round(shannon_test$p.value, 4))) +
  scale_fill_manual(values = c("AML" = "#4CAF50", "CHEM" = "#FF9800"))

# Simpson plot
p2 <- ggplot(alpha_div, aes(x = Group, y = Simpson, fill = Group)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_classic() +
  labs(title = "Simpson Index", 
       x = "", 
       y = "Simpson Diversity Index") +
  annotate("text", x = 1.5, y = max(alpha_div$Simpson), 
           label = paste0("P = ", round(simpson_test$p.value, 4))) +
  scale_fill_manual(values = c("AML" = "#4CAF50", "CHEM" = "#FF9800"))

gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()

# ==============================================================================
# 6. BETA DIVERSITY ANALYSIS
# ==============================================================================

# Calculate Bray-Curtis dissimilarity
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
bray_dist <- phyloseq::distance(ps_rel, method = "bray")

# Ordination (PCoA)
ord <- ordinate(ps_rel, method = "PCoA", distance = bray_dist)

# PERMANOVA test
metadata_for_adonis <- data.frame(sample_data(ps))
permanova_result <- adonis2(bray_dist ~ Group, 
                            data = metadata_for_adonis, 
                            permutations = 999)

print(permanova_result)
write.csv(as.data.frame(permanova_result), 
          file.path(output_path, "permanova_results.csv"))

# Figure 2: PCoA plot
pdf(file.path(output_path, "Figure_2_beta_diversity.pdf"), width = 8, height = 6)
plot_ordination(ps_rel, ord, color = "Group", shape = "Group") +
  geom_point(size = 4, alpha = 0.7) +
  theme_classic() +
  labs(title = "Beta Diversity Analysis (PCoA - Bray-Curtis)",
       subtitle = sprintf("PERMANOVA: R² = %.4f, F = %.2f, P = %.3f", 
                         permanova_result$R2[1], 
                         permanova_result$F[1], 
                         permanova_result$`Pr(>F)`[1])) +
  scale_color_manual(values = c("AML" = "#4CAF50", "CHEM" = "#FF9800")) +
  stat_ellipse(aes(group = Group), linetype = 2)
dev.off()

# ==============================================================================
# 7. TAXONOMIC COMPOSITION ANALYSIS
# ==============================================================================

# Agglomerate to genus level
ps_genus <- tax_glom(ps, taxrank = "Genus")
ps_genus_rel <- transform_sample_counts(ps_genus, function(x) x / sum(x))

# Get top 10 most abundant genera
top10_genera <- names(sort(taxa_sums(ps_genus_rel), decreasing = TRUE)[1:10])
ps_top10 <- prune_taxa(top10_genera, ps_genus_rel)

# Figure 3: Stacked bar plots per sample
pdf(file.path(output_path, "Figure_3_sample_composition.pdf"), width = 14, height = 6)
plot_bar(ps_top10, fill = "Genus") +
  facet_wrap(~Group, scales = "free_x") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        legend.position = "right") +
  labs(title = "Relative Abundance of Top 10 Genera per Sample",
       x = "Sample", y = "Relative Abundance") +
  scale_fill_brewer(palette = "Paired")
dev.off()

# Figure 4: Mean composition by group
genus_means <- ps_top10 %>%
  psmelt() %>%
  group_by(Group, Genus) %>%
  summarise(MeanAbundance = mean(Abundance), .groups = "drop")

pdf(file.path(output_path, "Figure_4_group_composition.pdf"), width = 10, height = 6)
ggplot(genus_means, aes(x = Group, y = MeanAbundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_classic() +
  labs(title = "Average Genus-Level Composition Across Groups",
       x = "Group", y = "Mean Relative Abundance") +
  scale_fill_brewer(palette = "Paired")
dev.off()

# ==============================================================================
# 8. DIFFERENTIAL ABUNDANCE ANALYSIS (DESeq2)
# ==============================================================================

# Convert to DESeq2 format
ps_genus_deseq <- phyloseq_to_deseq2(ps_genus, ~ Group)

# Run DESeq2
dds <- DESeq(ps_genus_deseq, test = "Wald", fitType = "parametric")

# Extract results
res <- results(dds, contrast = c("Group", "CHEM", "AML"), 
               cooksCutoff = FALSE)
res_df <- as.data.frame(res)
res_df$Genus <- tax_table(ps_genus)[rownames(res_df), "Genus"]
res_df$OTU <- rownames(res_df)

# Add significance category
res_df$Significance <- "Not Significant"
res_df$Significance[res_df$padj < 0.05 & res_df$log2FoldChange > 0] <- "CHEM-enriched"
res_df$Significance[res_df$padj < 0.05 & res_df$log2FoldChange < 0] <- "AML-enriched"

# Save results
write.csv(res_df, file.path(output_path, "deseq2_results.csv"))

# Figure 5: Volcano plot
pdf(file.path(output_path, "Figure_5_volcano_plot.pdf"), width = 10, height = 8)
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = Significance), alpha = 0.6, size = 3) +
  geom_text(data = subset(res_df, padj < 0.05), 
            aes(label = Genus), 
            vjust = -0.5, size = 3, check_overlap = TRUE) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  theme_classic() +
  labs(title = "Differential Abundance: AML vs CHEM",
       x = "Log2 Fold Change (CHEM vs AML)",
       y = "-Log10(Adjusted P-value)") +
  scale_color_manual(values = c("AML-enriched" = "#4CAF50", 
                                "CHEM-enriched" = "#FF9800",
                                "Not Significant" = "gray70"))
dev.off()

# ==============================================================================
# 9. FUNCTIONAL PREDICTION SETUP (PICRUSt2)
# ==============================================================================

# Export ASV table for PICRUSt2
# Note: PICRUSt2 runs in command line/conda environment

# Export sequence table
seqtab_export <- t(seqtab_nochim)
write.table(seqtab_export, 
            file.path(output_path, "asv_table_for_picrust2.txt"), 
            sep = "\t", quote = FALSE, col.names = NA)

# Export representative sequences
asv_seqs <- colnames(seqtab_nochim)
asv_headers <- paste0(">ASV", seq_along(asv_seqs))
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, file.path(output_path, "rep_seqs_for_picrust2.fasta"))

cat("\n=== PICRUSt2 Instructions ===\n")
cat("Run the following commands in your conda environment with PICRUSt2 installed:\n\n")
cat("conda activate picrust2\n")
cat("cd", output_path, "\n")
cat("picrust2_pipeline.py -s rep_seqs_for_picrust2.fasta \\\n")
cat("  -i asv_table_for_picrust2.txt \\\n")
cat("  -o picrust2_output \\\n")
cat("  -p 4 --verbose\n\n")

# ==============================================================================
# 10. FUNCTIONAL ANALYSIS (After running PICRUSt2)
# ==============================================================================

# This section assumes PICRUSt2 has been run and outputs are available
# Load PICRUSt2 outputs (uncomment when files are available)

# picrust_ko <- read_tsv(file.path(output_path, "picrust2_output/KO_metagenome_out/pred_metagenome_unstrat.tsv"))
# picrust_pathways <- read_tsv(file.path(output_path, "picrust2_output/pathways_out/path_abun_unstrat.tsv"))

# # Process and analyze functional data
# # Apply CLR transformation for PCA
# library(compositions)
# 
# ko_matrix <- as.matrix(picrust_ko[, -1])
# rownames(ko_matrix) <- picrust_ko$function.
# ko_clr <- clr(t(ko_matrix + 1))
# 
# # PCA on functional profiles
# pca_result <- prcomp(ko_clr, scale. = TRUE)
# 
# # Figure 6A: PCA of functional profiles
# pca_df <- as.data.frame(pca_result$x)
# pca_df$Group <- metadata$Group[match(rownames(pca_df), rownames(metadata))]
# 
# pdf(file.path(output_path, "Figure_6A_functional_PCA.pdf"), width = 8, height = 6)
# ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
#   geom_point(size = 4, alpha = 0.7) +
#   theme_classic() +
#   labs(title = "PCA of Predicted Functional Profiles",
#        x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
#        y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)")) +
#   scale_color_manual(values = c("AML" = "#4CAF50", "CHEM" = "#FF9800")) +
#   stat_ellipse()
# dev.off()

cat("\n=== Analysis Complete ===\n")
cat("All results saved to:", output_path, "\n")
