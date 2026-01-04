# Impact of Chemotherapy on Gut Microbiome in AML Patients

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This repository contains the complete computational pipeline for analyzing gut microbiome changes in acute myeloid leukemia (AML) patients before and after chemotherapy treatment. The analysis is based on 16S rRNA sequencing data from 54 stool samples (22 pre-chemotherapy, 32 post-chemotherapy).

## Citation

If you use this code, please cite:
 Dorra Rjaibi , Mohit Batra ,  Wee Wei Yee (2025). AML Gut Microbiome Analysis Pipeline.
GitHub: https://github.com/dorra28/AML-microbiome-chemotherapy/
DOI: 10.5281/zenodo.XXXXXXX
Original data source:
- Rashidi A, et al. (2022). NCBI SRA: SRP141394

## Key Findings

- **Alpha Diversity**: Significant decrease post-chemotherapy (Shannon: 5.66→5.17, P=0.006)
- **Beta Diversity**: PERMANOVA revealed 5.58% variance explained by treatment (P=0.001)
- **Taxonomic Shifts**: Expansion of opportunistic genera (Enterococcus, Streptococcus) post-treatment
- **Functional Changes**: Shift from biosynthetic pathways to stress response mechanisms

## Requirements

### Software
- R (v4.2.0 or higher)
- Conda (v4.14.0 or higher)
- Python (3.9 for PICRUSt2)

### R Packages
```r
install.packages("BiocManager")
BiocManager::install(c(
  "dada2",        # v1.18
  "phyloseq",     # v1.34
  "DESeq2",       # v1.38.0
  "ggplot2",      # v3.4.2
  "vegan",        # v2.5-7
  "pheatmap",     # v1.0.12
  "readr",        # v2.1.4
  "dplyr",
  "tidyr",
  "gridExtra"
))
```

### Python/Conda Environment
```bash
conda env create -f environment.yml
conda activate aml-microbiome
```

## Data Access

Raw sequencing data: [NCBI SRA SRP141394](https://www.ncbi.nlm.nih.gov/sra/?term=SRP141394)

Download using SRA Toolkit:
```bash
prefetch SRP141394
fastq-dump --split-files --gzip SRP141394
```

## Installation
```bash
git clone https://github.com/yourusername/AML-microbiome-chemotherapy.git
cd AML-microbiome-chemotherapy
conda env create -f environment.yml
```

## Usage

### Complete Pipeline
Run the full analysis:
```r
source("scripts/01_dada2_pipeline.R")
source("scripts/02_alpha_diversity.R")
source("scripts/03_beta_diversity.R")
source("scripts/04_taxonomic_composition.R")
source("scripts/05_differential_abundance.R")
```

### Functional Prediction
```bash
bash scripts/06_picrust2_analysis.sh
```
```r
source("scripts/07_functional_analysis.R")
```

## Workflow

1. **Quality Control & ASV Inference** (DADA2)
   - Quality filtering
   - Error learning
   - Denoising
   - Chimera removal

2. **Taxonomic Assignment** (SILVA v138)
   - Genus-level classification
   - Community composition analysis

3. **Diversity Analysis**
   - Alpha diversity (Shannon, Simpson indices)
   - Beta diversity (Bray-Curtis, PERMANOVA)

4. **Differential Abundance** (DESeq2)
   - Genus-level comparisons
   - Volcano plot visualization

5. **Functional Prediction** (PICRUSt2)
   - KEGG ortholog inference
   - Pathway analysis
   - Functional shifts

## Results

All figures and tables are saved in `results/`:
- **Figure 1**: Alpha diversity boxplots
- **Figure 2**: PCoA ordination (beta diversity)
- **Figure 3**: Per-sample taxonomic composition
- **Figure 4**: Mean composition by group
- **Figure 5**: Differential abundance volcano plot
- **Figure 6**: Functional profile PCA
- **Figure 7**: Pathway enrichment

## Contact

- **Author**: Dorra Rjaibi, Mohit Batra , Wee Wei Yee
- **Email**:  https://github.com/dorra28/

## Acknowledgments

- Original data: Rashidi et al. (2022)
- SILVA database maintainers
- DADA2, phyloseq, DESeq2, and PICRUSt2 developers

## References

1. Rashidi A, et al. (2022). [Original paper citation]
2. Callahan BJ, et al. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. Nature Methods.
3. McMurdie PJ, Holmes S (2013). phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data. PLoS ONE.
4. Love MI, et al. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology.
5. Douglas GM, et al. (2020). PICRUSt2 for prediction of metagenome functions. Nature Biotechnology.

##Data Directory

## Raw Data

Raw 16S rRNA sequencing data should be downloaded from NCBI SRA (SRP141394):
```bash
# Install SRA Toolkit
conda install -c bioconda sra-tools

# Download data
prefetch SRP141394
fastq-dump --split-files --gzip --outdir raw_fastq/ SRP141394
```

## Metadata

The `metadata.csv` file contains sample information:

| Column | Description |
|--------|-------------|
| SampleID | Unique sample identifier |
| Group | Treatment group (AML or CHEM) |
| Day | Day of collection (0-1 or 28-31) |
| Patient | Patient identifier |
| AntibioticExposure | Whether patient received antibiotics |

**Note**: Ensure sample IDs in metadata match FASTQ file names.

