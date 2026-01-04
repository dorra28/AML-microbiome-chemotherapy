#!/bin/bash
# ==============================================================================
# PICRUSt2 Functional Prediction Pipeline
# Script: 06_picrust2_analysis.sh
# ==============================================================================

set -e  # Exit on error

echo "=========================================="
echo "PICRUSt2 Functional Analysis Pipeline"
echo "=========================================="

# Define paths
RESULTS_DIR="results"
ASV_TABLE="${RESULTS_DIR}/asv_table_for_picrust2.txt"
REP_SEQS="${RESULTS_DIR}/rep_seqs_for_picrust2.fasta"
OUTPUT_DIR="${RESULTS_DIR}/picrust2_output"
THREADS=4

# Check if required files exist
if [ ! -f "$ASV_TABLE" ]; then
    echo "Error: ASV table not found at $ASV_TABLE"
    echo "Please run the DADA2 pipeline first (01_dada2_pipeline.R)"
    exit 1
fi

if [ ! -f "$REP_SEQS" ]; then
    echo "Error: Representative sequences not found at $REP_SEQS"
    echo "Please run the DADA2 pipeline first (01_dada2_pipeline.R)"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo ""
echo "Step 1: Running PICRUSt2 pipeline..."
echo "Input ASV table: $ASV_TABLE"
echo "Input sequences: $REP_SEQS"
echo "Output directory: $OUTPUT_DIR"
echo "Threads: $THREADS"
echo ""

# Run PICRUSt2 pipeline
picrust2_pipeline.py \
    -s "$REP_SEQS" \
    -i "$ASV_TABLE" \
    -o "$OUTPUT_DIR" \
    -p "$THREADS" \
    --stratified \
    --verbose

# Check if pipeline completed successfully
if [ $? -eq 0 ]; then
    echo ""
    echo "✓ PICRUSt2 pipeline completed successfully!"
else
    echo ""
    echo "✗ PICRUSt2 pipeline failed!"
    exit 1
fi

echo ""
echo "Step 2: Generated outputs:"
echo "  - KO predictions: ${OUTPUT_DIR}/KO_metagenome_out/"
echo "  - EC predictions: ${OUTPUT_DIR}/EC_metagenome_out/"
echo "  - Pathway predictions: ${OUTPUT_DIR}/pathways_out/"
echo "  - Pathway coverage: ${OUTPUT_DIR}/pathway_coverage/"

# Optional: Generate pathway coverage statistics
echo ""
echo "Step 3: Generating pathway coverage statistics..."

if [ -f "${OUTPUT_DIR}/pathways_out/path_abun_unstrat.tsv" ]; then
    echo "  - Unstratified pathway abundance table found"
    head -n 5 "${OUTPUT_DIR}/pathways_out/path_abun_unstrat.tsv"
fi

echo ""
echo "=========================================="
echo "PICRUSt2 Analysis Complete!"
echo "=========================================="
echo ""
echo "Next steps:"
echo "1. Run functional analysis in R: source('scripts/07_functional_analysis.R')"
echo "2. Results will be saved in: ${RESULTS_DIR}/figures/"
echo ""

# Create a summary file
SUMMARY_FILE="${OUTPUT_DIR}/analysis_summary.txt"
cat > "$SUMMARY_FILE" << EOF
PICRUSt2 Analysis Summary
========================
Date: $(date)
ASV Table: $ASV_TABLE
Representative Sequences: $REP_SEQS
Output Directory: $OUTPUT_DIR
Threads Used: $THREADS

Output Files:
- KO metagenome predictions (unstratified): KO_metagenome_out/pred_metagenome_unstrat.tsv
- EC metagenome predictions (unstratified): EC_metagenome_out/pred_metagenome_unstrat.tsv
- Pathway abundance (unstratified): pathways_out/path_abun_unstrat.tsv
- Pathway coverage: pathway_coverage/pathways_coverage.tsv

Status: Completed successfully
EOF

echo "Summary saved to: $SUMMARY_FILE"
