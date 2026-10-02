# Central path configuration for the complete Blood–MSG pipeline.
#
# CURRENT LAYOUT: all input files are directly under PROJECT_ROOT.
# Change only PROJECT_ROOT if the prefix changes.

PROJECT_ROOT <- "/data/salazarda/data/hackathon_nlm_trnascriptomics/inputs"

# Required input files
BLOOD_COUNTS <- file.path(INPUT_DIR, "RawCountFile_filtered_blood.txt")
MSG_COUNTS <- file.path(INPUT_DIR, "RawCountFile_filtered_msg.txt")
BLOOD_META <- file.path(INPUT_DIR, "sampletable_blood.txt")
MSG_META <- file.path(INPUT_DIR, "sampletable_msg.txt")
BLOOD_GROUPING <- file.path(INPUT_DIR, "sampletable_SjD_grouping_blood.txt")
MSG_GROUPING <- file.path(INPUT_DIR, "sampletable_SjD_grouping_msg.txt")

# Optional historical GCBC DESeq2 outputs
BLOOD_ORIGINAL_DE <- file.path(INPUT_DIR, "D_DE_neg-pos_blood.txt")
MSG_ORIGINAL_DE <- file.path(INPUT_DIR, "D_DE_neg-pos_msg.txt")

# Stage-specific output directories and hand-off files
PAIRING_OUT <- file.path(OUTPUT_DIR, "00_pairing")
DIFFERENTIAL_EXPRESSION_OUT <- file.path(OUTPUT_DIR, "01_differential_expression")
WGCNA_OUT <- file.path(OUTPUT_DIR, "02_wgcna_endotypes_projection")
ENRICHMENT_OUT <- file.path(OUTPUT_DIR, "03_functional_enrichment")

PAIRING_BUNDLE <- file.path(PAIRING_OUT, "pairing_bundle.rds")
PAIRED_MANIFEST <- file.path(PAIRING_OUT, "paired_manifest.tsv")

# Backward-compatible alias used internally by the Stage 01 script.
HACKATHON_OUT <- DIFFERENTIAL_EXPRESSION_OUT
