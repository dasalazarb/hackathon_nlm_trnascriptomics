# config_paths.R
#
# CURRENT LAYOUT: all input files are directly under PROJECT_ROOT.
# Change only PROJECT_ROOT if the prefix changes.

PROJECT_ROOT <- "/data/salazarda/data/hackathon_nlm_trnascriptomics/inputs"

# Required input files
BLOOD_COUNTS <- file.path(PROJECT_ROOT, "RawCountFile_filtered_blood.txt")
MSG_COUNTS   <- file.path(PROJECT_ROOT, "RawCountFile_filtered_msg.txt")

BLOOD_META <- file.path(PROJECT_ROOT, "sampletable_blood.txt")
MSG_META   <- file.path(PROJECT_ROOT, "sampletable_msg.txt")

BLOOD_GROUPING <- file.path(PROJECT_ROOT, "sampletable_SjD_grouping_blood.txt")
MSG_GROUPING   <- file.path(PROJECT_ROOT, "sampletable_SjD_grouping_msg.txt")

# Optional historical GCBC DESeq2 outputs
BLOOD_ORIGINAL_DE <- file.path(PROJECT_ROOT, "D_DE_neg-pos_blood.txt")
MSG_ORIGINAL_DE   <- file.path(PROJECT_ROOT, "D_DE_neg-pos_msg.txt")

# Outputs
OUTPUT_DIR <- file.path(PROJECT_ROOT, "output")
PAIRING_OUT <- file.path(OUTPUT_DIR, "00_pairing")
HACKATHON_OUT <- file.path(OUTPUT_DIR, "01_hackathon")

PAIRING_BUNDLE <- file.path(PAIRING_OUT, "pairing_bundle.rds")
