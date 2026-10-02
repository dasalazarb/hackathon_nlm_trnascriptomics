# Central path configuration for the complete Blood–MSG pipeline.
#
# Inputs and outputs intentionally use independent directories. This prevents
# generated results from being written under the immutable input directory.
# Every value can be overridden with an environment variable, for example:
#   PROJECT_ROOT=/data/project INPUT_DIR=/data/raw OUTPUT_DIR=/data/results \
#     Rscript src/00_validate_and_pair_samples.R

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", unset = "/mnt/file-systems")
INPUT_DIR <- Sys.getenv("INPUT_DIR", unset = file.path(PROJECT_ROOT, "inputs"))
OUTPUT_DIR <- Sys.getenv("OUTPUT_DIR", unset = file.path(PROJECT_ROOT, "output"))

# Fail early rather than contaminating the input tree with generated files.
clean_path <- function(path) sub("/+$", "", path.expand(path))
if (
  identical(clean_path(INPUT_DIR), clean_path(OUTPUT_DIR)) ||
  startsWith(clean_path(OUTPUT_DIR), paste0(clean_path(INPUT_DIR), "/"))
) {
  stop(
    "OUTPUT_DIR must be separate from, and not nested under, INPUT_DIR.\n",
    "INPUT_DIR: ", INPUT_DIR, "\nOUTPUT_DIR: ", OUTPUT_DIR
  )
}

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
