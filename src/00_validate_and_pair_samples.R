#!/usr/bin/env Rscript

# 00_validate_and_pair_samples.R
#
# No command-line paths are required.
# Configure PROJECT_ROOT, INPUT_DIR, and OUTPUT_DIR in config_paths.R.
#
# PURPOSE
# -------
# 1) Start from the already generated GCBC count matrices.
# 2) Validate Blood and MSG sample tables.
# 3) Read a SEPARATE sampletable_SjD_grouping.txt for Blood and MSG.
# 4) Pair Blood and MSG at the patient level.
# 5) Require concordant SjD / nonSjD / HV assignment across tissues.
# 6) Create pairing_bundle.rds for the hackathon analysis.
# 7) Run the minimal historical GCBC DESeq2 replication separately in Blood and MSG.
#
# IMPORTANT
# ---------
# sampletable.txt:
#   historical assay condition (neg / pos), used only for replication.
#
# sampletable_SjD_grouping.txt:
#   biological disease group (SjD / nonSjD / HV), used for the hackathon.
#
# The two grouping files are kept separate by tissue and checked against
# one another rather than merged blindly.

suppressPackageStartupMessages({
  library(DESeq2)
})

# -------------------------------------------------------------------------
# 0. Load central path configuration
# -------------------------------------------------------------------------

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0) return(getwd())
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE))
}

find_config <- function() {
  script_dir <- get_script_dir()
  candidates <- unique(c(
    file.path(script_dir, "config_paths.R"),
    file.path(dirname(script_dir), "config_paths.R"),
    file.path(getwd(), "config_paths.R"),
    file.path(dirname(getwd()), "config_paths.R")
  ))
  found <- candidates[file.exists(candidates)]
  if (length(found) == 0) {
    stop(
      "config_paths.R was not found. Checked:\n",
      paste0("  - ", candidates, collapse = "\n")
    )
  }
  found[1]
}

CONFIG_FILE <- find_config()
source(CONFIG_FILE)
cat("Using config:", CONFIG_FILE, "\n")
cat("INPUT_DIR:", INPUT_DIR, "\n")
cat("OUTPUT_DIR:", OUTPUT_DIR, "\n")

dir.create(PAIRING_OUT, recursive = TRUE, showWarnings = FALSE)

required_files <- c(
  BLOOD_COUNTS,
  BLOOD_META,
  BLOOD_GROUPING,
  MSG_COUNTS,
  MSG_META,
  MSG_GROUPING
)

missing_required <- required_files[!file.exists(required_files)]
if (length(missing_required) > 0) {
  stop(
    "Missing required input file(s):\n",
    paste0("  - ", missing_required, collapse = "\n"),
    "\n\nEdit INPUT_DIR in config_paths.R if the input location is incorrect."
  )
}

write_tsv <- function(x, path, row.names = FALSE) {
  write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = row.names,
    col.names = TRUE
  )
}

canonical_sample <- function(x) {
  x <- basename(as.character(x))
  x <- trimws(x)
  x <- sub("\\.star\\..*$", "", x)
  x <- sub("\\.ReadsPerGene.*$", "", x)
  x <- sub("\\.txt$", "", x)
  x
}

extract_patient_id <- function(sample_name) {
  x <- canonical_sample(sample_name)
  x <- sub("^[Bb]lood[_-]?", "", x)
  x <- sub("^[Mm][Ss][Gg][_-]?", "", x)
  x
}

read_counts <- function(path) {
  x <- read.delim(
    path,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  x <- as.matrix(x)
  suppressWarnings(storage.mode(x) <- "numeric")

  if (anyNA(x)) {
    stop("NA/non-numeric values found in count matrix: ", path)
  }

  if (any(x < 0)) {
    stop("Negative values found in count matrix: ", path)
  }

  if (any(abs(x - round(x)) > 1e-8)) {
    stop("Non-integer values found; expected raw counts: ", path)
  }

  storage.mode(x) <- "integer"
  colnames(x) <- canonical_sample(colnames(x))

  if (anyDuplicated(colnames(x))) {
    d <- unique(colnames(x)[duplicated(colnames(x))])
    stop(
      "Duplicated sample columns after canonicalization in ",
      path, ": ", paste(d, collapse = ", ")
    )
  }

  x
}

read_original_meta <- function(path, tissue) {
  m <- read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  required <- c("sampleName", "condition")
  missing <- setdiff(required, colnames(m))

  if (length(missing) > 0) {
    stop(
      "Missing column(s) in ", path, ": ",
      paste(missing, collapse = ", ")
    )
  }

  m$sampleName <- canonical_sample(m$sampleName)
  m$patient_id <- extract_patient_id(m$sampleName)
  m$tissue <- tissue
  m$assay_condition <- trimws(as.character(m$condition))

  if (anyDuplicated(m$sampleName)) {
    d <- unique(m$sampleName[duplicated(m$sampleName)])
    stop(
      "Duplicated sampleName in ", tissue, " sampletable: ",
      paste(d, collapse = ", ")
    )
  }

  if (anyDuplicated(m$patient_id)) {
    d <- unique(m$patient_id[duplicated(m$patient_id)])
    stop(
      "More than one ", tissue,
      " sample maps to the same patient_id: ",
      paste(d, collapse = ", ")
    )
  }

  m
}

read_tissue_grouping <- function(path, tissue) {
  g <- read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  required <- c("sampleName", "condition")
  missing <- setdiff(required, colnames(g))

  if (length(missing) > 0) {
    stop(
      tissue,
      " sampletable_SjD_grouping.txt must contain: ",
      paste(required, collapse = ", "),
      ". Missing: ",
      paste(missing, collapse = ", ")
    )
  }

  g$sampleName <- canonical_sample(g$sampleName)
  g$patient_id <- extract_patient_id(g$sampleName)
  g$disease_group <- trimws(as.character(g$condition))
  g$tissue <- tissue

  allowed <- c("SjD", "nonSjD", "HV")
  bad <- setdiff(unique(g$disease_group), allowed)

  if (length(bad) > 0) {
    stop(
      tissue,
      " grouping contains unexpected disease group(s): ",
      paste(bad, collapse = ", "),
      ". Expected only: SjD, nonSjD, HV"
    )
  }

  if (anyDuplicated(g$sampleName)) {
    d <- unique(g$sampleName[duplicated(g$sampleName)])
    stop(
      "Duplicated sampleName in ", tissue, " grouping file: ",
      paste(d, collapse = ", ")
    )
  }

  if (anyDuplicated(g$patient_id)) {
    dup_ids <- unique(g$patient_id[duplicated(g$patient_id)])

    for (pid in dup_ids) {
      vals <- unique(g$disease_group[g$patient_id == pid])

      if (length(vals) > 1) {
        stop(
          "Conflicting ", tissue,
          " disease groups for patient ", pid, ": ",
          paste(vals, collapse = ", ")
        )
      }
    }

    # If duplicates carry the same group, reduce to one row/patient.
    g <- g[!duplicated(g$patient_id), , drop = FALSE]
  }

  g
}

validate_and_reorder <- function(counts, meta, tissue) {
  missing_meta <- setdiff(colnames(counts), meta$sampleName)
  missing_counts <- setdiff(meta$sampleName, colnames(counts))

  if (length(missing_meta) > 0 || length(missing_counts) > 0) {
    report <- c(
      paste0("SAMPLE-ID MISMATCH: ", tissue),
      ""
    )

    if (length(missing_meta) > 0) {
      report <- c(
        report,
        "In counts but not sampletable.txt:",
        paste0("  ", missing_meta)
      )
    }

    if (length(missing_counts) > 0) {
      report <- c(
        report,
        "In sampletable.txt but not counts:",
        paste0("  ", missing_counts)
      )
    }

    writeLines(
      report,
      file.path(
        PAIRING_OUT,
        paste0(tolower(tissue), "_sample_id_mismatch.txt")
      )
    )

    stop("Count columns and sampletable.txt do not match for ", tissue)
  }

  meta <- meta[
    match(colnames(counts), meta$sampleName),
    ,
    drop = FALSE
  ]

  stopifnot(all(meta$sampleName == colnames(counts)))
  rownames(meta) <- meta$sampleName
  meta
}

attach_grouping_to_tissue <- function(meta, grouping, tissue) {
  z <- merge(
    meta,
    grouping[, c("patient_id", "disease_group"), drop = FALSE],
    by = "patient_id",
    all.x = TRUE,
    sort = FALSE
  )

  # Restore count/sampletable order after merge.
  z <- z[match(meta$patient_id, z$patient_id), , drop = FALSE]

  missing_group <- z[is.na(z$disease_group), , drop = FALSE]

  if (nrow(missing_group) > 0) {
    write_tsv(
      missing_group,
      file.path(
        PAIRING_OUT,
        paste0(tolower(tissue), "_samples_missing_disease_group.tsv")
      )
    )

    stop(
      tissue,
      ": ", nrow(missing_group),
      " sample(s) in sampletable.txt have no SjD/nonSjD/HV assignment."
    )
  }

  z
}

run_original_deseq2 <- function(counts, meta, tissue, out_path) {
  if (!all(c("neg", "pos") %in% unique(meta$assay_condition))) {
    warning(
      tissue,
      ": historical replication skipped because sampletable.txt does not contain both 'neg' and 'pos'. Found: ",
      paste(sort(unique(meta$assay_condition)), collapse = ", ")
    )
    return(NULL)
  }

  coldata <- meta
  rownames(coldata) <- coldata$sampleName
  coldata$condition <- factor(
    coldata$assay_condition,
    levels = c("neg", "pos")
  )

  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = coldata,
    design = ~ condition
  )

  # RawCountFile_filtered.txt already contains the historical GCBC filter.
  # Do not apply a second low-count filter in the replication gate.
  dds <- DESeq(dds, quiet = TRUE)

  res <- results(
    dds,
    contrast = c("condition", "neg", "pos"),
    alpha = 0.05
  )

  tab <- as.data.frame(res)
  tab$ensembl_gene_id <- rownames(tab)

  tab <- tab[, c(
    "ensembl_gene_id",
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj"
  )]

  tab <- tab[order(tab$padj, na.last = TRUE), ]

  write_tsv(tab, out_path)

  invisible(
    list(
      dds = dds,
      results = tab
    )
  )
}

compare_with_original <- function(new_tab, original_path, out_path) {
  if (is.null(new_tab)) return(NULL)

  if (!file.exists(original_path)) {
    return(NULL)
  }

  old <- read.delim(
    original_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (!"ensembl_gene_id" %in% colnames(old)) {
    possible <- grep(
      "ensembl|gene",
      colnames(old),
      ignore.case = TRUE,
      value = TRUE
    )

    if (length(possible) > 0) {
      colnames(old)[match(possible[1], colnames(old))] <- "ensembl_gene_id"
    } else {
      warning(
        "Could not identify ensembl_gene_id in original result: ",
        original_path,
        ". Concordance check skipped."
      )
      return(NULL)
    }
  }

  if (!"log2FoldChange" %in% colnames(old)) {
    warning(
      "Original result lacks log2FoldChange: ",
      original_path,
      ". Concordance check skipped."
    )
    return(NULL)
  }

  keep_old <- c("ensembl_gene_id", "log2FoldChange")
  if ("padj" %in% colnames(old)) {
    keep_old <- c(keep_old, "padj")
  }

  z <- merge(
    new_tab[, c("ensembl_gene_id", "log2FoldChange", "padj")],
    old[, keep_old, drop = FALSE],
    by = "ensembl_gene_id",
    suffixes = c("_new", "_old")
  )

  rho <- suppressWarnings(
    cor(
      z$log2FoldChange_new,
      z$log2FoldChange_old,
      method = "spearman",
      use = "complete.obs"
    )
  )

  sign_agreement <- mean(
    sign(z$log2FoldChange_new) ==
      sign(z$log2FoldChange_old),
    na.rm = TRUE
  )

  z$abs_log2FC_difference <- abs(
    z$log2FoldChange_new -
      z$log2FoldChange_old
  )

  write_tsv(z, out_path)

  list(
    n_genes_compared = nrow(z),
    spearman_log2FC = rho,
    sign_agreement = sign_agreement
  )
}

# -------------------------------------------------------------------------
# 1. Read required inputs
# -------------------------------------------------------------------------

blood_counts <- read_counts(BLOOD_COUNTS)
msg_counts   <- read_counts(MSG_COUNTS)

blood_meta <- read_original_meta(BLOOD_META, "Blood")
msg_meta   <- read_original_meta(MSG_META, "MSG")

blood_grouping <- read_tissue_grouping(
  BLOOD_GROUPING,
  "Blood"
)

msg_grouping <- read_tissue_grouping(
  MSG_GROUPING,
  "MSG"
)

blood_meta <- validate_and_reorder(
  blood_counts,
  blood_meta,
  "Blood"
)

msg_meta <- validate_and_reorder(
  msg_counts,
  msg_meta,
  "MSG"
)

blood_meta <- attach_grouping_to_tissue(
  blood_meta,
  blood_grouping,
  "Blood"
)

msg_meta <- attach_grouping_to_tissue(
  msg_meta,
  msg_grouping,
  "MSG"
)

# Save normalized views of the two grouping inputs for auditing.
write_tsv(
  blood_grouping,
  file.path(PAIRING_OUT, "blood_grouping_normalized.tsv")
)

write_tsv(
  msg_grouping,
  file.path(PAIRING_OUT, "msg_grouping_normalized.tsv")
)

# -------------------------------------------------------------------------
# 2. Build Blood-MSG patient pairing
# -------------------------------------------------------------------------

blood_manifest <- blood_meta[
  ,
  c(
    "patient_id",
    "sampleName",
    "assay_condition",
    "disease_group"
  ),
  drop = FALSE
]

colnames(blood_manifest) <- c(
  "patient_id",
  "blood_sample",
  "blood_assay_condition",
  "blood_disease_group"
)

msg_manifest <- msg_meta[
  ,
  c(
    "patient_id",
    "sampleName",
    "assay_condition",
    "disease_group"
  ),
  drop = FALSE
]

colnames(msg_manifest) <- c(
  "patient_id",
  "msg_sample",
  "msg_assay_condition",
  "msg_disease_group"
)

paired <- merge(
  blood_manifest,
  msg_manifest,
  by = "patient_id",
  all = FALSE,
  sort = TRUE
)

# -------------------------------------------------------------------------
# 3. Validate within-patient concordance
# -------------------------------------------------------------------------

paired$assay_condition_match <- (
  paired$blood_assay_condition ==
    paired$msg_assay_condition
)

assay_mismatch <- paired[
  !paired$assay_condition_match,
  ,
  drop = FALSE
]

if (nrow(assay_mismatch) > 0) {
  write_tsv(
    assay_mismatch,
    file.path(
      PAIRING_OUT,
      "assay_condition_mismatches.tsv"
    )
  )

  stop(
    "Found ",
    nrow(assay_mismatch),
    " Blood-MSG pair(s) with discordant historical neg/pos condition. ",
    "See assay_condition_mismatches.tsv"
  )
}

paired$disease_group_match <- (
  paired$blood_disease_group ==
    paired$msg_disease_group
)

disease_mismatch <- paired[
  !paired$disease_group_match,
  ,
  drop = FALSE
]

if (nrow(disease_mismatch) > 0) {
  write_tsv(
    disease_mismatch,
    file.path(
      PAIRING_OUT,
      "disease_group_mismatches.tsv"
    )
  )

  stop(
    "Found ",
    nrow(disease_mismatch),
    " Blood-MSG pair(s) with discordant SjD/nonSjD/HV classification. ",
    "See disease_group_mismatches.tsv"
  )
}

# Consolidated values after proving concordance.
paired$assay_condition <- paired$blood_assay_condition
paired$disease_group   <- paired$blood_disease_group

paired <- paired[
  order(paired$patient_id),
  ,
  drop = FALSE
]

rownames(paired) <- NULL

# Unmatched tissue samples.
unmatched_blood <- blood_manifest[
  !blood_manifest$patient_id %in% paired$patient_id,
  ,
  drop = FALSE
]

unmatched_msg <- msg_manifest[
  !msg_manifest$patient_id %in% paired$patient_id,
  ,
  drop = FALSE
]

write_tsv(
  unmatched_blood,
  file.path(PAIRING_OUT, "unmatched_blood.tsv")
)

write_tsv(
  unmatched_msg,
  file.path(PAIRING_OUT, "unmatched_msg.tsv")
)

# Detailed manifest preserves both tissue-specific source values.
write_tsv(
  paired,
  file.path(PAIRING_OUT, "paired_manifest_detailed.tsv")
)

# Compact analysis manifest.
paired_compact <- paired[
  ,
  c(
    "patient_id",
    "blood_sample",
    "msg_sample",
    "assay_condition",
    "disease_group"
  ),
  drop = FALSE
]

write_tsv(
  paired_compact,
  file.path(PAIRING_OUT, "paired_manifest.tsv")
)

# -------------------------------------------------------------------------
# 4. Build paired count matrices using shared genes
# -------------------------------------------------------------------------

shared_genes <- intersect(
  rownames(blood_counts),
  rownames(msg_counts)
)

if (length(shared_genes) == 0) {
  stop("Blood and MSG count matrices have zero genes in common.")
}

blood_paired <- blood_counts[
  shared_genes,
  paired_compact$blood_sample,
  drop = FALSE
]

msg_paired <- msg_counts[
  shared_genes,
  paired_compact$msg_sample,
  drop = FALSE
]

# Explicit paired sample names in downstream objects.
colnames(blood_paired) <- paste0(
  paired_compact$patient_id,
  "__Blood"
)

colnames(msg_paired) <- paste0(
  paired_compact$patient_id,
  "__MSG"
)

bundle <- list(
  manifest = paired_compact,
  manifest_detailed = paired,
  blood_counts = blood_paired,
  msg_counts = msg_paired,
  shared_genes = shared_genes,
  source = list(
    blood_counts = BLOOD_COUNTS,
    blood_sampletable = BLOOD_META,
    blood_grouping = BLOOD_GROUPING,
    msg_counts = MSG_COUNTS,
    msg_sampletable = MSG_META,
    msg_grouping = MSG_GROUPING
  )
)

saveRDS(
  bundle,
  PAIRING_BUNDLE
)

# -------------------------------------------------------------------------
# 5. Minimal historical replication gate
# -------------------------------------------------------------------------

blood_rep <- run_original_deseq2(
  blood_counts,
  blood_meta,
  "Blood",
  file.path(
    PAIRING_OUT,
    "replication_blood_deseq2.tsv"
  )
)

msg_rep <- run_original_deseq2(
  msg_counts,
  msg_meta,
  "MSG",
  file.path(
    PAIRING_OUT,
    "replication_msg_deseq2.tsv"
  )
)

blood_cmp <- compare_with_original(
  if (is.null(blood_rep)) NULL else blood_rep$results,
  BLOOD_ORIGINAL_DE,
  file.path(
    PAIRING_OUT,
    "replication_blood_vs_original.tsv"
  )
)

msg_cmp <- compare_with_original(
  if (is.null(msg_rep)) NULL else msg_rep$results,
  MSG_ORIGINAL_DE,
  file.path(
    PAIRING_OUT,
    "replication_msg_vs_original.tsv"
  )
)

# -------------------------------------------------------------------------
# 6. Summary
# -------------------------------------------------------------------------

summary_lines <- c(
  "PAIRING / REPLICATION SUMMARY",
  "=============================",
  paste("Project root:", PROJECT_ROOT),
  "",
  paste("Blood samples:", ncol(blood_counts)),
  paste("MSG samples:", ncol(msg_counts)),
  paste("Blood genes:", nrow(blood_counts)),
  paste("MSG genes:", nrow(msg_counts)),
  paste("Shared genes:", length(shared_genes)),
  paste("Matched Blood-MSG patients:", nrow(paired_compact)),
  paste("Unmatched Blood samples:", nrow(unmatched_blood)),
  paste("Unmatched MSG samples:", nrow(unmatched_msg)),
  "",
  "Disease-group counts among matched patients:",
  paste(
    capture.output(
      print(
        table(
          factor(
            paired_compact$disease_group,
            levels = c("SjD", "nonSjD", "HV")
          )
        )
      )
    ),
    collapse = "\n"
  ),
  "",
  "Historical assay-condition counts among matched patients:",
  paste(
    capture.output(
      print(table(paired_compact$assay_condition))
    ),
    collapse = "\n"
  ),
  "",
  "Cross-tissue validation:",
  "  Blood vs MSG assay_condition: concordant for all retained pairs",
  "  Blood vs MSG disease_group: concordant for all retained pairs",
  "",
  "Historical replication:",
  "  input = RawCountFile_filtered.txt",
  "  design = ~ condition",
  "  contrast = neg vs pos",
  "  second low-count filter = NONE"
)

if (!is.null(blood_cmp)) {
  summary_lines <- c(
    summary_lines,
    "",
    "Blood comparison with original GCBC result:",
    paste(
      "  genes compared:",
      blood_cmp$n_genes_compared
    ),
    paste(
      "  Spearman rho(log2FC):",
      round(blood_cmp$spearman_log2FC, 5)
    ),
    paste(
      "  sign agreement:",
      round(blood_cmp$sign_agreement, 5)
    )
  )
} else {
  summary_lines <- c(
    summary_lines,
    "",
    "Blood original DE table not found or comparison unavailable; skipped."
  )
}

if (!is.null(msg_cmp)) {
  summary_lines <- c(
    summary_lines,
    "",
    "MSG comparison with original GCBC result:",
    paste(
      "  genes compared:",
      msg_cmp$n_genes_compared
    ),
    paste(
      "  Spearman rho(log2FC):",
      round(msg_cmp$spearman_log2FC, 5)
    ),
    paste(
      "  sign agreement:",
      round(msg_cmp$sign_agreement, 5)
    )
  )
} else {
  summary_lines <- c(
    summary_lines,
    "",
    "MSG original DE table not found or comparison unavailable; skipped."
  )
}

writeLines(
  summary_lines,
  file.path(
    PAIRING_OUT,
    "pairing_summary.txt"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    PAIRING_OUT,
    "sessionInfo.txt"
  )
)

cat(
  paste(summary_lines, collapse = "\n"),
  "\n"
)
