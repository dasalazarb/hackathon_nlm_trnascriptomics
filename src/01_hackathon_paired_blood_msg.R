#!/usr/bin/env Rscript

# 01_hackathon_paired_blood_msg.R
#
# No command-line paths are required.
# Edit PROJECT_ROOT once in /filesystems/config_paths.R.
#
# INPUT
# -----
# output/00_pairing/pairing_bundle.rds
#
# BIOLOGICAL GROUP
# ----------------
# SjD / nonSjD / HV from the tissue-specific grouping files validated in 00_.
#
# MODEL
# -----
# limma-voom + duplicateCorrelation(patient_id)
# ~ disease_group * tissue
#
# Blood and MSG from the same patient are modeled as repeated measurements.
#
# CORE QUESTIONS
# --------------
# 1) SjD vs HV within Blood
# 2) SjD vs nonSjD within Blood
# 3) SjD vs HV within MSG
# 4) SjD vs nonSjD within MSG
# 5) Does the SjD-vs-HV effect differ between MSG and Blood?
# 6) Does the SjD-vs-nonSjD effect differ between MSG and Blood?

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(ggplot2)
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
    file.path(dirname(getwd()), "config_paths.R"),
    "/filesystems/config_paths.R"
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
cat("PROJECT_ROOT:", PROJECT_ROOT, "\n")

if (!file.exists(PAIRING_BUNDLE)) {
  stop(
    "Pairing bundle not found:\n  ",
    PAIRING_BUNDLE,
    "\nRun scripts/00_pair_and_validate.R first."
  )
}

dir.create(
  HACKATHON_OUT,
  recursive = TRUE,
  showWarnings = FALSE
)

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

bundle <- readRDS(PAIRING_BUNDLE)

required <- c(
  "manifest",
  "blood_counts",
  "msg_counts"
)

missing <- setdiff(
  required,
  names(bundle)
)

if (length(missing) > 0) {
  stop(
    "pairing_bundle.rds is missing: ",
    paste(missing, collapse = ", ")
  )
}

manifest <- bundle$manifest
blood <- bundle$blood_counts
msg   <- bundle$msg_counts

required_manifest <- c(
  "patient_id",
  "blood_sample",
  "msg_sample",
  "assay_condition",
  "disease_group"
)

missing_manifest <- setdiff(
  required_manifest,
  colnames(manifest)
)

if (length(missing_manifest) > 0) {
  stop(
    "Manifest is missing: ",
    paste(missing_manifest, collapse = ", "),
    ". Re-run 00_pair_and_validate.R."
  )
}

allowed_groups <- c(
  "HV",
  "nonSjD",
  "SjD"
)

bad_groups <- setdiff(
  unique(manifest$disease_group),
  allowed_groups
)

if (length(bad_groups) > 0) {
  stop(
    "Unexpected disease_group(s): ",
    paste(bad_groups, collapse = ", ")
  )
}

present_groups <- unique(
  as.character(manifest$disease_group)
)

if (!all(
  c("SjD", "nonSjD", "HV") %in% present_groups
)) {
  warning(
    "Not all three disease groups are represented among matched pairs. Present: ",
    paste(sort(present_groups), collapse = ", ")
  )
}

if (!identical(
  rownames(blood),
  rownames(msg)
)) {
  stop(
    "Blood and MSG gene rows are not identical in pairing_bundle.rds."
  )
}

if (
  ncol(blood) != nrow(manifest) ||
  ncol(msg) != nrow(manifest)
) {
  stop(
    "Expected exactly one Blood and one MSG sample per manifest row."
  )
}

# -------------------------------------------------------------------------
# 1. Build paired long-format sample metadata
# -------------------------------------------------------------------------

blood_meta <- data.frame(
  sample_id = colnames(blood),
  patient_id = manifest$patient_id,
  tissue = "Blood",
  disease_group = manifest$disease_group,
  assay_condition = manifest$assay_condition,
  stringsAsFactors = FALSE
)

msg_meta <- data.frame(
  sample_id = colnames(msg),
  patient_id = manifest$patient_id,
  tissue = "MSG",
  disease_group = manifest$disease_group,
  assay_condition = manifest$assay_condition,
  stringsAsFactors = FALSE
)

meta <- rbind(
  blood_meta,
  msg_meta
)

# Reference levels:
#   disease = HV
#   tissue  = Blood
meta$disease_group <- factor(
  meta$disease_group,
  levels = c(
    "HV",
    "nonSjD",
    "SjD"
  )
)

meta$tissue <- factor(
  meta$tissue,
  levels = c(
    "Blood",
    "MSG"
  )
)

meta$patient_id <- factor(
  meta$patient_id
)

combined_counts <- cbind(
  blood,
  msg
)

combined_counts <- combined_counts[
  ,
  meta$sample_id,
  drop = FALSE
]

stopifnot(
  all(
    colnames(combined_counts) ==
      meta$sample_id
  )
)

# -------------------------------------------------------------------------
# 2. Design-aware filter for the NEW analysis
# -------------------------------------------------------------------------

design <- model.matrix(
  ~ disease_group * tissue,
  data = meta
)

colnames(design) <- make.names(
  colnames(design)
)

write_tsv(
  cbind(
    sample_id = meta$sample_id,
    meta[
      ,
      c(
        "patient_id",
        "tissue",
        "disease_group",
        "assay_condition"
      )
    ],
    design
  ),
  file.path(
    HACKATHON_OUT,
    "model_matrix.tsv"
  )
)

y <- DGEList(
  counts = combined_counts
)

keep <- filterByExpr(
  y,
  design = design
)

y <- y[
  keep,
  ,
  keep.lib.sizes = FALSE
]

y <- calcNormFactors(
  y,
  method = "TMM"
)

# -------------------------------------------------------------------------
# 3. Repeated-measures paired model
# -------------------------------------------------------------------------

v0 <- voom(
  y,
  design,
  plot = FALSE
)

corfit <- duplicateCorrelation(
  v0,
  design,
  block = meta$patient_id
)

v <- voom(
  y,
  design,
  plot = FALSE,
  block = meta$patient_id,
  correlation = corfit$consensus.correlation
)

fit <- lmFit(
  v,
  design,
  block = meta$patient_id,
  correlation = corfit$consensus.correlation
)

# Expected after make.names():
# X.Intercept.
# disease_groupnonSjD
# disease_groupSjD
# tissueMSG
# disease_groupnonSjD.tissueMSG
# disease_groupSjD.tissueMSG

needed <- c(
  "disease_groupnonSjD",
  "disease_groupSjD",
  "tissueMSG",
  "disease_groupnonSjD.tissueMSG",
  "disease_groupSjD.tissueMSG"
)

if (!all(
  needed %in% colnames(design)
)) {
  stop(
    "Unexpected design columns: ",
    paste(
      colnames(design),
      collapse = ", "
    )
  )
}

cm <- makeContrasts(
  # Blood disease effects
  Blood_SjD_vs_HV =
    disease_groupSjD,

  Blood_nonSjD_vs_HV =
    disease_groupnonSjD,

  Blood_SjD_vs_nonSjD =
    disease_groupSjD -
    disease_groupnonSjD,

  # MSG disease effects
  MSG_SjD_vs_HV =
    disease_groupSjD +
    disease_groupSjD.tissueMSG,

  MSG_nonSjD_vs_HV =
    disease_groupnonSjD +
    disease_groupnonSjD.tissueMSG,

  MSG_SjD_vs_nonSjD =
    (
      disease_groupSjD -
      disease_groupnonSjD
    ) +
    (
      disease_groupSjD.tissueMSG -
      disease_groupnonSjD.tissueMSG
    ),

  # Formal tissue × disease interaction
  Interaction_SjD_vs_HV =
    disease_groupSjD.tissueMSG,

  Interaction_nonSjD_vs_HV =
    disease_groupnonSjD.tissueMSG,

  Interaction_SjD_vs_nonSjD =
    disease_groupSjD.tissueMSG -
    disease_groupnonSjD.tissueMSG,

  # Tissue effects within each disease group
  MSG_vs_Blood_HV =
    tissueMSG,

  MSG_vs_Blood_nonSjD =
    tissueMSG +
    disease_groupnonSjD.tissueMSG,

  MSG_vs_Blood_SjD =
    tissueMSG +
    disease_groupSjD.tissueMSG,

  levels = design
)

fit2 <- contrasts.fit(
  fit,
  cm
)

fit2 <- eBayes(
  fit2,
  trend = FALSE,
  robust = TRUE
)

extract_result <- function(
  fit_obj,
  coef_name
) {
  z <- topTable(
    fit_obj,
    coef = coef_name,
    number = Inf,
    sort.by = "P"
  )

  z$ensembl_gene_id <- rownames(z)

  z <- z[
    ,
    c(
      "ensembl_gene_id",
      "logFC",
      "AveExpr",
      "t",
      "P.Value",
      "adj.P.Val",
      "B"
    )
  ]

  z
}

results_list <- lapply(
  colnames(cm),
  function(coef_name) {
    extract_result(
      fit2,
      coef_name
    )
  }
)

names(results_list) <- colnames(cm)

for (nm in names(results_list)) {
  write_tsv(
    results_list[[nm]],
    file.path(
      HACKATHON_OUT,
      paste0(
        tolower(nm),
        ".tsv"
      )
    )
  )
}

# -------------------------------------------------------------------------
# 4. Blood-vs-MSG effect comparisons
# -------------------------------------------------------------------------

build_effect_comparison <- function(
  blood_name,
  msg_name,
  interaction_name,
  prefix
) {
  b <- results_list[[blood_name]]
  m <- results_list[[msg_name]]
  i <- results_list[[interaction_name]]

  cmp <- merge(
    b[
      ,
      c(
        "ensembl_gene_id",
        "logFC",
        "adj.P.Val"
      )
    ],
    m[
      ,
      c(
        "ensembl_gene_id",
        "logFC",
        "adj.P.Val"
      )
    ],
    by = "ensembl_gene_id",
    suffixes = c(
      "_blood",
      "_msg"
    )
  )

  cmp <- merge(
    cmp,
    i[
      ,
      c(
        "ensembl_gene_id",
        "logFC",
        "adj.P.Val"
      )
    ],
    by = "ensembl_gene_id",
    all.x = TRUE
  )

  colnames(cmp)[
    colnames(cmp) == "logFC"
  ] <- "interaction_logFC"

  colnames(cmp)[
    colnames(cmp) == "adj.P.Val"
  ] <- "interaction_FDR"

  cmp$same_direction <- (
    sign(cmp$logFC_blood) ==
      sign(cmp$logFC_msg)
  )

  cmp$class <- "no_strong_evidence"

  cmp$class[
    cmp$adj.P.Val_blood < 0.05 &
      cmp$adj.P.Val_msg < 0.05 &
      cmp$same_direction
  ] <- "shared_concordant"

  cmp$class[
    cmp$adj.P.Val_blood < 0.05 &
      cmp$adj.P.Val_msg >= 0.05
  ] <- "blood_only_candidate"

  cmp$class[
    cmp$adj.P.Val_blood >= 0.05 &
      cmp$adj.P.Val_msg < 0.05
  ] <- "msg_only_candidate"

  cmp$class[
    cmp$adj.P.Val_blood < 0.05 &
      cmp$adj.P.Val_msg < 0.05 &
      !cmp$same_direction
  ] <- "discordant"

  cmp$class[
    !is.na(cmp$interaction_FDR) &
      cmp$interaction_FDR < 0.05
  ] <- "significant_tissue_interaction"

  rho <- suppressWarnings(
    cor(
      cmp$logFC_blood,
      cmp$logFC_msg,
      method = "spearman",
      use = "complete.obs"
    )
  )

  write_tsv(
    cmp,
    file.path(
      HACKATHON_OUT,
      paste0(
        prefix,
        "_cross_tissue_effects.tsv"
      )
    )
  )

  p <- ggplot(
    cmp,
    aes(
      x = logFC_blood,
      y = logFC_msg
    )
  ) +
    geom_hline(
      yintercept = 0,
      linewidth = 0.3
    ) +
    geom_vline(
      xintercept = 0,
      linewidth = 0.3
    ) +
    geom_point(
      alpha = 0.35,
      size = 1
    ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = 2
    ) +
    labs(
      x = "Blood log2FC",
      y = "MSG log2FC",
      title = paste0(
        prefix,
        ": Blood vs MSG disease effect"
      ),
      subtitle = paste0(
        "Spearman rho = ",
        round(rho, 3)
      )
    ) +
    theme_minimal(
      base_size = 12
    )

  ggsave(
    file.path(
      HACKATHON_OUT,
      paste0(
        prefix,
        "_effect_scatter.png"
      )
    ),
    p,
    width = 7,
    height = 6,
    dpi = 300
  )

  list(
    data = cmp,
    rho = rho
  )
}

cmp_sjd_hv <- build_effect_comparison(
  "Blood_SjD_vs_HV",
  "MSG_SjD_vs_HV",
  "Interaction_SjD_vs_HV",
  "SjD_vs_HV"
)

cmp_sjd_non <- build_effect_comparison(
  "Blood_SjD_vs_nonSjD",
  "MSG_SjD_vs_nonSjD",
  "Interaction_SjD_vs_nonSjD",
  "SjD_vs_nonSjD"
)

# -------------------------------------------------------------------------
# 5. Joint PCA
# -------------------------------------------------------------------------

pca <- prcomp(
  t(v$E),
  scale. = FALSE
)

var_exp <- (
  pca$sdev^2 /
    sum(pca$sdev^2)
) * 100

pca_df <- data.frame(
  sample_id = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  patient_id = meta$patient_id,
  tissue = meta$tissue,
  disease_group = meta$disease_group,
  stringsAsFactors = FALSE
)

p <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    color = disease_group,
    shape = tissue
  )
) +
  geom_point(
    size = 2.5
  ) +
  labs(
    x = paste0(
      "PC1 (",
      round(var_exp[1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(var_exp[2], 1),
      "%)"
    ),
    title = "Matched Blood-MSG samples",
    color = "Disease group",
    shape = "Tissue"
  ) +
  theme_minimal(
    base_size = 12
  )

ggsave(
  file.path(
    HACKATHON_OUT,
    "paired_joint_PCA.png"
  ),
  p,
  width = 7,
  height = 6,
  dpi = 300
)

# -------------------------------------------------------------------------
# 6. Save model-ready expression for modules/WGCNA/endotypes
# -------------------------------------------------------------------------

saveRDS(
  list(
    E = v$E,
    weights = v$weights,
    meta = meta,
    design = design,
    contrasts = cm,
    consensus_correlation =
      corfit$consensus.correlation
  ),
  file.path(
    HACKATHON_OUT,
    "voom_expression.rds"
  )
)

# -------------------------------------------------------------------------
# 7. Summary
# -------------------------------------------------------------------------

n_sig <- function(x) {
  sum(
    x$adj.P.Val < 0.05,
    na.rm = TRUE
  )
}

summary_lines <- c(
  "HACKATHON PAIRED BLOOD-MSG SUMMARY",
  "=================================",
  paste(
    "Project root:",
    PROJECT_ROOT
  ),
  paste(
    "Matched patients:",
    nrow(manifest)
  ),
  paste(
    "Samples modeled:",
    nrow(meta)
  ),
  paste(
    "Genes entering combined matrix:",
    nrow(combined_counts)
  ),
  paste(
    "Genes after filterByExpr:",
    nrow(y)
  ),
  paste(
    "Estimated within-patient correlation:",
    round(
      corfit$consensus.correlation,
      4
    )
  ),
  "",
  "Matched patient counts:",
  paste(
    capture.output(
      print(
        table(
          factor(
            manifest$disease_group,
            levels = c(
              "SjD",
              "nonSjD",
              "HV"
            )
          )
        )
      )
    ),
    collapse = "\n"
  ),
  "",
  paste(
    "FDR<0.05 Blood SjD vs HV:",
    n_sig(
      results_list$Blood_SjD_vs_HV
    )
  ),
  paste(
    "FDR<0.05 MSG SjD vs HV:",
    n_sig(
      results_list$MSG_SjD_vs_HV
    )
  ),
  paste(
    "FDR<0.05 interaction SjD vs HV:",
    n_sig(
      results_list$Interaction_SjD_vs_HV
    )
  ),
  paste(
    "Blood-vs-MSG effect rho, SjD vs HV:",
    round(
      cmp_sjd_hv$rho,
      4
    )
  ),
  "",
  paste(
    "FDR<0.05 Blood SjD vs nonSjD:",
    n_sig(
      results_list$Blood_SjD_vs_nonSjD
    )
  ),
  paste(
    "FDR<0.05 MSG SjD vs nonSjD:",
    n_sig(
      results_list$MSG_SjD_vs_nonSjD
    )
  ),
  paste(
    "FDR<0.05 interaction SjD vs nonSjD:",
    n_sig(
      results_list$Interaction_SjD_vs_nonSjD
    )
  ),
  paste(
    "Blood-vs-MSG effect rho, SjD vs nonSjD:",
    round(
      cmp_sjd_non$rho,
      4
    )
  ),
  "",
  "Next stage:",
  "  pathway/module scoring -> WGCNA -> module projection -> patient endotypes."
)

writeLines(
  summary_lines,
  file.path(
    HACKATHON_OUT,
    "hackathon_summary.txt"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    HACKATHON_OUT,
    "sessionInfo.txt"
  )
)

cat(
  paste(
    summary_lines,
    collapse = "\n"
  ),
  "\n"
)
