#!/usr/bin/env Rscript

# ==============================================================================
# 02_05_day2_wgcna_endotypes_projection.R
#
# One-script Day 2 pipeline:
#   02) WGCNA discovery in MSG
#   03) module <-> disease associations
#   04) preliminary patient endotypes
#   05) MSG module projection into paired Blood
#
# Designed for the 97 paired Blood-MSG cohort:
#   47 SjD / 30 nonSjD / 20 HV
#
# Core principles:
#   - MSG is the discovery/reference tissue.
#   - WGCNA is fit ONLY in MSG.
#   - Endotypes are clustered from module eigengenes, not genes.
#   - Primary endotype analysis is within SjD (47 patients).
#   - Blood receives the fixed MSG module loadings; modules are NOT re-fit in Blood.
#   - Projection is paired by patient_id.
#
# Inputs expected from Stage 00 / original count files:
#   paired_manifest.tsv
#   RawCountFile_filtered_msg.txt
#   RawCountFile_filtered_blood.txt
#
# Run:
#   Rscript 02_05_day2_wgcna_endotypes_projection.R
#
# Optional environment variables:
#   PAIRED_MANIFEST=/path/to/paired_manifest.tsv
#   MSG_COUNTS=/path/to/RawCountFile_filtered_msg.txt
#   BLOOD_COUNTS=/path/to/RawCountFile_filtered_blood.txt
#   DAY2_OUT=/path/to/day2_results
#   WGCNA_THREADS=8
#
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(20260929)

# ------------------------------------------------------------------------------
# 0. CONFIGURATION
# ------------------------------------------------------------------------------

PAIRED_MANIFEST <- Sys.getenv("PAIRED_MANIFEST", "paired_manifest.tsv")
MSG_COUNTS      <- Sys.getenv("MSG_COUNTS", "RawCountFile_filtered_msg.txt")
BLOOD_COUNTS    <- Sys.getenv("BLOOD_COUNTS", "RawCountFile_filtered_blood.txt")
OUTDIR          <- Sys.getenv("DAY2_OUT", "day2_results")

# WGCNA settings: matched to the NIH-NLM endotypes-transcriptomics architecture.
SOFT_R2_TARGET  <- 0.90
SOFT_POWERS     <- c(1:10, seq(12, 20, by = 2))
SOFT_FALLBACK   <- 12       # WGCNA FAQ-style fallback for >40 samples, signed network
MIN_MODULE_SIZE <- 30
MERGE_CUT       <- 0.25
MAX_BLOCK_SIZE  <- 5000

# Endotype settings.
ENDOTYPE_K_MAX       <- 6
ENDOTYPE_MIN_CLUSTER <- 5

# Expected value from the completed Stage 01 context; only used as a QC note.
EXPECTED_FILTERED_GENES <- 22106L

# Threads.
n_threads <- suppressWarnings(as.integer(Sys.getenv(
  "WGCNA_THREADS",
  Sys.getenv("SLURM_CPUS_PER_TASK", "8")
)))
if (is.na(n_threads) || n_threads < 1) n_threads <- 1L

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. PACKAGE CHECK
# ------------------------------------------------------------------------------

required_pkgs <- c("WGCNA", "edgeR", "limma", "cluster")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  stop(
    paste0(
      "\nMissing required R packages: ", paste(missing_pkgs, collapse = ", "), "\n\n",
      "Install them before rerunning. Typical R commands:\n",
      "  if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')\n",
      "  BiocManager::install(c('edgeR','limma'))\n",
      "  install.packages(c('WGCNA','cluster'))\n"
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(WGCNA)
  library(edgeR)
  library(limma)
  library(cluster)
})

try(WGCNA::allowWGCNAThreads(nThreads = n_threads), silent = TRUE)

# ------------------------------------------------------------------------------
# 2. HELPERS
# ------------------------------------------------------------------------------

msg <- function(...) {
  cat(sprintf(...), "\n")
  flush.console()
}

wtsv <- function(x, file) {
  write.table(
    x,
    file = file.path(OUTDIR, file),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )
}

resolve_input <- function(path) {
  if (file.exists(path)) return(normalizePath(path))

  pat <- glob2rx(basename(path))
  hits <- list.files(".", pattern = pat, recursive = TRUE, full.names = TRUE)
  hits <- hits[basename(hits) == basename(path)]

  if (length(hits) == 1) return(normalizePath(hits))
  if (length(hits) > 1) {
    stop(
      "Multiple copies of ", basename(path), " were found:\n",
      paste0("  - ", hits, collapse = "\n"),
      "\nSet the corresponding environment variable to the exact file.",
      call. = FALSE
    )
  }

  stop(
    "Could not find input file: ", path,
    "\nSet PAIRED_MANIFEST, MSG_COUNTS, or BLOOD_COUNTS to the exact path.",
    call. = FALSE
  )
}

canon <- function(x) {
  x <- tolower(trimws(x))
  gsub("[^a-z0-9]+", "_", x)
}

find_col <- function(df, candidates, required = TRUE) {
  n <- names(df)
  cn <- canon(n)
  cc <- canon(candidates)

  hit <- match(cc, cn, nomatch = 0)
  hit <- hit[hit > 0]

  if (length(hit) > 0) return(n[hit[1]])

  if (required) {
    stop(
      "Could not identify a required column.\nCandidates: ",
      paste(candidates, collapse = ", "),
      "\nAvailable columns: ",
      paste(names(df), collapse = ", "),
      call. = FALSE
    )
  }
  NULL
}

normalize_group <- function(x) {
  y <- tolower(trimws(as.character(x)))
  y <- gsub("[[:space:]_-]+", "", y)

  out <- as.character(x)
  out[y %in% c("sjd", "sjogren", "sjogrens", "sjogrensdisease")] <- "SjD"
  out[y %in% c("nonsjd", "nonsjogren", "nonsjogrens", "negative")] <- "nonSjD"
  out[y %in% c("hv", "healthy", "healthyvolunteer", "healthyvolunteers", "control")] <- "HV"
  trimws(out)
}

normalize_tissue <- function(x) {
  y <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(y))
  out[grepl("blood|whole.?blood|pax", y)] <- "Blood"
  out[grepl("msg|minor.?salivary|salivary.?gland|gland", y)] <- "MSG"
  out
}

read_manifest_wide <- function(file) {
  m <- read.delim(
    file,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )

  patient_col <- find_col(
    m,
    c("patient_id", "patient", "subject_id", "subject", "participant_id", "participant")
  )

  # First try long format: one sample per row with tissue label.
  tissue_col <- find_col(m, c("tissue", "sample_tissue", "specimen_tissue"), required = FALSE)
  sample_col <- find_col(
    m,
    c("sample", "sample_id", "sample_name", "samplename", "rna_sample", "rnaseq_sample"),
    required = FALSE
  )
  group_col <- find_col(
    m,
    c("disease_group", "disease", "group", "sjd_group", "diagnosis", "classification"),
    required = FALSE
  )

  if (!is.null(tissue_col) && !is.null(sample_col)) {
    tt <- normalize_tissue(m[[tissue_col]])

    if (any(is.na(tt))) {
      bad <- unique(as.character(m[[tissue_col]])[is.na(tt)])
      stop(
        "Unrecognized tissue values in manifest: ",
        paste(bad, collapse = ", "),
        call. = FALSE
      )
    }

    if (is.null(group_col)) {
      stop("Long-format manifest detected, but no disease_group column was found.", call. = FALSE)
    }

    long <- data.frame(
      patient_id = as.character(m[[patient_col]]),
      sample = as.character(m[[sample_col]]),
      tissue = tt,
      disease_group = normalize_group(m[[group_col]]),
      stringsAsFactors = FALSE
    )

    if (anyDuplicated(long[, c("patient_id", "tissue")])) {
      stop("Manifest has duplicate patient_id x tissue rows.", call. = FALSE)
    }

    wb <- reshape(
      long,
      idvar = "patient_id",
      timevar = "tissue",
      direction = "wide"
    )

    need <- c("sample.Blood", "sample.MSG")
    if (!all(need %in% names(wb))) {
      stop("Could not reconstruct both Blood and MSG samples from long-format manifest.", call. = FALSE)
    }

    if ("disease_group.Blood" %in% names(wb) && "disease_group.MSG" %in% names(wb)) {
      if (any(wb$disease_group.Blood != wb$disease_group.MSG)) {
        stop("Disease-group mismatch between Blood and MSG in manifest.", call. = FALSE)
      }
      dg <- wb$disease_group.Blood
    } else {
      stop("Disease group could not be reconstructed from long-format manifest.", call. = FALSE)
    }

    out <- data.frame(
      patient_id = wb$patient_id,
      blood_sample = wb$sample.Blood,
      msg_sample = wb$sample.MSG,
      disease_group = dg,
      stringsAsFactors = FALSE
    )
  } else {
    # Wide format: one patient per row.
    blood_col <- find_col(
      m,
      c(
        "blood_sample", "blood_sample_id", "sample_blood", "samplename_blood",
        "blood_samplename", "sample_name_blood", "blood_sample_name", "blood"
      )
    )
    msg_col <- find_col(
      m,
      c(
        "msg_sample", "msg_sample_id", "sample_msg", "samplename_msg",
        "msg_samplename", "sample_name_msg", "msg_sample_name", "msg"
      )
    )

    if (!is.null(group_col)) {
      dg <- normalize_group(m[[group_col]])
    } else {
      gb <- find_col(
        m,
        c("disease_group_blood", "blood_disease_group", "group_blood"),
        required = FALSE
      )
      gm <- find_col(
        m,
        c("disease_group_msg", "msg_disease_group", "group_msg"),
        required = FALSE
      )

      if (is.null(gb) || is.null(gm)) {
        stop("No usable disease_group column found in paired_manifest.tsv.", call. = FALSE)
      }

      dg_b <- normalize_group(m[[gb]])
      dg_m <- normalize_group(m[[gm]])

      if (any(dg_b != dg_m)) {
        stop("Disease-group mismatch between Blood and MSG in wide manifest.", call. = FALSE)
      }
      dg <- dg_b
    }

    out <- data.frame(
      patient_id = as.character(m[[patient_col]]),
      blood_sample = as.character(m[[blood_col]]),
      msg_sample = as.character(m[[msg_col]]),
      disease_group = dg,
      stringsAsFactors = FALSE
    )
  }

  out$patient_id <- trimws(out$patient_id)
  out$blood_sample <- trimws(out$blood_sample)
  out$msg_sample <- trimws(out$msg_sample)
  out$disease_group <- normalize_group(out$disease_group)

  if (anyDuplicated(out$patient_id)) {
    stop("patient_id is not unique in the paired manifest.", call. = FALSE)
  }
  if (anyDuplicated(out$blood_sample)) {
    stop("blood_sample is not unique in the paired manifest.", call. = FALSE)
  }
  if (anyDuplicated(out$msg_sample)) {
    stop("msg_sample is not unique in the paired manifest.", call. = FALSE)
  }
  if (any(!out$disease_group %in% c("HV", "nonSjD", "SjD"))) {
    stop(
      "Unexpected disease_group values: ",
      paste(unique(out$disease_group[!out$disease_group %in% c("HV", "nonSjD", "SjD")]), collapse = ", "),
      call. = FALSE
    )
  }

  out
}

read_count_matrix <- function(file) {
  d <- read.delim(
    file,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )

  if (ncol(d) < 2) {
    stop("Count file has fewer than two columns: ", file, call. = FALSE)
  }

  gene <- as.character(d[[1]])
  if (anyNA(gene) || any(gene == "")) {
    stop("First column of count file must contain non-empty gene IDs: ", file, call. = FALSE)
  }

  mat <- as.matrix(d[, -1, drop = FALSE])
  suppressWarnings(storage.mode(mat) <- "numeric")

  if (anyNA(mat)) {
    stop(
      "Non-numeric or missing count values detected in: ", file,
      "\nExpected first column = gene ID, remaining columns = samples.",
      call. = FALSE
    )
  }

  rownames(mat) <- gene

  if (anyDuplicated(rownames(mat))) {
    msg("  Duplicated gene IDs found in %s; summing duplicate rows.", basename(file))
    mat <- rowsum(mat, group = rownames(mat), reorder = FALSE)
  }

  mat
}

check_samples <- function(required, available, tissue) {
  missing <- setdiff(required, available)
  if (length(missing) > 0) {
    stop(
      tissue, ": ", length(missing), " paired samples are missing from the count matrix.\n",
      "First missing IDs: ", paste(head(missing, 10), collapse = ", "),
      call. = FALSE
    )
  }
}

fallback_power <- function(n_samples) {
  # Signed-network fallback inspired by the WGCNA FAQ table.
  if (n_samples < 20) return(18)
  if (n_samples < 30) return(16)
  if (n_samples < 40) return(14)
  12
}

fit_module_loadings <- function(datExpr, module_colors, MEs) {
  modules <- setdiff(unique(module_colors), "grey")

  out <- lapply(setNames(modules, modules), function(mo) {
    g <- names(module_colors)[module_colors == mo]
    X <- datExpr[, g, drop = FALSE]

    mu <- colMeans(X)
    s <- apply(X, 2, sd)

    ok <- is.finite(s) & s > 0
    X <- X[, ok, drop = FALSE]
    mu <- mu[ok]
    s <- s[ok]
    g <- g[ok]

    Z <- sweep(X, 2, mu, "-")
    Z <- sweep(Z, 2, s, "/")

    sv <- svd(t(Z), nu = 1, nv = 1)
    wt <- sv$u[, 1] / sv$d[1]
    sc <- as.vector(Z %*% wt)

    me_col <- paste0("ME", mo)
    if (!me_col %in% colnames(MEs)) {
      stop("Missing eigengene column for module ", mo, call. = FALSE)
    }

    # Orient weights exactly as the WGCNA eigengene.
    if (cor(sc, MEs[, me_col], use = "pairwise.complete.obs") < 0) {
      wt <- -wt
      sc <- -sc
    }

    validation_r <- cor(sc, MEs[, me_col], use = "pairwise.complete.obs")
    validation_max_abs_diff <- max(abs(sc - MEs[, me_col]))

    data.frame(
      module = mo,
      gene = g,
      mean_MSG = as.numeric(mu),
      sd_MSG = as.numeric(s),
      weight = as.numeric(wt),
      validation_r = validation_r,
      validation_max_abs_diff = validation_max_abs_diff,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out)
}

module_disease_tests <- function(score_matrix, disease_group, score_family) {
  score_matrix <- as.matrix(score_matrix)
  disease_group <- factor(disease_group, levels = c("HV", "nonSjD", "SjD"))

  if (anyNA(disease_group)) {
    stop("Disease group contains NA or unexpected labels.", call. = FALSE)
  }

  informative <- apply(
    score_matrix,
    2,
    function(z) all(is.finite(z)) && is.finite(sd(z)) && sd(z) > 0
  )
  score_matrix <- score_matrix[, informative, drop = FALSE]

  if (ncol(score_matrix) == 0) {
    stop("No informative module scores were available for disease testing.", call. = FALSE)
  }

  # Standardize each module across patients. This makes contrast estimates SD units.
  Z <- scale(score_matrix)
  Z <- as.matrix(Z)

  # Omnibus three-group ANOVA.
  omnibus <- do.call(rbind, lapply(seq_len(ncol(Z)), function(j) {
    y <- Z[, j]
    fit <- lm(y ~ disease_group)
    a <- anova(fit)

    ss_between <- a$`Sum Sq`[1]
    ss_total <- sum(a$`Sum Sq`)

    data.frame(
      score_family = score_family,
      module = colnames(Z)[j],
      n = length(y),
      F = a$`F value`[1],
      p = a$`Pr(>F)`[1],
      eta2 = ss_between / ss_total,
      stringsAsFactors = FALSE
    )
  }))
  omnibus$fdr <- p.adjust(omnibus$p, method = "BH")
  omnibus <- omnibus[order(omnibus$p), ]

  # Limma pairwise contrasts.
  design <- model.matrix(~ 0 + disease_group)
  colnames(design) <- levels(disease_group)

  C <- matrix(
    0,
    nrow = ncol(design),
    ncol = 3,
    dimnames = list(
      colnames(design),
      c("SjD_vs_HV", "nonSjD_vs_HV", "SjD_vs_nonSjD")
    )
  )
  C["SjD", "SjD_vs_HV"] <- 1
  C["HV", "SjD_vs_HV"] <- -1

  C["nonSjD", "nonSjD_vs_HV"] <- 1
  C["HV", "nonSjD_vs_HV"] <- -1

  C["SjD", "SjD_vs_nonSjD"] <- 1
  C["nonSjD", "SjD_vs_nonSjD"] <- -1

  fit <- limma::lmFit(t(Z), design)
  fit <- limma::contrasts.fit(fit, C)
  fit <- limma::eBayes(fit)

  pairwise <- do.call(rbind, lapply(colnames(C), function(coef_name) {
    tt <- limma::topTable(
      fit,
      coef = coef_name,
      number = Inf,
      sort.by = "none",
      adjust.method = "BH"
    )

    data.frame(
      score_family = score_family,
      module = rownames(tt),
      contrast = coef_name,
      effect_SD = tt$logFC,
      t = tt$t,
      p = tt$P.Value,
      fdr = tt$adj.P.Val,
      stringsAsFactors = FALSE
    )
  }))
  rownames(pairwise) <- NULL

  list(z = Z, omnibus = omnibus, pairwise = pairwise)
}

cluster_endotypes <- function(score_matrix, patient_ids, prefix, min_cluster_size = 5, k_max = 6) {
  X <- as.matrix(score_matrix)

  # Re-standardize within the clustering subset.
  X <- scale(X)
  X <- X[, apply(X, 2, function(z) all(is.finite(z)) && sd(z) > 0), drop = FALSE]

  if (ncol(X) < 2) stop("Need at least two informative modules for clustering.", call. = FALSE)
  if (nrow(X) < 10) stop("Too few patients for endotype clustering.", call. = FALSE)

  d <- dist(X, method = "euclidean")
  hc <- hclust(d, method = "ward.D2")

  ks <- 2:min(k_max, nrow(X) - 1)

  sil <- data.frame(
    k = ks,
    average_silhouette = NA_real_,
    min_cluster_n = NA_integer_
  )

  for (i in seq_along(ks)) {
    k <- ks[i]
    cl <- cutree(hc, k = k)
    sil$min_cluster_n[i] <- min(table(cl))

    if (sil$min_cluster_n[i] >= min_cluster_size) {
      s <- cluster::silhouette(cl, d)
      sil$average_silhouette[i] <- mean(s[, "sil_width"])
    }
  }

  valid <- which(is.finite(sil$average_silhouette))

  # If the minimum-size rule is too strict, still return a result, but record it.
  if (length(valid) == 0) {
    for (i in seq_along(ks)) {
      k <- ks[i]
      cl <- cutree(hc, k = k)
      s <- cluster::silhouette(cl, d)
      sil$average_silhouette[i] <- mean(s[, "sil_width"])
    }
    valid <- seq_len(nrow(sil))
  }

  best_row <- valid[which.max(sil$average_silhouette[valid])]
  best_k <- sil$k[best_row]
  cl <- cutree(hc, k = best_k)

  # Deterministic human-readable labels ordered by centroid along PC1.
  pc <- prcomp(X, center = FALSE, scale. = FALSE)
  pc1 <- pc$x[, 1]
  centroid_pc1 <- tapply(pc1, cl, mean)
  ordered_clusters <- names(sort(centroid_pc1))
  label_map <- setNames(paste0(prefix, seq_along(ordered_clusters)), ordered_clusters)
  labels <- unname(label_map[as.character(cl)])

  assignments <- data.frame(
    patient_id = patient_ids,
    ward_cluster = cl,
    endotype = labels,
    stringsAsFactors = FALSE
  )

  list(
    X = X,
    dist = d,
    hc = hc,
    silhouette = sil,
    best_k = best_k,
    assignments = assignments,
    pc = pc
  )
}

safe_cor_test <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]

  if (length(x) < 4 || sd(x) == 0 || sd(y) == 0) {
    return(c(estimate = NA_real_, p = NA_real_, n = length(x)))
  }

  z <- suppressWarnings(
    cor.test(
      x,
      y,
      method = method,
      exact = if (method == "spearman") FALSE else NULL
    )
  )

  c(
    estimate = unname(z$estimate),
    p = z$p.value,
    n = length(x)
  )
}

project_loadings_target_z <- function(E_target, loadings) {
  # E_target: genes x patients.
  gene_sd <- apply(E_target, 1, sd)
  gene_mu <- rowMeans(E_target)

  usable <- is.finite(gene_sd) & gene_sd > 0
  Z <- E_target[usable, , drop = FALSE]
  Z <- sweep(Z, 1, gene_mu[usable], "-")
  Z <- sweep(Z, 1, gene_sd[usable], "/")

  split_ld <- split(loadings, loadings$module)

  coverage <- do.call(rbind, lapply(split_ld, function(d) {
    keep <- d$gene %in% rownames(Z)
    data.frame(
      module = d$module[1],
      genes_in_MSG_module = nrow(d),
      genes_used_in_Blood = sum(keep),
      share_used = sum(keep) / nrow(d),
      stringsAsFactors = FALSE
    )
  }))

  scores <- sapply(split_ld, function(d) {
    d <- d[d$gene %in% rownames(Z), , drop = FALSE]
    if (nrow(d) < 2) {
      return(rep(NA_real_, ncol(Z)))
    }
    colSums(Z[d$gene, , drop = FALSE] * d$weight)
  })

  scores <- as.matrix(scores)
  rownames(scores) <- colnames(E_target)

  list(scores = scores, coverage = coverage)
}

# ------------------------------------------------------------------------------
# 3. LOAD + VALIDATE THE 97 PAIRED PATIENTS
# ------------------------------------------------------------------------------

msg("=== DAY 2: 02-05 in one script ===")
msg("Threads: %d", n_threads)

manifest_file <- resolve_input(PAIRED_MANIFEST)
msg_counts_file <- resolve_input(MSG_COUNTS)
blood_counts_file <- resolve_input(BLOOD_COUNTS)

msg("Manifest: %s", manifest_file)
msg("MSG counts: %s", msg_counts_file)
msg("Blood counts: %s", blood_counts_file)

manifest <- read_manifest_wide(manifest_file)

msg("Paired patients in manifest: %d", nrow(manifest))
print(table(manifest$disease_group))

if (nrow(manifest) != 97) {
  warning("Expected 97 paired patients from Stage 00, but manifest has ", nrow(manifest), ".")
}

counts_msg <- read_count_matrix(msg_counts_file)
counts_blood <- read_count_matrix(blood_counts_file)

check_samples(manifest$msg_sample, colnames(counts_msg), "MSG")
check_samples(manifest$blood_sample, colnames(counts_blood), "Blood")

common_genes <- intersect(rownames(counts_msg), rownames(counts_blood))
if (length(common_genes) < 1000) {
  stop("Very few common genes between Blood and MSG count matrices.", call. = FALSE)
}

counts_msg <- counts_msg[common_genes, manifest$msg_sample, drop = FALSE]
counts_blood <- counts_blood[common_genes, manifest$blood_sample, drop = FALSE]

# Internal sample IDs guarantee uniqueness across tissues.
msg_internal <- paste0("MSG__", seq_len(nrow(manifest)))
blood_internal <- paste0("Blood__", seq_len(nrow(manifest)))

colnames(counts_msg) <- msg_internal
colnames(counts_blood) <- blood_internal

sample_meta <- rbind(
  data.frame(
    internal_sample = msg_internal,
    raw_sample = manifest$msg_sample,
    patient_id = manifest$patient_id,
    disease_group = manifest$disease_group,
    tissue = "MSG",
    stringsAsFactors = FALSE
  ),
  data.frame(
    internal_sample = blood_internal,
    raw_sample = manifest$blood_sample,
    patient_id = manifest$patient_id,
    disease_group = manifest$disease_group,
    tissue = "Blood",
    stringsAsFactors = FALSE
  )
)
rownames(sample_meta) <- sample_meta$internal_sample

counts_all <- cbind(counts_msg, counts_blood)
counts_all <- counts_all[, sample_meta$internal_sample, drop = FALSE]

# ------------------------------------------------------------------------------
# 4. JOINT FILTERING + NORMALIZATION
#    Recreates the Stage-01-like expression space without re-fitting DE.
# ------------------------------------------------------------------------------

sample_meta$disease_group <- factor(
  sample_meta$disease_group,
  levels = c("HV", "nonSjD", "SjD")
)
sample_meta$tissue <- factor(sample_meta$tissue, levels = c("Blood", "MSG"))

design <- model.matrix(~ disease_group * tissue, data = sample_meta)

y <- edgeR::DGEList(counts = round(counts_all))
keep <- edgeR::filterByExpr(y, design = design)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- edgeR::calcNormFactors(y, method = "TMM")

v <- limma::voom(y, design = design, plot = FALSE)
E <- v$E

msg("Common genes before Stage-01-style filtering: %d", length(common_genes))
msg("Genes after filterByExpr: %d", nrow(E))

if (nrow(E) != EXPECTED_FILTERED_GENES) {
  msg(
    "QC note: Stage 01 context reported %d genes; this script obtained %d.",
    EXPECTED_FILTERED_GENES,
    nrow(E)
  )
  msg("A small difference can occur if Stage 01 used slightly different filterByExpr inputs.")
}

# Paired patient expression matrices: genes x patients.
E_MSG <- E[, msg_internal, drop = FALSE]
E_BLOOD <- E[, blood_internal, drop = FALSE]
colnames(E_MSG) <- manifest$patient_id
colnames(E_BLOOD) <- manifest$patient_id

# ------------------------------------------------------------------------------
# 5. STAGE 02 — WGCNA IN MSG
# ------------------------------------------------------------------------------

msg("\n=== 02 WGCNA MSG ===")

datExpr <- t(E_MSG)  # patients x genes
rownames(datExpr) <- manifest$patient_id

gsg <- WGCNA::goodSamplesGenes(datExpr, verbose = 0)
msg("goodSamplesGenes: %d bad samples, %d bad genes",
    sum(!gsg$goodSamples), sum(!gsg$goodGenes))

datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]

if (sum(!gsg$goodSamples) > 0) {
  stop(
    "WGCNA excluded one or more MSG samples. For this paired hackathon pipeline, inspect them before continuing.",
    call. = FALSE
  )
}

# Soft-threshold selection.
msg("Estimating signed-network soft threshold...")

sft <- WGCNA::pickSoftThreshold(
  datExpr,
  powerVector = SOFT_POWERS,
  networkType = "signed",
  blockSize = min(MAX_BLOCK_SIZE, ncol(datExpr)),
  verbose = 3
)

fit <- sft$fitIndices
soft_table <- data.frame(
  power = fit$Power,
  SFT_R2 = fit$SFT.R.sq,
  slope = fit$slope,
  mean_connectivity = fit$mean.k.,
  median_connectivity = fit$median.k.,
  max_connectivity = fit$max.k.,
  stringsAsFactors = FALSE
)

candidate <- soft_table$power[
  is.finite(soft_table$SFT_R2) &
    soft_table$SFT_R2 >= SOFT_R2_TARGET &
    soft_table$slope < 0
]

if (length(candidate) > 0) {
  soft_power <- candidate[1]
  power_reason <- paste0("first power with SFT R2 >= ", SOFT_R2_TARGET, " and negative slope")
} else {
  soft_power <- fallback_power(nrow(datExpr))
  power_reason <- paste0("fallback because no tested power reached SFT R2 >= ", SOFT_R2_TARGET)
}

soft_table$selected <- soft_table$power == soft_power
wtsv(soft_table, "02_soft_threshold.tsv")

msg("Selected soft power: %d (%s)", soft_power, power_reason)

pdf(file.path(OUTDIR, "02_soft_threshold.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))

plot(
  soft_table$power,
  soft_table$SFT_R2,
  type = "b",
  xlab = "Soft threshold (power)",
  ylab = "Scale-free topology fit R^2",
  main = "MSG: scale-free fit"
)
abline(h = SOFT_R2_TARGET, lty = 2)
abline(v = soft_power, lty = 3)

plot(
  soft_table$power,
  soft_table$mean_connectivity,
  type = "b",
  xlab = "Soft threshold (power)",
  ylab = "Mean connectivity",
  main = "MSG: mean connectivity"
)
abline(v = soft_power, lty = 3)
dev.off()

# Module discovery.
msg("Running blockwiseModules on %d MSG patients x %d genes...", nrow(datExpr), ncol(datExpr))

net <- WGCNA::blockwiseModules(
  datExpr,
  power = soft_power,
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = MIN_MODULE_SIZE,
  mergeCutHeight = MERGE_CUT,
  maxBlockSize = MAX_BLOCK_SIZE,
  numericLabels = FALSE,
  pamRespectsDendro = FALSE,
  verbose = 3
)

module_colors <- net$colors
names(module_colors) <- colnames(datExpr)

module_sizes <- sort(table(module_colors), decreasing = TRUE)
module_sizes_df <- data.frame(
  module = names(module_sizes),
  n_genes = as.integer(module_sizes),
  is_grey = names(module_sizes) == "grey",
  stringsAsFactors = FALSE
)
wtsv(module_sizes_df, "02_module_sizes.tsv")

MEs <- WGCNA::orderMEs(net$MEs)
MEs <- MEs[, colnames(MEs) != "MEgrey", drop = FALSE]
rownames(MEs) <- rownames(datExpr)

msg(
  "WGCNA complete: %d biological modules; %d grey genes.",
  ncol(MEs),
  if ("grey" %in% names(module_sizes)) unname(module_sizes["grey"]) else 0
)

# Module membership / hub genes.
kME_mat <- cor(datExpr, MEs, use = "pairwise.complete.obs")

own_kME <- rep(NA_real_, length(module_colors))
for (i in seq_along(module_colors)) {
  me_col <- paste0("ME", module_colors[i])
  if (me_col %in% colnames(kME_mat)) {
    own_kME[i] <- kME_mat[i, me_col]
  }
}

module_genes <- data.frame(
  gene = names(module_colors),
  module = unname(module_colors),
  kME = own_kME,
  stringsAsFactors = FALSE
)
module_genes <- module_genes[order(module_genes$module, -module_genes$kME), ]
wtsv(module_genes, "02_module_genes.tsv")

hub_genes <- do.call(
  rbind,
  lapply(setdiff(unique(module_genes$module), "grey"), function(mo) {
    d <- module_genes[module_genes$module == mo, , drop = FALSE]
    d <- d[order(-d$kME), , drop = FALSE]
    head(d, 20)
  })
)
rownames(hub_genes) <- NULL
wtsv(hub_genes, "02_hub_genes_top20.tsv")

# Module eigengenes.
ME_out <- data.frame(
  patient_id = rownames(MEs),
  MEs,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
wtsv(ME_out, "02_MSG_module_eigengenes.tsv")

MEs_z <- scale(MEs)
MEz_out <- data.frame(
  patient_id = rownames(MEs_z),
  MEs_z,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
wtsv(MEz_out, "02_MSG_module_eigengenes_z.tsv")

# Dendrograms for all WGCNA blocks.
pdf(file.path(OUTDIR, "02_wgcna_dendrograms.pdf"), width = 12, height = 7)
for (b in seq_along(net$dendrograms)) {
  block_genes <- net$blockGenes[[b]]
  WGCNA::plotDendroAndColors(
    net$dendrograms[[b]],
    module_colors[block_genes],
    groupLabels = paste0("Block ", b),
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05,
    main = paste0("MSG WGCNA block ", b)
  )
}
dev.off()

# Fit explicit MSG PCA/SVD loadings for portable projection.
msg("Fitting portable MSG module loadings...")
loadings <- fit_module_loadings(datExpr, module_colors, MEs)

validation <- unique(
  loadings[, c("module", "validation_r", "validation_max_abs_diff")]
)
wtsv(validation, "02_module_loading_validation.tsv")
wtsv(loadings, "02_module_loadings.tsv")

if (any(validation$validation_r < 0.999999, na.rm = TRUE)) {
  warning("At least one module loading did not perfectly reproduce its WGCNA eigengene.")
}

# ------------------------------------------------------------------------------
# 6. STAGE 03 — MODULE <-> DISEASE ASSOCIATIONS
# ------------------------------------------------------------------------------

msg("\n=== 03 MODULE TRAITS ===")

meta_msg <- manifest[match(rownames(MEs), manifest$patient_id), , drop = FALSE]

# Use module names without the "ME" prefix for downstream tables.
score_MSG <- MEs
colnames(score_MSG) <- sub("^ME", "", colnames(score_MSG))

traits_MSG <- module_disease_tests(
  score_matrix = score_MSG,
  disease_group = meta_msg$disease_group,
  score_family = "MSG_WGCNA"
)

wtsv(traits_MSG$omnibus, "03_module_traits_omnibus.tsv")
wtsv(traits_MSG$pairwise, "03_module_traits_pairwise.tsv")

msg(
  "MSG modules associated with disease_group at omnibus FDR < 0.05: %d / %d",
  sum(traits_MSG$omnibus$fdr < 0.05, na.rm = TRUE),
  nrow(traits_MSG$omnibus)
)

# Quick presentation-ready plots for the six strongest omnibus modules.
top_modules <- head(traits_MSG$omnibus$module, 6)

pdf(file.path(OUTDIR, "03_top_module_by_disease_boxplots.pdf"), width = 10, height = 7)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
for (mo in top_modules) {
  yv <- traits_MSG$z[, mo]
  boxplot(
    yv ~ meta_msg$disease_group,
    xlab = "",
    ylab = "Module score (z)",
    main = mo
  )
  stripchart(
    yv ~ meta_msg$disease_group,
    vertical = TRUE,
    method = "jitter",
    add = TRUE,
    pch = 16,
    cex = 0.6
  )
}
dev.off()

# ------------------------------------------------------------------------------
# 7. STAGE 04 — PRELIMINARY PATIENT ENDOTYPES
# ------------------------------------------------------------------------------

msg("\n=== 04 PATIENT ENDOTYPES ===")
msg("Primary endotype discovery is within SjD to avoid diagnostic-group separation driving clusters.")

# Primary: SjD only.
idx_sjd <- meta_msg$disease_group == "SjD"
sjd_ids <- meta_msg$patient_id[idx_sjd]

end_sjd <- cluster_endotypes(
  score_matrix = score_MSG[idx_sjd, , drop = FALSE],
  patient_ids = sjd_ids,
  prefix = "E",
  min_cluster_size = ENDOTYPE_MIN_CLUSTER,
  k_max = ENDOTYPE_K_MAX
)

end_sjd$assignments$disease_group <- "SjD"
wtsv(end_sjd$assignments, "04_preliminary_endotypes_SjD.tsv")
wtsv(end_sjd$silhouette, "04_silhouette_SjD.tsv")

msg("Selected preliminary SjD endotype k: %d", end_sjd$best_k)
print(table(end_sjd$assignments$endotype))

# Secondary: all 97 patients, useful as a cohort-level map.
end_all <- cluster_endotypes(
  score_matrix = score_MSG,
  patient_ids = meta_msg$patient_id,
  prefix = "C",
  min_cluster_size = ENDOTYPE_MIN_CLUSTER,
  k_max = ENDOTYPE_K_MAX
)

end_all$assignments$disease_group <- meta_msg$disease_group[
  match(end_all$assignments$patient_id, meta_msg$patient_id)
]
wtsv(end_all$assignments, "04_patient_clusters_all97.tsv")
wtsv(end_all$silhouette, "04_silhouette_all97.tsv")

# Dendrogram + PCA plots.
pdf(file.path(OUTDIR, "04_endotypes_SjD_diagnostics.pdf"), width = 10, height = 7)
par(mfrow = c(2, 2))

plot(
  end_sjd$hc,
  labels = end_sjd$assignments$endotype,
  main = paste0("SjD Ward.D2 endotypes (k=", end_sjd$best_k, ")"),
  xlab = "",
  sub = ""
)

plot(
  end_sjd$silhouette$k,
  end_sjd$silhouette$average_silhouette,
  type = "b",
  xlab = "k",
  ylab = "Average silhouette",
  main = "Choice of k"
)
abline(v = end_sjd$best_k, lty = 2)

pc_sjd <- end_sjd$pc$x
grp_num <- as.integer(factor(end_sjd$assignments$endotype))
plot(
  pc_sjd[, 1],
  pc_sjd[, 2],
  pch = 19,
  col = grp_num,
  xlab = "PC1 of module space",
  ylab = "PC2 of module space",
  main = "SjD module-space PCA"
)
legend(
  "topright",
  legend = levels(factor(end_sjd$assignments$endotype)),
  col = seq_along(levels(factor(end_sjd$assignments$endotype))),
  pch = 19,
  cex = 0.8
)

plot(
  end_all$hc,
  labels = end_all$assignments$disease_group,
  main = paste0("All 97: Ward.D2 patient map (k=", end_all$best_k, ")"),
  xlab = "",
  sub = ""
)
dev.off()

# ------------------------------------------------------------------------------
# 8. STAGE 05 — PROJECT FIXED MSG MODULES INTO BLOOD
# ------------------------------------------------------------------------------

msg("\n=== 05 MSG -> BLOOD PROJECTION ===")
msg("Using fixed MSG module gene weights; modules are not re-fit in Blood.")
msg("For cross-tissue projection, each Blood gene is standardized within Blood before weights are applied.")

proj <- project_loadings_target_z(E_BLOOD, loadings)
blood_scores <- proj$scores
coverage <- proj$coverage

# Match module order with MSG eigengenes.
common_modules <- intersect(colnames(score_MSG), colnames(blood_scores))
score_MSG_common <- score_MSG[, common_modules, drop = FALSE]
blood_scores <- blood_scores[, common_modules, drop = FALSE]

# Ensure identical patient order.
blood_scores <- blood_scores[meta_msg$patient_id, , drop = FALSE]
stopifnot(identical(rownames(blood_scores), meta_msg$patient_id))

blood_out <- data.frame(
  patient_id = rownames(blood_scores),
  blood_scores,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
wtsv(blood_out, "05_Blood_projected_module_scores.tsv")
wtsv(coverage, "05_projection_gene_coverage.tsv")

# Optional structural check, analogous to the cross-study projection check:
# recompute a local Blood PC1 for each fixed MSG gene set and compare it with
# the fixed-weight projected score. These local Blood PCs are diagnostic only;
# they are NOT used to define the Blood scores.
projection_genes <- unique(loadings$gene)
projection_genes <- projection_genes[projection_genes %in% rownames(E_BLOOD)]
blood_gene_sd <- apply(E_BLOOD[projection_genes, , drop = FALSE], 1, sd)
projection_genes <- projection_genes[is.finite(blood_gene_sd) & blood_gene_sd > 0]

local_colors <- module_colors[projection_genes]
local_blood_MEs <- WGCNA::moduleEigengenes(
  t(E_BLOOD[projection_genes, , drop = FALSE]),
  colors = local_colors,
  excludeGrey = TRUE
)$eigengenes
local_blood_MEs <- WGCNA::orderMEs(local_blood_MEs)
rownames(local_blood_MEs) <- colnames(E_BLOOD)

local_check <- do.call(rbind, lapply(common_modules, function(mo) {
  me_col <- paste0("ME", mo)
  if (!me_col %in% colnames(local_blood_MEs)) {
    return(data.frame(
      module = mo,
      r_projected_vs_localBlood = NA_real_,
      abs_r_projected_vs_localBlood = NA_real_
    ))
  }

  rr <- cor(
    blood_scores[, mo],
    local_blood_MEs[rownames(blood_scores), me_col],
    use = "pairwise.complete.obs"
  )

  data.frame(
    module = mo,
    r_projected_vs_localBlood = rr,
    abs_r_projected_vs_localBlood = abs(rr)
  )
}))
wtsv(local_check, "05_projection_vs_local_Blood_PC1.tsv")

# Paired cross-tissue concordance.
projection_stats <- do.call(rbind, lapply(common_modules, function(mo) {
  x <- score_MSG_common[, mo]
  yb <- blood_scores[, mo]
  dg <- meta_msg$disease_group

  all_s <- safe_cor_test(x, yb, method = "spearman")

  sjd <- dg == "SjD"
  sjd_s <- safe_cor_test(x[sjd], yb[sjd], method = "spearman")

  # Correlation after removing disease-group mean differences.
  rx <- residuals(lm(x ~ dg))
  ry <- residuals(lm(yb ~ dg))
  resid_s <- safe_cor_test(rx, ry, method = "spearman")

  pear <- safe_cor_test(x, yb, method = "pearson")

  data.frame(
    module = mo,
    n_all = all_s["n"],
    rho_all = all_s["estimate"],
    p_all = all_s["p"],
    rho_SjD = sjd_s["estimate"],
    p_SjD = sjd_s["p"],
    rho_residualized_group = resid_s["estimate"],
    p_residualized_group = resid_s["p"],
    pearson_all = pear["estimate"],
    pearson_p_all = pear["p"],
    stringsAsFactors = FALSE
  )
}))

projection_stats$fdr_all <- p.adjust(projection_stats$p_all, method = "BH")
projection_stats$fdr_SjD <- p.adjust(projection_stats$p_SjD, method = "BH")
projection_stats$fdr_residualized_group <- p.adjust(
  projection_stats$p_residualized_group,
  method = "BH"
)
projection_stats$pearson_fdr_all <- p.adjust(
  projection_stats$pearson_p_all,
  method = "BH"
)

projection_stats <- merge(
  projection_stats,
  coverage,
  by = "module",
  all.x = TRUE,
  sort = FALSE
)
projection_stats <- merge(
  projection_stats,
  local_check,
  by = "module",
  all.x = TRUE,
  sort = FALSE
)
projection_stats <- projection_stats[
  order(projection_stats$fdr_all, -abs(projection_stats$rho_all)),
]

wtsv(projection_stats, "05_MSG_Blood_projection_concordance.tsv")

# Disease association of Blood-projected scores.
traits_Blood <- module_disease_tests(
  score_matrix = blood_scores,
  disease_group = meta_msg$disease_group,
  score_family = "Blood_projected_MSG_modules"
)

wtsv(traits_Blood$omnibus, "05_Blood_projected_module_traits_omnibus.tsv")
wtsv(traits_Blood$pairwise, "05_Blood_projected_module_traits_pairwise.tsv")

msg(
  "Projected modules with paired MSG-Blood concordance FDR < 0.05 (all 97): %d / %d",
  sum(projection_stats$fdr_all < 0.05, na.rm = TRUE),
  nrow(projection_stats)
)
msg(
  "Projected modules with concordance FDR < 0.05 within SjD: %d / %d",
  sum(projection_stats$fdr_SjD < 0.05, na.rm = TRUE),
  nrow(projection_stats)
)

# Scatter plots, one page per 9 modules.
pdf(file.path(OUTDIR, "05_MSG_vs_Blood_projection_scatterplots.pdf"), width = 10, height = 10)
par(mfrow = c(3, 3), mar = c(4, 4, 3, 1))
for (mo in projection_stats$module) {
  x <- score_MSG_common[, mo]
  yb <- blood_scores[, mo]
  st <- projection_stats[projection_stats$module == mo, ]

  plot(
    x,
    yb,
    pch = 19,
    cex = 0.7,
    xlab = "MSG module eigengene",
    ylab = "Blood projected score",
    main = sprintf("%s: rho=%.2f, FDR=%.3g", mo, st$rho_all, st$fdr_all)
  )
  abline(lm(yb ~ x), lty = 2)
}
dev.off()

# ------------------------------------------------------------------------------
# 9. READY-TO-MERGE PATIENT REPRESENTATION FOR CLINICAL DATA
# ------------------------------------------------------------------------------

msg("\n=== BUILDING PATIENT-LEVEL REPRESENTATION ===")

MSG_z <- scale(score_MSG_common)
Blood_z <- scale(blood_scores)

colnames(MSG_z) <- paste0("MSG_ME_", colnames(MSG_z))
colnames(Blood_z) <- paste0("BloodProj_", colnames(Blood_z))

patient_repr <- data.frame(
  patient_id = meta_msg$patient_id,
  disease_group = meta_msg$disease_group,
  MSG_z,
  Blood_z,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Add cohort cluster.
patient_repr$cohort_cluster <- end_all$assignments$endotype[
  match(patient_repr$patient_id, end_all$assignments$patient_id)
]

# Add SjD-only preliminary endotype; non-SjD patients remain NA.
patient_repr$preliminary_SjD_endotype <- NA_character_
idx <- match(patient_repr$patient_id, end_sjd$assignments$patient_id)
has <- !is.na(idx)
patient_repr$preliminary_SjD_endotype[has] <- end_sjd$assignments$endotype[idx[has]]

wtsv(patient_repr, "patient_module_representation.tsv")

# ------------------------------------------------------------------------------
# 10. SAVE COMPACT RDS BUNDLE + RUN SUMMARY
# ------------------------------------------------------------------------------

bundle <- list(
  manifest = manifest,
  normalization = list(
    n_common_genes = length(common_genes),
    n_filtered_genes = nrow(E),
    method = "edgeR filterByExpr + TMM + limma-voom jointly across paired Blood/MSG"
  ),
  wgcna = list(
    soft_power = soft_power,
    soft_power_reason = power_reason,
    soft_threshold_table = soft_table,
    module_colors = module_colors,
    module_sizes = module_sizes_df,
    module_genes = module_genes,
    hub_genes_top20 = hub_genes,
    MSG_module_eigengenes = score_MSG,
    module_loadings = loadings
  ),
  module_traits_MSG = traits_MSG,
  endotypes = list(
    SjD = list(
      best_k = end_sjd$best_k,
      silhouette = end_sjd$silhouette,
      assignments = end_sjd$assignments
    ),
    all97 = list(
      best_k = end_all$best_k,
      silhouette = end_all$silhouette,
      assignments = end_all$assignments
    )
  ),
  projection = list(
    Blood_scores = blood_scores,
    coverage = coverage,
    concordance = projection_stats,
    Blood_module_traits = traits_Blood
  ),
  patient_representation = patient_repr
)

saveRDS(bundle, file.path(OUTDIR, "day2_02_05_bundle.rds"), compress = TRUE)

summary_lines <- c(
  "DAY 2 COMPLETE: WGCNA -> module traits -> preliminary endotypes -> MSG-to-Blood projection",
  "",
  sprintf("Paired patients: %d", nrow(manifest)),
  sprintf(
    "Disease groups: SjD=%d, nonSjD=%d, HV=%d",
    sum(manifest$disease_group == "SjD"),
    sum(manifest$disease_group == "nonSjD"),
    sum(manifest$disease_group == "HV")
  ),
  sprintf("Genes after filtering: %d", nrow(E)),
  sprintf("WGCNA soft power: %d", soft_power),
  sprintf("Biological MSG modules: %d", ncol(score_MSG)),
  sprintf(
    "MSG disease-associated modules (omnibus FDR<0.05): %d",
    sum(traits_MSG$omnibus$fdr < 0.05, na.rm = TRUE)
  ),
  sprintf("Preliminary SjD endotypes selected by silhouette: k=%d", end_sjd$best_k),
  paste0(
    "SjD endotype sizes: ",
    paste(names(table(end_sjd$assignments$endotype)),
          as.integer(table(end_sjd$assignments$endotype)),
          sep = "=", collapse = ", ")
  ),
  sprintf(
    "MSG->Blood concordant modules, all patients (FDR<0.05): %d",
    sum(projection_stats$fdr_all < 0.05, na.rm = TRUE)
  ),
  sprintf(
    "MSG->Blood concordant modules within SjD (FDR<0.05): %d",
    sum(projection_stats$fdr_SjD < 0.05, na.rm = TRUE)
  ),
  "",
  "Most important downstream file for clinical integration:",
  "  patient_module_representation.tsv",
  "",
  "Interpretation note:",
  "  Endotype labels are preliminary/exploratory. Validate stability and clinical associations before treating them as final biological classes."
)

writeLines(summary_lines, file.path(OUTDIR, "DAY2_SUMMARY.txt"))

capture.output(sessionInfo(), file = file.path(OUTDIR, "sessionInfo.txt"))

msg("\n============================================================")
msg("DAY 2 PIPELINE FINISHED")
msg("Output directory: %s", normalizePath(OUTDIR))
msg("Next clinical-integration input: patient_module_representation.tsv")
msg("============================================================")
