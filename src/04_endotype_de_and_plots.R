#!/usr/bin/env Rscript

# ==============================================================================
# 04_endotype_de_and_plots.R
#
# Exploratory gene-level limma-voom contrast of SjD MSG endotypes (E1 - E2).
# Generates a genuine volcano plot (log2FC versus -log10(BH FDR)), per-patient
# expression heatmap, and annotated tables with HGNC symbols/full gene names.
#
# CRITICAL INTERPRETATION:
# E1/E2 were clustered using these same MSG transcriptomes. The differential
# expression p-values/FDR and visualizations are CONDITIONAL, EXPLORATORY
# characterization of discovered groups, NOT independent endotype validation.
# No biological/clinical covariates are adjusted for in this two-group contrast.
#
# Prerequisite: Stage 02 (and preferably Stage 03) successfully completed.
# Run from anywhere: Rscript src/04_endotype_de_and_plots.R
# Uses INPUT_DIR / OUTPUT_DIR configured by config_paths.R.
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(20261002)

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0L) return(getwd())
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE))
}

CONFIG_FILE <- file.path(dirname(get_script_dir()), "config_paths.R")
if (!file.exists(CONFIG_FILE)) stop("Cannot find config_paths.R: ", CONFIG_FILE)
source(CONFIG_FILE)

OUTDIR <- file.path(OUTPUT_DIR, "04_endotype_de")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

FDR_CUTOFF <- 0.05
ABS_LOG2FC_CUTOFF <- 1
TOP_GENES_PER_DIRECTION <- 15L
VOLCANO_LABELS_PER_DIRECTION <- 8L

required <- c("edgeR", "limma", "ggplot2", "org.Hs.eg.db", "AnnotationDbi")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop(
    "Missing R packages: ", paste(missing, collapse = ", "),
    "\nBiocManager::install(c('edgeR','limma','org.Hs.eg.db','AnnotationDbi')); ",
    "install.packages('ggplot2')",
    call. = FALSE
  )
}
suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

logmsg <- function(...) {
  cat(sprintf(...), "\n")
  flush.console()
}

read_tsv <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path, call. = FALSE)
  read.delim(
    path, sep = "\t", header = TRUE, check.names = FALSE,
    stringsAsFactors = FALSE, quote = "", comment.char = ""
  )
}

write_tsv <- function(x, name) {
  write.table(
    x, file.path(OUTDIR, name), sep = "\t",
    quote = FALSE, row.names = FALSE, col.names = TRUE, na = "NA"
  )
}

read_count_matrix <- function(path) {
  d <- read_tsv(path)
  if (ncol(d) < 3L) stop("Count matrix requires a gene column and samples.")
  ids <- as.character(d[[1]])
  if (anyNA(ids) || any(!nzchar(ids))) stop("Missing gene IDs in count matrix.")
  x <- as.matrix(d[, -1, drop = FALSE])
  suppressWarnings(storage.mode(x) <- "numeric")
  if (any(!is.finite(x)) || any(x < 0)) stop("Count matrix contains invalid counts.")
  if (any(abs(x - round(x)) > 1e-5)) stop("Raw gene counts must be integers.")
  rownames(x) <- ids
  if (anyDuplicated(rownames(x))) {
    logmsg("Summing duplicate gene-ID rows before fitting.")
    x <- rowsum(x, group = rownames(x), reorder = FALSE)
  }
  round(x)
}

# One annotation per original gene ID; Ensembl remains the stable row key.
# Only mapped official HGNC symbols appear as labels in the figures.
annotate_ids <- function(original_ids) {
  cleaned <- sub("\\.[0-9]+$", "", original_ids)
  sym <- setNames(rep(NA_character_, length(original_ids)), original_ids)
  full <- setNames(rep(NA_character_, length(original_ids)), original_ids)

  is_ens <- grepl("^ENSG[0-9]+$", cleaned)
  if (any(is_ens)) {
    valid_ens <- intersect(
      unique(cleaned[is_ens]),
      AnnotationDbi::keys(org.Hs.eg.db, keytype = "ENSEMBL")
    )
    if (length(valid_ens)) {
      mapped_symbol <- suppressMessages(AnnotationDbi::mapIds(
        org.Hs.eg.db, keys = valid_ens, column = "SYMBOL",
        keytype = "ENSEMBL", multiVals = "first"
      ))
      mapped_name <- suppressMessages(AnnotationDbi::mapIds(
        org.Hs.eg.db, keys = valid_ens, column = "GENENAME",
        keytype = "ENSEMBL", multiVals = "first"
      ))
      sym[is_ens] <- unname(mapped_symbol[cleaned[is_ens]])
      full[is_ens] <- unname(mapped_name[cleaned[is_ens]])
    }
  }

  # Allow matrices already keyed by official gene symbol.
  if (any(!is_ens)) {
    valid_sym <- intersect(
      unique(cleaned[!is_ens]),
      AnnotationDbi::keys(org.Hs.eg.db, keytype = "SYMBOL")
    )
    hit <- !is_ens & cleaned %in% valid_sym
    sym[hit] <- cleaned[hit]
    if (length(valid_sym)) {
      names_by_symbol <- suppressMessages(AnnotationDbi::mapIds(
        org.Hs.eg.db, keys = valid_sym, column = "GENENAME",
        keytype = "SYMBOL", multiVals = "first"
      ))
      full[hit] <- unname(names_by_symbol[cleaned[hit]])
    }
  }

  data.frame(
    ensembl_gene_id = original_ids,
    gene_symbol = unname(sym),
    gene_name = unname(full),
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------------------
# 1. SjD endotype manifest and original MSG counts
# ------------------------------------------------------------------------------

day2_bundle <- file.path(WGCNA_OUT, "day2_02_05_bundle.rds")
if (!file.exists(day2_bundle)) {
  stop("Stage 02 bundle not found: ", day2_bundle, call. = FALSE)
}
bundle <- readRDS(day2_bundle)
manifest <- bundle$manifest
if (!all(c("patient_id", "msg_sample", "disease_group") %in% names(manifest))) {
  stop("Stage 02 manifest requires patient_id, msg_sample, disease_group.")
}

endotypes <- read_tsv(file.path(WGCNA_OUT, "04_preliminary_endotypes_SjD.tsv"))
if (!all(c("patient_id", "endotype") %in% names(endotypes))) {
  stop("Endotype table requires patient_id and endotype.")
}
endotypes <- endotypes[
  !is.na(endotypes$endotype) & endotypes$endotype %in% c("E1", "E2"),
  c("patient_id", "endotype"), drop = FALSE
]
if (anyDuplicated(endotypes$patient_id)) stop("Duplicate patient IDs in endotype table.")

m <- merge(
  manifest[, c("patient_id", "msg_sample", "disease_group")],
  endotypes, by = "patient_id", all = FALSE, sort = FALSE
)
if (any(m$disease_group != "SjD")) {
  stop("E1/E2 analysis must include SjD participants only.")
}
if (anyDuplicated(m$patient_id) || anyDuplicated(m$msg_sample)) {
  stop("Manifest includes duplicate E1/E2 patients or MSG samples.")
}
m$endotype <- factor(m$endotype, levels = c("E1", "E2"))
m <- m[order(m$endotype, m$patient_id), , drop = FALSE]
if (any(table(m$endotype) < 3L)) {
  stop("At least three participants per endotype are required.")
}

counts <- read_count_matrix(MSG_COUNTS)
if (!all(m$msg_sample %in% colnames(counts))) {
  stop("MSG count matrix is missing samples: ",
       paste(setdiff(m$msg_sample, colnames(counts)), collapse = ", "))
}
counts <- counts[, m$msg_sample, drop = FALSE]
colnames(counts) <- m$patient_id
if (nrow(counts) < 50L) stop("Too few genes in count matrix.")

logmsg("SjD MSG participants: E1=%d; E2=%d; input genes=%d",
       sum(m$endotype == "E1"), sum(m$endotype == "E2"), nrow(counts))

# ------------------------------------------------------------------------------
# 2. Exploratory, unadjusted E1 versus E2 differential expression
# ------------------------------------------------------------------------------

design <- model.matrix(~ 0 + endotype, data = m)
colnames(design) <- c("E1", "E2")
rownames(design) <- as.character(m$patient_id)
stopifnot(identical(rownames(design), colnames(counts)))

y <- edgeR::DGEList(counts = counts)
keep <- edgeR::filterByExpr(y, design = design)
if (sum(keep) < 50L) stop("Fewer than 50 genes survived filterByExpr.")
y <- y[keep, , keep.lib.sizes = FALSE]
y <- edgeR::calcNormFactors(y, method = "TMM")

v <- limma::voom(y, design, plot = FALSE)
contrast <- limma::makeContrasts(E1_vs_E2 = E1 - E2, levels = design)
fit <- limma::lmFit(v, design)
fit <- limma::contrasts.fit(fit, contrast)
fit <- limma::eBayes(fit, trend = FALSE)

de <- limma::topTable(fit, coef = "E1_vs_E2", number = Inf, sort.by = "none")
de$ensembl_gene_id <- rownames(de)
de <- de[, c("ensembl_gene_id", "logFC", "AveExpr", "t",
             "P.Value", "adj.P.Val", "B"), drop = FALSE]
anno <- annotate_ids(de$ensembl_gene_id)
stopifnot(identical(de$ensembl_gene_id, anno$ensembl_gene_id))

de$gene_symbol <- anno$gene_symbol
de$gene_name <- anno$gene_name
idx_e1 <- m$endotype == "E1"
idx_e2 <- m$endotype == "E2"
de$mean_logCPM_E1 <- rowMeans(v$E[de$ensembl_gene_id, idx_e1, drop = FALSE])
de$mean_logCPM_E2 <- rowMeans(v$E[de$ensembl_gene_id, idx_e2, drop = FALSE])
de$delta_mean_logCPM_E1_minus_E2 <- de$mean_logCPM_E1 - de$mean_logCPM_E2
de$direction <- "not_significant"
de$direction[!is.na(de$adj.P.Val) & de$adj.P.Val < FDR_CUTOFF &
               de$logFC >= ABS_LOG2FC_CUTOFF] <- "E1_high"
de$direction[!is.na(de$adj.P.Val) & de$adj.P.Val < FDR_CUTOFF &
               de$logFC <= -ABS_LOG2FC_CUTOFF] <- "E2_high"

de <- de[order(de$adj.P.Val, -abs(de$logFC)), , drop = FALSE]
de <- de[, c("gene_symbol", "gene_name", "ensembl_gene_id", "logFC",
             "mean_logCPM_E1", "mean_logCPM_E2",
             "delta_mean_logCPM_E1_minus_E2", "AveExpr",
             "t", "P.Value", "adj.P.Val", "B", "direction")]
write_tsv(de, "06_E1_vs_E2_differential_expression.tsv")
write_tsv(de[order(-de$logFC), , drop = FALSE],
          "06_E1_vs_E2_gene_ranking_annotated.tsv")
write_tsv(head(de[order(-de$logFC), , drop = FALSE], 100),
          "06_E1_high_top100_genes_annotated.tsv")
write_tsv(head(de[order(de$logFC), , drop = FALSE], 100),
          "06_E2_high_top100_genes_annotated.tsv")
write_tsv(de[, c("ensembl_gene_id", "gene_symbol", "gene_name")],
          "06_gene_id_to_HGNC_map.tsv")

# ------------------------------------------------------------------------------
# 3. True volcano plot: x = log2FC, y = -log10(BH FDR)
# ------------------------------------------------------------------------------

de$minus_log10_FDR <- -log10(pmax(de$adj.P.Val, .Machine$double.xmin))
de$direction <- factor(
  de$direction, levels = c("not_significant", "E2_high", "E1_high")
)

labeled <- de[
  !is.na(de$gene_symbol) & nzchar(de$gene_symbol) &
    de$direction != "not_significant", , drop = FALSE
]
labeled <- labeled[order(labeled$adj.P.Val, -abs(labeled$logFC)), , drop = FALSE]
labeled <- labeled[!duplicated(labeled$gene_symbol), , drop = FALSE]
label_top <- do.call(rbind, lapply(c("E1_high", "E2_high"), function(direction) {
  head(labeled[as.character(labeled$direction) == direction, , drop = FALSE],
       VOLCANO_LABELS_PER_DIRECTION)
}))
if (is.null(label_top)) label_top <- labeled[FALSE, , drop = FALSE]

palette <- c(
  "not_significant" = "#B5BDC4",
  "E2_high" = "#C75867",
  "E1_high" = "#188C89"
)

p_volcano <- ggplot(de, aes(x = logFC, y = minus_log10_FDR)) +
  geom_point(aes(color = direction), alpha = 0.72, size = 1.2) +
  geom_vline(xintercept = c(-ABS_LOG2FC_CUTOFF, ABS_LOG2FC_CUTOFF),
             linetype = "dashed", linewidth = 0.35, color = "grey45") +
  geom_hline(yintercept = -log10(FDR_CUTOFF),
             linetype = "dashed", linewidth = 0.35, color = "grey45") +
  scale_color_manual(values = palette, drop = FALSE) +
  labs(
    title = "Exploratory E1 vs E2 differential expression in MSG",
    subtitle = sprintf("SjD only: E1 n=%d, E2 n=%d | TMM + limma-voom | log2FC=E1−E2",
                       sum(idx_e1), sum(idx_e2)),
    x = "log2 fold change (E1 − E2)",
    y = expression(-log[10]("BH-adjusted p-value")),
    color = "Gene status",
    caption = paste(
      "Dashed thresholds: BH FDR < 0.05 and |log2FC| >= 1.",
      "Clusters derived from these data: exploratory characterization, not independent validation."
    )
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom",
        plot.caption = element_text(size = 8.5))

if (nrow(label_top)) {
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p_volcano <- p_volcano +
      ggrepel::geom_text_repel(
        data = label_top, aes(label = gene_symbol),
        size = 3.2, max.overlaps = Inf, min.segment.length = 0,
        box.padding = 0.38, point.padding = 0.15,
        show.legend = FALSE, seed = 20261002
      )
  } else {
    p_volcano <- p_volcano +
      geom_text(
        data = label_top, aes(label = gene_symbol),
        size = 3, check_overlap = TRUE, vjust = -0.8,
        show.legend = FALSE
      )
  }
}

ggsave(file.path(OUTDIR, "06_E1_vs_E2_volcano.png"),
       p_volcano, width = 11, height = 7, dpi = 320)
ggsave(file.path(OUTDIR, "06_E1_vs_E2_volcano.pdf"),
       p_volcano, width = 11, height = 7)

# ------------------------------------------------------------------------------
# 4. Genuine PER-PATIENT heatmap of E1/E2-distinguishing genes
# ------------------------------------------------------------------------------

# Figures never display ENSG IDs. Unmapped genes stay in the complete TSV, but
# are not selected as heatmap row labels.
plot_de <- de[
  !is.na(de$gene_symbol) & nzchar(de$gene_symbol) &
    is.finite(de$logFC) & is.finite(de$adj.P.Val),
  , drop = FALSE
]
plot_de <- plot_de[order(plot_de$adj.P.Val, -abs(plot_de$logFC)), , drop = FALSE]
plot_de <- plot_de[!duplicated(plot_de$gene_symbol), , drop = FALSE]

choose_genes <- function(side) {
  if (side == "E1") {
    candidates <- plot_de[plot_de$logFC > 0, , drop = FALSE]
    significant <- candidates[candidates$direction == "E1_high", , drop = FALSE]
  } else {
    candidates <- plot_de[plot_de$logFC < 0, , drop = FALSE]
    significant <- candidates[candidates$direction == "E2_high", , drop = FALSE]
  }
  # Prefer genes passing BOTH thresholds; when fewer than 15 pass, supplement
  # by best FDR/absolute effect and flag this explicitly in the figure caption.
  significant <- head(significant, TOP_GENES_PER_DIRECTION)
  remainder <- candidates[
    !candidates$ensembl_gene_id %in% significant$ensembl_gene_id,
    , drop = FALSE
  ]
  rbind(
    significant,
    head(remainder, TOP_GENES_PER_DIRECTION - nrow(significant))
  )
}

h_e1 <- choose_genes("E1")
h_e2 <- choose_genes("E2")
heat_genes <- rbind(h_e2, h_e1)
if (nrow(heat_genes) < 2L) stop("Not enough HGNC-mapped genes for heatmap.")
stopifnot(!anyDuplicated(heat_genes$gene_symbol))

log_cpm <- edgeR::cpm(y, log = TRUE, prior.count = 2)
h <- log_cpm[heat_genes$ensembl_gene_id, , drop = FALSE]
rownames(h) <- heat_genes$gene_symbol
if (!identical(colnames(h), m$patient_id)) {
  stop("Heatmap sample order does not match endotype manifest.")
}

# Z-score ACROSS PATIENTS for each gene (not between group means).
h_z <- t(scale(t(h)))
h_z[!is.finite(h_z)] <- 0
h_z <- pmax(-2.5, pmin(2.5, h_z))

labels <- paste0(
  as.character(m$endotype), "_",
  ave(seq_len(nrow(m)), m$endotype, FUN = seq_along)
)
colnames(h_z) <- labels
write_tsv(
  data.frame(gene_symbol = rownames(h_z), h_z, check.names = FALSE),
  "06_E1_vs_E2_per_patient_heatmap_zscores.tsv"
)

long <- data.frame(
  gene_symbol = rep(rownames(h_z), times = ncol(h_z)),
  sample_label = rep(colnames(h_z), each = nrow(h_z)),
  z = as.vector(h_z),
  endotype = rep(as.character(m$endotype), each = nrow(h_z)),
  stringsAsFactors = FALSE
)
long$sample_label <- factor(long$sample_label, levels = labels)
long$gene_symbol <- factor(
  long$gene_symbol, levels = rev(rownames(h_z))
)
long$endotype <- factor(long$endotype, levels = c("E1", "E2"))
supplemented <- (nrow(h_e1[h_e1$direction == "E1_high", , drop = FALSE]) <
                 TOP_GENES_PER_DIRECTION) ||
                (nrow(h_e2[h_e2$direction == "E2_high", , drop = FALSE]) <
                 TOP_GENES_PER_DIRECTION)
selection_note <- if (supplemented) {
  "Top significant genes supplemented by descriptive best-ranked genes where fewer than 15 pass both thresholds."
} else {
  "Genes selected at BH FDR < 0.05 and |log2FC| >= 1."
}

p_heat <- ggplot(long, aes(x = sample_label, y = gene_symbol, fill = z)) +
  geom_tile() +
  facet_grid(. ~ endotype, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(
    low = "#315899", mid = "white", high = "#BE4B56",
    midpoint = 0, limits = c(-2.5, 2.5), name = "Gene-wise\nz-score"
  ) +
  labs(
    title = "E1- and E2-high genes across individual MSG samples",
    subtitle = sprintf("Top %d genes per direction | E1 n=%d; E2 n=%d",
                       TOP_GENES_PER_DIRECTION, sum(idx_e1), sum(idx_e2)),
    x = "Individual SjD MSG samples (anonymized labels)",
    y = "HGNC gene symbol",
    caption = paste(
      selection_note,
      "Colors show relative expression per gene, not log2 fold change.",
      "Endotype contrasts are exploratory and not independently validated.",
      sep = "\n"
    )
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold"),
    plot.caption = element_text(size = 8)
  )

ggsave(file.path(OUTDIR, "06_E1_vs_E2_per_patient_heatmap.png"),
       p_heat, width = 13, height = 9, dpi = 320)
ggsave(file.path(OUTDIR, "06_E1_vs_E2_per_patient_heatmap.pdf"),
       p_heat, width = 13, height = 9)

lines <- c(
  "STAGE 04 — EXPLORATORY E1 VS E2 GENE-LEVEL CHARACTERIZATION",
  sprintf("MSG samples: E1=%d; E2=%d", sum(idx_e1), sum(idx_e2)),
  sprintf("Genes after design-aware filterByExpr: %d", nrow(de)),
  sprintf("Genes mapped to HGNC symbols: %d", sum(!is.na(de$gene_symbol))),
  sprintf("E1-high genes: FDR<%.2f and log2FC>=%.1f: %d",
          FDR_CUTOFF, ABS_LOG2FC_CUTOFF, sum(de$direction == "E1_high")),
  sprintf("E2-high genes: FDR<%.2f and log2FC<=-%.1f: %d",
          FDR_CUTOFF, ABS_LOG2FC_CUTOFF, sum(de$direction == "E2_high")),
  "Contrast: E1 - E2; model includes endotype only (no covariate adjustment).",
  "Endotypes were derived from these same transcriptomes. P-values and FDR",
  "are conditional/exploratory; they are NOT independent validation.",
  "The volcano uses BH FDR on the y-axis; its labels use HGNC symbols.",
  "The heatmap uses logCPM z-scores across individual patients (one row per gene).",
  paste0("Heatmap selection note: ", selection_note)
)
writeLines(lines, file.path(OUTDIR, "04_ENDOTYPE_DE_SUMMARY.txt"))
capture.output(sessionInfo(), file = file.path(OUTDIR, "sessionInfo.txt"))

logmsg("Complete. Results: %s", OUTDIR)
logmsg("Volcano: 06_E1_vs_E2_volcano.png/.pdf")
logmsg("Heatmap: 06_E1_vs_E2_per_patient_heatmap.png/.pdf")
logmsg("Annotated DE: 06_E1_vs_E2_differential_expression.tsv")
