#!/usr/bin/env Rscript

# ==============================================================================
# 03_functional_enrichment.R
#
# Functional interpretation of:
#   A) all 23 MSG WGCNA modules
#   B) the disease-associated modules (expected n = 8)
#   C) the MSG -> Blood projected/concordant module (expected n = 1)
#   D) SjD E1 vs E2:
#        - module drivers
#        - descriptive gene-level GSEA
#
# IMPORTANT:
# E1/E2 were learned from the same molecular data. Therefore E1/E2 enrichment is
# descriptive characterization, NOT independent validation of the clusters.
#
# Run:
#   Rscript src/03_functional_enrichment.R
#
# Paths are obtained from config_paths.R.
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(20260929)

# ------------------------------------------------------------------------------
# 0. CONFIGURATION
# ------------------------------------------------------------------------------

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0) return(getwd())
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE))
}

CONFIG_FILE <- file.path(dirname(get_script_dir()), "config_paths.R")
if (!file.exists(CONFIG_FILE)) {
  stop("config_paths.R not found at: ", CONFIG_FILE, call. = FALSE)
}
source(CONFIG_FILE)
DAY2_DIR <- WGCNA_OUT
OUTDIR <- ENRICHMENT_OUT

cat("Using config:", CONFIG_FILE, "\n")
cat("INPUT_DIR:", INPUT_DIR, "\n")
cat("OUTPUT_DIR:", OUTPUT_DIR, "\n")

FDR_CUTOFF <- 0.05
TOP_TERMS_PER_MODULE <- 10
TOP_MODULES_DISEASE <- 3
TOP_ENDOTYPE_DRIVER_MODULES <- 5

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. PACKAGES
# ------------------------------------------------------------------------------

required <- c(
  "clusterProfiler",
  "org.Hs.eg.db",
  "ReactomePA",
  "AnnotationDbi",
  "msigdbr",
  "edgeR",
  "ggplot2"
)

missing <- required[
  !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing) > 0) {
  stop(
    "\nMissing packages: ",
    paste(missing, collapse = ", "),
    "\n\nInstall them before rerunning.\n",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(AnnotationDbi)
  library(msigdbr)
  library(edgeR)
  library(ggplot2)
})

# ------------------------------------------------------------------------------
# 2. HELPERS
# ------------------------------------------------------------------------------

msg <- function(...) {
  cat(sprintf(...), "\n")
  flush.console()
}

read_tsv <- function(file) {
  read.delim(
    file,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )
}

wtsv <- function(x, name) {
  write.table(
    x,
    file.path(OUTDIR, name),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )
}

need_file <- function(name) {
  f <- file.path(DAY2_DIR, name)
  if (!file.exists(f)) {
    stop("Missing required Day-2 output: ", f, call. = FALSE)
  }
  f
}

resolve_recursive <- function(path) {
  if (file.exists(path)) return(normalizePath(path))

  hits <- list.files(
    ".",
    pattern = paste0("^", basename(path), "$"),
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(hits) == 1) return(normalizePath(hits))

  if (length(hits) > 1) {
    stop(
      "Multiple files named ", basename(path), " found:\n",
      paste(hits, collapse = "\n"),
      "\nSet INPUT_DIR in config_paths.R to the exact input directory.",
      call. = FALSE
    )
  }

  stop(
    "Could not find ", path,
    ". Set INPUT_DIR in config_paths.R to the exact input directory.",
    call. = FALSE
  )
}

read_count_matrix <- function(file) {
  d <- read.delim(
    file,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )

  gene <- as.character(d[[1]])
  mat <- as.matrix(d[, -1, drop = FALSE])
  suppressWarnings(storage.mode(mat) <- "numeric")

  if (anyNA(mat)) {
    stop("Non-numeric values found in MSG count matrix.", call. = FALSE)
  }

  rownames(mat) <- gene

  if (anyDuplicated(rownames(mat))) {
    mat <- rowsum(mat, group = rownames(mat), reorder = FALSE)
  }

  mat
}

strip_ensembl_version <- function(x) {
  sub("\\.[0-9]+$", "", x)
}

detect_gene_keytype <- function(x) {
  y <- strip_ensembl_version(x)
  frac_ens <- mean(grepl("^ENSG[0-9]+", y))

  if (is.finite(frac_ens) && frac_ens > 0.5) {
    return("ENSEMBL")
  }

  "SYMBOL"
}

build_gene_map <- function(genes) {
  original <- unique(as.character(genes))
  cleaned <- strip_ensembl_version(original)
  keytype <- detect_gene_keytype(cleaned)

  msg("Detected gene ID type: %s", keytype)

  mapped <- suppressMessages(
    AnnotationDbi::select(
      org.Hs.eg.db,
      keys = unique(cleaned),
      keytype = keytype,
      columns = c("SYMBOL", "ENTREZID")
    )
  )

  names(mapped)[names(mapped) == keytype] <- "clean_gene"

  out <- merge(
    data.frame(
      gene = original,
      clean_gene = cleaned,
      stringsAsFactors = FALSE
    ),
    mapped,
    by = "clean_gene",
    all.x = TRUE
  )

  out <- out[
    !duplicated(out[, c("gene", "SYMBOL", "ENTREZID")]),
    ,
    drop = FALSE
  ]

  msg(
    "Mapped %.1f%% of WGCNA genes to ENTREZID",
    100 * mean(!is.na(out$ENTREZID))
  )

  out
}

msigdbr_safe <- function(collection, subcollection = NULL) {
  f <- msigdbr::msigdbr
  args <- list(species = "Homo sapiens")

  if ("collection" %in% names(formals(f))) {
    args$collection <- collection
  } else {
    args$category <- collection
  }

  if (!is.null(subcollection)) {
    if ("subcollection" %in% names(formals(f))) {
      args$subcollection <- subcollection
    } else {
      args$subcategory <- subcollection
    }
  }

  do.call(f, args)
}

extract_enrichment <- function(obj, module, database) {
  if (is.null(obj)) return(data.frame())

  d <- tryCatch(
    as.data.frame(obj),
    error = function(e) data.frame()
  )

  if (nrow(d) == 0) return(data.frame())

  d$module <- module
  d$database <- database

  keep <- intersect(
    c(
      "module", "database", "ID", "Description",
      "GeneRatio", "BgRatio",
      "pvalue", "p.adjust", "qvalue",
      "geneID", "Count"
    ),
    names(d)
  )

  d[, keep, drop = FALSE]
}

run_module_ora <- function(
    module,
    genes_entrez,
    genes_symbol,
    universe_entrez,
    universe_symbol,
    hallmark_term2gene) {

  out <- list()

  if (length(unique(genes_entrez)) >= 5) {
    ego <- suppressMessages(
      clusterProfiler::enrichGO(
        gene = unique(genes_entrez),
        universe = unique(universe_entrez),
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = "BP",
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        pAdjustMethod = "BH",
        minGSSize = 10,
        maxGSSize = 500,
        readable = TRUE
      )
    )

    er <- suppressMessages(
      ReactomePA::enrichPathway(
        gene = unique(genes_entrez),
        universe = unique(universe_entrez),
        organism = "human",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        minGSSize = 10,
        maxGSSize = 500,
        readable = TRUE
      )
    )

    out$GO_BP <- extract_enrichment(
      ego,
      module,
      "GO_BP"
    )

    out$Reactome <- extract_enrichment(
      er,
      module,
      "Reactome"
    )
  }

  if (length(unique(genes_symbol)) >= 5) {
    eh <- suppressMessages(
      clusterProfiler::enricher(
        gene = unique(genes_symbol),
        universe = unique(universe_symbol),
        TERM2GENE = hallmark_term2gene,
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        minGSSize = 5,
        maxGSSize = 500
      )
    )

    out$Hallmark <- extract_enrichment(
      eh,
      module,
      "Hallmark"
    )
  }

  do.call(rbind, out)
}

add_global_fdr <- function(d) {
  if (nrow(d) == 0) return(d)

  d$global_fdr <- NA_real_

  for (db in unique(d$database)) {
    idx <- which(d$database == db)
    d$global_fdr[idx] <- p.adjust(
      d$pvalue[idx],
      method = "BH"
    )
  }

  d
}

top_terms_by_module <- function(d, n = 10) {
  if (nrow(d) == 0) return(d)

  pieces <- split(
    d,
    interaction(d$module, d$database, drop = TRUE)
  )

  out <- lapply(pieces, function(x) {
    x <- x[order(x$p.adjust, x$pvalue), , drop = FALSE]
    head(x, n)
  })

  do.call(rbind, out)
}

collapse_rank <- function(df, id_col, stat_col) {
  d <- df[
    is.finite(df[[stat_col]]) &
      !is.na(df[[id_col]]) &
      df[[id_col]] != "",
    c(id_col, stat_col),
    drop = FALSE
  ]

  names(d) <- c("id", "stat")

  d <- d[
    order(abs(d$stat), decreasing = TRUE),
    ,
    drop = FALSE
  ]

  d <- d[!duplicated(d$id), , drop = FALSE]
  stats <- d$stat
  names(stats) <- d$id
  sort(stats, decreasing = TRUE)
}

extract_gsea <- function(obj, database) {
  if (is.null(obj)) return(data.frame())

  d <- tryCatch(
    as.data.frame(obj),
    error = function(e) data.frame()
  )

  if (nrow(d) == 0) return(d)

  d$database <- database
  d$direction <- ifelse(
    d$NES > 0,
    "E1_high",
    "E2_high"
  )

  keep <- intersect(
    c(
      "database", "ID", "Description",
      "setSize", "enrichmentScore", "NES",
      "pvalue", "p.adjust", "qvalue",
      "rank", "leading_edge", "core_enrichment",
      "direction"
    ),
    names(d)
  )

  d[, keep, drop = FALSE]
}

hedges_g <- function(x1, x2) {
  n1 <- length(x1)
  n2 <- length(x2)

  if (n1 < 2 || n2 < 2) return(NA_real_)

  s1 <- sd(x1)
  s2 <- sd(x2)

  sp <- sqrt(
    ((n1 - 1) * s1^2 + (n2 - 1) * s2^2) /
      (n1 + n2 - 2)
  )

  if (!is.finite(sp) || sp == 0) return(NA_real_)

  d <- (mean(x1) - mean(x2)) / sp
  J <- 1 - (3 / (4 * (n1 + n2) - 9))

  J * d
}

# ------------------------------------------------------------------------------
# 3. LOAD DAY-2 OUTPUTS
# ------------------------------------------------------------------------------

msg("=== Stage 06: module/endotype enrichment ===")

module_genes <- read_tsv(
  need_file("02_module_genes.tsv")
)

hub_genes <- read_tsv(
  need_file("02_hub_genes_top20.tsv")
)

traits_omnibus <- read_tsv(
  need_file("03_module_traits_omnibus.tsv")
)

traits_pairwise <- read_tsv(
  need_file("03_module_traits_pairwise.tsv")
)

endotypes <- read_tsv(
  need_file("04_preliminary_endotypes_SjD.tsv")
)

projection <- read_tsv(
  need_file("05_MSG_Blood_projection_concordance.tsv")
)

patient_repr <- read_tsv(
  need_file("patient_module_representation.tsv")
)

bundle <- readRDS(
  need_file("day2_02_05_bundle.rds")
)

all_modules <- sort(
  setdiff(
    unique(module_genes$module),
    "grey"
  )
)

disease_modules <- traits_omnibus$module[
  is.finite(traits_omnibus$fdr) &
    traits_omnibus$fdr < FDR_CUTOFF
]

disease_modules <- unique(disease_modules)

projection_modules <- projection$module[
  is.finite(projection$fdr_all) &
    projection$fdr_all < FDR_CUTOFF &
    is.finite(projection$fdr_SjD) &
    projection$fdr_SjD < FDR_CUTOFF
]

projection_modules <- unique(projection_modules)

if (length(projection_modules) == 0) {
  warning(
    "No module passed FDR<0.05 in both all-patient and SjD projection. ",
    "Using the module with the lowest fdr_all as the projection candidate."
  )

  projection_modules <- projection$module[
    which.min(projection$fdr_all)
  ]
}

msg("Biological WGCNA modules: %d", length(all_modules))
msg("Disease-associated modules FDR<0.05: %d", length(disease_modules))
msg(
  "MSG->Blood concordant module(s) FDR<0.05 in both analyses: %s",
  paste(projection_modules, collapse = ", ")
)

# Transparent top-3 disease modules:
top3_disease <- traits_omnibus[
  traits_omnibus$module %in% disease_modules,
  ,
  drop = FALSE
]

top3_disease <- top3_disease[
  order(top3_disease$fdr, top3_disease$p),
  ,
  drop = FALSE
]

top3_disease <- head(
  top3_disease,
  TOP_MODULES_DISEASE
)

wtsv(
  top3_disease,
  "06_top3_disease_modules.tsv"
)

# ------------------------------------------------------------------------------
# 4. MAP GENES
# ------------------------------------------------------------------------------

gene_map <- build_gene_map(
  module_genes$gene
)

mapped_modules <- merge(
  module_genes,
  gene_map,
  by = "gene",
  all.x = TRUE
)

universe_entrez <- unique(
  na.omit(mapped_modules$ENTREZID)
)

universe_symbol <- unique(
  na.omit(mapped_modules$SYMBOL)
)

hallmark <- msigdbr_safe("H")

hallmark_term2gene <- unique(
  hallmark[, c("gs_name", "gene_symbol")]
)

names(hallmark_term2gene) <- c(
  "term",
  "gene"
)

# ------------------------------------------------------------------------------
# 5. ORA — ALL 23 MODULES
# ------------------------------------------------------------------------------

msg("\nRunning ORA for all modules...")

all_ora <- lapply(all_modules, function(mo) {

  d <- mapped_modules[
    mapped_modules$module == mo,
    ,
    drop = FALSE
  ]

  genes_entrez <- unique(
    na.omit(d$ENTREZID)
  )

  genes_symbol <- unique(
    na.omit(d$SYMBOL)
  )

  msg(
    "  %s: %d genes (%d mapped Entrez)",
    mo,
    nrow(d),
    length(genes_entrez)
  )

  run_module_ora(
    module = mo,
    genes_entrez = genes_entrez,
    genes_symbol = genes_symbol,
    universe_entrez = universe_entrez,
    universe_symbol = universe_symbol,
    hallmark_term2gene = hallmark_term2gene
  )
})

all_ora <- do.call(rbind, all_ora)

if (is.null(all_ora)) {
  all_ora <- data.frame()
}

all_ora <- add_global_fdr(
  all_ora
)

wtsv(
  all_ora,
  "06_all23_modules_enrichment_full.tsv"
)

all_ora_top <- top_terms_by_module(
  all_ora,
  TOP_TERMS_PER_MODULE
)

wtsv(
  all_ora_top,
  "06_all23_modules_enrichment_top10.tsv"
)

# ------------------------------------------------------------------------------
# 6. ENRICHMENT — 8 DISEASE-ASSOCIATED MODULES
# ------------------------------------------------------------------------------

disease_ora <- all_ora[
  all_ora$module %in% disease_modules,
  ,
  drop = FALSE
]

wtsv(
  disease_ora,
  "06_disease_associated_modules_enrichment_full.tsv"
)

disease_ora_top <- top_terms_by_module(
  disease_ora,
  TOP_TERMS_PER_MODULE
)

wtsv(
  disease_ora_top,
  "06_disease_associated_modules_enrichment_top10.tsv"
)

# Top-3 module enrichment
top3_names <- top3_disease$module

top3_ora <- all_ora[
  all_ora$module %in% top3_names,
  ,
  drop = FALSE
]

wtsv(
  top3_ora,
  "06_top3_disease_modules_enrichment.tsv"
)

# ------------------------------------------------------------------------------
# 7. ENRICHMENT — MSG -> BLOOD PROJECTABLE MODULE
# ------------------------------------------------------------------------------

projectable_ora <- all_ora[
  all_ora$module %in% projection_modules,
  ,
  drop = FALSE
]

wtsv(
  projectable_ora,
  "06_projectable_module_enrichment.tsv"
)

projectable_hubs <- hub_genes[
  hub_genes$module %in% projection_modules,
  ,
  drop = FALSE
]

wtsv(
  projectable_hubs,
  "06_projectable_module_hub_genes.tsv"
)

projectable_stats <- projection[
  projection$module %in% projection_modules,
  ,
  drop = FALSE
]

wtsv(
  projectable_stats,
  "06_projectable_module_projection_stats.tsv"
)

# ------------------------------------------------------------------------------
# 8. MODULE-LEVEL DRIVERS OF E1 VS E2
# ------------------------------------------------------------------------------

msg("\nCharacterizing E1 vs E2 module-space differences...")

if (!"preliminary_SjD_endotype" %in% names(patient_repr)) {
  stop(
    "patient_module_representation.tsv lacks preliminary_SjD_endotype.",
    call. = FALSE
  )
}

msg_cols <- grep(
  "^MSG_ME_",
  names(patient_repr),
  value = TRUE
)

sjd <- patient_repr[
  patient_repr$disease_group == "SjD" &
    patient_repr$preliminary_SjD_endotype %in% c("E1", "E2"),
  ,
  drop = FALSE
]

if (nrow(sjd) == 0) {
  stop("No SjD E1/E2 patients found.", call. = FALSE)
}

msg(
  "Endotype counts: %s",
  paste(
    names(table(sjd$preliminary_SjD_endotype)),
    table(sjd$preliminary_SjD_endotype),
    sep = "=",
    collapse = ", "
  )
)

driver_table <- do.call(
  rbind,
  lapply(msg_cols, function(col) {

    mo <- sub("^MSG_ME_", "", col)

    e1 <- sjd[
      sjd$preliminary_SjD_endotype == "E1",
      col
    ]

    e2 <- sjd[
      sjd$preliminary_SjD_endotype == "E2",
      col
    ]

    data.frame(
      module = mo,
      n_E1 = length(e1),
      n_E2 = length(e2),
      mean_E1 = mean(e1, na.rm = TRUE),
      mean_E2 = mean(e2, na.rm = TRUE),
      delta_E1_minus_E2 = mean(e1, na.rm = TRUE) -
        mean(e2, na.rm = TRUE),
      hedges_g_descriptive = hedges_g(
        e1[is.finite(e1)],
        e2[is.finite(e2)]
      ),
      stringsAsFactors = FALSE
    )
  })
)

driver_table$abs_delta <- abs(
  driver_table$delta_E1_minus_E2
)

driver_table <- driver_table[
  order(
    -driver_table$abs_delta
  ),
  ,
  drop = FALSE
]

wtsv(
  driver_table,
  "06_E1_vs_E2_module_drivers.tsv"
)

top_driver_modules <- head(
  driver_table$module,
  TOP_ENDOTYPE_DRIVER_MODULES
)

driver_ora <- all_ora[
  all_ora$module %in% top_driver_modules,
  ,
  drop = FALSE
]

wtsv(
  driver_ora,
  "06_E1_vs_E2_top_driver_modules_enrichment.tsv"
)

# Heatmap-like matrix for external plotting if desired
driver_matrix <- sjd[
  ,
  c(
    "patient_id",
    "preliminary_SjD_endotype",
    msg_cols
  ),
  drop = FALSE
]

ord <- order(
  driver_matrix$preliminary_SjD_endotype
)

driver_matrix <- driver_matrix[
  ord,
  ,
  drop = FALSE
]

wtsv(
  driver_matrix,
  "06_E1_E2_module_matrix.tsv"
)

# ------------------------------------------------------------------------------
# 9. DESCRIPTIVE GENE-LEVEL E1 VS E2 RANK
# ------------------------------------------------------------------------------

msg("\nBuilding descriptive E1 vs E2 gene ranking from MSG counts...")

counts_file <- resolve_recursive(
  MSG_COUNTS
)

counts_msg <- read_count_matrix(
  counts_file
)

manifest <- bundle$manifest

if (
  !"patient_id" %in% names(manifest) ||
  !"msg_sample" %in% names(manifest)
) {
  stop(
    "day2 bundle manifest does not contain patient_id and msg_sample.",
    call. = FALSE
  )
}

endo_map <- endotypes[
  endotypes$endotype %in% c("E1", "E2"),
  c("patient_id", "endotype"),
  drop = FALSE
]

m <- merge(
  manifest[, c("patient_id", "msg_sample")],
  endo_map,
  by = "patient_id",
  all = FALSE
)

missing_samples <- setdiff(
  m$msg_sample,
  colnames(counts_msg)
)

if (length(missing_samples) > 0) {
  stop(
    "MSG samples missing from count matrix: ",
    paste(head(missing_samples, 10), collapse = ", "),
    call. = FALSE
  )
}

counts_e <- counts_msg[
  ,
  m$msg_sample,
  drop = FALSE
]

colnames(counts_e) <- m$patient_id

# Restrict to genes that entered WGCNA/background
genes_wgcna <- intersect(
  rownames(counts_e),
  unique(module_genes$gene)
)

counts_e <- counts_e[
  genes_wgcna,
  ,
  drop = FALSE
]

# Filter only obvious unusable rows, then TMM logCPM.
keep_gene <- rowSums(counts_e > 5) >= 2
counts_e <- counts_e[
  keep_gene,
  ,
  drop = FALSE
]

dge <- edgeR::DGEList(
  counts = round(counts_e)
)

dge <- edgeR::calcNormFactors(
  dge,
  method = "TMM"
)

logcpm <- edgeR::cpm(
  dge,
  log = TRUE,
  prior.count = 2
)

grp <- m$endotype[
  match(
    colnames(logcpm),
    m$patient_id
  )
]

e1_idx <- grp == "E1"
e2_idx <- grp == "E2"

mean_E1 <- rowMeans(
  logcpm[, e1_idx, drop = FALSE]
)

mean_E2 <- rowMeans(
  logcpm[, e2_idx, drop = FALSE]
)

rank_gene <- mean_E1 - mean_E2

rank_table <- data.frame(
  gene = names(rank_gene),
  mean_logCPM_E1 = mean_E1,
  mean_logCPM_E2 = mean_E2,
  delta_E1_minus_E2 = rank_gene,
  stringsAsFactors = FALSE
)

rank_table <- rank_table[
  order(
    rank_table$delta_E1_minus_E2,
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]

wtsv(
  rank_table,
  "06_E1_vs_E2_gene_ranking.tsv"
)

wtsv(
  head(rank_table, 100),
  "06_E1_high_top100_genes.tsv"
)

wtsv(
  head(
    rank_table[
      order(
        rank_table$delta_E1_minus_E2,
        decreasing = FALSE
      ),
      ,
      drop = FALSE
    ],
    100
  ),
  "06_E2_high_top100_genes.tsv"
)

# ------------------------------------------------------------------------------
# 10. MAP E1/E2 RANK TO SYMBOL AND ENTREZ
# ------------------------------------------------------------------------------

rank_map <- merge(
  rank_table[, c("gene", "delta_E1_minus_E2")],
  gene_map,
  by = "gene",
  all.x = TRUE
)

rank_symbol <- collapse_rank(
  rank_map,
  id_col = "SYMBOL",
  stat_col = "delta_E1_minus_E2"
)

rank_entrez <- collapse_rank(
  rank_map,
  id_col = "ENTREZID",
  stat_col = "delta_E1_minus_E2"
)

# ------------------------------------------------------------------------------
# 11. GSEA — HALLMARK
# ------------------------------------------------------------------------------

msg("Running Hallmark GSEA for E1 vs E2...")

hallmark_gsea <- suppressMessages(
  clusterProfiler::GSEA(
    geneList = rank_symbol,
    TERM2GENE = hallmark_term2gene,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    verbose = FALSE,
    seed = TRUE
  )
)

gsea_hallmark_df <- extract_gsea(
  hallmark_gsea,
  "Hallmark"
)

wtsv(
  gsea_hallmark_df,
  "06_E1_vs_E2_GSEA_Hallmark.tsv"
)

# ------------------------------------------------------------------------------
# 12. GSEA — GO BIOLOGICAL PROCESS
# ------------------------------------------------------------------------------

msg("Running GO-BP GSEA for E1 vs E2...")

go_gsea <- suppressMessages(
  clusterProfiler::gseGO(
    geneList = rank_entrez,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    verbose = FALSE
  )
)

gsea_go_df <- extract_gsea(
  go_gsea,
  "GO_BP"
)

wtsv(
  gsea_go_df,
  "06_E1_vs_E2_GSEA_GO_BP.tsv"
)

# ------------------------------------------------------------------------------
# 13. GSEA — REACTOME
# ------------------------------------------------------------------------------

msg("Running Reactome GSEA for E1 vs E2...")

reactome_gsea <- suppressMessages(
  ReactomePA::gsePathway(
    geneList = rank_entrez,
    organism = "human",
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    verbose = FALSE
  )
)

gsea_reactome_df <- extract_gsea(
  reactome_gsea,
  "Reactome"
)

wtsv(
  gsea_reactome_df,
  "06_E1_vs_E2_GSEA_Reactome.tsv"
)

all_gsea <- rbind(
  gsea_hallmark_df,
  gsea_go_df,
  gsea_reactome_df
)

if (nrow(all_gsea) > 0) {
  all_gsea <- all_gsea[
    order(
      all_gsea$p.adjust,
      -abs(all_gsea$NES)
    ),
    ,
    drop = FALSE
  ]
}

wtsv(
  all_gsea,
  "06_E1_vs_E2_GSEA_all_databases.tsv"
)

# ------------------------------------------------------------------------------
# 14. BUILD MASTER MODULE SUMMARY
# ------------------------------------------------------------------------------

msg("\nBuilding master module summary...")

best_pairwise <- do.call(
  rbind,
  lapply(
    split(
      traits_pairwise,
      traits_pairwise$module
    ),
    function(d) {

      d <- d[
        order(
          d$fdr,
          -abs(d$effect_SD)
        ),
        ,
        drop = FALSE
      ]

      head(d, 1)
    }
  )
)

best_pairwise <- best_pairwise[
  ,
  c(
    "module",
    "contrast",
    "effect_SD",
    "fdr"
  ),
  drop = FALSE
]

names(best_pairwise) <- c(
  "module",
  "best_disease_contrast",
  "best_disease_effect_SD",
  "best_disease_contrast_fdr"
)

master <- data.frame(
  module = all_modules,
  stringsAsFactors = FALSE
)

master <- merge(
  master,
  traits_omnibus[
    ,
    c(
      "module",
      "F",
      "p",
      "fdr",
      "eta2"
    ),
    drop = FALSE
  ],
  by = "module",
  all.x = TRUE
)

names(master)[
  names(master) == "fdr"
] <- "disease_omnibus_fdr"

names(master)[
  names(master) == "p"
] <- "disease_omnibus_p"

master <- merge(
  master,
  best_pairwise,
  by = "module",
  all.x = TRUE
)

master <- merge(
  master,
  projection[
    ,
    intersect(
      c(
        "module",
        "rho_all",
        "fdr_all",
        "rho_SjD",
        "fdr_SjD",
        "rho_residualized_group",
        "fdr_residualized_group",
        "genes_in_MSG_module",
        "genes_used_in_Blood",
        "share_used"
      ),
      names(projection)
    ),
    drop = FALSE
  ],
  by = "module",
  all.x = TRUE
)

master <- merge(
  master,
  driver_table[
    ,
    c(
      "module",
      "mean_E1",
      "mean_E2",
      "delta_E1_minus_E2",
      "hedges_g_descriptive"
    ),
    drop = FALSE
  ],
  by = "module",
  all.x = TRUE
)

master$disease_associated_FDR05 <- master$module %in%
  disease_modules

master$projectable_MSG_to_Blood <- master$module %in%
  projection_modules

master$top3_disease_module <- master$module %in%
  top3_names

master$top5_E1E2_driver <- master$module %in%
  top_driver_modules

# add top enriched term per database
for (db in c("Hallmark", "GO_BP", "Reactome")) {

  dd <- all_ora[
    all_ora$database == db,
    ,
    drop = FALSE
  ]

  if (nrow(dd) == 0) next

  best <- do.call(
    rbind,
    lapply(
      split(dd, dd$module),
      function(x) {
        x <- x[
          order(
            x$p.adjust,
            x$pvalue
          ),
          ,
          drop = FALSE
        ]
        head(x, 1)
      }
    )
  )

  add <- best[
    ,
    c(
      "module",
      "Description",
      "p.adjust"
    ),
    drop = FALSE
  ]

  names(add) <- c(
    "module",
    paste0("top_", db, "_term"),
    paste0("top_", db, "_fdr")
  )

  master <- merge(
    master,
    add,
    by = "module",
    all.x = TRUE
  )
}

master <- master[
  order(
    master$disease_omnibus_fdr
  ),
  ,
  drop = FALSE
]

wtsv(
  master,
  "06_MASTER_module_summary.tsv"
)

# ------------------------------------------------------------------------------
# 15. SIMPLE FIGURES
# ------------------------------------------------------------------------------

# A) Disease-associated modules: Hallmark top terms.
plot_disease <- disease_ora[
  disease_ora$database == "Hallmark" &
    is.finite(disease_ora$p.adjust),
  ,
  drop = FALSE
]

if (nrow(plot_disease) > 0) {

  plot_disease <- do.call(
    rbind,
    lapply(
      split(plot_disease, plot_disease$module),
      function(d) {
        d <- d[
          order(
            d$p.adjust
          ),
          ,
          drop = FALSE
        ]
        head(d, 5)
      }
    )
  )

  plot_disease$minuslog10FDR <- -log10(
    pmax(
      plot_disease$p.adjust,
      1e-300
    )
  )

  p <- ggplot(
    plot_disease,
    aes(
      x = module,
      y = reorder(
        Description,
        minuslog10FDR
      ),
      size = minuslog10FDR
    )
  ) +
    geom_point() +
    labs(
      title = "Disease-associated MSG modules: Hallmark enrichment",
      x = "WGCNA module",
      y = "Hallmark pathway",
      size = "-log10(FDR)"
    ) +
    theme_bw(base_size = 10)

  ggsave(
    file.path(
      OUTDIR,
      "06_disease_modules_Hallmark_dotplot.pdf"
    ),
    p,
    width = 10,
    height = 8
  )
}

# B) E1/E2 Hallmark GSEA.
if (nrow(gsea_hallmark_df) > 0) {

  gd <- gsea_hallmark_df[
    order(
      gsea_hallmark_df$p.adjust,
      -abs(gsea_hallmark_df$NES)
    ),
    ,
    drop = FALSE
  ]

  gd <- head(
    gd,
    20
  )

  p2 <- ggplot(
    gd,
    aes(
      x = NES,
      y = reorder(
        Description,
        NES
      ),
      size = -log10(
        pmax(
          p.adjust,
          1e-300
        )
      )
    )
  ) +
    geom_point() +
    geom_vline(
      xintercept = 0,
      linetype = 2
    ) +
    labs(
      title = "E1 vs E2 descriptive Hallmark GSEA",
      subtitle = "NES > 0 = E1-high; NES < 0 = E2-high",
      x = "Normalized Enrichment Score",
      y = "Hallmark pathway",
      size = "-log10(FDR)"
    ) +
    theme_bw(base_size = 10)

  ggsave(
    file.path(
      OUTDIR,
      "06_E1_vs_E2_Hallmark_GSEA.pdf"
    ),
    p2,
    width = 10,
    height = 8
  )
}

# ------------------------------------------------------------------------------
# 16. KEY FINDINGS TEXT
# ------------------------------------------------------------------------------

sig_all <- all_ora[
  is.finite(all_ora$p.adjust) &
    all_ora$p.adjust < FDR_CUTOFF,
  ,
  drop = FALSE
]

sig_gsea <- all_gsea[
  is.finite(all_gsea$p.adjust) &
    all_gsea$p.adjust < FDR_CUTOFF,
  ,
  drop = FALSE
]

lines <- c(
  "STAGE 06 — MODULE AND ENDOTYPE ENRICHMENT",
  "",
  sprintf(
    "Biological MSG modules evaluated: %d",
    length(all_modules)
  ),
  sprintf(
    "Disease-associated modules at FDR<0.05: %d",
    length(disease_modules)
  ),
  paste0(
    "Disease-associated modules: ",
    paste(disease_modules, collapse = ", ")
  ),
  "",
  paste0(
    "Top 3 disease modules by omnibus FDR: ",
    paste(top3_names, collapse = ", ")
  ),
  "",
  paste0(
    "MSG->Blood projectable/concordant module(s): ",
    paste(projection_modules, collapse = ", ")
  ),
  "",
  paste0(
    "Top E1/E2 module drivers by absolute mean difference: ",
    paste(top_driver_modules, collapse = ", ")
  ),
  "",
  sprintf(
    "Significant module-pathway ORA results (within-module FDR<0.05): %d",
    nrow(sig_all)
  ),
  sprintf(
    "Significant E1/E2 GSEA pathways (FDR<0.05): %d",
    nrow(sig_gsea)
  ),
  "",
  "Interpretation caveat:",
  "E1/E2 were learned from these molecular data. Their module differences and GSEA",
  "are descriptive characterization of the discovered clusters, not independent",
  "validation of endotypes.",
  "",
  "Most useful integrated output:",
  "06_MASTER_module_summary.tsv"
)

writeLines(
  lines,
  file.path(
    OUTDIR,
    "KEY_FINDINGS.txt"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    OUTDIR,
    "sessionInfo.txt"
  )
)

msg("\n============================================================")
msg("STAGE 06 COMPLETE")
msg("Output: %s", normalizePath(OUTDIR))
msg("Main summary: 06_MASTER_module_summary.tsv")
msg("============================================================")
