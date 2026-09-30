# Paired Blood–MSG transcriptomics pipeline

R pipeline for analyzing paired peripheral blood (**Blood**) and minor salivary gland (**MSG**) samples from patients with Sjögren's disease (**SjD**), non-Sjögren controls (**nonSjD**), and healthy volunteers (**HV**).

## Analysis goals

The pipeline addresses four related questions:

1. **Are Blood and MSG samples paired correctly?** It validates count matrices, metadata, clinical groups, and patient-level matching.
2. **Which genes change by disease group and tissue?** It fits a paired `limma-voom` model and estimates SjD–HV, SjD–nonSjD, and tissue-by-disease interaction contrasts.
3. **Which transcriptomic programs define disease and its endotypes?** It discovers WGCNA modules in MSG, characterizes SjD endotypes, and projects the modules into Blood without retraining them.
4. **Which biological functions do those modules represent?** It performs functional enrichment and GSEA for disease-associated modules, projectable modules, and E1–E2 differences.

> E1/E2 endotypes are discovered and characterized using the same molecular data. Their enrichment results are descriptive and require independent validation.

## Repository structure

```text
.
├── README.md
├── config_paths.R                          # single path configuration file
└── src/
    ├── 00_validate_and_pair_samples.R      # validation, pairing, historical replication
    ├── 01_paired_differential_expression.R # paired differential expression
    ├── 02_wgcna_endotypes_and_blood_projection.R
    │                                       # WGCNA, endotypes, and projection
    └── 03_functional_enrichment.R          # enrichment and GSEA
```

The numeric prefixes represent the actual execution order. Do not run scripts 01–03 before the preceding stage has completed successfully.

## Data and output organization

**Inputs and outputs are deliberately separated.** The recommended layout is:

```text
/mnt/file-systems/
├── inputs/                                  # original data; never written to
│   ├── RawCountFile_filtered_blood.txt
│   ├── RawCountFile_filtered_msg.txt
│   ├── sampletable_blood.txt
│   ├── sampletable_msg.txt
│   ├── sampletable_SjD_grouping_blood.txt
│   ├── sampletable_SjD_grouping_msg.txt
│   ├── D_DE_neg-pos_blood.txt               # optional
│   └── D_DE_neg-pos_msg.txt                 # optional
└── output/                                   # created automatically
    ├── 00_pairing/
    ├── 01_differential_expression/
    ├── 02_wgcna_endotypes_projection/
    └── 03_functional_enrichment/
```

Count files must be tab-delimited matrices with genes in the first column and samples in the remaining columns. Every `sampletable_*.txt` file must contain `sampleName` and `condition` columns:

- `sampletable_blood.txt` and `sampletable_msg.txt`: historical `neg/pos` condition, used for the DESeq2 replication.
- `sampletable_SjD_grouping_*.txt`: biological `SjD/nonSjD/HV` group, used for the main analysis.

Sample names are normalized, and patient identifiers are derived by removing the tissue prefix. Stage 00 stops on duplicates, invalid counts, missing samples, or disease-group disagreement between tissues.

## Path configuration

Review `config_paths.R` before running the pipeline. Its defaults are:

```r
PROJECT_ROOT <- "/mnt/file-systems"
INPUT_DIR <- file.path(PROJECT_ROOT, "inputs")
OUTPUT_DIR <- file.path(PROJECT_ROOT, "output")
```

`INPUT_DIR` and `OUTPUT_DIR` are independent. Changing the input location does **not** place results inside that directory. Paths may also be configured without editing the file:

```bash
export PROJECT_ROOT=/path/to/project
export INPUT_DIR=/path/to/original_data
export OUTPUT_DIR=/path/to/results
```

Environment variables take precedence over defaults. All four scripts load the same `config_paths.R`; stage-specific path arguments are not required.

The configuration fails immediately if `OUTPUT_DIR` is equal to or nested below `INPUT_DIR`. This protects original data from generated pipeline files.

## Dependencies

- R 4.2 or later is recommended.
- Bioconductor: `DESeq2`, `edgeR`, `limma`, `clusterProfiler`, `org.Hs.eg.db`, `ReactomePA`, and `AnnotationDbi`.
- CRAN: `ggplot2`, `WGCNA`, `cluster`, and `msigdbr`.

Suggested installation from R:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "DESeq2", "edgeR", "limma", "clusterProfiler",
  "org.Hs.eg.db", "ReactomePA", "AnnotationDbi"
))
install.packages(c("ggplot2", "WGCNA", "cluster", "msigdbr"))
```

WGCNA can use multiple CPU cores. Set the thread count with:

```bash
export WGCNA_THREADS=8
```

## Sequential execution

Run the pipeline from the repository root. `set -euo pipefail` makes the command block stop as soon as a stage fails:

```bash
cd /workspace/hackathon_nlm_trnascriptomics
set -euo pipefail
Rscript src/00_validate_and_pair_samples.R
Rscript src/01_paired_differential_expression.R
Rscript src/02_wgcna_endotypes_and_blood_projection.R
Rscript src/03_functional_enrichment.R
```

### Stage 00 — validate and pair samples

```bash
Rscript src/00_validate_and_pair_samples.R
cat /mnt/file-systems/output/00_pairing/pairing_summary.txt
```

This stage validates the six required inputs, creates `paired_manifest.tsv` and `pairing_bundle.rds`, and attempts to reproduce the historical Blood and MSG analyses separately. The `D_DE_neg-pos_*.txt` files are optional: when present, they are compared with the replicated results, but they are not required for pairing.

**Checkpoint before continuing:** review the number of matched patients, unmatched samples, SjD/nonSjD/HV distribution, and condition concordance.

### Stage 01 — paired differential expression

```bash
Rscript src/01_paired_differential_expression.R
cat /mnt/file-systems/output/01_differential_expression/hackathon_summary.txt
```

This stage consumes `pairing_bundle.rds`, applies `filterByExpr`, TMM normalization, and `voom`, and models each patient's two samples as repeated measurements with `duplicateCorrelation(patient_id)`. It produces contrast tables, QC figures, and `voom_expression.rds`.

### Stage 02 — WGCNA, endotypes, and Blood projection

```bash
WGCNA_THREADS=8 Rscript src/02_wgcna_endotypes_and_blood_projection.R
cat /mnt/file-systems/output/02_wgcna_endotypes_projection/DAY2_SUMMARY.txt
```

This stage uses `paired_manifest.tsv` and the original count matrices. The network and modules are learned **only in MSG**; endotypes are clustered from module eigengenes, and the fixed MSG loadings are projected into the paired Blood samples. The main hand-off file for the next stage is `day2_02_05_bundle.rds`.

### Stage 03 — functional enrichment

```bash
Rscript src/03_functional_enrichment.R
```

This stage consumes Stage 02 results and generates GO, Reactome, and Hallmark enrichment tables; hub genes; disease-associated modules; E1–E2 drivers; GSEA results; and PDF figures. The integrated summary is:

```text
/mnt/file-systems/output/03_functional_enrichment/06_MASTER_module_summary.tsv
```

Internal output names beginning with `02`–`06` are retained for traceability to the original analytical steps.

## Validation and troubleshooting

Before starting a long run, verify the configured paths and required files:

```bash
Rscript -e 'source("config_paths.R"); cat("INPUT:", INPUT_DIR, "\nOUTPUT:", OUTPUT_DIR, "\n"); stopifnot(dir.exists(INPUT_DIR)); print(file.exists(c(BLOOD_COUNTS, MSG_COUNTS, BLOOD_META, MSG_META, BLOOD_GROUPING, MSG_GROUPING)))'
```

Common issues:

- **`Missing required input file(s)`**: verify `INPUT_DIR`, exact filenames, and capitalization.
- **`pairing_bundle.rds` is absent**: Stage 00 did not finish; inspect the console error and `pairing_summary.txt` when available.
- **Disease groups disagree between tissues**: correct the tissue-specific grouping files. The pipeline does not select one group automatically.
- **Packages are missing**: install them into the library used by the same R version that runs `Rscript`.
- **WGCNA runs out of memory**: use a node with more RAM. Reducing threads does not necessarily reduce the size of the main in-memory objects.
- **Restarting after failure**: every stage writes to its own directory. You may remove only the failed stage directory and rerun that stage, but any input change requires regenerating all stages from Stage 00 onward.

## Reproducibility

The WGCNA and enrichment scripts set a random seed, and Stage 02 records session information in its output. Preserve `config_paths.R`, the repository commit, the original inputs, and the complete output directory together to document a run.
