# Stage 04: exploratory E1 vs E2 differential expression

Run after Stage 02 (and optionally after Stage 03):

```bash
Rscript src/04_endotype_de_and_plots.R
```

Uses the existing `config_paths.R`, `RawCountFile_filtered_msg.txt`, `02_wgcna_endotypes_projection/day2_02_05_bundle.rds`, and `04_preliminary_endotypes_SjD.tsv`. Outputs are written separately to `OUTPUT_DIR/04_endotype_de/`, preserving earlier stages.

**Methods:** SjD MSG only; design-aware `filterByExpr`, TMM, `limma::voom`, E1−E2 model, empirical Bayes, BH-FDR. Gene labels are mapped from Ensembl IDs (version suffix removed) via `org.Hs.eg.db` to HGNC symbols and full gene names. Tables retain original IDs for traceability. Unmapped genes remain in the full DE table but are not displayed as ENSG labels in plots. A gene is classified as E1-high or E2-high if BH FDR < 0.05 and absolute log2FC >= 1.

**Outputs:** `06_E1_vs_E2_differential_expression.tsv`, `06_E1_vs_E2_gene_ranking_annotated.tsv`, `06_E1_high_top100_genes_annotated.tsv`, `06_E2_high_top100_genes_annotated.tsv`, `06_gene_id_to_HGNC_map.tsv`, `06_E1_vs_E2_volcano.png/.pdf`, `06_E1_vs_E2_per_patient_heatmap.png/.pdf`, `06_E1_vs_E2_per_patient_heatmap_zscores.tsv`, `04_ENDOTYPE_DE_SUMMARY.txt`.

The volcano x-axis is E1−E2 log2FC and y-axis −log10(BH FDR); gene labels use HGNC symbols. The heatmap is **per patient**, with gene-wise z-scores of normalized MSG expression across the SjD E1/E2 samples and patient IDs replaced by within-endotype labels. If there are fewer than 15 significant mapped genes per direction, the plot supplements with best-ranked descriptive genes and documents this in its caption.

**Scientific limitation:** E1/E2 were discovered from the same transcriptomic data. This is exploratory, conditional molecular characterization **not independent validation**. Group-only DE does not adjust clinical or technical covariates; investigate confounding and validate findings in a held-out/external dataset.

**Dependencies:** existing Stage 03 packages plus `limma` (already required in Stages 01/02); `ggrepel` is optional for improved volcano label spacing.

**Smoke checks after running:**

```bash
test -s "$OUTPUT_DIR/04_endotype_de/06_E1_vs_E2_differential_expression.tsv"
test -s "$OUTPUT_DIR/04_endotype_de/06_E1_vs_E2_volcano.png"
test -s "$OUTPUT_DIR/04_endotype_de/06_E1_vs_E2_per_patient_heatmap.png"
cat "$OUTPUT_DIR/04_endotype_de/04_ENDOTYPE_DE_SUMMARY.txt"
```
