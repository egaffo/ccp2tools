# ccp2tools (development version)

# ccp2tools 0.2.0

## New features

* `to_summarized_experiment()` converts `combine_ccp2_runs()` results into a
  `RangedSummarizedExperiment` at circRNA, gene (`assay = "gene"`) or
  transcript (`assay = "trx"`) level.
  - circRNA level: assays `counts` (backspliced reads) and `linCounts`
    (linearly spliced reads on the backsplice ends); host genes in
    `rowData()` as `host_gene_id` and `host_gene_name`.
  - gene/transcript level: assays `counts`, `abundance` and `length` from
    tximport; `metadata()$countsFromAbundance` records how counts were derived.
* `attach_coldata()` adds extra sample annotation to `colData()`.
* New converters to differential-expression frameworks:
  - `se_to_dds()` for DESeq2 (same logic as `DESeqDataSetFromTximport()`);
  - `se_to_dgelist()` for edgeR (length offsets as in the tximport vignette);
  - `se_to_voom_input()` for limma-voom.
  `se_to_dgelist()` and `se_to_voom_input()` carry `rowData()` (e.g. circRNA
  host genes) into `$genes`, so it appears in `topTags()` / `topTable()`.
  DESeq2 and edgeR are suggested, not required.

## Bug fixes

* `to_summarized_experiment()` works with a single sample; circRNA and
  linear count matrices now have consistent row names.
* Missing linear counts are read as zero; circRNAs or samples absent from the
  linear-count table are reported with a warning.
* Fixed an error with unstranded data.
* Fixed an error when different runs used an uneven number of circRNA
  detection methods.