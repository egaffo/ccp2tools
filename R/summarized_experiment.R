## ---- helpers ---------------------------------------------------------------

#' Turn an identifier column into row names
#'
#' Converts `x` to a `data.frame`, uses `id_col` as row names and drops it.
#' Avoids `data.frame(x, row.names = id_col)`, which takes `id_col` as a
#' literal row name when `x` has a single row.
#'
#' @param x A `data.frame`-like object (e.g. a `data.table`).
#' @param id_col Name of the column holding unique row identifiers.
#'
#' @returns A `data.frame` without `id_col`, with its values as row names.
#' @noRd
.df_with_rownames <- function(x, id_col) {
  x <- as.data.frame(x, check.names = FALSE)
  if (!id_col %in% names(x))
    stop("column '", id_col, "' not found", call. = FALSE)
  if (anyDuplicated(x[[id_col]]))
    stop("duplicated values in '", id_col, "'", call. = FALSE)
  rownames(x) <- as.character(x[[id_col]])
  x[[id_col]] <- NULL
  x
}

#' Subset and reorder sample metadata to match count-matrix columns
#'
#' @param sample_metadata A `data.frame` with sample IDs as row names.
#' @param samples Character vector of sample IDs, in count-matrix column order.
#' @param cols Optional character vector of columns to keep; columns not in
#'   `sample_metadata` are silently skipped. `NULL` keeps all columns.
#'
#' @returns `sample_metadata` with rows in the order of `samples`.
#' @noRd
.align_coldata <- function(sample_metadata, samples, cols = NULL) {
  miss <- setdiff(samples, rownames(sample_metadata))
  if (length(miss))
    stop("samples without read_stats: ",
         paste(miss, collapse = ", "),
         call. = FALSE)
  cols <- if (is.null(cols))
    names(sample_metadata)
  else
    intersect(cols, names(sample_metadata))
  sample_metadata[samples, cols, drop = FALSE]
}

#' Collapse unique non-NA values into one string
#'
#' @param x An atomic vector.
#'
#' @returns A length-1 character: unique non-`NA` values of `x` joined by
#'   `"|"`, or `NA_character_` if all values are `NA`.
#' @noRd
.collapse_na <- function(x) {
  x <- unique(x[!is.na(x)])
  if (length(x)) paste(x, collapse = "|") else NA_character_
}

#' Build the circRNA-level SummarizedExperiment
#'
#' Missing linear counts in the source table are read as zero. circRNAs or
#' samples absent from the linear table get `NA` in `linCounts`, with a
#' warning.
#'
#' @param combined_prjs Output of [combine_ccp2_runs()].
#' @param sample_metadata Sample `data.frame` from `.df_with_rownames()`.
#'
#' @returns A `RangedSummarizedExperiment` with assays `counts` (backspliced
#'   reads) and `linCounts` (linearly spliced reads on the backsplice ends),
#'   and host genes in `rowData()` as `host_gene_id` and `host_gene_name`
#'   (multiple hosts joined by `"|"`; `NA` if intergenic or unannotated).
#' @importFrom data.table .SD
#' @noRd
.se_circrna <- function(combined_prjs, sample_metadata) {
  bks <- as.matrix(.df_with_rownames(combined_prjs$circ_read_count_mt, "circ_id"))
  circ_ids <- rownames(bks)
  samples  <- colnames(bks)
  lin <- as.matrix(.df_with_rownames(combined_prjs$lin_read_count_mt, "circ_id"))
  lin[is.na(lin)] <- 0  # assumes quantifier omits zero linear counts
  
  ri <- match(circ_ids, rownames(lin))
  ci <- match(samples, colnames(lin))
  if (anyNA(ri))
    warning(sum(is.na(ri)), " circRNAs absent from linear counts (set to NA)",
            call. = FALSE)
  if (anyNA(ci))
    warning(sum(is.na(ci)), " samples absent from linear counts (set to NA)",
            call. = FALSE)
  
  lin <- lin[ri, ci, drop = FALSE]
  dimnames(lin) <- list(circ_ids, samples)
  
  ## one host-gene line per circRNA; intergenic or unannotated -> NA
  hosts <- unique(data.table::as.data.table(combined_prjs$circ_gene_anno)[
    , c("circ_id", "gene_id", "gene_name")])
  hosts <- hosts[, lapply(.SD, .collapse_na), by = "circ_id"]
  h <- match(circ_ids, hosts$circ_id)
  
  ## TODO: check circ_id format (strand?) parses as intended
  rr <- GenomicRanges::GRanges(circ_ids)
  names(rr) <- circ_ids
  S4Vectors::mcols(rr) <- S4Vectors::DataFrame(host_gene_id   = hosts$gene_id[h],
                                               host_gene_name = hosts$gene_name[h])
  
  SummarizedExperiment::SummarizedExperiment(
    assays    = list(counts = bks, linCounts = lin),
    rowRanges = rr,
    colData   = .align_coldata(sample_metadata, samples),
    metadata  = list(
      circ_gene_anno = combined_prjs$circ_gene_anno,
      quantification = "backsplice"
    )
  )
}

#' Build a gene- or transcript-level SummarizedExperiment from a tximport list
#'
#' @param txi A tximport-style list with `counts`, `abundance`, `length` and
#'   `countsFromAbundance`.
#' @param anno A `GRanges` annotation with a `type` metadata column.
#' @param feature Annotation `type` to keep (`"gene"` or `"transcript"`).
#' @param id_col Metadata column of `anno` matching `rownames(txi$counts)`.
#' @param sample_metadata Sample `data.frame` from `.df_with_rownames()`.
#'
#' @returns A `RangedSummarizedExperiment` with assays `counts`, `abundance`
#'   and `length`, and `countsFromAbundance` stored in `metadata()`.
#' @noRd
.se_txi <- function(txi,
                    anno,
                    feature,
                    id_col,
                    sample_metadata) {
  ids     <- rownames(txi$counts)
  samples <- colnames(txi$counts)
  
  rr  <- anno[anno$type == feature]
  idx <- match(ids, S4Vectors::mcols(rr)[[id_col]])
  if (anyNA(idx))
    stop(
      sum(is.na(idx)),
      " ",
      feature,
      "s lack annotation, e.g. ",
      paste(utils::head(ids[is.na(idx)], 3), collapse = ", "),
      call. = FALSE
    )
  rr <- rr[idx]
  names(rr) <- ids
  
  cfa <- if (is.null(txi$countsFromAbundance))
    "no"
  else
    txi$countsFromAbundance
  if (cfa == "no" && (anyNA(txi$length) || any(txi$length <= 0)))
    warning("some lengths are NA or <= 0: DESeq2 will refuse them as avgTxLength",
            call. = FALSE)
  
  SummarizedExperiment::SummarizedExperiment(
    assays    = list(
      counts    = txi$counts,
      abundance = txi$abundance,
      length    = txi$length
    ),
    rowRanges = rr,
    colData   = .align_coldata(
      sample_metadata,
      samples,
      c(
        "Raw reads",
        "Clean reads",
        "Dropped",
        "Linearly unmapped",
        "Linearly mapped"
      )
    ),
    metadata  = list(
      countsFromAbundance = cfa,
      quantification = paste0("tximport_", feature)
    )
  )
}

#' Check that a converter's suggested package is installed
#'
#' @param pkg Package name.
#'
#' @returns `TRUE`, invisibly; stops with an install hint if `pkg` is missing.
#' @noRd
.require_suggested <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop(
      "package '",
      pkg,
      "' is required for this function; install it with ",
      "BiocManager::install(\"",
      pkg,
      "\")",
      call. = FALSE
    )
  invisible(TRUE)
}

#' Row annotation as a data.frame for `edgeR::DGEList(genes = )`
#'
#' Keeps only atomic `rowData()` columns, since list columns (e.g.
#' `CharacterList`) do not convert to a `data.frame`.
#'
#' @param se A `SummarizedExperiment`.
#'
#' @returns A `data.frame` with the original column names, or `NULL` if no
#'   atomic column is left.
#' @noRd
.rowdata_df <- function(se) {
  rd <- SummarizedExperiment::rowData(se)
  keep <- vapply(rd, is.atomic, logical(1))
  if (any(keep))
    as.data.frame(rd[, keep, drop = FALSE], optional = TRUE)
  else
    NULL
}

#' Add sample annotation columns to an existing column-data table
#'
#' Appends the columns of `extra_coldata` to `se_coldata`, matching rows by
#' row name (sample ID). Fails if any sample of `se_coldata` is missing from
#' `extra_coldata`, or if a column would be duplicated. Extra rows in
#' `extra_coldata` are ignored.
#'
#' @param se_coldata A `data.frame` or [S4Vectors::DataFrame] with sample IDs
#'   as row names, typically `SummarizedExperiment::colData(se)`.
#' @param extra_coldata A `data.frame` with sample IDs as row names. It must
#'   contain every row name of `se_coldata` and none of its column names.
#'
#' @returns An object of the same class as `se_coldata`, with the columns of
#'   `extra_coldata` appended and rows in the original order.
#' @export
#'
#' @examples
#' cd <- S4Vectors::DataFrame(`Raw reads` = c(1e6, 2e6), check.names = FALSE,
#'                            row.names = c("s1", "s2"))
#' extra <- data.frame(condition = c("B", "A"), row.names = c("s2", "s1"))
#' attach_coldata(cd, extra)
attach_coldata <- function(se_coldata, extra_coldata) {
  missing_rows <- setdiff(rownames(se_coldata), rownames(extra_coldata))
  if (length(missing_rows) > 0)
    stop("extra_coldata lacks these samples: ",
         paste(missing_rows, collapse = ", "),
         call. = FALSE)
  
  dup_cols <- intersect(colnames(se_coldata), colnames(extra_coldata))
  if (length(dup_cols) > 0)
    stop("coldata already has these columns: ",
         paste(dup_cols, collapse = ", "),
         call. = FALSE)
  
  cbind(se_coldata, extra_coldata[rownames(se_coldata), , drop = FALSE])
}

## ---- main ------------------------------------------------------------------

#' Convert combined CirComPara2 runs into a SummarizedExperiment
#'
#' Builds a tool-neutral \link[SummarizedExperiment:RangedSummarizedExperiment-class]{RangedSummarizedExperiment}
#' at circRNA, gene or transcript level from the output of
#' [combine_ccp2_runs()]. Counts are stored as quantified (not rounded); use
#' [se_to_dds()], [se_to_dgelist()] or [se_to_voom_input()] to prepare the
#' object for a specific differential-expression framework.
#'
#' @param combined_prjs A named list, as returned by [combine_ccp2_runs()].
#'   Must include `read_stats` and, depending on `assay`:
#'   * `"circrna"`: `circ_read_count_mt`, `lin_read_count_mt`, `circ_gene_anno`;
#'   * `"gene"`, `"trx"`: `lin_xpr` (with `txi` and `tx2gene`) and `gene_anno`.
#' @param assay Expression level to extract: `"circrna"`, `"gene"` or `"trx"`
#'   (transcript).
#' @param coldata Optional `data.frame` of extra sample annotation, with sample
#'   IDs as row names; added with [attach_coldata()].
#'
#' @returns A \link[SummarizedExperiment:RangedSummarizedExperiment-class]{RangedSummarizedExperiment}:
#'   * `"circrna"`: assays `counts` (backspliced reads) and `linCounts`
#'     (linearly spliced reads on the backsplice ends; missing source values
#'     read as zero); host genes in `rowData()` as `host_gene_id` and
#'     `host_gene_name` (multiple hosts joined by `"|"`; `NA` if intergenic or
#'     unannotated); the full circRNA-gene annotation in
#'     `metadata()$circ_gene_anno`.
#'   * `"gene"`, `"trx"`: assays `counts`, `abundance` (TPM) and `length`
#'     (average transcript length); `metadata()$countsFromAbundance` records
#'     how counts were derived.
#'
#'   Sample read statistics are in `colData()`, in count-matrix column order.
#'
#' @seealso [se_to_dds()], [se_to_dgelist()], [se_to_voom_input()]
#' @export
#'
#' @examples
#' \dontrun{
#' prjs <- c("/home/user/ccp2_nf/sample1",
#'           "/home/user/ccp2_nf/sample2")
#' combined_prjs <- combine_ccp2_runs(prjs)
#'
#' circrna_se <- to_summarized_experiment(combined_prjs)
#' gene_se    <- to_summarized_experiment(combined_prjs, assay = "gene")
#' trx_se     <- to_summarized_experiment(combined_prjs, assay = "trx")
#'
#' dds <- se_to_dds(gene_se, design = ~ 1)
#' }
to_summarized_experiment <- function(combined_prjs,
                                     assay = c("circrna", "gene", "trx"),
                                     coldata = NULL) {
  assay <- match.arg(assay)
  
  required <- c("read_stats", switch(
    assay,
    circrna = c(
      "circ_read_count_mt",
      "lin_read_count_mt",
      "circ_gene_anno"
    ),
    c("lin_xpr", "gene_anno")
  ))
  missing_entries <- setdiff(required, names(combined_prjs))
  if (length(missing_entries))
    stop("combined_prjs lacks: ",
         paste(missing_entries, collapse = ", "),
         call. = FALSE)
  
  sample_metadata <- .df_with_rownames(combined_prjs$read_stats, "Sample")
  lx <- combined_prjs$lin_xpr
  
  se <- switch(
    assay,
    circrna = .se_circrna(combined_prjs, sample_metadata),
    gene    = .se_txi(
      tximport::summarizeToGene(lx$txi, lx$tx2gene),
      combined_prjs$gene_anno,
      "gene",
      "gene_id",
      sample_metadata
    ),
    trx     = .se_txi(
      lx$txi,
      combined_prjs$gene_anno,
      "transcript",
      "transcript_id",
      sample_metadata
    )
  )
  
  if (!is.null(coldata))
    SummarizedExperiment::colData(se) <-
    attach_coldata(SummarizedExperiment::colData(se), coldata)
  se
}

## ---- converters --------------------------------------------------------------

#' Prepare a SummarizedExperiment for DESeq2
#'
#' Creates a `DESeqDataSet` following the same logic as
#' `DESeq2::DESeqDataSetFromTximport()`. For tximport-derived objects (those
#' with `metadata(se)$countsFromAbundance`), counts are rounded and, when
#' `countsFromAbundance == "no"`, the `length` assay is renamed `avgTxLength`
#' so that DESeq2 uses it as a gene-specific normalization offset. Other
#' objects (e.g. circRNA-level) are passed to `DESeq2::DESeqDataSet()`
#' unchanged.
#'
#' @param se A `SummarizedExperiment` from [to_summarized_experiment()], with
#'   `counts` as first assay.
#' @param design A formula or model matrix, passed to `DESeq2::DESeqDataSet()`.
#' @param ... Further arguments passed to `DESeq2::DESeqDataSet()`.
#'
#' @returns A `DESeqDataSet`; `rowData(se)` (e.g. host genes for circRNAs) is
#'   kept and returned by `DESeq2::results()` via `mcols()`.
#' @seealso [se_to_dgelist()], [se_to_voom_input()]
#' @export
#'
#' @examples
#' if (requireNamespace("DESeq2", quietly = TRUE)) {
#'   set.seed(1)
#'   cts <- matrix(stats::rnbinom(60, mu = 100, size = 5) + 0.4, nrow = 10,
#'                 dimnames = list(paste0("g", 1:10), paste0("s", 1:6)))
#'   len <- matrix(stats::runif(60, 500, 2000), nrow = 10,
#'                 dimnames = dimnames(cts))
#'   se <- SummarizedExperiment::SummarizedExperiment(
#'     assays = list(counts = cts, length = len),
#'     colData = S4Vectors::DataFrame(condition = gl(2, 3, labels = c("A", "B")),
#'                                    row.names = colnames(cts)),
#'     metadata = list(countsFromAbundance = "no"))
#'
#'   dds <- se_to_dds(se, design = ~ condition)
#'   SummarizedExperiment::assayNames(dds)
#' }
se_to_dds <- function(se, design, ...) {
  .require_suggested("DESeq2")
  cfa <- S4Vectors::metadata(se)$countsFromAbundance
  if (!is.null(cfa)) {
    cfa <- match.arg(cfa,
                     c("no", "scaledTPM", "lengthScaledTPM", "dtuScaledTPM"))
    SummarizedExperiment::assay(se, "counts") <-
      round(SummarizedExperiment::assay(se, "counts"))
    an <- SummarizedExperiment::assayNames(se)
    if (cfa == "no" && "length" %in% an) {
      len <- SummarizedExperiment::assay(se, "length")
      if (anyNA(len) || any(len <= 0))
        stop("NA or non-positive lengths in 'length' assay", call. = FALSE)
      SummarizedExperiment::assayNames(se)[an == "length"] <- "avgTxLength"
    } else {
      message("countsFromAbundance = '",
              cfa,
              "': using counts only, no length offset")
    }
  }
  DESeq2::DESeqDataSet(se, design = design, ...)
}

#' Prepare a SummarizedExperiment for edgeR
#'
#' Creates a `DGEList`. When `metadata(se)$countsFromAbundance == "no"` and a
#' `length` assay is present, average transcript lengths are turned into an
#' offset matrix following the tximport vignette recipe
#' (`edgeR::scaleOffset()`); otherwise TMM normalization factors are computed
#' with `edgeR::calcNormFactors()`.
#'
#' Filter (e.g. with `edgeR::filterByExpr()`) after calling this function:
#' subsetting a `DGEList` also subsets its offset matrix.
#'
#' @param se A `SummarizedExperiment` from [to_summarized_experiment()].
#' @param group Optional vector or factor of sample groups, passed to
#'   `edgeR::DGEList()`.
#'
#' @returns A `DGEList`, with `colData(se)` in `$samples` and the atomic
#'   `rowData(se)` columns in `$genes` (e.g. `host_gene_id` and
#'   `host_gene_name` for circRNAs), carried through to `edgeR::topTags()`.
#' @seealso [se_to_dds()], [se_to_voom_input()]
#' @export
#'
#' @examples
#' if (requireNamespace("edgeR", quietly = TRUE)) {
#'   set.seed(1)
#'   cts <- matrix(stats::rnbinom(60, mu = 100, size = 5) + 0.4, nrow = 10,
#'                 dimnames = list(paste0("g", 1:10), paste0("s", 1:6)))
#'   len <- matrix(stats::runif(60, 500, 2000), nrow = 10,
#'                 dimnames = dimnames(cts))
#'   se <- SummarizedExperiment::SummarizedExperiment(
#'     assays = list(counts = cts, length = len),
#'     colData = S4Vectors::DataFrame(condition = gl(2, 3, labels = c("A", "B")),
#'                                    row.names = colnames(cts)),
#'     metadata = list(countsFromAbundance = "no"))
#'
#'   y <- se_to_dgelist(se, group = se$condition)
#'   dim(y$offset)
#' }
se_to_dgelist <- function(se, group = NULL) {
  .require_suggested("edgeR")
  cts <- SummarizedExperiment::assay(se, "counts")
  cfa <- S4Vectors::metadata(se)$countsFromAbundance
  
  ## row annotation -> y$genes, carried through to topTags()
  y <- edgeR::DGEList(
    cts,
    group   = group,
    genes   = .rowdata_df(se),
    samples = as.data.frame(SummarizedExperiment::colData(se), optional = TRUE)
  )
  if (identical(cfa, "no") &&
      "length" %in% SummarizedExperiment::assayNames(se)) {
    len <- SummarizedExperiment::assay(se, "length")
    if (anyNA(len) || any(len <= 0))
      stop("NA or non-positive lengths in 'length' assay", call. = FALSE)
    len <- len / exp(rowMeans(log(len)))
    norm_cts <- cts / len
    eff_lib  <- edgeR::calcNormFactors(norm_cts) * colSums(norm_cts)
    y <- edgeR::scaleOffset(y, log(sweep(len, 2, eff_lib, "*")))
  } else {
    y <- edgeR::calcNormFactors(y)
  }
  y
}

#' Prepare a SummarizedExperiment for limma-voom
#'
#' Creates a TMM-normalized `DGEList` ready for `limma::voom()`. voom cannot
#' use a length-offset matrix, so tximport-derived gene- or transcript-level
#' counts should already carry the length correction: import them with
#' `countsFromAbundance = "lengthScaledTPM"` (or `"scaledTPM"`,
#' `"dtuScaledTPM"`). A warning is issued otherwise.
#'
#' @inheritParams se_to_dgelist
#'
#' @returns A `DGEList` with TMM normalization factors, `colData(se)` in
#'   `$samples` and the atomic `rowData(se)` columns in `$genes` (e.g.
#'   `host_gene_id` and `host_gene_name` for circRNAs), carried through
#'   `limma::voom()` and `limma::lmFit()` to `limma::topTable()`.
#' @seealso [se_to_dds()], [se_to_dgelist()]
#' @export
#'
#' @examples
#' if (requireNamespace("edgeR", quietly = TRUE)) {
#'   set.seed(1)
#'   cts <- matrix(stats::rnbinom(60, mu = 100, size = 5), nrow = 10,
#'                 dimnames = list(paste0("g", 1:10), paste0("s", 1:6)))
#'   se <- SummarizedExperiment::SummarizedExperiment(
#'     assays = list(counts = cts),
#'     colData = S4Vectors::DataFrame(condition = gl(2, 3, labels = c("A", "B")),
#'                                    row.names = colnames(cts)),
#'     metadata = list(countsFromAbundance = "lengthScaledTPM"))
#'
#'   y <- se_to_voom_input(se, group = se$condition)
#'   y$samples
#' }
se_to_voom_input <- function(se, group = NULL) {
  .require_suggested("edgeR")
  cfa <- S4Vectors::metadata(se)$countsFromAbundance
  if (!is.null(cfa) &&
      !cfa %in% c("lengthScaledTPM", "scaledTPM", "dtuScaledTPM"))
    warning(
      "countsFromAbundance = '",
      cfa,
      "': voom ignores transcript-length bias; ",
      "re-import with countsFromAbundance = 'lengthScaledTPM'",
      call. = FALSE
    )
  
  ## row annotation -> y$genes, carried through voom/lmFit to topTable()
  y <- edgeR::DGEList(
    SummarizedExperiment::assay(se, "counts"),
    group   = group,
    genes   = .rowdata_df(se),
    samples = as.data.frame(SummarizedExperiment::colData(se), optional = TRUE)
  )
  edgeR::calcNormFactors(y)
}