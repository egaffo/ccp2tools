#' Convert the result of combine_ccp2_runs() into a SummarizedExperiment object
#'
#' @param combined_prjs a named list, as resulting from the
#'  \link[combine_ccp2_runs]{ccp2tools::combine_ccp2_runs()} function. Must include the following entries:
#'  circ_read_count_mt, lin_read_count_mt, circ_gene_anno, lin_xpr, read_stats,
#'  gene_anno.
#' @param assays select which RNA assay to extract among "circrna", "gene", or
#'  "trx" for circrna-, gene-, or transcript-level expression analysis.
#'
#' @returns a \link{SummarizedExperiment} object for circrna assay; a
#'  \link{RangedSummarizedExperiment} if gene or transcript assay
#'
#' @import data.table tximport GenomicRanges SummarizedExperiment
#'
#' @export
#'
#' @examples \dontrun{
#'
#' library(ccp2tools)
#'
#' prjs <- c("/home/user/ccp2_nf/sample1",
#'           "/home/user/ccp2_nf/sample2")
#'
#' combined_prjs <- combine_ccp2_runs(prjs)
#'
#' circrna_se <- to_summarized_experiment(combined_prjs = combined_prjs)
#'
#' gene_se <- to_summarized_experiment(combined_prjs = combined_prjs,
#'                                     assay = "gene")
#'
#' trx_se <- to_summarized_experiment(combined_prjs = combined_prjs,
#'                                     assay = "trx")
#' }
#'
to_summarized_experiment <-
  function(combined_prjs,
           assay = c("circrna", "gene", "trx")) {

  assay <- match.arg(assay)

  sample_metadata <- data.frame(combined_prjs$read_stats,
                                row.names = "Sample",
                                check.names = F)

  # ---- circRNA-level ----
  if (assay == "circrna") {

    ## get one annotation line per circRNA id
    circrna_hosts <-
      unique(combined_prjs$circ_gene_anno[, .(circ_id, gene_id, gene_name)])[, lapply(.SD, function(x) {
        paste0(x, collapse = "|")
      }), by = circ_id]

    intergenic_circs <-
      combined_prjs$circ_read_count_mt$circ_id[!combined_prjs$circ_read_count_mt$circ_id %in% circrna_hosts$circ_id]

    circrna_hosts <- merge(
      circrna_hosts,
      data.table::data.table(circ_id = intergenic_circs),
      by = "circ_id",
      all = T
    )

    ## the backspliced read counts
    count_data <- data.frame(
      combined_prjs$circ_read_count_mt,
      row.names = "circ_id",
      check.names = F
    )

    ## linearly-spliced read counts on the backsplices
    lin_count_data <- data.frame(
      combined_prjs$lin_read_count_mt,
      row.names = "circ_id",
      check.names = F
    )

    ## gene annotation for the circRNAs
    gene_metadata <- data.frame(circrna_hosts,
                                row.names = "circ_id",
                                check.names = F)

    se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = count_data, linCounts = lin_count_data),
      colData = sample_metadata,
      rowData = gene_metadata
    )
  }

  # ---- gene-level ----
  if (assay == "gene") {

    gnx <- tximport::summarizeToGene(combined_prjs$lin_xpr$txi,
                                     combined_prjs$lin_xpr$tx2gene)

    rowRanges <- combined_prjs$gene_anno[combined_prjs$gene_anno$type == "gene", ]
    rowRanges <- rowRanges[rowRanges$gene_id %in% rownames(gnx$counts), ]
    rowRanges_order <- match(rownames(gnx$counts), mcols(rowRanges)$gene_id)

    se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(
        counts = gnx$counts,
        abundance = gnx$abundance,
        length = gnx$length
      ),
      colData = sample_metadata[, c("Raw reads",
                                    "Clean reads",
                                    "Dropped",
                                    "Linearly unmapped",
                                    "Linearly mapped")],
      rowRanges = rowRanges[rowRanges_order, ]
    )
  }

  # ---- transcript-level ----
  if (assay == "trx") {
    rowRanges <- combined_prjs$gene_anno[combined_prjs$gene_anno$type == "transcript", ]
    rowRanges <- rowRanges[rowRanges$transcript_id %in% rownames(combined_prjs$lin_xpr$txi$counts), ]
    rowRanges_order <- match(rownames(combined_prjs$lin_xpr$txi$counts),
                             mcols(rowRanges)$transcript_id)

    se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(
        counts = combined_prjs$lin_xpr$txi$counts,
        abundance = combined_prjs$lin_xpr$txi$abundance,
        length = combined_prjs$lin_xpr$txi$length
      ),
      colData = sample_metadata[, c("Raw reads",
                                    "Clean reads",
                                    "Dropped",
                                    "Linearly unmapped",
                                    "Linearly mapped")],
      # rowData = data.frame(combined_prjs$lin_xpr$tx2gene, row.names = "t_name"),
      rowRanges = rowRanges[rowRanges_order, ]
    )

  }

  se
}
