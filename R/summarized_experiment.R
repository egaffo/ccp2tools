#' Convert the result of [combine_ccp2_runs] into a SummarizedExperiment object
#' @param combined_prjs
#' A named list, as returned by the [combine_ccp2_runs] function.
#' Must include the following entries:
#' \itemize{
#'   \item \code{circ_read_count_mt}
#'   \item \code{lin_read_count_mt}
#'   \item \code{circ_gene_anno}
#'   \item \code{lin_xpr}
#'   \item \code{read_stats}
#'   \item \code{gene_anno}
#' }
#' @param assay select which RNA assay to extract among "circrna", "gene", or
#'  "trx" for circrna-, gene-, or transcript-level expression analysis.
#' @param coldata a data.frame with additional sample data to include into the
#'  \link{SummarizedExperiment}
#'
#' @returns a \link{RangedSummarizedExperiment} object
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
           assay = c("circrna", "gene", "trx"),
           coldata = NULL) {
    assay <- match.arg(assay)

    sample_metadata <- data.frame(combined_prjs$read_stats,
                                  row.names = "Sample",
                                  check.names = F)

    # ---- circRNA-level ----
    if (assay == "circrna") {
      ## gene annotation for the circRNAs
      ## get one gene annotation line per circRNA id
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

      se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = count_data, linCounts = lin_count_data),
        colData = sample_metadata,
        # rowData = gene_metadata
        metadata = list(circrna_host_gene = circrna_hosts,
                        circ_gene_anno = combined_prjs$circ_gene_anno),
        rowRanges = GRanges(rownames(count_data))
      )
    }

    # ---- gene-level ----
    if (assay == "gene") {
      gnx <- tximport::summarizeToGene(combined_prjs$lin_xpr$txi,
                                       combined_prjs$lin_xpr$tx2gene)

      rowRanges <- combined_prjs$gene_anno[combined_prjs$gene_anno$type == "gene", ]
      rowRanges <- rowRanges[rowRanges$gene_id %in% rownames(gnx$counts), ]
      rowRanges_order <- match(rownames(gnx$counts), GenomicRanges::mcols(rowRanges)$gene_id)

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
        rowRanges = rowRanges[rowRanges_order, ],
        metadata = list(countsFromAbundance = gnx$countsFromAbundance)
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
        rowRanges = rowRanges[rowRanges_order, ],
        metadata = list(countsFromAbundance = combined_prjs$lin_xpr$txi$countsFromAbundance)
      )

    }

    if (!is.null(coldata)) {
      colData(se) <- attach_coldata(colData(se), coldata)
    }

    return(se)
  }


#' Attach columns to existing data.frame and checks whether row and column names
#' are consistent.
#'
#' @param se_coldata the original data.frame to extend
#' @param extra_coldata a data.frame with the additional columns to add. Row
#' and column names must be consistent with the original data.frame, i.e., the
#' extra_coldata must include all row names and none colnames of the original
#' data.frame
#'
#' @returns a data.frame with the extra columns attached
#' @export
#'
#' @examples \dontrun{
#' extra_cols <-
#'   data.frame(condition = rep_len("A", rownames(colData(se))),
#'              row.names = rownames(colData(se)))
#' colData(se) <- attach_coldata(colData(se), extra_cols)
#' }
attach_coldata <- function(se_coldata, extra_coldata) {
  ## check all samples are included into the extra coldata
  if (!all(rownames(se_coldata) %in% rownames(extra_coldata))) {
    stop(
      "The extra_coldata does not include all row names of the original data.frame. ",
      paste("Missing", rownames(se_coldata)[!rownames(se_coldata) %in% rownames(extra_coldata)])
    )
  }

  ## check colnames already present in the SummarizedExperiment object
  duplicated_cols <- which(colnames(se_coldata) %in% colnames(extra_coldata))
  if (length(duplicated_cols > 0)) {
    stop("coldata already has the following columns: ",
         colnames(se_coldata)[duplicated_cols])
  }

  ## attach the extra columns
  cbind(se_coldata, extra_coldata[rownames(se_coldata), ])

}
