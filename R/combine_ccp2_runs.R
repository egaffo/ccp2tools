#' Combine the results of multiple CirComPara2 analyses
#'
#' This function will combine multiple CirComPara2 analyses into single output
#' result files. Specifically, it combines the circRNA expression of all 
#' samples from the different projects into one expression matrix reporting the 
#' backsplice junction read counts. Moreover, it also combines the read counts 
#' of the reads linearly spliced on the backsplice junctions. For the circRNAs 
#' not detected in all the projects it will calculate the linearly spliced read 
#' counts by parsing the linear alignment files.
#'
#' @param files a character indicating the path of the parent directory of the
#' CirComPara2 runs to be combined.
#' \code{files} child directories will be scan to search for the CirComPara2 
#' results to merge (i.e. the circRNA backsplice read counts and linear read 
#' counts).
#' \code{files} might also be a file listing either (1) directories,
#' (2) bks.counts.union.csv files, or (3) a mix of directories and files.
#' @param min_methods the minimum number of circRNA-detection methods a circRNA
#' must be identified by. 2 by default.
#' @param min_reads the minimum number of backspliced reads a circRNA must have
#' to be kept. 2 by default.
#' @param merge_lin logical indicating whether to merge also the files
#' reporting the number of reads linearly mapped into the backsplices.
#' FALSE by default.
#' @param recycle_existing_lincount a logical indicating whether to use
#' previously computed backsplice linear alignment read count files and only
#' compute the possible missed backsplices that have been detected only in other
#' samples. FALSE by default. (this option is not currently implemented yet).
#' @param is_paired_end logical indicating if the linear alignments are paired-end.
#' TRUE by default.
#' @param cpus integer indicating the number of parallel threads to use. 1 by
#' default.
#'
#' @return a list of two elements: (1) the matrix of the merged samples'
#' backspliced read counts (\code{ccp_counts_dt}), and (2) the matrix of the
#' merged samples' backsplice linear read counts (\code{lin_bks_counts}). The
#' latter is \code{NA} when \code{merge_lin = FALSE}.
#' @import data.table
#' @import Rsubread
#' @export
#'
#' @examples
combine_ccp2_runs <-
  function(files,
           min_methods = 2,
           min_reads = 2,
           merge_lin = FALSE,
           recycle_existing_lincount = FALSE,
           is_paired_end = TRUE,
           cpus = 1) {

    # files <- "/sharedfs01/circrna/zicirc/analysis/GSE159225/batch"
    # require(data.table)
    # require(Rsubread)

    if (dir.exists(files)) {
      #' if files is a directory, search all sub directories for the CCP2
      #' results to merge

      ccp_count_files <-
        dir(path = files,
            pattern = "bks.counts.union.csv",
            full.names = TRUE,
            recursive = TRUE)
    }else {
      ## 'files' might be one file listing either
      ## 1) the CCP2 project directories to merge,
      ## 2) the 'bks.counts.union.csv' files to merge, or
      ## 3) a mix of directories and files

      files <- readLines(files)

      ccp_count_files <-
        sapply(files, function(f) {
          ifelse(file.exists(f), ## if not a file, assume it is a directory
                 f,
                 dir(path = files_dir,
                     pattern = "bks.counts.union.csv",
                     full.names = TRUE,
                     recursive = TRUE))
        })
    }

    message("Combining BJR counts from ", length(ccp_count_files),
            "projects...")
    ccp_counts <- data.table::rbindlist(sapply(ccp_count_files,
                                               data.table::fread,
                                               simplify = FALSE),
                                        use.names = TRUE)

    ccp_counts[, circ_id := paste0(chr, ":", start, "-", end)]
    message("Detected ", nrow(ccp_counts), " circRNAs")

    message("Removing circRNAs detected with <= ",
            min_reads,
            " BJRs by <= ",
            min_methods,
            " circRNA detection methods, in all samples...")
    reliable_circ_ids <-
      ccp_counts[n_methods >= min_methods &
                   read.count >= min_reads,
                 unique(circ_id)]
    message(length(reliable_circ_ids), " reliable circRNAs were kept.")

    ccp_counts_dt <-
      dcast(data = ccp_counts[circ_id %in% reliable_circ_ids],
            formula = circ_id ~ sample_id,
            value.var = "read.count",
            fill = 0)

    lin_bks_counts <- NA

    if (merge_lin) {
      ## merge also the linear counts
      ## N.B. this might require to compute the linear counts for circRNAs not
      ## detected in the sample (but detected in other samples)
      message("Merging the linearly spliced read counts...")

      if (recycle_existing_lincount) {
        ## recycle the existing lincount files and compute the lincounts
        ## for the circRNAs expressed only in other samples
        ## TODO
        message("Use of pre-computed linearly spliced read counts ",
                "not yet implemented. Skipping.")
      }else {
        ## count the linearly spliced reads on the backsplice ends

        ## retrieve HISAT2 BAM files from the CCP2 project directories
        ccp_count_file_postfix <- file.path("circular_expression",
                                            "circrna_analyze",
                                            "counts",
                                            "bks.counts.union.csv")

        ccp2_runs_basedirs <- unique(sub(ccp_count_file_postfix,
                                         "",
                                         ccp_count_files))
        message("Counting the reads linearly spliced on the backsplice ",
                "junctions from ",
                length(ccp2_runs_basedirs),
                " projects...")

        bam_files <-
          unlist(sapply(ccp2_runs_basedirs,
                        FUN = function(f) {
                          dir(path = f,
                              pattern = "_hisat2.bam$",
                              full.names = TRUE,
                              recursive = TRUE) },
                        simplify = TRUE,
                        USE.NAMES = FALSE))

        names(bam_files) <-
          sapply(bam_files,
                 function(f)sub( "_hisat2.bam", "", basename(f) ),
                 USE.NAMES = FALSE)
        message("Parsing ", length(bam_files), " linear alignment files...")

        ## make SAF data.frame of start/end of circRNAs
        ## TODO: strandedness
        sn_unique_circ <-
          melt(ccp_counts_dt[, tstrsplit(circ_id, ":|-", type.convert = FALSE),
                             by = circ_id],
               id.vars = c("circ_id", "V1"),
               variable.name = "featureType",
               value.name = "Start")

        sn_unique_circ[featureType == "V2", featureType := "Start"]
        sn_unique_circ[featureType == "V3", featureType := "Stop"]
        sn_unique_circ[, Start := as.integer(Start) + 1][, `:=`(End = Start,
                                                                Strand = "+")]

        circ_ids_saf <- as.data.frame(sn_unique_circ[, .(GeneID = circ_id,
                                                         Chr = V1,
                                                         Start,
                                                         End,
                                                         Strand,
                                                         featureType)])

        ## count spliced reads on circRNA junctions gene expression
        ## using the featureCounts function of the Rsubread package
        fc_linear_backsplices <-
          featureCounts(files = bam_files,
                        annot.ext = circ_ids_saf,
                        splitOnly = TRUE,
                        nonSplitOnly = FALSE,
                        useMetaFeatures = TRUE,
                        allowMultiOverlap = TRUE,
                        countMultiMappingReads = TRUE,
                        countReadPairs = TRUE,
                        countChimericFragments = FALSE,
                        largestOverlap = TRUE,
                        primaryOnly = TRUE,
                        strandSpecific = 0,
                        isPairedEnd = is_paired_end,
                        nthreads = cpus)

        lin_bks_counts <-
          data.table(fc_linear_backsplices$counts,
                     keep.rownames = "circ_id")
        colnames(lin_bks_counts) <- sub("_hisat2.bam", "",
                                        colnames(lin_bks_counts))
      }
    }

    list(circ_read_count_mt = ccp_counts_dt,
         lin_read_count_mt = lin_bks_counts)
  }
