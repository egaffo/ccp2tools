#' Combine results of different CirComPara2 runs
#'
#' @param files a character indicating the path of the parent directory of the
#' CirComPara2 runs to be combined.
#' _files_ child directories will be scan to search for the CirComPara2 results
#' to merge (i.e. the circRNA backsplice read counts and linear read counts).
#' _files_ might also be a file listing either (1) directories,
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
#' @param isPairedEnd logical indicating if the linear alignments are paired-end.
#' TRUE by default.
#' @param cpus integer indicating the number of parallel threads to use. 1 by
#' default.
#'
#' @return a list of two elements: (1) the matrix of the merged samples'
#' backspliced read counts ( _ccp\_counts\_dt_ ), and (2) the matrix of the
#' merged samples' backsplice linear read counts ( _lin\_bks\_counts_ ). The
#' latter is _NA_ when _merge\_lin = FALSE_.
#' @import data.table
#' @import Rsubread
#' @export
#'
#' @examples
combine_ccp2_runs <-
  function(files,
           min_methods = 2,
           min_reads = 2,
           merge_lin = F,
           recycle_existing_lincount = F,
           isPairedEnd = T,
           cpus = 1) {

    if(dir.exists(files)) {
      ## if files is a directory, search all subdirectories for the CCP2 results
      ## to merge

      ccp_count_files <-
        dir(path = files,
            pattern = "bks.counts.union.csv",
            full.names = T,
            recursive = T)
    }else{
      ## files might be a file listing either
      ## 1) directories to merge,
      ## 2) bks.counts.union.csv files to merge, or
      ## 3) a mix of directories and files

      files <- readLines(files)

      ccp_count_files <-
        sapply(files, function(f) {
          ifelse(file.exists(f), ## if not a file, assume it is a directory
                 f,
                 dir(path = files_dir,
                     pattern = "bks.counts.union.csv",
                     full.names = T,
                     recursive = T))
        })
    }

    ccp_counts <- rbindlist(sapply(ccp_count_files,
                                   fread,
                                   simplify = F),
                            use.names = T)

    ccp_counts[, circ_id := paste0(chr, ":", start, "-", end)]

    reliable_circ_ids <-
      ccp_counts[n_methods >= min_methods &
                   read.count >= min_reads,
                 unique(circ_id)]

    ccp_counts_dt <-
      dcast(data = ccp_counts[circ_id %in% reliable_circ_ids],
            formula = circ_id ~ sample_id,
            value.var = "read.count",
            fill = 0)

    lin_bks_counts <- NA

    if(merge_lin) {
      ## merge also the linear counts
      ## N.B. this might require to compute the linear counts for circRNAs not
      ## detected in the sample (but detected in other samples)

      if(recycle_existing_lincount) {
        ## recycle the existing lincount files and compute the lincounts
        ## for the circRNAs expressed only in other samples
        ## TODO
      }else{
        ## count the linearly spliced reads on the backsplice ends

        ## retrieve HISAT2 BAM files from the CCP2 project directories
        ccp_count_file_postfix <- file.path("circular_expression",
                                            "circrna_analyze",
                                            "counts",
                                            "bks.counts.union.csv")

        ccp2_runs_basedirs <- unique(sub(ccp_count_file_postfix,
                                         "",
                                         ccp_count_files))

        bam_files <-
          unlist(sapply(ccp2_runs_basedirs,
                        FUN = function(f){dir(path = f,
                                              pattern = "_hisat2.bam$",
                                              full.names = T,
                                              recursive = T)},
                        simplify = T,
                        USE.NAMES = F))

        names(bam_files) <-
          sapply(bam_files,
                 function(f)sub("_hisat2.bam", "", basename(f)),
                 USE.NAMES = F)

        ## make SAF data.frame of start/end of circRNAs
        sn_unique_circ <-
          melt(ccp_counts_dt[, tstrsplit(circ_id, ":|-", type.convert = F),
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
        ## usign the featureCounts function of the Rsubread package
        fc_linear_backsplices <-
          featureCounts(files = bam_files,
                        annot.ext = circ_ids_saf,
                        splitOnly = T,
                        nonSplitOnly = F,
                        useMetaFeatures = T,
                        allowMultiOverlap = T,
                        countMultiMappingReads = T,
                        countReadPairs = T,
                        countChimericFragments = F,
                        largestOverlap = T,
                        primaryOnly = T,
                        strandSpecific = 0,
                        isPairedEnd = isPairedEnd,
                        nthreads = cpus)

        lin_bks_counts <-
          data.table(fc_linear_backsplices$counts,
                     keep.rownames = "circ_id")
        colnames(lin_bks_counts) <- sub("_hisat2.bam", "", colnames(lin_bks_counts))
      }
    }

    list(circ_read_count_mt = ccp_counts_dt,
         lin_read_count_mt = lin_bks_counts)

  }

