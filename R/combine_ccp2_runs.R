#' Try to resolve strandedness ambiguities of the detected circRNAs
#'
#' @param strand_pattern a character in the form of 
#' "strand_reads_nMethods|strand_reads_nMethods". E.g. "+_30_6|-_26_1"
#' @param circ_methods (optional) the methods corresponding to the strand 
#' pattern, separated by an "@" character. It helps to solve some ambiguities.
#' E.g: "CIRCexplorer2_bwa|CIRCexplorer2_star|CIRCexplorer2_tophat@dcc"
#'
#' @return "+", "-", or "." when the ambiguity cannot be resolved
#' @export
#'
#' @examples guess_strand("+_30_6|-_26_1") # "+"
#' guess_strand("+_30_6|-_36_1") # "-"
#' guess_strand("+_30_6|-_30_1") # "+"
#' guess_strand("+_30_1|-_30_1") # "."
#' guess_strand("+_30_6") # "+"
guess_strand <- function(strand_pattern, circ_methods = NULL) {
  
  pp <- strsplit(strand_pattern, "_|\\|")[[1]]
  ## if there was no ambiguous pattern, just return the strand
  if(length(pp) < 6) return(pp[1])
  ## return the strand with the largest read count
  if(as.integer(pp[2]) > as.integer(pp[5])) return(pp[1])
  if(as.integer(pp[2]) < as.integer(pp[5])) return(pp[4])
  ## if the number of reads is the same, then return the strand with the 
  ## largest consensus among methods, "." if the ambiguity was not resolved
  if(as.integer(pp[2]) == as.integer(pp[5])) {
    if(as.integer(pp[3]) > as.integer(pp[6])) return(pp[1])
    if(as.integer(pp[3]) < as.integer(pp[6])) return(pp[4])
    if(as.integer(pp[3]) == as.integer(pp[6])) {
      if(!is.null(circ_methods)) {
        ## check methods to resolve ambiguities
        mm <- strsplit(circ_methods, "@")[[1]]
        ## do not select the strand given by dcc as it showed many discordant 
        ## strand in previous tests 
        ss <- c(pp[1], pp[4])[which(!grepl("dcc", mm))][1]
        ifelse(is.na(ss),
               return("."),
               return(ss))
      }else{
        return(".")
      }
    }
  }
}

#' Get circRNA BJR count expression matrix
#' 
#' Compute the backsplice junction read (BJR) counts according to the 
#' CirComPara2's framework.  
#' This function collect the results of one or more CirComPara2 runs and merges 
#' the results.
#'
#' @param files the path of a CirComPara2 run, or a text file listing the paths 
#' of multiple CirComPara2 run to merge.
#' @param is_list_file set TRUE if "files" is a text file listing the projects 
#' to merge
#'
#' @return a data.table with the following columns: 
#' "sample_id"
#' "chr", "start", "end", "strand"
#' "read.count": the BJR count
#' "n_methods": the number of the circRNA-detection tools that commonly 
#' identified the circRNA
#' "circ_methods": the names of the circRNA-detection tools that commonly 
#' identified the circRNA in a bar-separated list alphanumerically sorted. 
#' E.g: "CIRCexplorer2_bwa|ciri_out|dcc"
#' "circ_id": the circRNA identifier composed as chr:start-end[:strand]
#' "strand_pattern": a string representing the number of reads and methods 
#' according to their strand(s). See also \link[ccp2tools]{guess_strand}.
#' @export
#'
#' @examples
merge_ccp_counts <- function(files) {
  
  ccp_count_files <-
    sapply(files, function(f) {
      ifelse(R.utils::isFile(f), ## if not a file, assume it is a directory
             f,
             dir(path = f,
                 pattern = "bks.counts.union.csv",
                 full.names = TRUE,
                 recursive = TRUE))
    })
  
  message("Combining BJR counts from ", length(ccp_count_files),
          " projects...")
  ccp_counts <- data.table::rbindlist(sapply(ccp_count_files,
                                             data.table::fread,
                                             simplify = FALSE),
                                      use.names = TRUE)
  
  ccp_counts[, circ_id := paste0(chr, ":", start, "-", end)]
  
  ## check strandedness inconsistencies and fix 
  ## check if any circ_id has two strands
  ambig_strnd_candidates <- 
    ccp_counts[, .N, 
               by = .(sample_id, circ_id, 
                      strand)][, .N, 
                               by = .(sample_id, 
                                      circ_id)][N > 1, unique(circ_id)]
  
  if (length(ambig_strnd_candidates) > 1) {
    message("Found ",
            length(ambig_strnd_candidates),
            " circRNAs with ambiguous strand assignment.")
    
    message("Trying to fix ambiguous strand circRNAs...")
    circ_reads <- 
      data.table::rbindlist(sapply(file.path(dirname(ccp_count_files),
                                             "bks.counts.collected_reads.csv"),
                                   data.table::fread,
                                   simplify = F))[, circ_id := paste0(chr, ":",
                                                                      start,
                                                                      "-", end),
                                                  by = .(chr, start,
                                                         end)][circ_id %in%
                                                                 ambig_strnd_candidates]
    
    ## count read fragments irrespective of the alignment strand
    uns_frags_count <- 
      circ_reads[, .N, by = .(sample_id, circ_id, 
                              read_id)][, .(frag.count = .N), 
                                        by = .(sample_id, circ_id)]
    
    ## compare the actual number of fragments with the sum of the number of 
    ## reads from different strands.
    fixed_strand <- 
      merge(ccp_counts[circ_id %in% ambig_strnd_candidates,
                       .(read.count = sum(read.count),
                         n_methods = length(unique(unlist(strsplit(paste0(circ_methods, 
                                                                          collapse = "|"), 
                                                                   "\\|")))),
                         circ_methods = paste0(circ_methods, collapse = "@"),
                         strand_pattern = paste0(strand, "_", read.count, "_",
                                                 n_methods, collapse = "|")),
                       by = .(sample_id, circ_id)],
            uns_frags_count, 
            by = c("sample_id", 
                   "circ_id"))[, .(strand = guess_strand(strand_pattern, 
                                                         circ_methods),
                                   strand_pattern, 
                                   circ_methods = paste0(sort(unique(unlist(strsplit(circ_methods, 
                                                                                     "\\|")))),
                                                         collapse = "|")), 
                               by = .(read.count = frag.count,
                                      sample_id, 
                                      circ_id, 
                                      n_methods)]
    
    message("Resolved ambiguous strand of ", 
            length(ambig_strnd_candidates) - 
              length(fixed_strand[strand == ".", unique(circ_id)]),
            " circRNAs.")
    
    fixed_strand[, c("chr", "start", "end") := data.table::tstrsplit(circ_id, ":|-"), 
                 by = circ_id]
    
    ccp_counts <- 
      data.table::rbindlist(list(ccp_counts[!circ_id %in% ambig_strnd_candidates],
                                 fixed_strand), 
                            use.names = T, fill = T)
  }
  message("Detected ", length(unique(ccp_counts$circ_id)), " circRNAs")
  
  ccp_counts
}

#' Count the reads linearly spliced on the backsplice junctions
#'
#' @param files the CirComPara2 project directories
#' @param circ_ids the list of backsplices in the form of chr:start-end[:strand] 
#' @param is_stranded should the strand be considered? (default: TRUE)
#' @param is_paired_end were the reads mapped as paired-end? (default: TRUE)
#' @param cpus number of CPUs for parallel execution (default: 1)
#'
#' @return a matrix with circRNAs as rows and samples as columns and counts in 
#' the cells
#' 
#' @importFrom Rsubread featureCounts
#' @import data.table
#' @export
#'
#' @examples
compute_lin_bks_counts <- function(files, circ_ids, 
                                   is_stranded = TRUE, 
                                   is_paired_end = TRUE, 
                                   cpus = 1) {
  
  circ_ids <- data.table::data.table(circ_id = circ_ids)
  
  ## retrieve HISAT2 BAM files from the CCP2 project directories
  message("Counting the reads linearly spliced on the backsplice ",
          "junctions from ",
          length(files),
          " projects...")
  
  bam_files <-
    unlist(sapply(files,
                  function(f) {
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
  id_fileds <- c("Chr", "Start", "Stop")
  if (is_stranded) id_fileds <- c(id_fileds, "Strand")
  
  circ_ids[, c(id_fileds) := data.table::tstrsplit(circ_id,
                                                   ":|-",
                                                   type.convert = FALSE),
           by = circ_id]
  if (is_stranded) circ_ids[Strand == "", Strand := "-"]
  
  sn_unique_circ <-
    data.table::melt(circ_ids,
                     id.vars = c("circ_id", "Chr", "Strand"),
                     variable.name = "featureType",
                     value.name = "Start")
  
  ## convert from 0-based into 1-based
  sn_unique_circ[, Start := as.integer(Start) + 1][, `:=`(End = Start)]
  if (!is_stranded) sn_unique_circ[, Strand := "+"]
  
  ## compose SAF input
  circ_ids_saf <- as.data.frame(sn_unique_circ[, .(GeneID = circ_id,
                                                   Chr,
                                                   Start,
                                                   End,
                                                   Strand,
                                                   featureType)])
  
  ## count spliced reads on circRNA junctions gene expression
  ## using the featureCounts function of the Rsubread package
  fc_linear_backsplices <-
    Rsubread::featureCounts(files = bam_files,
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
                            strandSpecific = ifelse(is_stranded, 1, 0),
                            isPairedEnd = is_paired_end,
                            nthreads = cpus)
  
  lin_bks_counts <-
    data.table::data.table(fc_linear_backsplices$counts,
                           keep.rownames = "circ_id")
  colnames(lin_bks_counts) <- sub("_hisat2.bam", "",
                                  colnames(lin_bks_counts))
  
  lin_bks_counts
}

#' Merge the counts of the reads linearly spliced on the backsplice junctions
#'
#' @param files the paths of CirComPara2 runs to merge.
#'
#' @return
#' @export
#'
#' @examples
merge_lin_bks_counts <- function(files) {
  
  ## find the ccp_bks_linexp.csv files
  ccp_bks_linexp_files <-
    unlist(sapply(files,
                  function(f) {
                    dir(path = f,
                        pattern = "ccp_bks_linexp.csv",
                        full.names = TRUE,
                        recursive = TRUE) },
                  simplify = TRUE,
                  USE.NAMES = FALSE))
  
  message("Combining the reads linearly spliced on the backsplice ",
          "junctions from ",
          length(files),
          " projects...")
  
  ccp_bks_linexp <- data.table::rbindlist(sapply(ccp_bks_linexp_files,
                                                 data.table::fread,
                                                 simplify = FALSE),
                                          use.names = TRUE)
  
  ## return a matrix
  data.table::dcast(data = ccp_bks_linexp, 
                    formula = circ_id ~ sample_id, 
                    value.var = "lin.reads", 
                    fill = NA)
}

#' Find circRNA host-gene(s)
#'
#' @param circ_ids a list of circRNA identifiers in the form of 
#' chr:start-end[:strand]
#' @param gtf_file an Ensemble GTF gene annotation file
#'
#' @return a data.table of the genes overlapping each circRNA, one row for each 
#' overlap.
#' 
#' @import GenomicRanges 
#' @importFrom rtracklayer import
#' @export
#'
#' @examples 
#' \dontrun{
#' get_circrna_host_genes("10:1000676-1000868:+", 
#'                        "Homo_sapiens.GRCh38.108.gtf")} 
get_circrna_host_genes <- function(circ_ids, gtf_file) {
  
  circ_dt <- 
    data.table::data.table(circ_id = circ_ids)[, data.table::tstrsplit(circ_id, 
                                                                       ":|-"), 
                                               by = circ_id][]
  cols <- c("circ_id", "chr", "start", "end")
  
  stranded <- FALSE
  ## is stranded?
  if (ncol(circ_dt) > 4) {
    ## fix "-" strand
    circ_dt[V4 == "", V4 := "-"]
    cols <- c(cols, "strand")
    stranded <- TRUE
  }
  colnames(circ_dt) <- cols
  
  circ_dt[, `:=`(start = as.integer(start),
                 end = as.integer(end))]
  
  circ_gr <- 
    GenomicRanges::makeGRangesFromDataFrame(df = circ_dt, 
                                            keep.extra.columns = TRUE,
                                            ignore.strand = !stranded, 
                                            starts.in.df.are.0based = TRUE)
  
  message("Importing gene annotations from ", gtf_file)
  gtf_gr <- rtracklayer::import(gtf_file)
  
  message("Finding host-gene for ", length(circ_ids), " circRNAs...")
  hits <- GenomicRanges::findOverlaps(query = circ_gr, 
                                      subject = gtf_gr[gtf_gr$type == "gene"],
                                      type = "any",
                                      select = "all",
                                      ignore.strand = !stranded)
  
  hits_df <- as.data.frame(hits)
  hits_df$circ_id <- circ_gr$circ_id[hits_df$queryHits]
  
  data.table::data.table(cbind(hits_df, 
                               data.frame(gtf_gr[gtf_gr$type == 
                                                   "gene"])[hits_df$subjectHits, 
                                                   ]))
}

#' Merge multiple CirComPara2 runs to get linear transcript/gene expression
#' 
#' @param prj_paths 
#' @param ... additional parameters used in tximport::tximport()
#'
#' @return a named list of two elements: 
#'  1) txi: a tximport object at transcript level
#'  2) tx2gene: the transcript/gene ids table used in the txi object
#' 
#' @importFrom tximport tximport
#' @export
#'
#' @examples \dontrun{
#' prjs <- c("/home/user/CirComPara2/batch1",
#'           "/home/user/CirComPara2/batch2")
#' 
#' ## get the linear transcripts' expression
#' trx_ls <- merge_lin_counts(prjs)
#' 
#' trx <- trx_ls$txi
#' 
#' ## get the linear gene expression#' 
#' gnx <- tximport::summarizeToGene(trx_ls$txi, trx_ls$tx2gene)
#' }
merge_lin_counts <- function(prj_paths, ...) {
  
  t_data_files <- 
    unlist(sapply(prj_paths, function(f) {
      dir(path = f,
          pattern = "t_data.ctab",
          full.names = TRUE,
          recursive = TRUE)
    }))
  names(t_data_files) <-
    gsub(file.path(".*", "samples", "([^/]+)", 
                   "processings", "stringtie", "ballgown_ctabs", 
                   "t_data.ctab"),
         "\\1",
         t_data_files)
  
  tx2gene <- 
    data.table::fread(c(t_data_files)[1], 
                      data.table = F)[, c("t_name", "gene_id", "gene_name", 
                                          "chr", "start", "end", "strand", 
                                          "num_exons", "length")]
  
  txi <- tximport::tximport(files = t_data_files,
                            importer = data.table::fread,
                            type = "stringtie",
                            tx2gene = tx2gene,
                            txIdCol = "t_name",
                            geneIdCol = "gene_id",
                            txOut = T,
                            ...)
  
  list(txi = txi, tx2gene = tx2gene)
}

#' Combine the results of multiple CirComPara2 analyses
#'
#' This function will merge the otputs of multiple CirComPara2 analyses into
#' single output result files. Specifically, it combines the circRNA expression
#' of all samples from the different projects into one expression matrix
#' reporting the backsplice junction read counts. Moreover, it also combines the
#' read counts of the reads linearly spliced on the backsplice junctions. For
#' the circRNAs not detected in all the projects it will calculate the linearly
#' spliced read counts by parsing the linear alignment files.
#'
#' @param files a character indicating the path of the parent directory of the
#'   CirComPara2 runs to be combined. \code{files} child directories will be
#'   scan to search for the CirComPara2 results to merge (i.e. the circRNA
#'   backsplice read counts and linear read counts). \code{files} might also be
#'   a file listing either (1) directories, (2) bks.counts.union.csv files, or
#'   (3) a mix of directories and files.
#' @param merge_circs \code{TRUE} to get merged backsplice read count matrix and
#'   gene annotation of circRNAs (default: \code{TRUE})
#' @param merge_lin_bks \code{TRUE} indicates to merge the files reporting the
#'   number of reads linearly mapped into the backsplices (default: \code{TRUE})
#' @param merge_lin \code{TRUE} to get linear transcript/gene expression merged
#'   from the projects (default: \code{TRUE})
#' @param min_methods the minimum number of circRNA-detection methods a circRNA
#'   must be identified by. 2 by default.
#' @param min_reads the minimum number of backspliced reads a circRNA must have
#'   to be kept. 2 by default.
#' @param recycle_existing_lincount a logical indicating whether to use
#'   previously computed backsplice linear alignment read count files and only
#'   compute the possible missed backsplices that have been detected only in
#'   other samples. FALSE by default. (this option is not currently implemented
#'   yet).
#' @param is_paired_end logical indicating if the linear alignments are
#'   paired-end (default: \code{TRUE})
#' @param cpus integer indicating the number of parallel threads to use. 1 by
#'   default.
#' @param gtf_file the path of the Ensemble GTF gene annotation file;
#'   \code{auto} will use the GTF path from the vars.py file; \code{NULL} will
#'   skip this step.
#' @param is_stranded TRUE if circRNA strandedness has to be considered
#' @param ... additional parameters that will be passed to the tximport function
#'   for linear gene expression
#' @return a list of four elements: (1) \code{ccp_counts_dt}: the matrix of the
#'   merged samples'backspliced read counts; (2) \code{lin_bks_counts}: the
#'   matrix of the merged samples' backsplice linear read counts if
#'   \code{merge_lin_bks = FALSE}, \code{NA} otherwise; (3)
#'   \code{circ_gene_anno}: the circRNA host-gene annotation; and (4)
#'   \code{lin_xpr}: the linear transcript expression as read counts, if
#'   \code{merge_lin = TRUE}, \code{NULL} otherwise.
#' @import data.table Rsubread tximport
#' @export
#' 
#' @examples \dontrun{
#' prjs <- c("/home/user/circompara2_batch1",
#'            "/home/user/circompara2_batch2")
#'            
#' ## merge the projects           
#' combined_prjs <- combine_ccp2_runs(prjs)
#' 
#' ## circRNA BJR matrix
#' circrna_raw_counts <- 
#'   data.frame(combined_prjs$circ_read_count_mt, 
#'              row.names = "circ_id")
#'              
#' ## compute the circular to linear expression proportions (a.k.a. CLPs)
#' clps <- 
#'   data.table::dcast(merge(data.table::melt(combined_prjs$circ_read_count_mt, 
#'                                            id.vars = "circ_id", 
#'                                            variable.name = "sample_id", 
#'                                            value.name = "BJR"),
#'                           data.table::melt(combined_prjs$lin_read_count_mt, 
#'                                            id.vars = "circ_id", 
#'                                            variable.name = "sample_id", 
#'                                            value.name = "lin.reads"),
#'                            by = c("sample_id", 
#'                                   "circ_id"))[, CLP := BJR/(BJR + lin.reads)],
#'                     formula = circ_id ~ sample_id, value.var = "CLP")
#' 
#' ## get circRNA host-gene(s)
#' circrna_hosts <- 
#'   unique(combined_prjs$circ_gene_anno[, .(circ_id, gene_id, 
#'                                           gene_name)])[, lapply(.SD, 
#'                                                         function(x){paste0(x, 
#'                                                            collapse = "|")}), 
#'                                                        by = circ_id]
#' ## prepare for differential gene expression analysis with DESeq2
#' gnx <- tximport::summarizeToGene(combined_prjs$lin_xpr$txi, 
#'                                  combined_prjs$lin_xpr$tx2gene)
#' sampleTable <- data.frame(sample_id = colnames(combined_prjs$lin_xpr$txi$counts))
#' dds <- DESeq2::DESeqDataSetFromTximport(gnx, sampleTable, ~1)
#' }
combine_ccp2_runs <-
  function(files,
           merge_circs = TRUE,
           merge_lin_bks = TRUE,
           merge_lin = TRUE,
           min_methods = 2,
           min_reads = 2,
           recycle_existing_lincount = TRUE,
           is_paired_end = TRUE,
           is_stranded = TRUE,
           cpus = 1,
           gtf_file = "auto",
           ...) {
    
    is_list_file <- FALSE ## TODO: determine automatically
    if (is_list_file) {
      ## 'files' might be one file listing either
      ## 1) the CCP2 project directories to merge,
      ## 2) the 'bks.counts.union.csv' files to merge, or
      ## 3) a mix of directories and files
      
      files <- readLines(files)
    }
    
    ccp_counts_dt <- NA
    circ_gene_anno <- NA
    ## -------- merge backsplice junction read counts -------- ##
    if (merge_circs) {
      
      ccp_counts <- merge_ccp_counts(files)
      
      if (is_stranded) {
        ## add strand to the circRNA identifier
        ccp_counts[, circ_id := paste0(circ_id, ":", strand), 
                   by = .(circ_id, strand)]
      }
      
      ## Select relible circRNAs
      message("Removing circRNAs detected with <= ",
              min_reads,
              " BJRs and by <= ",
              min_methods,
              " circRNA detection methods, in all samples...")
      reliable_circ_ids <-
        ccp_counts[n_methods >= min_methods &
                     read.count >= min_reads,
                   unique(circ_id)]
      message(length(reliable_circ_ids), " reliable circRNAs were kept.")
      
      ## make the reliable circRNA expression matrix
      ccp_counts_dt <-
        data.table::dcast(data = ccp_counts[circ_id %in% reliable_circ_ids],
                          formula = circ_id ~ sample_id,
                          value.var = "read.count",
                          fill = 0)
      
      ## -------- compute circRNA host-gene annotation -------- ##
      if (!is.null(gtf_file)) {
        if (gtf_file == "auto") {
          gtf_file <- 
            gsub(" |\"|\'", "", 
                 strsplit(grep("ANNOTATION", 
                               readLines(file.path(files[1], 
                                                   "vars.py")), 
                               value = TRUE), "=")[[1]][2])
          message("Using gene annotation from file ", gtf_file)
          # TODO: check file.exists(gtf_file)
        }
        
        circ_gene_anno <- get_circrna_host_genes(reliable_circ_ids, gtf_file)
      }
    }
    
    
    ## -------- merge the linearly spliced reads on the BJ -------- ##
    lin_bks_counts <- NA
    if (merge_lin_bks) {
      ## merge also the linear counts
      ## N.B. this might require to compute the linear counts for circRNAs not
      ## detected in the sample (but detected in other samples)
      message("Merging the linearly spliced read counts...")
      
      if (recycle_existing_lincount) {
        ## recycle the existing lincount files and compute the lincounts
        ## for the circRNAs expressed only in other samples
        message("Use of pre-computed linearly spliced read counts")
        
        lin_bks_counts <- merge_lin_bks_counts(files)
        
        # ## here we have two choices: 
        # ##   1) just overlook the linear bks reads of the circRNAs not detected
        # ##      in some samples (if circRNA expression is 0, then the CLR and 
        # ##      CLP will be 0)
        # ##   2) calculate the lin bks reads of the missed circRNAs in their
        # ##      samples
        # compute_missing_linbks <- FALSE
        # if (compute_missing_linbks) {
        #   
        #   ## get missed circRNAs in samples
        #   missed_bks <- 
        #     data.table::melt(data = lin_bks_counts,
        #                      na.rm = FALSE,
        #                      id.vars = "circ_id",
        #                      variable.name = "sample_id",
        #                      value.name = "lin.reads")[is.na(lin.reads)]
        #   
        #   if (nrow(missed_bks) > 0) {
        #     ## TODO
        #   }
        # }
        
        ## update circ_ids
        if (is_stranded) {
          
          if (merge_circs) {
            circid_map <- data.table::data.table(circ_id = ccp_counts_dt$circ_id)
            circid_map[, circ_id_s := circ_id][, circ_id := sub(":.$", "", 
                                                                circ_id_s)]
            
            # TODO: check duplicate IDs
            # dups <- circid_map[, .N, by = circ_id][N > 1, circ_id]
            # circid_map[, .N, by = circ_id_s][N > 1]
            
            lin_bks_counts <- 
              merge(circid_map,
                    lin_bks_counts,
                    all.x = TRUE, all.y = FALSE,
                    by = "circ_id")[, circ_id := circ_id_s][, circ_id_s := NULL][] 
          } else {
            warning("Merging the existing linearly-spliced-read-aligned-to-backsplice ", 
                    "counts accounting the strandedness is only available when",
                    " merge_circs = TRUE. Strand information will not be ", 
                    "considered.")
          }
        }
      }else {
        
        ## count the linearly spliced reads on the backsplice ends
        ## this might take some time to run...
        lin_bks_counts <- compute_lin_bks_counts(files,
                                                 ccp_counts_dt$circ_id, 
                                                 is_stranded,
                                                 is_paired_end,
                                                 cpus)
      }
    }
    
    ## -------- merge the linear transcript/gene read counts -------- ##
    lin_xpr <- NA
    if (merge_lin) {
      lin_xpr <- merge_lin_counts(files, ...)
    }
    
    ## -------- END -------- ##
    message("Done merging CirComPara2 projects!")
    list(circ_read_count_mt = ccp_counts_dt,
         lin_read_count_mt = lin_bks_counts,
         circ_gene_anno = circ_gene_anno,
         lin_xpr = lin_xpr)
    
  }
