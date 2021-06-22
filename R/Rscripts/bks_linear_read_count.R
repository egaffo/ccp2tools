#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Rsubread))

option_list <- list(
  make_option(c("-c", "--circrna_ids_file"),
              action = "store",
              type = "character",
              help = paste0("Either the reliable_circRNAs.csv file or a ",
                            "file with a list of circRNA identifiers ",
                            "formatted as chr:start-end. ",
                            "N.B.: circRNA coordinates must be BED format like")),
  make_option(c("-b", "--bam_files_list_file"),
              action = "store",
              type = "character",
              default = "",
              help = paste0("The list of the BAM files to process. ",
                            "A text file with one line per BAM file")),
  make_option(c("-t", "--cpus"),
              action = "store",
              type = "integer",
              default = 1,
              help = "The number of threads for parallel computing"),
  make_option(c("-o", "--out_dir"),
              action = "store",
              type = "character",
              default = "./",
              help = paste0("Output directory. The bks_lin_spld_read_count.csv ",
                            "and sn_circ_id_split.gtf output file will be ",
                            "saved in this directory. ",
                            "The bks_lin_spld_read_count.csv file will be a ",
                            "table with circ_ids as rows and samples as columns",
                            "with cells reporting the count of the linearly ",
                            "spliced reads at the backsplicing ends")),
  make_option(c("-p", "--isPE"),
              action = "store_true",
              type = "logical",
              default = TRUE,
              help = paste0("Treat reads as paired end: count read fragments ",
                            "instead of single reads. The default is to count ",
                            "the fragments as to comply with the CirComPara2 ",
                            "default counting approach")),
  make_option(c("-n", "--tidy_sample_names"),
              action = "store_true",
              type = "logical",
              default = TRUE,
              help = paste0("Whether to process the BAM file name nad path to ",
                            "obtain simple sample identifiers. BAM file names ",
                            "should be in the path/to/SAMPLENAME_histat2.bam ",
                            "format; then, sample name will be the file path ",
                            "basename amended from the '_hisat2.bam' suffix")
  )
)

parser <- OptionParser(usage = paste0("%prog -c reliable_circRNAs.csv ",
                                      "-b bam_files.txt -t 12 ",
                                      "-o lin_bks_counts"),
                       option_list = option_list)
arguments <- parse_args(parser, positional_arguments = F)

circrna_ids_file <- arguments$circrna_ids_file
bam_files_list_file <- arguments$bam_files_list_file
cpus <- arguments$cpus
out_dir <- arguments$out_dir
tidy_sample_names <- arguments$tidy_sample_names

dir.create(out_dir, recursive = T, showWarnings = F)

bam_files <- readLines(bam_files_list_file)

in_col_num <- dim(fread(circrna_ids_file, nrows = 1, header = F))[2]

if (in_col_num == 1) {
  ## input is a list
  circrna_ids <- fread(circrna_ids_file, header = F, col.names = "circ_id")
}
if (in_col_num > 1) {
  ## input is a table (the reliable circRNA expression matrix)
  circrna_ids <- fread(circrna_ids_file, select = "circ_id")
}

circrna_ids[, c("chr", "start", "stop") := tstrsplit(circ_id, ":|-",
                                                     type.convert = T)]

## compose the single nucleotide backsplice representation
sn_unique_circ <-
  melt(unique(circrna_ids[, .(circ_id,
                              chr,
                              start = as.integer(start + 1),
                              stop = as.integer(stop))]),
       id.vars = c("chr", "circ_id"))

sn_unique_circ_gtf <-
  sn_unique_circ[, .(chr, ".", variable, start = value, stop = value,
                     ".", "+", ".",
                     V9 = paste0("gene_id \"", circ_id, "\";"))]

sn_unique_circ_gtf_file <-
  file.path(out_dir, "sn_circ_id_split.gtf")

fwrite(x = sn_unique_circ_gtf,
       file = sn_unique_circ_gtf_file,
       col.names = F, row.names = F, sep = "\t", quote = F)

## spliced reads on circRNA junctions gene expression
fc_linear_backsplices <-
  featureCounts(files = bam_files,
                annot.ext = sn_unique_circ_gtf_file,
                isGTFAnnotationFile = T,
                GTF.featureType = c("start", "stop"),
                splitOnly = T,
                nonSplitOnly = F,
                GTF.attrType = "gene_id",
                useMetaFeatures = T,
                allowMultiOverlap = T,
                countMultiMappingReads = T,
                largestOverlap = T,
                primaryOnly = T,
                strandSpecific = 0,
                isPairedEnd = isPE,
                read2pos = NULL,
                countChimericFragments = F,
                requireBothEndsMapped = ifelse(isPE, T, F),
                nthreads = cpus)

lin_bks_counts <-
  data.table(fc_linear_backsplices$counts,
             keep.rownames = "circ_id")

if (tidy_sample_names) {
  colnames(lin_bks_counts) <-
    sub("\\.hisat2\\.bam", "", basename(colnames(lin_bks_counts)))
}

lin_bks_counts_file <- file.path(out_dir, "bks_lin_spld_read_count.csv")

fwrite(x = lin_bks_counts,
       file = lin_bks_counts_file,
       col.names = F,
       row.names = F,
       sep = "\t",
       quote = F)
