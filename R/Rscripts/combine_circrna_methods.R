#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ccp2tools))

option_list <- list(
  make_option(c("-i", "--file_names"),
              action = "store",
              type = "character",
              help = paste0("",
                            "")),
  make_option(c("-f", "--filter"),
              action = "store_true",
              type = "logical",
              default = TRUE,
              help = paste0("Filter circRNAs according to the number of reads ",
                            "detection methods")),
  make_option(c("-r", "--min_reads"),
              action = "store",
              type = "integer",
              default = 2,
              help = paste0("The minimum read count a circRNA must have to ",
                            "be included in the results when the _filter_  ",
                            "option is enabled")),
  make_option(c("-m", "--min_methods"),
              action = "store",
              type = "integer",
              default = 2,
              help = paste0("The minimum number of methods a circRNA must be ",
                            "detected by to be included in the results when ",
                            "the _filter_ option is enabled")),
  make_option(c("-t", "--tidy_sample_names"),
              action = "store_true",
              type = "logical",
              default = T,
              help = paste0("Simplify the sample names. If FALSE then sample ",
                            "names will be the full file path.")),
  make_option(c("-n", "--sample_names"),
              action = "store",
              type = "character",
              default = NA,
              help = paste0("A text file with samples names for the input ",
                            "files. Order must be the same as the input files"))
)

parser <- OptionParser(usage = paste0("%prog -i ",
                                      "-o "),
                       option_list = option_list)
arguments <- parse_args(parser, positional_arguments = F)

file_names <- arguments$file_names
out_dir <- arguments$out_dir
filter <- arguments$filter
min_reads <- arguments$min_reads
min_methods <- arguments$min_methods

dir.create(out_dir, recursive = T, showWarnings = F)

if (tidy_sample_names) {
  ## write and apply a function to infer the sample name from the file path
  sample_names <- sub(".*/samples/([^/]+/.*)", "\\1", file_names)
} else {
  sample_names <- file_names
}

if (!is.na(arguments$sample_names)) {
  sample_names <- readLines(arguments$sample_names)
}

names(file_names) <- sample_names

combined_method_circrnas <-
  combine_circrna_methods(file_names, filter = filter,
                          min_reads = min_reads, min_methods = min_methods)

circrna_expression_matrix_file <-
  file.path(out_dir, "circrna_expression_matrix.csv")

fwrite(x = combined_method_circrnas$wide_tab,
       file = circrna_expression_matrix_file,
       quote = F, sep = "\t", row.names = F)

circrna_expression_table_file <-
  file.path(out_dir, "circrna_expression_table.csv")

fwrite(x = combined_method_circrnas$long_tab,
       file = circrna_expression_table_file,
       quote = F, sep = "\t", row.names = F)
