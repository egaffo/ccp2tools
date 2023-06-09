% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccp2tools.R
\name{get_bks_reads_ciri}
\alias{get_bks_reads_ciri}
\title{Obtain the backspliced reads for the circRNAs predicted by CIRI2}
\usage{
get_bks_reads_ciri(
  input = "ciri.out",
  reads_bed = "circular.reads.bed.gz",
  read_list = "bks.reads"
)
}
\arguments{
\item{input}{CIRI2 output file}

\item{reads_bed}{Output filename of BED file with read IDs for each
backsplice}

\item{read_list}{Output filename of all read ID list}
}
\description{
Obtain the backspliced reads for the circRNAs predicted by CIRI2
}
