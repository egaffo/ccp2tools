% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccp2tools.R
\name{get_circompara_counts}
\alias{get_circompara_counts}
\title{Compute CirComPara2 circRNA expression}
\usage{
get_circompara_counts(
  input = "circular.reads.bed.gz.txt",
  min_methods = 2,
  output_prefix = "bks_counts.",
  stranded = FALSE,
  circrnas_gtf = NA,
  xprtypes = "UN"
)
}
\arguments{
\item{input}{File with the list of circular.reads.bed.gz files to merge,
each sample and method}

\item{min_methods}{Minimum number of circRNA detection methods required to
keep a circRNA}

\item{output_prefix}{A prefix for output file names.
Default outputs: bks.counts.intersect.csv, bks.counts.union.csv,
bks.counts.union.intersected.csv}

\item{stranded}{Set TRUE if circRNA strand has to be considered}

\item{circrnas_gtf}{A circrnas.gtf file or a text file listing circrnas.gtf
file paths to merge}

\item{xprtypes}{An underscore separated list of the strategy(ies) that will
be used to combine and report the expression estimates. The options available
are:
UN = combine and count all the unique backsplice junction read fragments
(BJRs) from any circRNA detection method;
IN = count only the BJRs commonly identified by at least 'min_methods' methods
IU = at least 1 BJRs is required to be commonly identified by at least
'min_methods' methods, then count all other BJRs
MD = compute the median of the read count reported by each method that
detected the circRNA
ME = (not yet implemented) compute the mean of the read count reported by
each method that detected the circRNA}
}
\description{
Compute CirComPara2 circRNA expression
}
