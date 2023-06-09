% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccp2tools.R
\name{get_bks_reads_ce2_star}
\alias{get_bks_reads_ce2_star}
\title{Obtain the backspliced reads for the circRNAs predicted by CIRCexplorer2 on
STAR mappings}
\usage{
get_bks_reads_ce2_star(
  chimreads,
  circrnas,
  output = "circular.reads.bed.gz",
  range = 10
)
}
\arguments{
\item{chimreads}{Chimeric.out.junction BED transformed by
chimout_junc_to_bed.py}

\item{circrnas}{circularRNA_known.txt as output by CIRCexplorer2 annotate,
or back_spliced_junction.bed as output by CIRCexplorer2 parse}

\item{output}{The circRNA read IDs for each circRNA in compressed BED}

\item{range}{Number of basepairs tolerated in realigning circRNAs from
CIRCexplorer2 annotate}
}
\description{
Obtain the backspliced reads for the circRNAs predicted by CIRCexplorer2 on
STAR mappings
}
