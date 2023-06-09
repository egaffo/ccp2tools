% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccp2tools.R
\name{collect_read_stats}
\alias{collect_read_stats}
\title{Title}
\usage{
collect_read_stats(
  read_stats_collect.file,
  circrna.reads.stats.file = NA,
  results.dir = "./",
  linear.mapper = "hisat2"
)
}
\arguments{
\item{read_stats_collect.file}{Read statistic files generated by parsing
FASTQC and (optionally) the linear aligner log files}

\item{results.dir}{The output directory path}

\item{linear.mapper}{The linear read aligner applied}
}
\description{
Title
}
