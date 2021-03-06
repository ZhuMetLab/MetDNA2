% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MetDNA2Plots.R
\name{plotVolcano}
\alias{plotVolcano}
\title{plotVolcano}
\usage{
plotVolcano(
  raw_data,
  p_cutoff = 0.05,
  fc_cutoff = 1.25,
  log_tran = c("log2", "log10")
)
}
\arguments{
\item{raw_data}{1st column: feature name; 2nd column: p_values; 3rd column: fold_changes}

\item{p_cutoff}{the p cutoff of significant; All p-values would be -log10 transformed in the plot. Default: 0.05}

\item{fc_cutoff}{the fold change cutoff of significant: Default: 1.25}

\item{log_tran}{log_transformation method for fold change, "log2" or "log10". Default: log2}
}
\description{
plotVolcano
}
\examples{
data('stat_result')
plotVolcano(raw_data = stat_result, p_cutoff = 0.05, fc_cutoff = 1.25, log_tran = 'log2')
}
\author{
Zhiwei Zhou
}
