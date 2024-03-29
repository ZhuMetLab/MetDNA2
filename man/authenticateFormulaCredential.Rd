% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AnnotationCredentialFormula.R
\name{authenticateFormulaCredential}
\alias{authenticateFormulaCredential}
\title{authenticateFormulaCredential}
\usage{
authenticateFormulaCredential(
  list_peak_group,
  ms2_data,
  polarity = c("positive", "negative"),
  platform = c("windows", "linux"),
  path_dir = ".",
  thread = 4,
  dir_GenForm = "G:/software/GeneForm",
  ppm = 10,
  acc = 15,
  elements = "CHNOPS",
  num_formula_candidate = 3
)
}
\arguments{
\item{polarity}{'positive', 'negative'}

\item{path_dir}{Default: '.'}

\item{dir_GenForm}{directory path of GenForm}

\item{ppm}{m/z tolerance of MS1; Default: 10 ppm}

\item{acc}{m/z tolerance of MS2; Default: 15 ppm}

\item{elements}{CHNOPS}

\item{num_formula_candidate}{3}
}
\description{
predict formula for peak group
}
\author{
Zhiwei Zhou
}
