% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RecursiveAnnotationMRN.R
\name{changeTags}
\alias{changeTags}
\title{changeTags: Chages tags (matrix or data.frame) to S4 class peakInfo}
\usage{
changeTags(
  tags,
  ms2,
  weight_mz = 0.25,
  weight_rt = 0.25,
  weight_ccs = 0,
  weight_dp = 0.5,
  mz_tol = 25,
  rt_tol = 30,
  dp_tol = 0.5,
  ccs_tol = 4,
  ...
)
}
\arguments{
\item{tags}{The tags information.}

\item{ms2}{The ms2.}

\item{weight_mz}{mz weight.}

\item{weight_rt}{RT weight}

\item{weight_dp}{Dop product weight}

\item{mz_tol}{The mz tolerance.}

\item{rt_tol}{The RT tolerance.}

\item{dp_tol}{The dp tolerance.}

\item{...}{other parameters.}

\item{adduct_table}{The adduct table.}
}
\value{
Tags2 data.
}
\description{
changeTages is used to change tags information (matrix or data frame) to peakInfo (S4 class).
}
\author{
Zhiwei Zhou, Xiaotao Shen
}
