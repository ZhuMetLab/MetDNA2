% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RecursiveAnnotationInitialSeed.R
\name{convertSpecAnnotationClass2Table}
\alias{convertSpecAnnotationClass2Table}
\title{convertSpecAnnotationClass2Table}
\usage{
convertSpecAnnotationClass2Table(
  ms1_data,
  result_annotation,
  lib_meta,
  instrument = c("SciexTripleTOF", "AgilentQTOF", "BrukerQTOF", "ThermoOrbitrap",
    "ThermoExploris", "WatersQTOF", "WatersTWIMMS", "AgilentDTIMMS", "BrukerTIMS"),
  is_rt_score = TRUE,
  is_ccs_score = FALSE,
  is_msms_score = TRUE,
  rt_cutoff = 30,
  dp_cutoff = 0.8,
  matched_frag_cutoff = 1,
  ...
)
}
\description{
convertSpecAnnotationClass2Table
}
