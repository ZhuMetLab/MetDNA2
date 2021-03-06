% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AnnotationCredentialPeakGroup.R
\name{authenticatePeakGroupCredential}
\alias{authenticatePeakGroupCredential}
\title{authenticatePeakGroupCredential}
\usage{
authenticatePeakGroupCredential(
  list_peak_group,
  ms2_data,
  path_dir = ".",
  polarity = c("positive", "negative"),
  tol_mz = 25,
  isotope_int_ratio_check = TRUE,
  isotope_int_ratio_cutoff = 500,
  is_ms2_check = TRUE,
  ms2_score_cutoff = -1,
  cutoff_ssc = 0.3,
  cutoff_ssc_int = 3000,
  is_rule_limitation = TRUE,
  cutoff_topN = 5,
  is_plot_pseudo_MS1 = TRUE,
  thread = 4,
  type_order = c("level1", "level2", "level3"),
  ...
)
}
\arguments{
\item{path_dir}{'.'}

\item{polarity}{ionzation polarity, 'positive', 'negative'; Default: 'positive'}

\item{tol_mz}{mz tolerance. Default: 10 ppm}

\item{isotope_int_ratio_check}{whether check isotope intensity; Default: TRUE}

\item{isotope_int_ratio_cutoff}{isotope intensity ratio cutoff; Default: 500%}

\item{is_ms2_check}{whether compare ms2 of NL with base peak; Default: TRUE}

\item{ms2_score_cutoff}{Default: -1; # -1 represent not filter}

\item{is_plot_pseudo_MS1}{whether output the pseudo MS1 spectrum.}
}
\description{
semi-targeted feature clustering (isotope, adduct, neutral loss, in-source fragments)
}
\examples{
load(system.file("tempdata", "list_peak_group_200805.RData", package="MetDNA2"))
load(system.file("tempdata", "raw_msms_200805.RData", package="MetDNA2"))
authenticatePeakGroupCredential(list_peak_group = list_peak_group,
                                ms2_data = raw_msms,
                                path_dir = 'I:/00_projects/03_MetDNA2/00_data/20200805_targeted_annotation_evaluation_200STD/00_demo_data',
                                polarity = 'negative',
                                tol_mz = 10,
                                isotope_int_ratio_check = TRUE,
                                isotope_int_ratio_cutoff = 500,
                                is_ms2_check = TRUE,
                                ms2_score_cutoff = -1,
                                is_plot_pseudo_MS1 = TRUE)
}
\author{
Zhiwei Zhou
}
