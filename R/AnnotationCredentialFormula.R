################################################################################
# authenticateFormulaCredential ------------------------------------------------
#' @title authenticateFormulaCredential
#' @author Zhiwei Zhou
#' @description predict formula for peak group
#' @param list_peak_group
#' @param ms2_data
#' @param polarity 'positive', 'negative'
#' @param path_dir Default: '.'
#' @param dir_GenForm directory path of GenForm
#' @param ppm m/z tolerance of MS1; Default: 10 ppm
#' @param acc m/z tolerance of MS2; Default: 15 ppm
#' @param elements CHNOPS
#' @param num_formula_candidate 3
#' @export
#' @examples

# load('./inst/tempdata/list_peak_group_annotation_concised_pg2pg_200805.RData')
# load('./inst/tempdata/raw_msms_200805.RData')
# list_peak_group <- list_peak_group_annotation_concised
# ms2_data <- raw_msms
# path_dir <- 'I:/00_projects/03_MetDNA2/00_data/20200810_conbine_Genform_to_feature_credential/demo_formula'
# polarity <- 'negative'
# authenticateFormulaCredential(list_peak_group = list_peak_group_annotation_concised,
#                               ms2_data = raw_msms,
#                               polarity = 'negative',
#                               path_dir = 'I:/00_projects/03_MetDNA2/00_data/20200810_conbine_Genform_to_feature_credential/demo_formula',
#                               dir_GenForm = 'H:/software/GeneForm',
#                               ppm = 10,
#                               acc = 15,
#                               elements = "CHNOPS",
#                               num_formula_candidate = 5)

setGeneric(name = 'authenticateFormulaCredential',
           def = function(
             list_peak_group,
             ms2_data,
             polarity = c('positive', 'negative'),
             platform = c('windows', 'linux'),
             path_dir = '.',
             thread = 4,
             dir_GenForm = 'G:/software/GeneForm',
             ppm = 10, # ms1 tolerance
             acc = 15, # ms2 tolerance
             elements = "CHNOPS",
             num_formula_candidate = 3 # return top 3 candidate

           ){
             polarity <- match.arg(polarity)
             # platform <- match.arg(platform)

             if (!('GenForm.exe' %in% list.files(dir_GenForm)) & !('GenForm' %in% list.files(dir_GenForm))) {
               stop('GenForm progrem is not found, please check the path\n')
             }

             # change rela path to absolute path, avoid obj missing in runGenForm
             if (nchar(path_dir) == 1) {
               path_dir <- getwd()
             }

             cat('Generate files for formula prediction...\n')
             progress <- mapProgress(n = length(list_peak_group))
             list_peak_group_formula <- purrr::map(seq_along(list_peak_group), function(i){
               mapProgressPrint(progress)
               # cat(i, ' ')
               temp_peak_group <- list_peak_group[[i]]

               idx <- which(names(ms2_data) == temp_peak_group@base_peak_name)

               if (length(idx) == 0) {
                 temp_ms2_spec <- NULL
               } else {
                 temp_ms2_spec <- ms2_data[[idx]]$spec
               }

               result <- generateFiles4FormulaPrediction(peak_group = temp_peak_group,
                                                         ms2_spec = temp_ms2_spec,
                                                         polarity = polarity,
                                                         path_dir = path_dir)

               return(result)
             })

             names(list_peak_group_formula) <- names(list_peak_group)

             dir.create(file.path(path_dir, '03_annotation_credential', "00_intermediate_data"),
                        showWarnings = FALSE,
                        recursive = TRUE)
             save(list_peak_group_formula,
                  file = file.path(path_dir,
                                   '03_annotation_credential',
                                   "00_intermediate_data",
                                   'list_peak_group_formula_empty.RData'),
                  version = 2, compress = 'gzip')



             cat('\n\n'); cat('Predict formula for each peak group...\n')
             progress <- mapProgress(n = length(list_peak_group_formula))
             purrr::walk(seq_along(list_peak_group_formula), function(i){
               # cat(i, ' ')
               mapProgressPrint(progress)
               # browser()

               peak_group <- list_peak_group_formula[[i]]

               temp_path <- file.path(path_dir,
                                      '03_annotation_credential',
                                      '02_formula_prediction',
                                      paste(peak_group@base_peak_name,
                                            peak_group@base_peak_adduct,
                                            sep = '_'))

               result <- try(runGenForm(ms_file = file.path(temp_path, 'ms1_file.txt'),
                                        msms_file = file.path(temp_path, 'ms2_file.txt'),
                                        ion = ifelse(polarity == 'positive', '+H', '-H'),
                                        dir_GenForm = dir_GenForm,
                                        dir_result = temp_path,
                                        platform = platform,
                                        ppm = ppm,
                                        acc = acc,
                                        elements = elements),
                             silent = TRUE)

             })

             cat('\n\n');cat('Rank predicted formula...\n')
             progress <- mapProgress(n = length(list_peak_group_formula))
             list_peak_group_formula <- purrr::map(seq_along(list_peak_group_formula), function(i){
               mapProgressPrint(progress)
               # cat(i, ' ')
               peak_group <- list_peak_group_formula[[i]]

               temp_path <- file.path(path_dir,
                                      '03_annotation_credential',
                                      '02_formula_prediction',
                                      paste(peak_group@base_peak_name,
                                            peak_group@base_peak_adduct,
                                            sep = '_'))

               file_genform <- 'ms2_file_GenForm_sorted.txt'
               if(!file.exists(file.path(temp_path, file_genform))) {
                 return(peak_group)
               }

               result <- rankGenFormFormula(file_genform = file.path(temp_path, file_genform))

               # export ranked formula list
               readr::write_tsv(result,
                                path = file.path(temp_path, 'formula_ranked_lib_formula_final.txt'))

               peak_group@pred_formula_list <- result
               return(peak_group)
             })

             # data('lib_formula', envir = environment())
             #
             # temp_fun <- function(i, path_dir, list_peak_group_formula, rankGenFormFormula, lib_formula, `%>%`) {
             #   # require(packageName())
             #   # require('AnnotationCredential')
             #
             #   # data('lib_formula', package = 'MetDNA2')
             #   peak_group <- list_peak_group_formula[[i]]
             #
             #   temp_path <- file.path(path_dir,
             #                          'annot_credential',
             #                          'formula_prediction',
             #                          paste(peak_group@base_peak_name,
             #                                peak_group@base_peak_adduct,
             #                                sep = '_'))
             #
             #   file_genform <- 'ms2_file_GenForm_sorted.txt'
             #   if(!file.exists(file.path(temp_path, file_genform))) {
             #     return(peak_group)
             #   }
             #
             #   result <- rankGenFormFormula(file_genform = file.path(temp_path, file_genform))
             #
             #   # export ranked formula list
             #   readr::write_tsv(result,
             #                    path = file.path(temp_path, 'formula_ranked_lib_formula_final.txt'))
             #
             #   peak_group@pred_formula_list <- result
             #   return(peak_group)
             # }
             #
             # list_peak_group_formula <- BiocParallel::bplapply(X = seq_along(list_peak_group_formula),
             #                                                   FUN = temp_fun,
             #                                                   BPPARAM = BiocParallel::SnowParam(workers = thread,
             #                                                                                     progressbar = TRUE),
             #                                                   path_dir = path_dir,
             #                                                   list_peak_group_formula = list_peak_group_formula,
             #                                                   rankGenFormFormula = rankGenFormFormula,
             #                                                   lib_formula = lib_formula,
             #                                                   `%>%` = `%>%`)

             names(list_peak_group_formula) <- names(list_peak_group)

             dir.create(file.path(path_dir, '03_annotation_credential', "00_intermediate_data"),
                        showWarnings = FALSE,
                        recursive = TRUE)
             save(list_peak_group_formula,
                  file = file.path(path_dir,
                                   '03_annotation_credential',
                                   "00_intermediate_data",
                                   'list_peak_group_formula.RData'),
                  version = 2, compress = 'gzip')

             cat('\n\n');cat('Generate formula prediction result...\n')
             result_table_formula <- generateFormulaResult(list_peak_group_formula =
                                                             list_peak_group_formula,
                                                           num_formula_candidate =
                                                             num_formula_candidate)

             dir.create(file.path(path_dir, '03_annotation_credential'),
                        showWarnings = FALSE,
                        recursive = TRUE)

             readr::write_csv(result_table_formula,
                              path = file.path(path_dir,
                                               '03_annotation_credential',
                                               'formula_prediction_result.csv'))

           })




################################################################################
# generateFilesFormula4Prediction ----------------------------------------------
#' @title generateFiles4FormulaPrediction
#' @author Zhiwei Zhou
#' @description Generate files (ms1_file & ms2_file) for formula prediction
#' @return a object of PeakGroupFormula
#' @param peak_group
#' @param ms2_spec
#' @param polarity c('positive', 'negative')
#' @param path_dir '.'
#' @export
#' @examples
#' load('./inst/tempdata/list_peak_group_annotation_concised_pg2pg_200805.RData')
#' load('./inst/tempdata/raw_msms_200805.RData')
#' peak_group <- list_peak_group_annotation_concised$`M482T929_[M-H]-`
#' polarity <- 'negative'
#' ms2_spec <- raw_msms$M482T929$spec
#' path_dir <- 'I:/00_projects/03_MetDNA2/00_data/20200810_conbine_Genform_to_feature_credential/demo_formula'
#' test <- generateFiles4FormulaPrediction(peak_group = peak_group,
#'                                         ms2_spec = ms2_spec,
#'                                         polarity = 'negative',
#'                                         path_dir = path_dir)

# load('./inst/tempdata/list_peak_group_annotation_concised_pg2pg_200805.RData')
# load('./inst/tempdata/raw_msms_200805.RData')
# peak_group <- list_peak_group_annotation_concised$`M482T929_[M-H]-`
# polarity <- 'negative'
# ms2_spec <- raw_msms$M482T929$spec
# path_dir <- 'I:/00_projects/03_MetDNA2/00_data/20200810_conbine_Genform_to_feature_credential/demo_formula'
# test <- generateFiles4FormulaPrediction(peak_group = peak_group,
#                                         ms2_spec = ms2_spec,
#                                         polarity = 'negative',
#                                         path_dir = path_dir)

setGeneric(name = 'generateFiles4FormulaPrediction',
           def = function(
             peak_group,
             ms2_spec,
             polarity = c('positive', 'negative'),
             path_dir = '.'
           ){
             dir.create(file.path(path_dir,  '03_annotation_credential', '02_formula_prediction'),
                        showWarnings = FALSE, recursive = TRUE)

             temp_adduct <- ifelse(polarity == 'positive', '[M+H]+', '[M-H]-')

             mz_precursor <- convertMz2Adduct(base_mz = peak_group@base_peak_mz,
                                              base_adduct = peak_group@base_peak_adduct,
                                              adduct_list = temp_adduct) %>%
               dplyr::pull(mz)

             ms1_file <- tibble::tibble(mz = as.numeric(mz_precursor),
                                        intensity = 1)

             dir.create(file.path(path_dir,
                                  '03_annotation_credential',
                                  '02_formula_prediction',
                                  paste(peak_group@base_peak_name,
                                        peak_group@base_peak_adduct,
                                        sep = '_')),
                        showWarnings = FALSE, recursive = TRUE)

             readr::write_delim(x = ms1_file,
                                delim = ' ',
                                path = file.path(path_dir,
                                                 '03_annotation_credential',
                                                 '02_formula_prediction',
                                                 paste(peak_group@base_peak_name,
                                                       peak_group@base_peak_adduct,
                                                       sep = '_'),
                                                 'ms1_file.txt'),
                                col_names = FALSE)

             if (length(ms2_spec) > 0) {
               target_msms <- purifyMs2(spec = ms2_spec,
                                         is_include_precursor = TRUE,
                                         is_deisotope = FALSE,
                                         mz_range_ms2 = c(0, mz_precursor+0.01),
                                         int_ms2_min_abs = 30,
                                         int_ms2_min_relative = 0.01,
                                         mz_precursor = peak_group@base_peak_mz,
                                         ppm_precursor_filter = 10) %>%
                 as.data.frame()


               if (length(target_msms) > 0) {
                 ms2_file <- target_msms %>% tibble::as_tibble()

                 readr::write_delim(x = ms2_file,
                                    delim = ' ',
                                    path = file.path(path_dir,
                                                     '03_annotation_credential',
                                                     '02_formula_prediction',
                                                     paste(peak_group@base_peak_name,
                                                           peak_group@base_peak_adduct,
                                                           sep = '_'),
                                                     'ms2_file.txt'),
                                    col_names = FALSE)
               }

             }

             result <- new('PeakGroupFormula',
                           base_peak_name = peak_group@base_peak_name,
                           base_peak_mz = peak_group@base_peak_mz,
                           base_peak_rt = peak_group@base_peak_rt,
                           base_peak_ccs = peak_group@base_peak_ccs,
                           base_peak_adduct = peak_group@base_peak_adduct)

             return(result)

           })



setClass(Class = "PeakGroupFormula",
         slots = list(base_peak_name = "character",
                      base_peak_mz = "numeric",
                      base_peak_rt = "numeric",
                      base_peak_ccs = 'numeric',
                      base_peak_adduct = 'character',
                      pred_formula_list = 'data.frame')
)


setMethod(f = "show",
          signature = "PeakGroupFormula",
          definition = function(object){
            # cat("-----------Meta information------------\n")
            cat("Base peak name:", object@base_peak_name, "\n")
            cat("m/z:", object@base_peak_mz, "\n")
            cat("RT:", object@base_peak_rt, "\n")
            cat("CCS:", object@base_peak_ccs, "\n")
            cat("Adduct:", object@base_peak_adduct, "\n")
            cat("Top 3 formula:",
                paste(object@pred_formula_list$formula[1:3], collapse = ';'), "\n")
            cat("No. of predicted formula:",
                nrow(object@pred_formula_list), "\n")
          }
)




################################################################################
# runGenForm -------------------------------------------------------------------
# modified from GenFormR [https://github.com/schymane/GenFormR]

#' @title runGenForm
#' @author Emma Schymanski (R wrapper, \code{\link{emma.schymanski@@uni.lu}}) and Markus Meringer (GenForm); Zhiwei Zhou (Modified for MetDNA2)
#' @description This is a wrapper function for GenForm, requiring MS, MS/MS, ion settings
#' and directories. This function is modified from GenFormR
#' @param ms_file File name of a plain text, two column file (e.g. \code{MS.txt}) containing
#' space-separated mz and intensity values of the MS1 spectrum (full scan). The most intense
#' peak is taken as "m" unless this is defined with \code{"mz"}.
#' @param msms_file A plain text, two column file (e.g. \code{MSMS.txt}) containing
#' space-separated mz and intensity values of the tandem (MS2) spectrum, i.e. the fragments.
#' While GenForm can be run without MSMS, this function requires its use at this stage.
#' @param ion Ion setting to use. Options: \code{c("-e","+e","-H","+H","+Na")}.
#' @param dir_GenForm Directory containing \code{GenForm.exe}.
#' @param dir_result Directory to save the GenForm results. Existing files will be overwritten.
#' @param ppm Accuracy setting for MS1 (default \code{5} ppm).
#' @param acc Accuracy setting for MS2 (default \code{15} ppm).
#' @param elements Selection of elements to consider. Default \code{"CHNOPS"}, implemented
#' are \code{"CHBrClFINOPSSi"}.
#' @param ff A "fuzzy formula" input to offer more control beyond \code{elements}. This should be
#' in format \code{C0-6H0-20O1-3}.
#' @param dbe Default \code{TRUE} writes double bond equivalents to output; \code{FALSE} suppresses
#' this output.
#' @param oei Default \code{TRUE} allows odd electron ions for explaining MS/MS peaks, \code{FALSE}
#' means only even electron ions are considered.
#' @param mz If mz is >0, this is the mass GenForm will use to calculate the formulas, considering the \code{ion} setting. Default is the most intense peak in \code{ms_file}.
#' @param exist Default \code{TRUE} activates the existance filter (allowing only formulae for which at least one structural formula exists). \code{FALSE} is not recommended for most use cases.
#' @param an_loss Activate the annotation of MSMS fragments with subformulas and losses.
#' @param an_loss_limit Define the cut-off limit for annotating subformulas and losses. For many peaks and formulas, this becomes time consuming and produces massive files of limited use.
#' @param sort Default \code{FALSE} will produce default output. If \code{TRUE}, generates an additional
#' sorted file sorted by combined match value.
#' @param cleanMSMS Default \code{FALSE}. If \code{TRUE}, produces an additional cleaned MSMS file containing only peaks with subformula annotation.
#' @return Several files containing various GenForm outputs.
#' @references GenForm on Source Forge: https://sourceforge.net/projects/genform/
#' @export
#' @examples
#' runGenForm(ms_file = 'I:/00_projects/03_MetDNA2/00_data/20200810_conbine_Genform_to_feature_credential/example/SinapinicAcidMs.txt',
#'            msms_file = 'I:/00_projects/03_MetDNA2/00_data/20200810_conbine_Genform_to_feature_credential/example/SinapinicAcidMsMs.txt',
#'            ion = c('+H'),
#'            dir_GenForm = 'H:/software/GeneForm',
#'            dir_result = 'I:/00_projects/03_MetDNA2/00_data/20200810_conbine_Genform_to_feature_credential/example/result2',
#'            ppm = 10,
#'            elements = "CHNOPS",
#'            ff = "",
#'            dbe = TRUE,
#'            oei = TRUE,
#'            mz = 0,
#'            exist = TRUE,
#'            an_loss = TRUE,
#'            an_loss_limit = 500,
#'            sort = TRUE,
#'            cleanMSMS = FALSE
#'            )


# runGenForm(ms_file = 'I:/00_projects/03_MetDNA2/00_data/20200810_conbine_Genform_to_feature_credential/example/SinapinicAcidMs.txt',
#            msms_file = 'I:/00_projects/03_MetDNA2/00_data/20200810_conbine_Genform_to_feature_credential/example/SinapinicAcidMsMs.txt',
#            ion = c('+H'),
#            dir_GenForm = 'H:/software/GeneForm',
#            dir_result = 'I:/00_projects/03_MetDNA2/00_data/20200810_conbine_Genform_to_feature_credential/example/result2',
#            ppm = 10,
#            elements = "CHNOPS",
#            ff = "",
#            dbe = TRUE,
#            oei = TRUE,
#            mz = 0,
#            exist = TRUE,
#            an_loss = TRUE,
#            an_loss_limit = 500,
#            sort = TRUE,
#            cleanMSMS = FALSE
#            )

# runGenForm(ms_file = 'I:/00_projects/03_MetDNA2/00_data/20200810_conbine_Genform_to_feature_credential/example/SinapinicAcidMs_Na.txt',
#            msms_file = 'I:/00_projects/03_MetDNA2/00_data/20200810_conbine_Genform_to_feature_credential/example/SinapinicAcidMsMs.txt',
#            ion = c('+Na'),
#            dir_GenForm = 'H:/software/GeneForm',
#            dir_result = 'I:/00_projects/03_MetDNA2/00_data/20200810_conbine_Genform_to_feature_credential/example/result_Na2',
#            ppm = 10,
#            elements = "CHNOPS",
#            ff = "",
#            dbe = TRUE,
#            oei = TRUE,
#            mz = 0,
#            exist = TRUE,
#            an_loss = TRUE,
#            an_loss_limit = 500,
#            sort = TRUE,
#            cleanMSMS = FALSE
#            )

# ms_file = 'I:/00_projects/03_MetDNA2/00_data/20200810_conbine_Genform_to_feature_credential/example/SinapinicAcidMs.txt'
# msms_file = 'I:/00_projects/03_MetDNA2/00_data/20200810_conbine_Genform_to_feature_credential/example/SinapinicAcidMsMs.txt'
# dir_result = 'I:/00_projects/03_MetDNA2/00_data/20200810_conbine_Genform_to_feature_credential/example/result'
#
# ms_file = '/home/zhouzw/Data_processing/20201019_formula_prediction_test_ubuntu/SinapinicAcidMs.txt'
# msms_file = '/home/zhouzw/Data_processing/20201019_formula_prediction_test_ubuntu/SinapinicAcidMsMs.txt'
# ion = c('+H')
# dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/bin/Debug'
# # dir_GenForm = system.file("extdata", package="MetDNA2I")
# dir_result = '/home/zhouzw/Data_processing/20201019_formula_prediction_test_ubuntu/result'
# ppm=10
# acc=25
# elements="CHNOPS"
# ff=""
# dbe=TRUE
# oei=TRUE
# mz=0
# exist=TRUE
# an_loss=TRUE
# an_loss_limit=500
# sort=TRUE
# cleanMSMS=FALSE

setGeneric(name = 'runGenForm',
           def = function(
             ms_file,
             msms_file,
             ion,
             dir_GenForm,
             dir_result,
             ppm=5,
             acc=15,
             elements="CHNOPS",
             ff="",
             dbe=TRUE,
             oei=TRUE,
             mz=0,
             exist=TRUE,
             an_loss=TRUE,
             an_loss_limit=500,
             sort=TRUE,
             cleanMSMS=FALSE,
             platform = c('windows', 'linux')
           ){
             curr_dir <- getwd()

             # check files exist
             if(!file.exists(ms_file)) {
               stop(paste(ms_file, "doesn't exist, please try again"),sep=" ")
             }
             if(!file.exists(msms_file)) {
               stop(paste(msms_file, "doesn't exist, please try again"),sep=" ")
             }
             if(!file.exists(dir_GenForm)) {
               stop(paste(dir_GenForm, "doesn't exist, please try again"),sep=" ")
             }
             if(!file.exists(dir_result)) {
               dir.create(dir_result)
             }
             possible_ions <- c("+H","-H","-e","+e","+Na")
             if(!(ion %in% possible_ions)) {
               stop(paste("Incorrect ion state defined, try one of", paste(possible_ions,collapse=", ")))
             }


             ms_file <- ms_file
             msms_file <- msms_file
             msms_file_name <- basename(msms_file)
             out_file <- paste(dir_result,"/",sub(".txt","_GenForm.txt",msms_file_name),sep="")
             log_file <- paste(dir_result,"/",sub(".txt","_GenForm_log.txt",msms_file_name),sep="")
             other_options <- ""
             if(dbe) other_options <- paste(other_options, "dbe")
             if(oei) other_options <- paste(other_options, "oei")
             if(exist) other_options <- paste(other_options, "exist")
             #   if(!grepl("none",sort)) other_options <- paste(other_options," sort=",sort,sep="")
             # check if mz is way too high
             ms <- read.table(ms_file)

             # build the command
             if (nchar(ff)>0 && mz>0) {
               GenFormBaseCommand <- paste("GenForm ms=", ms_file,
                                           " msms=", msms_file,
                                           " ion=", ion,
                                           " ppm=", ppm,
                                           " acc=", acc,
                                           " m=", mz,
                                           " ff=",ff,
                                           other_options,
                                           sep="")
             } else if (nchar(ff)>0 && mz==0) {
               GenFormBaseCommand <- paste("GenForm ms=", ms_file,
                                           " msms=", msms_file,
                                           " ion=", ion,
                                           " ppm=", ppm,
                                           " acc=", acc,
                                           " ff=", ff,
                                           other_options,
                                           sep="")
             } else if (nchar(ff)==0 && mz>0) {
               GenFormBaseCommand <- paste("GenForm ms=", ms_file,
                                           " msms=", msms_file,
                                           " ion=", ion,
                                           " ppm=", ppm,
                                           " acc=", acc,
                                           " m=", mz,
                                           " el=", elements,
                                           other_options,
                                           sep="")
             } else {
               GenFormBaseCommand <- paste("GenForm ms=", ms_file,
                                           " msms=", msms_file,
                                           " ion=", ion,
                                           " ppm=", ppm,
                                           " acc=", acc,
                                           " el=", elements,
                                           other_options,
                                           sep="")
             }

             GenFormCommand <- paste(GenFormBaseCommand,
                                     " out=", out_file,
                                     sep="")

             setwd(dir_GenForm)
             switch (platform,
               'linux' = {
                 system(paste0('./', GenFormCommand))
                 # system.file("extdata", "GenForm", ackage="AllCcsLib")
               },
               'windows' = {
                 system(GenFormCommand, show.output.on.console = FALSE)
               }
             )

             write(GenFormCommand, log_file, append=TRUE)

             # test to see whether to run analyze/loss
             mz_to_test <- ms[which.max(ms[,2]), 1]
             if(an_loss && (mz_to_test < an_loss_limit)) {
               out_file_an <- paste(dir_result, "/",
                                    sub(".txt","_GenForm_an.txt", msms_file_name),
                                    sep="")
               out_file_loss <- paste(dir_result, "/",
                                      sub(".txt","_GenForm_loss.txt", msms_file_name),
                                      sep="")
               GenFormCommand_an <- paste(GenFormBaseCommand, " analyze out=",
                                          out_file_an, sep="")
               GenFormCommand_loss <- paste(GenFormBaseCommand, " analyze loss out=",
                                            out_file_loss, sep="")
               system(GenFormCommand_an, show.output.on.console = FALSE)
               system(GenFormCommand_loss, show.output.on.console = FALSE)

               write(GenFormCommand_an, log_file, append=TRUE)
               write(GenFormCommand_loss, log_file, append=TRUE)

             }

             # else {
             #   print("No analyze/loss output performed, either mz is too large or option is off.")
             # }

             #create a sorted output file too
             if (sort) {
               sorted_out_file <- file.path(dir_result,
                                            sub(".txt","_GenForm_sorted.txt",
                                                msms_file_name))
               output <- read.table(out_file)
               o <- order(output$V6, decreasing=TRUE)
               ordered_out <- output[o,]
               write.table(ordered_out, sorted_out_file,
                           quote=FALSE, col.names=FALSE,
                           row.names=FALSE, sep="\t")
             }
             if (cleanMSMS) {
               cleanedMSMS_file <- file.path(dir_result,
                                             sub(".txt","_GenForm_cleanMSMS.txt",
                                                 msms_file_name))

               GenFormCommand_clean <- paste(GenFormBaseCommand,
                                             " oclean=",
                                             cleanedMSMS_file, " out", sep="")
               system(GenFormCommand_clean, show.output.on.console = FALSE)
               write(GenFormCommand_clean, log_file,append=TRUE)
             }


             setwd(curr_dir)

           })



#' @title annotateCleanMSMS
#' @description This takes the output from \code{\link{RunGenForm}} and creates a peaklist with the annotations for further use.
#' @usage annotateCleanMSMS(cleanMSMS, GenForm_an_out, GenForm_loss_out, lowestPpm=TRUE, round=TRUE, removeAnnoSpace=FALSE)
#' @param cleanMSMS File name for Cleaned MSMS text file output from \code{\link{RunGenForm}}
#' @param GenForm_an_out File name for annotated MSMS text file output from \code{\link{RunGenForm}}
#' @param GenForm_loss_out File name for annotated "loss" MSMS text file output from \code{\link{RunGenForm}}
#' @param lowestPpm If \code{TRUE}, choses the lowest ppm subformula for assignment of formula
#' to fragments. If \code{FALSE}, takes the first formula.
#' @param round If \code{TRUE}, rounds the output to 4 (mz) and 1 (Int) decimal places for aesthetics.
#' If \code{FALSE}, no rounding is done.
#' @param removeAnnoSpace If \code{TRUE}, removes the space in the \code{"no loss"} annotation.
#' @return Returns annotated peaks for further use
#' @author Emma Schymanski (\code{\link{emma.schymanski@@uni.lu}})
#' @export
#' @examples
#'
annotateCleanMSMS <- function(cleanMSMS, GenForm_an_out, GenForm_loss_out, lowestPpm=TRUE,
                              round=TRUE, removeAnnoSpace=FALSe) {
  peaks <- read.table(cleanMSMS)
  colnames(peaks) <- c("mz", "Int")
  peaks$Int <- 999*peaks$Int
  subform_an <- read.table(GenForm_an_out, header=F, sep="\t", skip=1)
  colnames(subform_an) <- c("mz", "subform_an", "ppm")
  subform_loss <- read.table(GenForm_loss_out, header=F, sep="\t", skip=1)
  colnames(subform_loss) <- c("mz", "subform_loss", "ppm")
  subform_an[,4] <- subform_loss$subform_loss
  subform_an[,5] <- 0
  colnames(subform_an) <- c("mz", "subform_an", "ppm", "subform_loss", "Int")
  # replace missing mzs with numbers, add intensity
  for (i in 1:length(subform_an$mz)) {
    if (is.na(subform_an$mz[i])) {
      subform_an$mz[i] <- subform_an$mz[(i-1)]
      #subform_loss$mz[i] <- subform_loss$mz[(i-1)]
    }
    subform_an$Int[i] <- peaks$Int[which(peaks$mz==subform_an$mz[i])]
  }
  # now deal with duplicate entries ...
  peaks$subform_an <- ""
  peaks$subform_loss <- ""
  peaks$ppm <- ""
  peaks$n_formulas <- 0
  for (n in 1:length(peaks$mz)) {
    subform_an_entries <- subform_an[which(subform_an$mz==peaks$mz[n]),]
    if (length(subform_an_entries$mz)==1) {
      peaks$subform_an[n] <- as.character(subform_an_entries$subform_an[1])
      peaks$subform_loss[n] <- as.character(subform_an_entries$subform_loss[1])
      peaks$ppm[n] <- subform_an_entries$ppm[1]
      peaks$n_formulas[n] <- 1
    } else if (length(subform_an_entries$mz)>1) {
      # choose either first entry or smallest ppm value
      if (lowestPpm) {
        ppm_ind <- which(abs(subform_an_entries$ppm)==min(abs(subform_an_entries$ppm)))
        peaks$subform_an[n] <- as.character(subform_an_entries$subform_an[ppm_ind[1]])
        peaks$subform_loss[n] <- as.character(subform_an_entries$subform_loss[ppm_ind[1]])
        peaks$ppm[n] <- subform_an_entries$ppm[ppm_ind[1]]
        peaks$n_formulas[n] <- length(subform_an_entries$mz)
      } else { #take first entry
        peaks$subform_an[n] <- as.character(subform_an_entries$subform_an[1])
        peaks$subform_loss[n] <- as.character(subform_an_entries$subform_loss[1])
        peaks$ppm[n] <- subform_an_entries$ppm[1]
        peaks$n_formulas[n] <- length(subform_an_entries$mz)
      }
    }
  }
  # adjust the rounding ... to make output look nicer
  if (round) {
    peaks$Int <- round(peaks$Int,1)
    peaks$mz <- round(peaks$mz,4)
  }

  return(peaks)

}

# Fix Whitespace Error in GenForm Loss Annotation
#' @title fixNoLossSpace
#' @description This takes a filename with e.g. the output from \code{\link{annotateCleanMSMS}} and removes a space in the \code{"no loss"} annotation to allow for easier handling.
#' @usage fixNoLossSpace(anno_MS_file, replace=TRUE)
#' @param anno_MS_file File name for annotated cleaned MSMS text file output from \code{\link{annotateCleanMSMS}}
#' @param replace If \code{TRUE}, replace the input file. If \code{FALSE}, saved under the new name
#' \code{"*_ed.txt"}.
#' @return File name of the "fixed" file.
#' @author Emma Schymanski (\code{\link{emma.schymanski@@uni.lu}})
#' @seealso \code{\link{annotateCleanMSMS}}
#' @export
#' @examples

fixNoLossSpace <- function(anno_MS_file, replace=TRUE) {
  fileConnection <- file(anno_MS_file)
  anno_MS_lines <- readLines(fileConnection)
  close(fileConnection)
  anno_MS_lines[length(anno_MS_lines)] <- sub("no loss", "noLoss",
                                              anno_MS_lines[length(anno_MS_lines)], fixed=TRUE)
  if (replace) {
    file_name <- anno_MS_file
  } else {
    file_name <- sub(".txt", "_ed.txt",anno_MS_file, fixed=TRUE)
  }
  write.table(anno_MS_lines, file_name, row.names=F, col.names = F, quote=F)
  return(file_name)
}




################################################################################
# rankGenFormFormula -----------------------------------------------------------
#' @title rankGenFormFormula
#' @description rank predicted formula with known formula DB
#' @author Zhiwei Zhou
#' @param file_genform the file name of GenForm result
#' @export
#' @examples
#' rankGenFormFormula(file_genform = 'I:/00_projects/03_MetDNA2/00_data/20200810_conbine_Genform_to_feature_credential/example/result/SinapinicAcidMsMs_GenForm_sorted.txt')
# rankGenFormFormula(file_genform = 'I:/00_projects/03_MetDNA2/00_data/20200810_conbine_Genform_to_feature_credential/example/result/SinapinicAcidMsMs_GenForm_sorted.txt')

setGeneric(name = 'rankGenFormFormula',
           def = function(
             file_genform = 'ms2_file_GenForm_sorted.txt'
           ){

             # browser()
             options(readr.num_columns = 0)
             result_genform <- readr::read_tsv(file_genform,
                                               col_names = FALSE)

             colnames(result_genform) <- c('formula',
                                           'num_frag_msms',
                                           'ms1_error',
                                           'score_ms1',
                                           'score_ms2',
                                           'score_combined')

             result_genform <- result_genform %>%
               dplyr::rowwise() %>%
               dplyr::mutate(score_formula_db = searchFormulaDB(target_formula = formula)) %>%
               dplyr::mutate(score_inte = score_combined + score_formula_db) %>%
               dplyr::arrange(dplyr::desc(score_inte)) %>%
               tibble::as_tibble()

             # lapply(result_genform$formula, function(x){
             #   cat(x, ' ')
             #   searchFormulaDB(target_formula = x)
             # })

             return(result_genform)
           }
)


#' @title searchFormulaDB
#' @description search formula in known compound formula DB to plus additional score for ranking. The MsfinderFormulaDB-VS10.efd (MSFINDER DB) was used as known compound DB.
#' @author Zhiwei Zhou
#' @return additional score
#' @param target_formula
#' @param score_bio_db Default: 0.3
#' @param score_pred_db Default: 0.15
#' @export
#' @examples
#' searchFormulaDB(target_formula = 'C6H12O6')

# searchFormulaDB(target_formula = 'C6H12O6')

setGeneric(name = 'searchFormulaDB',
           def = function(target_formula,
                          # lib_formula,
                          score_bio_db = 0.3,
                          score_pred_db = 0.15,
                          ...){

             result <- lib_formula %>%
               dplyr::filter(formula == target_formula)

             if (nrow(result) == 0) {
               return(0)
             }

             num_records <- result %>% dplyr::pull(records)
             is_mine_db <- result %>% dplyr::pull(mine_included)
             is_emrn_db <- result %>% dplyr::pull(emrn_included)

             if (num_records > 2) {
               return(score_bio_db)
             } else if(all(is.na(is_mine_db), is.na(is_emrn_db))) {
               return(score_bio_db)
             } else {
               return(score_pred_db)
             }

           }
)

################################################################################
# generateFormulaResult --------------------------------------------------------
#' @title generateFormulaResult
#' @author Zhiwei Zhou
#' @param list_peak_group_formula
#' @param num_formula_candidate top N formula reserved. Default: 5
#' @export
#' @examples
#' load('./inst/tempdata/list_peak_group_formula.RData')
#' test <- generateFormulaResult(list_peak_group_formula, num_formula_candidate = 3)

# load('./inst/tempdata/list_peak_group_formula.RData')
# test <- generateFormulaResult(list_peak_group_formula, num_formula_candidate = 3)

setGeneric(name = 'generateFormulaResult',
           function(
             list_peak_group_formula,
             num_formula_candidate = 5
           ){
             result_table <- purrr::map(seq_along(list_peak_group_formula), function(i){
               # cat(i, ' ')
               temp_data <- list_peak_group_formula[[i]]

               if (nrow(temp_data@pred_formula_list) == 0) {
                 return(NULL)
               }

               temp_formula_list <- temp_data@pred_formula_list %>%
                 dplyr::arrange(desc(score_inte)) %>%
                 dplyr::slice(seq(num_formula_candidate))

               result <- tibble::tibble(name = temp_data@base_peak_name,
                                        mz = temp_data@base_peak_mz,
                                        rt = temp_data@base_peak_rt,
                                        ccs = temp_data@base_peak_ccs,
                                        adduct = temp_data@base_peak_adduct,
                                        formula = paste(temp_formula_list$formula,
                                                        collapse = ';'),
                                        score = paste(temp_formula_list$score_inte,
                                                      collapse = ';'))

               return(result)
             })

             result_table <- result_table %>% dplyr::bind_rows()
             return(result_table)

           })
