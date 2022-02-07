################################################################################
# authenticateAnnotationCredential ---------------------------------------------

#' @title authenticateAnnotationCredential
#' @description authenticate annotations with semi-targeted feature clustering and formula prediction.
#' @author Zhiwei Zhou
#' @param peak_table_file 'peak_table.csv'
#' @param annotation_initial_file 'annotation_initial.csv'
#' @param ms2_data_file 'ms2_data.RData'
#' @param path_dir path of working dictory. Default: '.'
#' @param thread number of thread. Default: 4
#' @param tol_rt rt tolerance for peak group generation, Default: 3s
#' @param cutoff_ssc sample-sample correlation for peak group generation, Default: 0.3
#' @param polarity ionzation polarity, 'positive', 'negative'; Default: 'positive'
#' @param tol_mz mz tolerance. Default: 10 ppm
#' @param isotope_int_ratio_check whether check isotope intensity; Default: TRUE
#' @param isotope_int_ratio_cutoff isotope intensity ratio cutoff; Default: 500%
#' @param is_ms2_check whether compare ms2 of NL with base peak; Default: TRUE
#' @param ms2_score_cutoff Default: -1; # -1 represent not filter
#' @param is_plot_pseudo_MS1 whether output the pseudo MS1 spectrum.
#' @param dir_GenForm directory path of GenForm
#' @param is_pred_formula_all whether predict formula for all peak group (TRUE) or concised peak groups (FALSE). Default: FALSE
#' @param ppm m/z tolerance of MS1; Default: 10 ppm
#' @param acc m/z tolerance of MS2; Default: 15 ppm
#' @param elements CHNOPS
#' @param num_formula_candidate 3
#' @importFrom magrittr '%>%' '%$%'
#' @export
#' @examples
#' peak_table_file <- system.file("extdata", "peak_table.csv", package="MetDNA2")
#' annotation_initial_file <- system.file("extdata", "annotation_initial.csv", package="MetDNA2")
#' ms2_data <- system.file("tempdata", "ms2_data.RData", package="MetDNA2")
#' authenticateAnnotationCredential(peak_table_file = "peak_table.csv",
#'                                  annotation_initial_file = "annotation_initial.csv",
#'                                  ms2_data = "ms2_data.RData",
#'                                  path_dir = 'I:/00_projects/03_MetDNA2/00_data/20200811_annotationCredential_workflow/demo_data2',
#'                                  polarity = 'negative')




# peak_table_file <- system.file("extdata", "peak_table.csv", package="MetDNA2")
# annotation_initial_file <- system.file("extdata", "annotation_initial.csv", package="MetDNA2")
# ms2_data <- system.file("tempdata", "ms2_data.RData", package="MetDNA2")
# authenticateAnnotationCredential(peak_table_file = "peak_table.csv",
#                                  annotation_initial_file = "annotation_initial.csv",
#                                  ms2_data_file = "ms2_data.RData",
#                                  path_dir = 'I:/00_projects/03_MetDNA2/00_data/20200811_annotationCredential_workflow/demo_data2',
#                                  polarity = 'negative')


# peak_table_file <- 'peak_table.csv'
# annotation_initial_file <- 'annotation_initial.csv'
# ms2_data_file <- 'ms2_data.RData'
# path_dir <- '/home/zhouzw/Data_processing/20210319_different_instrument_test/QE_mice_liver/SNCE_rerun/'
# # path_dir <- '/home/zhouzw/Data_processing/20210315_tims_data_processing_metdna/'
# ## path_dir <- '/home/zhouzw/Data_processing/20201019_fly_aging_pos'
# ## polarity <- 'negative'
# polarity <- 'negative'
# tol_mz = 25
# thread = 4
# isotope_int_ratio_check = TRUE # para annotateIsotope
# isotope_int_ratio_cutoff = 500
# is_ms2_check = TRUE # compare ms2 similarity of neutral loss
# ms2_score_cutoff = -1 # -1 represent not filter
# is_plot_pseudo_MS1 = TRUE
# cutoff_ssc = 0.3
# cutoff_ssc_int = 3000
# is_rule_limitation = TRUE
# cutoff_topN = 5
# is_plot_pseudo_MS1 = TRUE
# dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release'
# ppm = 20 # ms1 tolerance
# acc = 20 # ms2 tolerance
# elements = "CHNOPS"
# num_formula_candidate = 3 # return top 3 candidate
# is_pred_formula_all = FALSE
# platform = 'linux'
# type_order = c('level1', 'level2', 'level3')
# tol_rt = 3


setGeneric(name = 'authenticateAnnotationCredential',
           def = function(
             peak_table_file,
             annotation_initial_file,
             ms2_data_file = 'ms2_data.RData',
             path_dir = '.',
             thread = 4,
             tol_rt = 3, # pg generation
             # cutoff_ssc = 0.3, # pg generation

             # feature clustering
             polarity = c('positive', 'negative'),
             tol_mz = 25,
             isotope_int_ratio_check = TRUE, # para annotateIsotope
             isotope_int_ratio_cutoff = 500,
             is_ms2_check = TRUE, # compare ms2 similarity of neutral loss
             ms2_score_cutoff = -1, # -1 represent not filter
             cutoff_ssc = 0.3,
             cutoff_ssc_int = 3000,
             is_rule_limitation = TRUE,
             cutoff_topN = 5,
             is_plot_pseudo_MS1 = TRUE,
             type_order = c('level1', 'level2', 'level3'),

             # formula prediction
             is_pred_formula_all = FALSE, # predict formula for all peak group
             dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/bin/Debug',
             ppm = 10, # ms1 tolerance
             acc = 15, # ms2 tolerance
             elements = "CHNOPS",
             num_formula_candidate = 3, # return top 3 candidate
             platform = c('windows', 'linux')

           ){
             polarity <- match.arg(polarity)
             # browser()

             list_para <- list(
               # 'package_verison' = packageVersion('AnnotationCredential'),
                               'peak_table_file' = peak_table_file,
                               'annotation_initial_file' = paste(annotation_initial_file, collapse = ';'),
                               'ms2_data_file' = ms2_data_file,
                               'path_dir' = path_dir,
                               'thread' = thread,
                               'polarity' = polarity,
                               'tol_mz' = tol_mz,
                               'isotope_int_ratio_check' = isotope_int_ratio_check,
                               'isotope_int_ratio_cutoff' = isotope_int_ratio_check,
                               'is_ms2_check' = is_ms2_check,
                               'ms2_score_cutoff' = ms2_score_cutoff,
                               'cutoff_ssc' = cutoff_ssc,
                               'cutoff_ssc_int' = cutoff_ssc_int,
                               'is_rule_limitation' = is_rule_limitation,
                               'cutoff_topN' = cutoff_topN,
                               'is_plot_pseudo_MS1' = is_plot_pseudo_MS1,
                               'dir_GenForm' = dir_GenForm,
                               'is_pred_formula_all' = is_pred_formula_all,
                               'ppm' = ppm,
                               'acc' = acc,
                               'elements' = elements,
                               'num_formula_candidate' = num_formula_candidate) %>%
               tibble::as_tibble() %>%
               dplyr::mutate_all(as.character) %>%
               tidyr::pivot_longer(cols = dplyr::everything(), names_to = 'para')

             dir.create(file.path(path_dir, '03_annotation_credential'),
                        showWarnings = FALSE,
                        recursive = TRUE)

             readr::write_tsv(list_para, path = file.path(path_dir, '03_annotation_credential', 'para_list.txt'))

             # 01 load data ----------------------------------------------------
             cat('=============================================================\n')
             # cat('AnnotationCredential (Verison', as.character(packageVersion('AnnotationCredential')), ')\n')
             cat('Start run AnnotationCredential... \n')
             cat('=============================================================\n\n')

             options(readr.num_columns = 0)
             peak_table <- readr::read_csv(file.path(path_dir, '03_annotation_credential', peak_table_file)) %>%
               dplyr::distinct(name, .keep_all = TRUE) %>%
               dplyr::mutate(name = as.character(name),
                             ccs = as.numeric(ccs))

             options(readr.num_columns = 0)
             annotation_initial <- readr::read_csv(file.path(path_dir, '03_annotation_credential', annotation_initial_file)) %>%
               dplyr::mutate(name = as.character(name),
                             ccs = as.numeric(ccs))

             # load ms2
             load(file.path(path_dir, '03_annotation_credential', ms2_data_file))

             dir.create(file.path(path_dir, '03_annotation_credential', "00_intermediate_data"),
                        showWarnings = FALSE,
                        recursive = TRUE)

             save(annotation_initial,
                  file = file.path(path_dir,
                                   '03_annotation_credential',
                                   "00_intermediate_data",
                                   'annotation_initial.RData'),
                  version = 2)

             # 02 peak grouping ------------------------------------------------
             message('Divide peak table and annotation table into peak group...\n')

             if ('list_peak_group.RData' %in%
                 list.files(file.path(path_dir,
                                      '03_annotation_credential',
                                      "00_intermediate_data"))) {
               cat('Skip peak grouping.\n')
               load(file.path(path_dir,
                              '03_annotation_credential',
                              "00_intermediate_data",
                              'list_peak_group.RData'))

             } else {
               # peak_table_annotated <- annotation_initial %>%
               #   dplyr::distinct(name, adduct, .keep_all = TRUE) %>%
               #   dplyr::select(name, adduct)
               #
               # # add seed tag in PeakGroup
               # if ('type' %in% colnames(annotation_initial)){
               #   seed_peak_list <- annotation_initial %>%
               #     dplyr::filter(type == 'seed') %>%
               #     dplyr::distinct(name, adduct, .keep_all = TRUE) %>%
               #     dplyr::select(name, adduct, type)
               #
               #   peak_table_annotated <- peak_table_annotated %>%
               #     dplyr::left_join(seed_peak_list, by = c('name', 'adduct')) %>%
               #     tidyr::replace_na(list(type = 'none_seed'))
               # } else {
               #   peak_table_annotated <- peak_table_annotated %>%
               #     dplyr::mutate(type = 'none_seed')
               # }

               # add seed tag in PeakGroup
               if ('type' %in% colnames(annotation_initial)){
                 peak_table_annotated <- annotation_initial %>%
                   dplyr::arrange(match(type, type_order)) %>%
                   dplyr::distinct(name, adduct, .keep_all = TRUE) %>%
                   dplyr::select(name, adduct, type)
               } else {
                 peak_table_annotated <- annotation_initial %>%
                   dplyr::distinct(name, adduct, .keep_all = TRUE) %>%
                   dplyr::mutate(type = 'level3') %>%
                   dplyr::select(name, adduct, type)
               }

               progress <- mapProgress(n = length(peak_table_annotated$name))

               list_peak_group <- purrr::map(seq_along(peak_table_annotated$name), function(i){

                 mapProgressPrint(progress)
                 temp_name <- peak_table_annotated$name[i]
                 # temp_rt <- peak_table_annotated$rt[i]
                 # temp_ccs <- peak_table_annotated
                 temp_adduct <- peak_table_annotated$adduct[i]
                 temp_seed <- peak_table_annotated$type[i]

                 result <- generatePeakGroup(peak_table = peak_table,
                                             base_peak_name = temp_name,
                                             base_peak_adduct = temp_adduct,
                                             base_peak_seed = temp_seed,
                                             tol_rt = tol_rt,
                                             cutoff_ssc = 0) # not apply ssc in peak group generation
                 return(result)

               })

               names(list_peak_group) <- paste(peak_table_annotated$name,
                                               peak_table_annotated$adduct,
                                               sep = '_')

               dir.create(file.path(path_dir, '03_annotation_credential', "00_intermediate_data"),
                          showWarnings = FALSE,
                          recursive = TRUE)

               save(list_peak_group,
                    file = file.path(path_dir,
                                     '03_annotation_credential',
                                     "00_intermediate_data",
                                     'list_peak_group.RData'),
                    version = 2)

             }

             # 03 semi-targeted feature clustering -----------------------------

             # message('Start to separate peak group...\n')
             cat('\n\n'); message('Start credential with semi-targeted feature clustering ...\n')

             if ('list_peak_group_annotation_concised.RData' %in%
                 list.files(file.path(path_dir,
                                      '03_annotation_credential',
                                      "00_intermediate_data"))) {
               cat('Skip semi-targeted feature clustering.\n')
               load(file.path(path_dir,
                              '03_annotation_credential',
                              "00_intermediate_data",
                              'list_peak_group_annotation_concised.RData'))

             } else {

               # semi-targeted feature clustering
               authenticatePeakGroupCredential(list_peak_group = list_peak_group,
                                               ms2_data = raw_msms,
                                               path_dir = path_dir,
                                               thread = thread,
                                               polarity = polarity,
                                               tol_mz = tol_mz,
                                               isotope_int_ratio_check = isotope_int_ratio_check,
                                               isotope_int_ratio_cutoff = isotope_int_ratio_cutoff,
                                               is_ms2_check = is_ms2_check,
                                               ms2_score_cutoff = ms2_score_cutoff,
                                               cutoff_ssc = cutoff_ssc,
                                               cutoff_ssc_int = cutoff_ssc_int,
                                               is_rule_limitation = is_rule_limitation,
                                               cutoff_topN = cutoff_topN,
                                               is_plot_pseudo_MS1 = is_plot_pseudo_MS1,
                                               type_order = type_order)


               load(file.path(path_dir,
                              '03_annotation_credential',
                              "00_intermediate_data",
                              'list_peak_group_annotation_concised.RData'))
             }

             # 04 formula prediction -------------------------------------------
             cat('\n\n'); message('Start credential with formula prediction...\n')
             if ('list_peak_group_formula.RData' %in%
                 list.files(file.path(path_dir,
                                      '03_annotation_credential',
                                      "00_intermediate_data"))) {
               cat('Skip formula prediction.\n')
               load(file.path(path_dir,
                              '03_annotation_credential',
                              "00_intermediate_data",
                              'list_peak_group_formula.RData'))

             } else {
               temp_path <- getwd()
               if (is_pred_formula_all) {
                 cat('Note: predict formula for all annotated peak groups\n')
                 temp_list_peak_group <- list_peak_group
               } else {
                 cat('Note: predict formula for concised peak groups\n')
                 temp_list_peak_group <- list_peak_group_annotation_concised
               }

               authenticateFormulaCredential(list_peak_group = temp_list_peak_group,
                                             ms2_data = raw_msms,
                                             polarity = polarity,
                                             path_dir = path_dir,
                                             thread = thread,
                                             dir_GenForm = dir_GenForm,
                                             ppm = ppm,
                                             acc = acc,
                                             platform = platform,
                                             elements = elements,
                                             num_formula_candidate = num_formula_candidate)

               setwd(temp_path)
             }

             message('Start refine annotations with credential...\n')
             load(file.path(path_dir,
                            '03_annotation_credential',
                            "00_intermediate_data",
                            'list_peak_group_formula.RData'))

             refineAnnotation(annotation_initial,
                              list_peak_group_annotation_concised,
                              list_peak_group_formula,
                              path_dir = path_dir)

             cat('\n'); cat('Annotation credential is done\n')


           })


################################################################################
# convertInputData4AnnotationCredential ----------------------------------------
#' @title convertAnnotationTable2InitialId
#' @description convert Annotation Table (MetAnalyzer) to initial ID table
#' @author Zhiwei Zhou
#' @param peak_table_file Default: 'idresults.csv'
#' @param ms2_file 'spectra.msp'
#' @param sample_info_file 'sample.info.csv'
#' @param lib library of IDs. 'zhulib', 'zhumetlib', 'kegg'. Default: 'zhumetlib'
#' @param path_dir Default: '.'
#' @param match_method 'forward', 'reverse'; Default: 'reverse'
#' @param polarity 'positive', 'negative'
#' @param tool 'MetAnalyzer1', 'MetAnalyzer2', 'MetDNA'; Default:'MetAnalyzer1'
#' @export
#' @examples




# convertAnnotationTable2InitialId(peak_table_file = 'data.csv',
#                                  sample_info_file =  'sample.info.csv',
#                                  path_dir = 'H:/00_projects/03_MetDNA2/00_data/20200922_metdna2_development/pos_201012_metdna2/',
#                                  lib = 'kegg',
#                                  polarity = 'positive',
#                                  tool = 'MetDNA2')



setGeneric(name = 'convertAnnotationTable2InitialId',
           def = function(
             peak_table_file = 'idresults.csv',
             ms2_file = 'spectra.msp',
             sample_info_file = 'sample.info.csv',
             # sample_name,
             path_dir = '.',
             lib = c('zhulib', 'zhumetlib', 'kegg'),
             match_method = c('forward', 'reverse'),
             polarity = c('positive', 'negative'),
             tool = c('MetAnalyzer1', 'MetAnalyzer2', 'MetDNA', 'MetDNA2'),
             # test_assi_confidence_initial_seed = FALSE,
             # dp_cutoff = 0.8,
             ...
           ){
             lib <- match.arg(lib)
             match_method <- match.arg(match_method)
             polarity <- match.arg(polarity)
             tool <- match.arg(tool)

             # if (missing(sample) & missing(sample_info_file)) {
             #   stop('Please input ')
             # }

             if (tool == 'MetDNA') {
               lib <- 'kegg'
             }

             switch (tool,
                     'MetAnalyzer1' = {
                       prepareMetAnalyzer(peak_table_file = peak_table_file,
                                          ms2_file = ms2_file,
                                          sample_info_file = sample_info_file,
                                          path_dir = path_dir,
                                          lib = lib,
                                          match_method = match_method,
                                          polarity = polarity,
                                          verison = '1',
                                          ...)
                     },
                     'MetAnalyzer2' = {
                       prepareMetAnalyzer(peak_table_file = peak_table_file,
                                          ms2_file = ms2_file,
                                          sample_info_file = sample_info_file,
                                          path_dir = path_dir,
                                          lib = lib,
                                          match_method = match_method,
                                          polarity = polarity,
                                          verison = '2',
                                          ...)
                     },
                     'MetDNA' = {
                       prepareMetDNA(peak_table_file = peak_table_file,
                                     sample_info_file = sample_info_file,
                                     path_dir = path_dir,
                                     polarity = polarity)
                     },
                     'MetDNA2' = {
                       prepareMetDNA2(peak_table_file = peak_table_file,
                                      sample_info_file = sample_info_file,
                                      path_dir = path_dir,
                                      polarity = polarity,
                                      # test_assi_confidence_initial_seed = test_assi_confidence_initial_seed,
                                      ...)
                     }
             )

             cat('\n');cat('Convert MSP to RData\n')
             convertMs2RData(ms2_file = 'ms2_data.msp',
                             path_dir = file.path(path_dir, '03_annotation_credential'))



           })





# convertMs2RData --------------------------------------------------------------

#' @title convertMs2RData
#' @author Zhiwei Zhou
#' @param ms2_file
#' @param path_dir '.'
#' @export
#' @examples
#' convertMs2RData(ms2_file = 'spectra.msp',
#'                 path_dir = 'I:/00_projects/03_MetDNA2/00_data/20200811_annotationCredential_workflow/demo_data')

# convertMs2RData(ms2_file = 'spectra.msp',
#                 path_dir = 'I:/00_projects/03_MetDNA2/00_data/20200811_annotationCredential_workflow/demo_data')
#
# raw_msms <- readMSP("/home/zhouzw/Data_processing/20210522_debug/neg_mode/03_annotation_credential/ms2_data.msp",
#                     mode = 'all',
#                     source = 'Other')

# ms2_file <- 'spectra.msp'
setGeneric(name = 'convertMs2RData',
           def = function(
             ms2_file,
             path_dir = '.'
           ){
             raw_msms <- readMSP(file.path(path_dir, 'ms2_data.msp'),
                                 mode = 'all',
                                 source = 'Other')
             temp <- sapply(raw_msms, function(x){
               x[[1]]$NAME
             })
             names(raw_msms) <- temp

             save(raw_msms,
                  file = file.path(path_dir, 'ms2_data.RData'),
                  version = 2)
           })

# prepareFuns ------------------------------------------------------------------
#    prepareMetAnalyzer ---------------------------------------------------------
#' @title prepareMetAnalyzer
#' @author Zhiwei Zhou
#' @description generate files (MetAnalyzer result) for annotation credential

# path_dir <- 'H:/00_projects/03_MetDNA2/00_data/20200916_annotation_credential_modification/02_MetAnalyzer1'
# prepareMetAnalyzer(path_dir = path_dir,
#                    lib = 'zhumetlib',
#                    match_method = 'reverse',
#                    polarity = 'positive',
#                    verison = '1')

# path_dir <- 'H:/00_projects/03_MetDNA2/00_data/20200916_annotation_credential_modification/01_MetAnalyzer2/200std_pos'
# peak_table_file <- 'result_MSMSmatch_zhulib.csv'
# ms2_file = 'spectra_ext0.msp'
# match_method <- 'forward'
# polarity <- 'positive'
# verison <- '2'

# prepareMetAnalyzer(peak_table_file = 'result_MSMSmatch_zhulib.csv',
#                    ms2_file = 'spectra_ext0.msp',
#                    path_dir = path_dir,
#                    lib = 'zhumetlib',
#                    match_method = 'reverse',
#                    polarity = 'positive',
#                    verison = '2')

setGeneric(name = 'prepareMetAnalyzer',
           def = function(
             peak_table_file = 'idresults.csv',
             ms2_file = 'spectra.msp',
             sample_info_file = 'sample.info.csv',
             # sample_name,
             path_dir = '.',
             lib = c('zhulib', 'zhumetlib'),
             match_method = c('forward', 'reverse'),
             polarity = c('positive', 'negative'),
             verison = c('1', '2'),
             ...
           ){

             lib <- match.arg(lib)
             match_method <- match.arg(match_method)
             polarity <- match.arg(polarity)
             verison <- match.arg(verison)

             options(readr.num_columns = 0)
             temp_peak_table <- readr::read_csv(file.path(path_dir, peak_table_file))

             if (!all(c(peak_table_file, ms2_file) %in% list.files(path_dir))) {
               stop(peak_table_file, ' and ', ms2_file,
                    ' is not found, please check the path\n')
             }


             # if (missing(sample_name)) {
             #
             # }

             if (!(sample_info_file %in% list.files(path_dir))) {
               stop(sample_info_file, 'not found, please check the path')
             }

             options(readr.num_columns = 0)
             temp_sample_info <- readr::read_csv(file.path(path_dir, sample_info_file))
             sample_name <- temp_sample_info$sample.name

             if (!all(sample_name %in% names(temp_peak_table))) {
               stop('Some samples not existed in peak table, please check!\n')
             }


             cat('Export peak_table for credential...\n')
             # generate peak table for annotation credential
             switch (verison,
                     '1' = {
                       peak_table <- temp_peak_table %>%
                         dplyr::select(name, mzmed, rtmed, sample_name) %>%
                         dplyr::rename(mz = mzmed,
                                       rt = rtmed)

                       readr::write_csv(peak_table,
                                        path = file.path(path_dir, 'peak_table.csv'))
                     },
                     '2' = {
                       peak_table <- temp_peak_table %>%
                         dplyr::select(name, mz, rt, sample_name)

                       readr::write_csv(peak_table,
                                        path = file.path(path_dir, 'peak_table.csv'))
                     }
             )



             cat('Export annotation_initial for credential...\n')
             # generate annotation initial table for annotation credential
             switch (lib,
                     'zhumetlib' = {
                       # lib_cpd <- cpd_zhulib
                       lib_cpd <- zhuMetlib$meta$compound
                     },
                     'zhulib' = {
                       # lib_cpd <- cpd_zhulib
                       lib_cpd <- zhuMetlib$meta$compound
                     }
             )

             switch (verison,
                     '1' = {
                       if (match_method == 'reverse') {

                         result_annotation_initial <- temp_peak_table %>%
                           dplyr::select(name, mzmed, rtmed,
                                         dplyr::starts_with('hits.reverse')) %>%
                           dplyr::rename(mz = mzmed,
                                         rt = rtmed,
                                         temp_id = dplyr::starts_with('hits.reverse')) %>%
                           dplyr::filter(!is.na(temp_id)) %>%
                           tidyr::separate_rows(temp_id, sep = ';') %>%
                           dplyr::mutate(id = stringr::str_extract(temp_id,
                                                                   pattern = 'labid\\{.+\\}'),
                                         adduct = stringr::str_extract(temp_id,
                                                                       pattern = 'adduct\\{.+\\}name')) %>%
                           dplyr::mutate(id = stringr::str_replace(id,
                                                                   pattern = 'labid\\{',
                                                                   replacement = ''),
                                         adduct = stringr::str_replace(adduct,
                                                                       pattern = 'adduct\\{',
                                                                       replacement = '')) %>%
                           dplyr::mutate(id = stringr::str_replace(id,
                                                                   pattern = '\\}',
                                                                   replacement = ''),
                                         adduct = stringr::str_replace(adduct,
                                                                       pattern = '\\}name',
                                                                       replacement = '')) %>%
                           dplyr::select(name, mz, rt, id, adduct)
                       } else {
                         result_annotation_initial <- temp_peak_table %>%
                           dplyr::select(name, mzmed, rtmed,
                                         dplyr::starts_with('hits.forward')) %>%
                           dplyr::rename(mz = mzmed,
                                         rt = rtmed,
                                         temp_id = dplyr::starts_with('hits.forward')) %>%
                           dplyr::filter(!is.na(temp_id)) %>%
                           tidyr::separate_rows(temp_id, sep = ';') %>%
                           dplyr::mutate(id = stringr::str_extract(temp_id,
                                                                   pattern = 'labid\\{.+\\}'),
                                         adduct = stringr::str_extract(temp_id,
                                                                       pattern = 'adduct\\{.+\\}name')) %>%
                           dplyr::mutate(id = stringr::str_replace(id,
                                                                   pattern = 'labid\\{',
                                                                   replacement = ''),
                                         adduct = stringr::str_replace(adduct,
                                                                       pattern = 'adduct\\{',
                                                                       replacement = '')) %>%
                           dplyr::mutate(id = stringr::str_replace(id,
                                                                   pattern = '\\}',
                                                                   replacement = ''),
                                         adduct = stringr::str_replace(adduct,
                                                                       pattern = '\\}name',
                                                                       replacement = '')) %>%
                           dplyr::select(name, mz, rt, id, adduct)
                       }

                       # modify adduct format
                       if (polarity == 'positive') {
                         result_annotation_initial <- result_annotation_initial %>%
                           dplyr::mutate(adduct = dplyr::case_when(
                             adduct == '(M+H)+' ~ '[M+H]+',
                             adduct == '(M+H-H2O)+' ~ '[M-H2O+H]+',
                             adduct == "(M+H-2H2O)+" ~ "[M-2H2O+H]+",
                             adduct == '(M+NH4)+' ~ '[M+NH4]+',
                             adduct == '(M+Na)+' ~ '[M+Na]+',
                             adduct == "(M-H+2Na)+" ~ "[M-H+2Na]+",
                             adduct == "(M-2H+3Na)+" ~ "[M-2H+3Na]+",
                             adduct == "(M+K)+" ~ "[M+K]+",
                             adduct == "(M-H+2K)+" ~ "[M-H+2K]+",
                             adduct == "(M+CH3CN+H)+" ~ "[M+CH3CN+H]+",
                             adduct == "(M+CH3CN+Na)+" ~ "[M+CH3CN+Na]+",
                             adduct == "(2M+H)+" ~ "[2M+H]+",
                             adduct == "(2M+NH4)+" ~ "[2M+NH4]+",
                             adduct == "(2M+Na)+" ~ "[2M+Na]+",
                             adduct == "(2M+K)+" ~ "[2M+K]+",
                             adduct == "(M+CH3COO+2H)+" ~ "[M+CH3COO+2H]+",
                             adduct == "(M+HCOO+2H)+" ~ "[M+HCOO+2H]+",

                             # other possible adducts
                             adduct == 'M+' ~ 'M+'
                           ))
                       } else {
                         result_annotation_initial <- result_annotation_initial %>%
                           dplyr::mutate(adduct = dplyr::case_when(
                             adduct == "M-" ~ "M-",
                             adduct == '(M-H)-' ~ '[M-H]-',
                             adduct == '(M-H2O-H)-' ~ '[M-H2O-H]-',
                             adduct == '(M+Na-2H)-' ~ '[M+Na-2H]-',
                             adduct == '(M+K-2H)-' ~ '[M+K-2H]-',
                             adduct == '(M+NH4-2H)-' ~ '[M+NH4-2H]-',
                             adduct == '(2M-H)-' ~ '[2M-H]-',
                             adduct == '(M+CH3COO)-' ~ '[M+CH3COO]-',
                             adduct == '(M+F)-' ~ '[M+F]-',

                             # other possible adducts
                             adduct == '(M+HCOO)-' ~ '[M+HCOO]-'

                           ))
                       }

                       # add formula and cpd_name
                       temp_formula <- match(result_annotation_initial$id, lib_cpd$id) %>% lib_cpd$formula[.]
                       temp_cpd_name <- match(result_annotation_initial$id, lib_cpd$id) %>% lib_cpd$name[.]

                       result_annotation_initial <- result_annotation_initial %>%
                         dplyr::mutate(formula = temp_formula,
                                       cpd_name = temp_cpd_name) %>%
                         dplyr::select(name:adduct, formula, cpd_name, dplyr::everything())

                       readr::write_csv(result_annotation_initial,
                                        path = file.path(path_dir, 'annotation_initial.csv'))


                     },
                     '2' = {

                       if (match_method == 'reverse') {

                         result_annotation_initial <- temp_peak_table %>%
                           dplyr::select(name, mz, rt,
                                         dplyr::starts_with('hits.reverse')) %>%
                           dplyr::rename(temp_id = dplyr::starts_with('hits.reverse')) %>%
                           dplyr::filter(!is.na(temp_id)) %>%
                           tidyr::separate_rows(temp_id, sep = ';') %>%
                           dplyr::mutate(id = stringr::str_extract(temp_id,
                                                                   pattern = 'labid\\{.+\\}'),
                                         adduct = stringr::str_extract(temp_id,
                                                                       pattern = 'adduct\\{.+\\}name')) %>%
                           dplyr::mutate(id = stringr::str_replace(id,
                                                                   pattern = 'labid\\{',
                                                                   replacement = ''),
                                         adduct = stringr::str_replace(adduct,
                                                                       pattern = 'adduct\\{',
                                                                       replacement = '')) %>%
                           dplyr::mutate(id = stringr::str_replace(id,
                                                                   pattern = '\\}',
                                                                   replacement = ''),
                                         adduct = stringr::str_replace(adduct,
                                                                       pattern = '\\}name',
                                                                       replacement = '')) %>%
                           dplyr::select(name, mz, rt, id, adduct)
                       } else {
                         result_annotation_initial <- temp_peak_table %>%
                           dplyr::select(name, mz, rt,
                                         dplyr::starts_with('hits.forward')) %>%
                           dplyr::rename(temp_id = dplyr::starts_with('hits.forward')) %>%
                           dplyr::filter(!is.na(temp_id)) %>%
                           tidyr::separate_rows(temp_id, sep = ';') %>%
                           dplyr::mutate(id = stringr::str_extract(temp_id,
                                                                   pattern = 'labid\\{.+\\}'),
                                         adduct = stringr::str_extract(temp_id,
                                                                       pattern = 'adduct\\{.+\\}name')) %>%
                           dplyr::mutate(id = stringr::str_replace(id,
                                                                   pattern = 'labid\\{',
                                                                   replacement = ''),
                                         adduct = stringr::str_replace(adduct,
                                                                       pattern = 'adduct\\{',
                                                                       replacement = '')) %>%
                           dplyr::mutate(id = stringr::str_replace(id,
                                                                   pattern = '\\}',
                                                                   replacement = ''),
                                         adduct = stringr::str_replace(adduct,
                                                                       pattern = '\\}name',
                                                                       replacement = '')) %>%
                           dplyr::select(name, mz, rt, id, adduct)
                       }

                       # modify adduct format
                       if (polarity == 'positive') {
                         result_annotation_initial <- result_annotation_initial %>%
                           dplyr::mutate(adduct = dplyr::case_when(
                             adduct == '[M+H]+' ~ '[M+H]+',
                             adduct == '[M+H-H2O]+' ~ '[M-H2O+H]+',
                             adduct == "[M+H-2H2O]+" ~ "[M-2H2O+H]+",
                             adduct == '[M+NH4]+' ~ '[M+NH4]+',
                             adduct == '[M+Na]+' ~ '[M+Na]+',
                             adduct == "[M+2Na-H]+" ~ "[M-H+2Na]+",
                             adduct == "[M+3Na-2H]+" ~ "[M-2H+3Na]+",
                             adduct == "[M+K]+" ~ "[M+K]+",
                             adduct == "[M+2K-H]+" ~ "[M-H+2K]+",
                             adduct == "[M+CH3CN+H]+" ~ "[M+CH3CN+H]+",
                             adduct == "[M+CH3CN+Na]+" ~ "[M+CH3CN+Na]+",
                             adduct == "[M+CH3COO+2H]+" ~ "[M+CH3COO+2H]+",

                             # other possible adducts
                             adduct == "[2M+H]+" ~ "[2M+H]+",
                             adduct == "[2M+NH4]+" ~ "[2M+NH4]+",
                             adduct == "[2M+Na]+" ~ "[2M+Na]+",
                             adduct == "[2M+K]+" ~ "[2M+K]+",
                             adduct == "[M+HCOO+2H]+" ~ "[M+HCOO+2H]+",
                             adduct == 'M+' ~ 'M+'
                           )) %>%
                           dplyr::filter(!is.na(adduct))
                       } else {
                         result_annotation_initial <- result_annotation_initial %>%
                           dplyr::mutate(adduct = dplyr::case_when(
                             adduct == '[M-H]-' ~ '[M-H]-',
                             adduct == '[M+Cl]-' ~ '[M+Cl]-',
                             adduct == '[M-2H+Na]-' ~ '[M+Na-2H]-',
                             adduct == '[M-2H+K]-' ~ '[M+K-2H]-',
                             adduct == '[M-H+CH3CN]-' ~ '[M+CH3CN-H]-',
                             adduct == '[M+Cl+NH3]-' ~ '[M+NH3+Cl]-',
                             adduct == '[M-H+NH3]-' ~ '[M+NH4-2H]-',
                             adduct == '[M-H-H2O]-' ~ '[M-H2O-H]-',
                             adduct == '[M+CH3COO]-' ~ '[M+CH3COO]-',
                             adduct == '[2M-H]-' ~ '[2M-H]-',


                             # other possible adducts
                             adduct == "M-" ~ "M-",
                             adduct == '[M+Na-2H]-' ~ '[M+Na-2H]-',
                             adduct == '[M-H-H2O]-' ~ '[M-H2O-H]-',
                             adduct == '[M-H2O-H]-' ~ '[M-H2O-H]-',
                             adduct == '[M+K-2H]-' ~ '[M+K-2H]-',
                             adduct == '[M+F]-' ~ '[M+F]-'
                           )) %>%
                           dplyr::filter(!is.na(adduct))
                       }

                       # add formula and cpd_name
                       temp_formula <- match(result_annotation_initial$id, lib_cpd$id) %>% lib_cpd$formula[.]
                       temp_cpd_name <- match(result_annotation_initial$id, lib_cpd$id) %>% lib_cpd$name[.]

                       result_annotation_initial <- result_annotation_initial %>%
                         dplyr::mutate(formula = temp_formula,
                                       cpd_name = temp_cpd_name) %>%
                         dplyr::select(name:adduct, formula, cpd_name, dplyr::everything())

                       readr::write_csv(result_annotation_initial,
                                        path = file.path(path_dir, 'annotation_initial.csv'))

                     }
             )


             cat('Export MS/MS spectra for credential...\n')

             switch (verison,
                     '1' = {
                       raw_msms <- readMSP(file.path(path_dir, ms2_file), mode = 'all')

                       progress <- mapProgress(n = length(raw_msms))
                       purrr::walk(seq_along(raw_msms), function(i){
                         mapProgressPrint(progress)
                         temp_name <- raw_msms[[i]]$info$NAME
                         temp_mz <- raw_msms[[i]]$info$PRECURSORMZ
                         temp_spec <- raw_msms[[i]]$spec

                         generateMSP(file_name = file.path(path_dir, 'ms2_data.msp'),
                                     cmp_name = temp_name,
                                     precusormz = temp_mz,
                                     spec = temp_spec)

                       })
                     },
                     '2' = {
                       raw_msms <- readMSP(file.path(path_dir, ms2_file), mode = 'all')

                       progress <- mapProgress(n = length(raw_msms))
                       purrr::walk(seq_along(raw_msms), function(i){
                         mapProgressPrint(progress)
                         temp_name <- raw_msms[[i]]$info$NAME
                         temp_mz <- raw_msms[[i]]$info$PRECURSORMZ
                         temp_rt <- raw_msms[[i]]$info$RETENTIONTIME
                         temp_spec <- raw_msms[[i]]$spec

                         generateMSP(file_name = file.path(path_dir, 'ms2_data.msp'),
                                     cmp_name = temp_name,
                                     precusormz = temp_mz,
                                     spec = temp_spec)

                       })

                     }
             )



           })

#    prepareMetDNA -------------------------------------------------------------
#' @title prepareMetDNA
#' @author Zhiwei Zhou
#' @description generate files for annotation credential

# path_dir <- 'H:/00_projects/03_MetDNA2/00_data/20200922_metdna2_development/'
# prepareMetDNA(path_dir = path_dir, polarity = 'positive')

setGeneric(name = 'prepareMetDNA',
           def = function(peak_table_file = 'data.csv',
                          sample_info_file = 'sample.info.csv',
                          path_dir = '.',
                          polarity = c('positive', 'negative')){
             polarity <- match.arg(polarity)

             if (!('MRN_annotation_result' %in% list.files(path_dir))) {
               stop('MRN_annotation_result directory is not found, please check the path\n')
             }

             if (!all(c(peak_table_file, sample_info_file) %in% list.files(path_dir))) {
               stop(peak_table_file, ' and ', sample_info_file,
                    ' is not found, please check the path\n')
             }

             load(file.path(path_dir, 'MRN_annotation_result/intermediate_data/tags.result'))
             load(file.path(path_dir, 'MRN_annotation_result/intermediate_data/tags2.after.redundancy.remove'))

             cat('Export peak_table for credential...\n')
             options(readr.num_columns = 0)
             temp_peak_table <- readr::read_csv(file.path(path_dir, peak_table_file))
             temp_sample_info <- readr::read_csv(file.path(path_dir, sample_info_file))
             peak_table <- temp_peak_table %>%
               dplyr::select('name', 'mz', 'rt', temp_sample_info$sample.name)
             readr::write_csv(peak_table, path = file.path(path_dir, 'peak_table.csv'))

             # dir.create(file.path(path_dir, 'cre'))
             cat('Export annotation_initial for credential...\n')

             tags_result <- tags.result %>%
               readr::type_convert() %>%
               tibble::as_tibble()

             switch(polarity,
                    'positive' = {
                      annotation_initial <- tags_result %>%
                        dplyr::filter(type %in% c("metAnnotation", "adductAnnotation", "seed")) %>%
                        dplyr::select(name, mz, rt, to, adduct, Formula) %>%
                        dplyr::rename(formula = Formula,
                                      id = to) %>%
                        # pull(adduct) %>% unique()
                        dplyr::mutate(adduct = dplyr::case_when(
                          adduct == 'M+H' ~ '[M+H]+',
                          adduct == 'M+H-H2O' ~ '[M-H2O+H]+',
                          adduct == "M+H-2H2O" ~ "[M-2H2O+H]+",
                          adduct == 'M+NH4' ~ '[M+NH4]+',
                          adduct == 'M+Na' ~ '[M+Na]+',
                          adduct == "M+2Na-H" ~ "[M-H+2Na]+",
                          adduct == "M+3Na-2H" ~ "[M-2H+3Na]+",
                          adduct == "M+K" ~ "[M+K]+",
                          adduct == 'M+2K-H' ~ '[M-H+2K]+',
                          adduct == "M-2H+3K" ~ "[M-2H+3K]+",
                          adduct == 'M+CH3CN+H' ~ '[M+CH3CN+H]+',
                          adduct == 'M+CH3CN+H' ~ '[M+CH3CN+H]+',
                          adduct == "M+CH3CN+Na" ~ "[M+CH3CN+Na]+",
                          adduct == "M+CH3COO+2H" ~ "[M+CH3COO+2H]+",
                          adduct == 'M+Na+HCOOH' ~ '[M+HCOO+H+Na]+',
                          adduct == 'M+K+HCOOH' ~ '[M+HCOO+H+K]+',
                          adduct == 'M+H+HCOOH' ~ '[M+HCOO+2H]+',
                          adduct == '2M+H' ~ '[2M+H]+',
                          adduct == "2M+K" ~ '[2M+K]+',
                          # other possible
                          adduct == 'M+' ~ '[M]+',
                          adduct == '2M+Na' ~ '[2M+Na]+',
                          adduct == '2M+NH4' ~ '[2M+NH4]+',
                          adduct == "M-H+2Na" ~ '[M-H+2Na]+',
                          adduct == 'M-2H+3Na' ~ '[M-2H+3Na]+',
                          adduct == 'M-H+2K' ~ '[M-H+2K]+'
                        )) %>%
                        dplyr::filter(!is.na(adduct))

                    },
                    'negative' = {
                      annotation_initial <- tags_result %>%
                        dplyr::filter(type %in% c("metAnnotation", "adductAnnotation", "seed")) %>%
                        dplyr::select(name, mz, rt, to, adduct, Formula) %>%
                        dplyr::rename(formula = Formula,
                                      id = to) %>%
                        # pull(adduct) %>% unique()
                        dplyr::mutate(adduct = dplyr::case_when(
                          adduct == 'M-H' ~ '[M-H]-',
                          adduct == 'M+Cl' ~ '[M+Cl]-',
                          adduct == 'M-2H+Na' ~ '[M+Na-2H]-',
                          adduct == 'M-2H+K' ~ '[M+K-2H]-',
                          adduct == 'M-H+CH3CN' ~ '[M+CH3CN-H]-',
                          adduct == 'M+Cl+NH3' ~ '[M+NH3+Cl]-',
                          adduct == 'M-H+NH3' ~ '[M+NH4-2H]-',
                          adduct == 'M+F' ~ '[M+F]-',
                          adduct == '2M-H' ~ '[2M-H]-',
                          adduct == 'M+CH3COO' ~ '[M+CH3COO]-',

                          # other possible adducts
                          adduct == "M-" ~ "[M]-",
                          adduct == 'M+Na-2H' ~ '[M+Na-2H]-',
                          adduct == 'M-H-H2O' ~ '[M-H2O-H]-',
                          adduct == 'M-H2O-H' ~ '[M-H2O-H]-',
                          adduct == 'M+K-2H' ~ '[M+K-2H]-',

                        )) %>%
                        dplyr::filter(!is.na(adduct))
                    }
             )

             idx <- match(annotation_initial$id, name_formula_kegg$id)
             annotation_initial <- annotation_initial %>%
               dplyr::mutate(cpd_name = name_formula_kegg$name[idx])

             readr::write_csv(annotation_initial, path = file.path(path_dir, 'annotation_initial.csv'))

             cat('Export MS/MS spectra for credential...\n')
             progress <- mapProgress(n = length(tags2.after.redundancy.remove))
             purrr::walk(seq_along(tags2.after.redundancy.remove), function(i){
               mapProgressPrint(progress = progress)
               temp_data <- tags2.after.redundancy.remove[[i]]

               if (nrow(temp_data@ms2) > 0) {
                 temp_name <- temp_data@name
                 temp_mz <- temp_data@mz
                 temp_rt <- temp_data@rt

                 temp_spec <- temp_data@ms2

                 generateMSP(file_name = file.path(path_dir, 'ms2_data.msp'),
                             cmp_name = temp_name,
                             precusormz = temp_mz,
                             rt = temp_rt,
                             spec = temp_spec)
               }

             })


           })



#    prepareMetDNA2 -------------------------------------------------------------
#' @title prepareMetDNA2
#' @author Zhiwei Zhou
#' @description generate files for annotation credential

# path_dir <- 'H:/00_projects/03_MetDNA2/00_data/20200922_metdna2_development/pos_201012_metdna2/'
# prepareMetDNA2(path_dir = path_dir, polarity = 'positive')

# path_dir <- '/home/zhouzw/Data_processing/20201130_modify_metdna2_aging_fly/'
# path_dir <- '/home/zhouzw/Data_processing/20210315_tims_data_processing_metdna'
# direction <- 'reverse'
# polarity <- 'positive'

setGeneric(name = 'prepareMetDNA2',
           def = function(peak_table_file = 'data.csv',
                          sample_info_file = 'sample.info.csv',
                          direction = c('reverse', 'forward'),
                          path_dir = '.',
                          polarity = c('positive', 'negative'),
                          # test_assi_confidence_initial_seed = FALSE,
                          # instrument = c("SciexTripleTOF", "AgilentQTOF", "BrukerQTOF", "ThermoOrbitrap", "IMMS"),
                          ...){
             polarity <- match.arg(polarity)
             direction <- match.arg(direction)

             if (!('02_result_MRN_annotation' %in% list.files(path_dir))) {
               stop('MRN_annotation_result directory is not found, please check the path\n')
             }

             if (!all(c(peak_table_file, sample_info_file) %in% list.files(path_dir))) {
               stop(peak_table_file, ' and ', sample_info_file,
                    ' is not found, please check the path\n')
             }

             dir.create(file.path(path_dir, '03_annotation_credential'),
                        showWarnings = FALSE,
                        recursive = TRUE)
             path_output <- file.path(path_dir, '03_annotation_credential')

             # load(file.path(path_dir, 'MS2_match_result/intermediate_data/result_annotation'))
             # load(file.path(path_dir, 'MRN_annotation_result/intermediate_data/id_result_redun_rm'))
             # load(file.path(path_dir, 'MRN_annotation_result/intermediate_data/tags2_after_redundancy_remove'))

             load(file.path(path_dir, '03_annotation_credential/00_intermediate_data/id_merge_before_credential.RData'))
             load(file.path(path_dir, '02_result_MRN_annotation/00_intermediate_data/tags2_after_redundancy_remove'))

             cat('Export peak_table for credential...\n')
             # if (instrument != 'IMMS') {
             #   options(readr.num_columns = 0)
             #   temp_peak_table <- readr::read_csv(file.path(path_dir, peak_table_file))
             #   temp_sample_info <- readr::read_csv(file.path(path_dir, sample_info_file))
             #   peak_table <- temp_peak_table %>%
             #     dplyr::select('name', 'mz', 'rt', temp_sample_info$sample.name)
             #   readr::write_csv(peak_table, path = file.path(path_dir, 'peak_table.csv'))
             # } else {
             #   temp_sample_info <- readr::read_csv(file.path(path_dir, sample_info_file))
             #   load(file.path(path_dir, '01_result_initial_seed_annotation/00_intermediate_data/ms1_data'))
             #   peak_table <- ms1_data$info %>%
             #     dplyr::bind_cols(ms1_data$subject) %>%
             #     dplyr::select('name', 'mz', 'rt', temp_sample_info$sample.name)
             #   readr::write_csv(peak_table, path = file.path(path_dir, 'peak_table.csv'))
             # }

             temp_sample_info <- readr::read_csv(file.path(path_dir, sample_info_file))
             load(file.path(path_dir, '01_result_initial_seed_annotation/00_intermediate_data/ms1_data'))
             peak_table <- ms1_data$info %>%
               dplyr::bind_cols(ms1_data$subject) %>%
               dplyr::select('name', 'mz', 'rt', 'ccs', temp_sample_info$sample.name)
             readr::write_csv(peak_table, path = file.path(path_output, 'peak_table.csv'))


             # dir.create(file.path(path_dir, 'cre'))
             cat('Export annotation_initial for credential...\n')

             # initial seed annotation
             # result_initial_seed <- lapply(seq_along(result_annotation), function(i){
             #   x <- result_annotation[[i]]
             #
             #   if (nrow(x@annotation_result) > 0) {
             #     result <- x@annotation_result %>%
             #       dplyr::arrange(desc(msms_score_forward),
             #                      desc(msms_score_reverse),
             #                      desc(ccs_score),
             #                      desc(rt_score)) %>%
             #       dplyr::mutate(feature = x@peak_info$name,
             #                     mz = x@peak_info$mz,
             #                     rt = x@peak_info$rt) %>%
             #       dplyr::select(feature, mz, rt, dplyr::everything())
             #     return(result)
             #   } else {
             #     return(NULL)
             #   }
             # }) %>% dplyr::bind_rows()



             # switch(direction,
             #        'forward' = {
             #          annotation_initial_seed <- result_initial_seed %>%
             #            dplyr::filter(msms_score_forward >= dp_cutoff) %>%
             #            dplyr::select(feature:rt, id, adduct, formula, name) %>%
             #            dplyr::rename(name = feature,
             #                          cpd_name = name) %>%
             #            dplyr::mutate(type = 'seed')
             #        },
             #        'reverse' = {
             #          annotation_initial_seed <- result_initial_seed %>%
             #            dplyr::filter(msms_score_reverse >= dp_cutoff) %>%
             #            dplyr::select(feature:rt, id, adduct, formula, name) %>%
             #            dplyr::rename(name = feature,
             #                          cpd_name = name) %>%
             #            dplyr::mutate(type = 'seed')
             #        })
             #
             #
             #
             #
             # # recursive annotation
             # tags_result <- id_result_redun_rm %>%
             #   readr::type_convert() %>%
             #   tibble::as_tibble()
             #
             # switch(polarity,
             #        'positive' = {
             #          annotation_initial_recursive <- tags_result %>%
             #            dplyr::filter(type %in% c("metAnnotation", "adductAnnotation", "seed")) %>%
             #            dplyr::select(name, mz, rt, to, adduct, formula, type) %>%
             #            dplyr::rename(id = to) %>%
             #            # pull(adduct) %>% unique()
             #            dplyr::mutate(adduct = dplyr::case_when(
             #              adduct == 'M+' ~ 'M+',
             #              adduct == 'M+H' ~ '[M+H]+',
             #              adduct == 'M+NH4' ~ '[M+NH4]+',
             #              adduct == 'M+Na' ~ '[M+Na]+',
             #              adduct == "M-H+2Na" ~ "[M-H+2Na]+",
             #              adduct == 'M-2H+3Na' ~ '[M-2H+3Na]+',
             #              adduct == "M+K" ~ "[M+K]+",
             #              adduct == 'M-H+2K' ~ '[M-H+2K]+',
             #              adduct == "M-2H+3K" ~ "[M-2H+3K]+",
             #              adduct == 'M+CH3CN+H' ~ '[M+CH3CN+H]+',
             #              adduct == "M+CH3CN+Na" ~ "[M+CH3CN+Na]+",
             #              adduct == '2M+H' ~ '[2M+H]+',
             #              adduct == '2M+NH4' ~ '[2M+NH4]+',
             #              adduct == "2M+K" ~ '[2M+K]+',
             #              adduct == '2M+Na' ~ '[2M+Na]+',
             #              adduct == "M+CH3COO+2H" ~ "[M+CH3COO+2H]+",
             #              adduct == 'M+HCOO+2H' ~ '[M+HCOO+2H]+',
             #              adduct == 'M+HCOO+H+K' ~ '[M+HCOO+H+K]+',
             #              adduct == 'M+HCOO+H+Na' ~ '[M+HCOO+H+Na]+',
             #              adduct == 'M-H2O+H' ~ '[M-H2O+H]+',
             #              adduct == "M-2H2O+H" ~ "[M-2H2O+H]+"
             #
             #            )) %>%
             #            dplyr::filter(!is.na(adduct))
             #
             #        },
             #        'negative' = {
             #          annotation_initial_recursive <- tags_result %>%
             #            dplyr::filter(type %in% c("metAnnotation", "adductAnnotation", "seed")) %>%
             #            dplyr::select(name, mz, rt, to, adduct, formula, type) %>%
             #            dplyr::rename(id = to) %>%
             #            # pull(adduct) %>% unique()
             #            dplyr::mutate(adduct = dplyr::case_when(
             #              adduct == "M-" ~ "M-",
             #              adduct == 'M-H' ~ '[M-H]-',
             #              adduct == 'M+Na-2H' ~ '[M+Na-2H]-',
             #              adduct == 'M+K-2H' ~ '[M+K-2H]-',
             #              adduct == 'M+NH4-2H' ~ '[M+NH4-2H]-',
             #              adduct == '2M-H' ~ '[2M-H]-',
             #              adduct == 'M+CH3COO' ~ '[M+CH3COO]-',
             #              adduct == 'M+F' ~ '[M+F]-',
             #              adduct == 'M+HCOO' ~ '[M+HCOO]-',
             #              adduct == "2M+Na-2H" ~ "[2M+Na-2H]-",
             #              adduct == "3M+Na-2H" ~ "[3M+Na-2H]-",
             #              adduct == "M+2Na-3H" ~ "[M+2Na-3H]-",
             #              adduct == "3M-H" ~ "[3M-H]-",
             #              adduct == 'M+Cl' ~ '[M+Cl]-',
             #              adduct == 'M+CH3CN-H' ~ '[M+CH3CN-H]-',
             #              adduct == 'M+NH3+Cl' ~ '[M+NH3+Cl]-',
             #              adduct == "M-2H" ~ "[M-2H]2-",
             #              adduct == 'M-H2O-H' ~ '[M-H2O-H]-'
             #
             #            )) %>%
             #            dplyr::filter(!is.na(adduct))
             #        }
             # )
             #
             # idx <- match(annotation_initial_recursive$id, cpd_emrn$id)
             # annotation_initial_recursive <- annotation_initial_recursive %>%
             #   dplyr::mutate(cpd_name = cpd_emrn$name[idx])
             #
             # annotation_initial <- annotation_initial_seed %>%
             #   dplyr::bind_rows(annotation_initial_recursive)

             # For adductAnnotation, and isotopeAnnotation that directly annotated from initial seed,
             #   assign confidence 'level2' in annotation credential because they have higher accuracy than others
             # if (test_assi_confidence_initial_seed) {
             #   idx <- which(id_merge$round == 2 & (id_merge$recursive_type %in% c('adductAnnotation')))
             #   id_merge$confidence_level[idx] <- 'level2'
             # }

             annotation_initial <- id_merge %>%
               dplyr::filter(!is.na(id)) %>%
               dplyr::filter((confidence_level %in% c('level1', 'level2')) |
                               (recursive_type %in% c("metAnnotation", "adductAnnotation", "seed"))) %>%
               dplyr::select(peak_name, mz, rt, ccs, id, adduct, formula, confidence_level, name) %>%
               dplyr::rename(name = peak_name,
                             type = confidence_level,
                             cpd_name = name)

             readr::write_csv(annotation_initial, path = file.path(path_output, 'annotation_initial.csv'))

             cat('Export MS/MS spectra for credential...\n')

             # remove existed MSP files to avoid replicated extend
             if ('ms2_data.msp' %in% list.files(path_output)) {
               file.remove(file.path(path_output, 'ms2_data.msp'))
             }

             progress <- mapProgress(n = length(tags2_after_redundancy_remove))
             purrr::walk(seq_along(tags2_after_redundancy_remove), function(i){
               mapProgressPrint(progress = progress)
               temp_data <- tags2_after_redundancy_remove[[i]]

               if (nrow(temp_data@ms2) > 0) {
                 temp_name <- temp_data@name
                 temp_mz <- temp_data@mz
                 temp_rt <- temp_data@rt

                 temp_spec <- temp_data@ms2

                 generateMSP(file_name = file.path(path_output, 'ms2_data.msp'),
                             cmp_name = temp_name,
                             precusormz = temp_mz,
                             rt = temp_rt,
                             polarity = polarity,
                             spec = temp_spec)
               }

             })


           })





################################################################################
# refineAnnotation -------------------------------------------------------------

#' @title refineAnnotation
#' @author Zhiwei Zhou
#' @description refine annotation credential with semi-targeted clustering and formula prediction
#' @param annotation_initial
#' @param list_peak_group_annotation_concised
#' @param list_peak_group_formula
#' @param path_dir Default: '.'
#' @export
#' @examples

# annotation_initial <- readr::read_csv(file.path(path_dir, 'annotation'))
# load(file.path(path_dir,
#                '03_annotation_credential',
#                '00_intermediate_data',
#                'annotation_initial.RData'))

setGeneric(name = 'refineAnnotation',
           function(annotation_initial,
                    list_peak_group_annotation_concised,
                    list_peak_group_formula,
                    path_dir = '.'){

             # message('Refine annotations with semi-targeted feature clustering and formula prediction\n')
             # browser()

             annotation_credential <- refineAnnotationFeatureClustering(annotation_initial = annotation_initial,
                                                                        list_peak_group_annotation_concised = list_peak_group_annotation_concised)

             annotation_credentia2 <- refineAnnotationFormula(annotation_initial = annotation_credential,
                                                              list_peak_group_formula = list_peak_group_formula,
                                                              path_dir = path_dir)

             annotation_credential_long <- annotation_credentia2 %>%
               dplyr::select(name:cpd_name,
                             peak_group_scale, score_formula, rank_formula,
                             dplyr::everything())

             dir.create(file.path(path_dir, '03_annotation_credential'),
                        showWarnings = FALSE, recursive = TRUE)

             readr::write_csv(annotation_credential_long,
                              path = file.path(path_dir, '03_annotation_credential',
                                               'annontation_credential_long.csv'))


             annotation_credential_wide <- annotation_credential_long %>%
               dplyr::group_by(name) %>%
               dplyr::summarise(id = paste(id, collapse = ';'),
                                adduct = paste(adduct, collapse = ';'),
                                formula = paste(formula, collapse = ';'),
                                cpd_name = paste(cpd_name, collapse = ';'),
                                peak_group_scale = paste(peak_group_scale, collapse = ';'),
                                score_formula = paste(score_formula, collapse = ';'),
                                rank_formula = paste(rank_formula, collapse = ';'))

             readr::write_csv(annotation_credential_wide,
                              path = file.path(path_dir, '03_annotation_credential',
                                               'annontation_credential_wide.csv'))

             # export type and reason for all filtered annotations
             load(file.path(path_dir, '03_annotation_credential', "00_intermediate_data", 'record_rule_filter_peak_group.RData'))
             load(file.path(path_dir, '03_annotation_credential', "00_intermediate_data", 'record_conflict_peak_group.RData'))
             load(file.path(path_dir, '03_annotation_credential', "00_intermediate_data", 'record_overlap_peak_group.RData'))
             load(file.path(path_dir, '03_annotation_credential', "00_intermediate_data", 'record_formula_filter.RData'))

             annotation_filter <- annotation_initial %>%
               dplyr::mutate(peak_group = paste(name, adduct, sep = '_'),
                             peak_group_formula = paste(name, adduct, sep = '_') %>% paste(formula, sep = '@'),
                             temp = paste(name, id, sep = '_'),
                             type = NA) %>%
               dplyr::filter(!(temp %in% paste(annotation_credential_long$name, annotation_credential_long$id, sep = '_')))

             idx <- which(!(match(annotation_filter$peak_group, record_rule_filter_peak_group) %>% is.na()))
             if (length(idx) > 0) {
               annotation_filter$type[idx] <- 'type0'
             }

             idx <- which(!(match(annotation_filter$peak_group, record_conflict_peak_group) %>% is.na()))
             if (length(idx) > 0) {
               annotation_filter$type[idx] <- 'type1'
             }

             record_overlap_peak_group <- record_overlap_peak_group %>% unlist() %>% unique()
             temp <- record_overlap_peak_group %>%
               stringr::str_split('@') %>%
               do.call(rbind, .) %>%
               tibble::as_tibble() %>%
               dplyr::rename(peak_group = V1,
                             source = V2)
             annotation_filter <- annotation_filter %>% dplyr::left_join(temp, by = c('peak_group'))
             annotation_filter$type[!is.na(annotation_filter$source)] <- 'type2'


             idx <- which(!(match(annotation_filter$peak_group_formula, record_formula_filter) %>% is.na()))
             if (length(idx) > 0) {
               annotation_filter$type[idx] <- 'type3'
             }

             annotation_filter <- annotation_filter %>% dplyr::select(-c('peak_group', 'peak_group_formula', 'temp'))


             readr::write_csv(annotation_filter,
                              path = file.path(path_dir, '03_annotation_credential',
                                               'annontation_credential_filter.csv'))

           })









#' @title refineAnnotationFeatureClustering
#' @author Zhiwei Zhou
#' @param annotation_initial
#' @param list_peak_group_annotation_concised
#' @export

# load('./inst/tempdata/list_peak_group_annotation_concised_pg2pg_200805.RData')
#
# annotation_initial <- readr::read_csv('I:/00_projects/03_MetDNA2/00_data/20200805_targeted_annotation_evaluation_200STD/01_targeted_annotation/200STD_neg_reverse/initial_id_result.csv')
#
# test <- refineAnnotation(annotation_initial = annotation_initial,
#                          list_peak_group_annotation_concised = list_peak_group_annotation_concised)

setGeneric(name = 'refineAnnotationFeatureClustering',
           def = function(
             annotation_initial,
             list_peak_group_annotation_concised
           ){
             concised_table <- purrr::map(seq_along(list_peak_group_annotation_concised), function(i){
               result <- tibble::tibble(peak_name_concised = list_peak_group_annotation_concised[[i]]@base_peak_name,
                                        adduct_concised = list_peak_group_annotation_concised[[i]]@base_peak_adduct,
                                        peak_group_scale = nrow(list_peak_group_annotation_concised[[i]]@peak_list_annotated))

             }) %>%
               dplyr::bind_rows() %>%
               dplyr::mutate(label = paste(peak_name_concised, adduct_concised, sep = '_'))

             annotation_concise <- annotation_initial %>%
               dplyr::mutate(label = paste(name, adduct, sep = '_')) %>%
               dplyr::left_join(concised_table, by = 'label') %>%
               dplyr::filter(!is.na(peak_name_concised)) %>%
               dplyr::select(-c('peak_name_concised', 'adduct_concised', 'label'))

             num_peak_initial <- unique(annotation_initial$name) %>% length()
             num_id_initial <- unique(annotation_initial$id) %>% length()
             num_annotation_initial <- nrow(annotation_initial)

             num_peak_concise <- unique(annotation_concise$name) %>% length()
             num_id_concise <- unique(annotation_concise$id) %>% length()
             num_annotation_concise <- nrow(annotation_concise)

             num_stat <- tibble::tibble('No. peak' = c(num_peak_initial,
                                                       num_peak_concise),
                                        'No. metabolite' = c(num_id_initial,
                                                             num_id_concise),
                                        'No. annotation (peak-met)' = c(num_annotation_initial,
                                                                        num_annotation_concise))

             cat('Summary of ID refinement\n')
             rownames(num_stat) <- c('Initial ID', 'Credential ID')
             print(knitr::kable(num_stat))

             return(annotation_concise)

           })


#' @title refineAnnotationFormula
#' @author Zhiwei Zhou
#' @param annotation_initial
#' @param list_peak_group_formula
#' @param num_formula_candidate
#' @export
#' @examples

# load('./inst/tempdata/list_peak_group_formula.RData')
#
# annotation_initial <- readr::read_csv('I:/00_projects/03_MetDNA2/00_data/20200805_targeted_annotation_evaluation_200STD/01_targeted_annotation/200STD_neg_reverse/initial_id_result.csv')

# test <- refineAnnotation(annotation_initial = annotation_initial,
#                          list_peak_group_annotation_concised = list_peak_group_annotation_concised)

setGeneric(name = 'refineAnnotationFormula',
           def = function(
             annotation_initial,
             list_peak_group_formula,
             num_formula_candidate = 3,
             path_dir,
             ...
           ){


             formula_list <- purrr::map(seq_along(list_peak_group_formula), function(i){
               # cat(i, ' ')
               temp_data <- list_peak_group_formula[[i]]


               if (nrow(temp_data@pred_formula_list) > 0) {
                 temp_formula_list <- temp_data@pred_formula_list %>%
                   dplyr::arrange(dplyr::desc(score_inte)) %>%
                   dplyr::slice(seq(num_formula_candidate)) %>%
                   dplyr::mutate(rank_formula = seq(dplyr::n()),
                                 name = temp_data@base_peak_name,
                                 adduct = temp_data@base_peak_adduct) %>%
                   dplyr::select(name, adduct, dplyr::everything())

                 return(temp_formula_list)
               } else {
                 return(NULL)
               }

               # if (nrow(temp_formula_list) > 0) {
               #   temp_formula_list <- temp_formula_list %>%
               #     dplyr::mutate(rank_formula = seq(dplyr::n()),
               #                   name = temp_data@base_peak_name,
               #                   adduct = temp_data@base_peak_adduct) %>%
               #     dplyr::select(name, adduct, dplyr::everything())
               # } else {
               #   return(NULL)
               # }


             }) %>% dplyr::bind_rows()

             # progress <- mapProgress(n = nrow(annotation_initial))
             annotation_formula <- purrr::map(seq_along(annotation_initial$name), function(i){
               # mapProgressPrint(progress)
               # cat(i, ' ')
               temp_peak_name <- annotation_initial$name[i]
               temp_adduct <- annotation_initial$adduct[i]
               temp_annotation_formula <- annotation_initial$formula[i]

               result <- formula_list %>%
                 dplyr::filter(name == temp_peak_name &
                                 adduct == temp_adduct &
                                 formula == temp_annotation_formula) %>%
                 dplyr::select(score_inte, rank_formula) %>%
                 dplyr::rename(score_formula = score_inte)

               if (nrow(result) == 0) {
                 result <- tibble::tibble(score_formula = 0,
                                          rank_formula = 0)
               }

               result
             }) %>% dplyr::bind_rows()

             annotation_formula_filter <- annotation_initial %>%
               dplyr::bind_cols(annotation_formula) %>%
               dplyr::filter(rank_formula > 0)

             # export record_formula_filter
             temp <- sapply(list_peak_group_formula, function(x){paste(x@base_peak_name, x@base_peak_adduct, sep = '_')})
             record_formula_filter <- annotation_initial %>%
               dplyr::bind_cols(annotation_formula) %>%
               dplyr::filter(rank_formula == 0) %>%
               dplyr::mutate(peak_group = paste(name, adduct, sep = '_')) %>%
               dplyr::filter(peak_group %in% temp) %>%
               dplyr::mutate(record = paste(peak_group, formula, sep = '@')) %>%
               dplyr::pull(record) %>%
               unique()

             save(record_formula_filter,
                  file = file.path(path_dir,
                                   '03_annotation_credential',
                                   "00_intermediate_data",
                                   'record_formula_filter.RData'))

             num_peak_initial <- unique(annotation_initial$name) %>% length()
             num_id_initial <- unique(annotation_initial$id) %>% length()
             num_annotation_initial <- nrow(annotation_initial)

             num_peak_formula_filter <- unique(annotation_formula_filter$name) %>% length()
             num_id_formula_filter <- unique(annotation_formula_filter$id) %>% length()
             num_annotation_formula_filter <- nrow(annotation_formula_filter)

             num_stat <- tibble::tibble('No. peak' = c(num_peak_initial,
                                                       num_peak_formula_filter),
                                        'No. metabolite' = c(num_id_initial,
                                                             num_id_formula_filter),
                                        'No. annotation (peak-met)' = c(num_annotation_initial,
                                                                        num_annotation_formula_filter))

             cat('\n\n'); cat('Summary of ID refinement\n')
             rownames(num_stat) <- c('Initial ID', 'Credential formula')
             print(knitr::kable(num_stat))

             return(annotation_formula_filter)

           })
