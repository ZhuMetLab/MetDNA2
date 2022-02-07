################################################################################
# annotateInitialSeed ----------------------------------------------------------

#' @title annotateInitialSeed
#' @description Annotate peak using MS/MS spectra in in-house database.
#' @author Zhiwei Zhou, Xiaotao Shen, Yandong Yin
#' @param ms1_file The name of ms1 peak table. Column 1 is "name", Column 2 is "mz" and column is "rt".
#' @param ms2_file Default: NULL
# #' @param sample_info_file Default: sample.info.csv
#' @param ms2_type "mgf", "mzXML", "msp", "cef"
#' @param metdna_version 'version1', 'version2'. Default: "version1"
#' @param path Default: '.'
#' @param instrument The instrument you used to acquire data. "SciexTripleTOF", "AgilentQTOF", "BrukerQTOF", "ThermoOrbitrap", "ThermoExploris", "AgilentDTIMMS", "BrukerTIMS", 'WatersQTOF','WatersTWIMMS'. Default: "SciexTripleTOF"
#' @param lib Default: 'zhuMetLib'
#' @param column "hilic", "rp"
#' @param ce "10", "20", "30", "35,15", "40", "50". Default: "30"
#' @param method_lc 'Amide12min', 'Amide23min', 'Other', 'MetlinRP'. Default: 'Amide12min'
#' @param excluded_adduct adduct list for exclusion. Default: NULL
#' @param is_rt_calibration Default: FALSE
#' @param mz_tol ms1 match. Default: 25 ppm
#' @param pf_rt_tol penalty-free rt range. Default: 0 s
#' @param tolerance_rt_range maxmium rt tolerance. Default: 30 s
#' @param pf_ccs_tol penalty-free rt range. Default: 0 %
#' @param tolerance_ccs_range maxmium rt tolerance. Default: 2 %
#' @param is_filter Whether filter candidates whose rt/ccs exceed the maximum tolerance. Default: FALSE
#' @param is_rt_score logical vector. Default: TRUE
#' @param is_ccs_score logical vector. Default: FALSE
#' @param is_ms2_score logical vector. Default: TRUE
#' @param is_include_precursor Default: TRUE
#' @param int_ms2_min_abs Default: 30
#' @param int_ms2_min_relative Default: 0.01
#' @param mz_tol_combine_ms1_ms2 Default: 25 ppm
#' @param rt_tol_combine_ms1_ms2 Default: 10 s
#' @param ccs_tol_combine_ms1_ms2 Only appliable for IM-MS data. Default: NULL
#' @param mz_tol_ms2 Default: 35 ppm
#' @param dp_cutoff The cutoff of dot product. Default is 0.8.
#' @param matched_frag_cutoff The cutoff of matched fragment number. Default: 1
#' @param direction Reserve which direct dot product. 'reverse', 'forward'. Default: 'reverse'
#' @param scoring_approach spectral match approach, including 'dp', 'bonanza', 'hybrid', 'gnps'. Default: 'dp'.
#' @param is_plot_ms2 Output MS2 match plot or not. Default: TRUE
#' @return Return the annotation result.
#' @export


# ms1_file = "data.csv"
# ms2_file = NULL
# sample_info_file = "sample.info.csv"
# metdna_version = 'version2'
# ms2_type = 'mgf'
# # # ms2_type = 'msp'
# # path = "/home/zhouzw/Data_processing/20201102_metdna2_pos_development/"
# path = '/home/zhouzw/Data_processing/20210323_metdna2_scoring_evaluation/para_1_dp/'
#
# # para of loadDB
# polarity = 'positive'
# instrument = "SciexTripleTOF"
# lib = 'zhuMetLib'
# column = 'hilic'
# ce = '30'
# # method_lc = 'Amide23min'
# method_lc = 'Amide23min'
# excluded_adduct = NULL
# is_rt_calibration = TRUE
#
# # para of matchMs1WithSpecLib
# mz_tol = 25
# pf_rt_range = 0
# tolerance_rt_range = 30
# pf_ccs_range = 0
# tolerance_ccs_range = 4
# is_filter = TRUE
# is_rt_score = TRUE
# is_ccs_score = TRUE
# is_msms_score = TRUE
#
# # parameters of matchMs2WithSpecLib
# is_include_precursor = TRUE
# is_deisotope = FALSE
# int_ms2_min_abs = 50
# int_ms2_min_relative = 0.01
# ppm_precursor_filter = 20
# mz_range_ms2 = c(0, 1700)
#
# mz_tol_combine_ms1_ms2 = 25 # ppm
# rt_tol_combine_ms1_ms2 = 10 # s
# ccs_tol_combine_ms1_ms2 = NULL # %
# mz_tol_ms2 = 35
# dp_cutoff = 0.8
# matched_frag_cutoff = 1
# direction = 'reverse'
# scoring_approach = 'dp'
#
# is_plot_ms2 <- TRUE

# annotateInitialSeed(ms1_file = "data.csv",
#                     ms2_file = NULL,
#                     metdna_version = 'version2',
#                     ms2_type = 'mgf',
#                     path = "/home/zhouzw/Data_processing/20210330_debug/",
#
#                     polarity = 'positive',
#                     instrument = "SciexTripleTOF",
#                     lib = 'zhuMetLib',
#                     column = 'hilic',
#                     ce = '30',
#                     method_lc = 'Amide23min',
#                     excluded_adduct = NULL,
#                     is_rt_calibration = FALSE,
#
#                     mz_tol = 25,
#                     pf_rt_range = 0,
#                     tolerance_rt_range = 30,
#                     pf_ccs_range = 0,
#                     tolerance_ccs_range = 2,
#                     is_filter = TRUE,
#                     is_rt_score = TRUE,
#                     is_ccs_score = FALSE,
#                     is_ms2_score = TRUE,
#
#                     is_include_precursor = TRUE,
#                     int_ms2_min_abs = 50,
#                     int_ms2_min_relative = 0.01,
#                     mz_tol_combine_ms1_ms2 = 25, # ppm
#                     rt_tol_combine_ms1_ms2 = 10, # s
#                     ccs_tol_combine_ms1_ms2 = NULL, # %
#                     mz_tol_ms2 = 35,
#                     dp_cutoff = 0.8,
#                     matched_frag_cutoff = 1,
#                     direction = 'reverse',
#                     scoring_approach = 'dp',
#                     # direction = 'reverse',
#
#                     is_plot_ms2 <- TRUE
#                     )

setGeneric(name = "annotateInitialSeed",
           def = function(ms1_file = "data.csv",
                          ms2_file = NULL,
                          # sample_info_file = "sample.info.csv",
                          metdna_version = c('version2', 'version1'),
                          ms2_type = c("mgf", "mzXML", "msp", "cef"),
                          path = ".",

                          # parameters of loadDB
                          polarity = c("positive", "negative"),
                          instrument = c("SciexTripleTOF", "AgilentQTOF", "BrukerQTOF", "ThermoOrbitrap", "ThermoExploris",
                                         'WatersQTOF', 'WatersTWIMMS', "AgilentDTIMMS", "BrukerTIMS"),
                          lib = c('zhuMetLib', 'zhuMetLib_orbitrap', 'fiehnHilicLib', 'zhuRPLib'),
                          column = c("hilic", "rp"),
                          ce = c("30", "10", "20", "35,15", "40", "50",
                                 'NCE10', 'NCE20', 'NCE30', 'NCE40', 'NCE50',
                                 'SCE20_30_40%', "SNCE20_30_40%"),
                          method_lc = c('Amide12min', 'Amide23min', 'Other', 'MetlinRP', 'zhulabRP'),
                          excluded_adduct = NULL,
                          is_rt_calibration = FALSE,

                          # parameters of matchMs1WithSpecLib
                          mz_tol = 25,
                          mz_ppm_thr = 400,
                          pf_rt_range = 0,
                          tolerance_rt_range = 30,
                          pf_ccs_range = 0,
                          tolerance_ccs_range = 2,
                          is_filter = TRUE,
                          is_rt_score = TRUE,
                          is_ccs_score = FALSE,
                          is_ms2_score = TRUE,

                          # parameters of matchMs2WithSpecLib
                          is_include_precursor = TRUE,
                          # is_deisotope = FALSE,
                          int_ms2_min_abs = 50,
                          int_ms2_min_relative = 0.01,
                          # ppm_precursor_filter = 20,
                          # mz_range_ms2 = c(0, 1700),
                          mz_tol_combine_ms1_ms2 = 25, # ppm
                          rt_tol_combine_ms1_ms2 = 10, # s
                          ccs_tol_combine_ms1_ms2 = NULL, # %
                          mz_tol_ms2 = 35,
                          dp_cutoff = 0.8,
                          matched_frag_cutoff = 1,
                          direction = c('reverse', 'forward'),
                          scoring_approach = c('dp', 'bonanza', 'hybrid', 'gnps'),
                          test_adduct_version = c('version2', 'version1'),
                          test_evaluation = c('No', '200STD', '46STD'),
                          test_force_filtering_rt = NULL, # numeric para, if applied, it will remove annotations in final result larger than this value (s)

                          is_plot_ms2 = TRUE
           ){

             polarity <- match.arg(polarity)
             column <- match.arg(column)
             ms2_type <- match.arg(ms2_type)
             instrument <- match.arg(instrument)
             ce <- match.arg(ce)
             direction <- match.arg(direction)
             test_adduct_version <- match.arg(test_adduct_version)
             test_evaluation <- match.arg(test_evaluation)

             # browser()
             # path_output <- file.path(path, "MS2_match_result")
             path_output <- file.path(path, "01_result_initial_seed_annotation")
             dir.create(file.path(path_output, "00_intermediate_data"), showWarnings = FALSE, recursive = TRUE)

             # check ms1_file and ms2_file
             file_list <- list.files(path)
             if (all(file_list != ms1_file)) {stop("There is not ms1_file in you directory.")}

             # if (length(grep(ms2_type, file_list)) == 0) {stop("There are not ms2_file in you directory.")}

             if (length(ms2_file) == 0) {
               ms2_file <- grep(paste0('\\.(', ms2_type, ')$'), file_list, value = TRUE)
             }

             # change parameters according to version
             switch(metdna_version,
                    'version1' = {
                      if (polarity == 'positive') {
                        is_include_precursor = FALSE
                      } else {
                        is_include_precursor = TRUE
                      }

                      mz_ppm_thr <- 0
                      int_ms2_min_abs <- 30
                      direction <- 'reverse'
                      method_lc <- 'Other'
                      is_check_cation <- FALSE
                    },
                    'version2' = {
                      # mz_ppm_thr <- 400
                      is_check_cation <- TRUE
                    })


             # cat("\n", "--------------------------------------------------------------\n")
             cat('Load database ...\n')

             # # metdna_version1: keep same with MetDNA v1.21
             # switch (metdna_version,
             #         'version1' = {
             #           if (polarity == 'positive') {
             #             if (column == 'hilic') {
             #               adduct_list <- c("M+", "[M+H]+", "[M+NH4]+", "[M-H2O+H]+",  "[M-2H2O+H]+",
             #                                "[M+Na]+", "[M-H+2Na]+", "[M-2H+3Na]+", "[M+K]+", "[M-H+2K]+",
             #                                "[M-2H+3K]+", "[M+CH3CN+H]+", "[M+CH3CN+Na]+", "[2M+H]+", "[2M+NH4]+",
             #                                "[2M+Na]+", "[2M+K]+", "[M+CH3COO+2H]+")
             #             } else {
             #               adduct_list <- c("M+", "[M+H]+", "[M+NH4]+", "[M-H2O+H]+",  "[M-2H2O+H]+",
             #                                "[M+Na]+", "[M-H+2Na]+", "[M-2H+3Na]+", "[M+K]+", "[M-H+2K]+",
             #                                "[M-2H+3K]+", "[M+CH3CN+H]+", "[M+CH3CN+Na]+", "[2M+H]+", "[2M+NH4]+",
             #                                "[2M+Na]+", "[2M+K]+", "[M+CH3COO+2H]+")
             #             }
             #           } else {
             #             if (column == 'hilic') {
             #               adduct_list <- c("M-", "[M-H]-", "[M-H2O-H]-", "[M+Na-2H]-", "[M+K-2H]-",
             #                                "[M+NH4-2H]-", "[2M-H]-", "[M+CH3COO]-")
             #             } else {
             #               adduct_list <- c("M-", "[M-H]-", "[M-H2O-H]-", "[M+Na-2H]-", "[M+K-2H]-",
             #                                "[M+NH4-2H]-", "[2M-H]-", "[M+F]-")
             #             }
             #           }
             #
             #         },
             #         'version2' = {
             #           if (polarity == 'positive') {
             #             adduct_list <- lib_adduct_nl$positive %>% dplyr::filter(annotation == 'Yes') %>% dplyr::pull(adduct)
             #           } else {
             #             adduct_list <- lib_adduct_nl$negative %>% dplyr::filter(annotation == 'Yes' | adduct %in% c('[M-2H]-')) %>% dplyr::pull(adduct)
             #           }
             #
             #           if (length(excluded_adduct) > 0) {
             #             adduct_list <- adduct_list[!which(adduct_list %in% excluded_adduct)]
             #           }
             #         }
             # )


             # metdna_version1: keep same with MetDNA v1.21
             switch (test_adduct_version,
                     'version1' = {
                       if (polarity == 'positive') {
                         if (column == 'hilic') {
                           adduct_list <- c("[M]+", "[M+H]+", "[M+NH4]+", "[M-H2O+H]+",  "[M-2H2O+H]+",
                                            "[M+Na]+", "[M-H+2Na]+", "[M-2H+3Na]+", "[M+K]+", "[M-H+2K]+",
                                            "[M-2H+3K]+", "[M+CH3CN+H]+", "[M+CH3CN+Na]+", "[2M+H]+", "[2M+NH4]+",
                                            "[2M+Na]+", "[2M+K]+", "[M+CH3COO+2H]+")
                         } else {
                           adduct_list <- c("[M]+", "[M+H]+", "[M+NH4]+", "[M-H2O+H]+",  "[M-2H2O+H]+",
                                            "[M+Na]+", "[M-H+2Na]+", "[M-2H+3Na]+", "[M+K]+", "[M-H+2K]+",
                                            "[M-2H+3K]+", "[M+CH3CN+H]+", "[M+CH3CN+Na]+", "[2M+H]+", "[2M+NH4]+",
                                            "[2M+Na]+", "[2M+K]+", "[M+CH3COO+2H]+")
                         }
                       } else {
                         if (column == 'hilic') {
                           adduct_list <- c("[M]-", "[M-H]-", "[M-H2O-H]-", "[M+Na-2H]-", "[M+K-2H]-",
                                            "[M+NH4-2H]-", "[2M-H]-", "[M+CH3COO]-")
                         } else {
                           adduct_list <- c("[M]-", "[M-H]-", "[M-H2O-H]-", "[M+Na-2H]-", "[M+K-2H]-",
                                            "[M+NH4-2H]-", "[2M-H]-", "[M+F]-")
                         }
                       }

                     },
                     'version2' = {
                       if (polarity == 'positive') {
                         adduct_list <- lib_adduct_nl$positive %>% dplyr::filter(annotation == 'Yes' | adduct %in% c('[M]+')) %>% dplyr::pull(adduct)
                       } else {
                         adduct_list <- lib_adduct_nl$negative %>% dplyr::filter(annotation == 'Yes' | adduct %in% c('[M-2H]-')) %>% dplyr::pull(adduct)
                       }

                       if (length(excluded_adduct) > 0) {
                         adduct_list <- adduct_list[which(!(adduct_list %in% excluded_adduct))]
                       }
                     }
             )

             lib_db <- loadSpecDB(lib = lib,
                                  instrument = instrument,
                                  column = column,
                                  method_lc = method_lc,
                                  ce = ce,
                                  polarity = polarity,
                                  adduct_list = adduct_list,
                                  is_rt_calibration = is_rt_calibration,
                                  path = path)

             lib_meta <- lib_db$lib_meta
             lib_spec <- lib_db$lib_spec

             if (test_evaluation == '200STD') {
               data("cpd_200stdExd", envir = environment())
               id_200std <- cpd_200stdExd %>% dplyr::filter(source == 'Precursor') %>% dplyr::pull(id)
               idx <- which(names(lib_spec) %in% id_200std)
               lib_meta <- lib_meta %>% dplyr::filter(id %in% id_200std)
               lib_spec <- lib_spec[idx]

               rm(list = c('id_200std', 'idx', 'cpd_200stdExd'));gc()
             }

             if (test_evaluation == '46STD') {
               data("cpd_46stdExd", envir = environment())
               id_46std <- cpd_46stdExd %>% dplyr::filter(source == 'Precursor') %>% dplyr::pull(id)
               idx <- which(names(lib_spec) %in% id_46std)
               lib_meta <- lib_meta %>% dplyr::filter(id %in% id_46std)
               lib_spec <- lib_spec[idx]

               rm(list = c('id_46std', 'idx', 'cpd_46stdExd'));gc()
             }

             save(lib_meta,
                  file = file.path(path_output, "00_intermediate_data", 'lib_meta'),
                  compress = 'gzip',
                  version = 2)

             save(lib_spec,
                  file = file.path(path_output, "00_intermediate_data", 'lib_spec'),
                  compress = 'gzip',
                  version = 2)

             # mz, rt, ccs match ...............................................
             ms1_data <- readMs1(filename = ms1_file,
                                 instrument = instrument,
                                 path = path)

             save(ms1_data,
                  file = file.path(path_output, "00_intermediate_data", 'ms1_data'),
                  compress = 'gzip',
                  version = 2)

             ms1_result <- matchMs1WithSpecLib(ms1_data = ms1_data,
                                               lib_meta = lib_meta,
                                               mz_tol = mz_tol,
                                               mz_ppm_thr = mz_ppm_thr,
                                               pf_rt_range = pf_rt_range,
                                               tolerance_rt_range = tolerance_rt_range,
                                               pf_ccs_range = pf_ccs_range, # penalty-free CCS range-percentage
                                               tolerance_ccs_range = tolerance_ccs_range, # percentage
                                               is_filter = is_filter,
                                               is_rt_score = is_rt_score,
                                               is_ccs_score = is_ccs_score,
                                               is_check_cation = is_check_cation)

             save(ms1_result,
                  file = file.path(path_output, "00_intermediate_data", 'ms1_result'),
                  compress = 'gzip',
                  version = 2)


             if (is_ms2_score) {
               cat("\n")

               if ('ms2' %in% list.files(file.path(path_output, "00_intermediate_data"))) {
                 cat('Read, purify ms2 spec, and combine with ms1 & ms2...\n')
                 cat('Note: load the existed ms2\n\n')
                 load(file.path(path_output, "00_intermediate_data", 'ms2'))
               } else {
                 cat('\n');cat('Read, purify ms2 spec, and combine with ms1 & ms2...\n')
                 # browser()

                 ms2_data <- readMs2(ms2_file = file.path(path, ms2_file),
                                     ms2_type = ms2_type)

                 ms2_data <- inteMs2(ms2_data = ms2_data,
                                     ms2_type = ms2_type,
                                     metdna_version = metdna_version,
                                     instrument = instrument,
                                     is_include_precursor = is_include_precursor,
                                     is_deisotope = FALSE,
                                     int_ms2_min_abs = int_ms2_min_abs,
                                     int_ms2_min_relative = int_ms2_min_relative,
                                     mz_range_ms2 = NULL)

                 ms2 <- combineMs1Ms2(ms1_data = ms1_data,
                                      ms2_data = ms2_data,
                                      ms2_type = ms2_type,
                                      mz_tol_combine_ms1_ms2 = mz_tol_combine_ms1_ms2,
                                      rt_tol_combine_ms1_ms2 = rt_tol_combine_ms1_ms2,
                                      ccs_tol_combine_ms1_ms2 = ccs_tol_combine_ms1_ms2,
                                      path = file.path(path_output, "00_intermediate_data"))

                 remove(list = c('ms2_data'))
                 gc()
               }

               # browser()

               if (!('result_annotation' %in% list.files(file.path(path_output, "00_intermediate_data")))) {
                 cat('\n');cat('start ms2 matching...\n')
                 # ms2 order and name matched with ms1
                 ms2_name <- unname(unlist(lapply(ms2, function(x) x[[1]][1,])))
                 ms1_name <- ms1_data$info$name
                 temp_idx <- match(ms2_name, ms1_name)
                 exp_spec <- vector(mode = "list", length = nrow(ms1_data$info))
                 names(exp_spec) <- ms1_name
                 exp_spec[temp_idx] <- ms2

                 exp_spec <- lapply(exp_spec, function(x){
                   if(is.null(x)) return(x)
                   x[[2]]
                 })

                 result_annotation <- matchMs2WithSpecLib(ms1_result = ms1_result,
                                                          exp_spec = exp_spec,
                                                          lib_meta = lib_meta,
                                                          lib_spec = lib_spec,
                                                          metdna_version = metdna_version,
                                                          mz_tol_ms2 = mz_tol_ms2,
                                                          dp_cutoff = dp_cutoff,
                                                          matched_frag_cutoff = matched_frag_cutoff,
                                                          direction = direction,
                                                          scoring_approach = scoring_approach,
                                                          path = path_output,
                                                          is_include_precursor = is_include_precursor)

                 save(result_annotation,
                      file = file.path(path_output, '00_intermediate_data', 'result_annotation'),
                      version = 2)
               } else {
                 cat('start ms2 matching...\n')
                 cat('Note: load the existed result_annotation\n')

                 load(file.path(path_output, '00_intermediate_data', 'result_annotation'))
               }

             } else {
               result_annotation <- ms1_result
             }

             # enforce filter absolute RT error less than this cutoff in final result table
             # modified in 20220115
             if (length(test_force_filtering_rt) > 0) {
               result_annotation <- lapply(result_annotation, function(x){
                 if (nrow(x@annotation_result) > 0) {
                   x@annotation_result <- x@annotation_result %>% dplyr::filter(rt_error <= test_force_filtering_rt | is.na(rt_error))

                   return(x)
                 } else {
                   return(x)
                 }
               })

               save(result_annotation,
                    file = file.path(path_output, '00_intermediate_data', 'result_annotation'),
                    version = 2)
             }

             table_annotation <- convertSpecAnnotationClass2Table(ms1_data = ms1_data,
                                                                  result_annotation = result_annotation,
                                                                  lib_meta = lib_meta,
                                                                  instrument = instrument,
                                                                  is_rt_score = is_rt_score,
                                                                  is_ccs_score = is_ccs_score,
                                                                  is_msms_score = is_msms_score,
                                                                  dp_cutoff = dp_cutoff)

             readr::write_csv(table_annotation,
                              file.path(path_output, "ms2_match_annotation_result.csv"))

             if (is_plot_ms2 & is_ms2_score) {
               load(file.path(path_output, '00_intermediate_data', 'ms2_result'))

               cat('\n')
               # table_annotation <- readr::read_csv(file.path(path_output, "ms2_match_annotation_result.csv"))

               if (scoring_approach == 'dp') {
                 purrr::walk(c('forward', 'reverse'), function(x){
                   cat('Plot', x, 'spec match...\n')

                   if (x == 'forward') {
                     temp_result <- table_annotation %>%
                       dplyr::filter(!is.na(id_forward_summary)) %>%
                       dplyr::rename(id = id_forward_summary)
                   } else {
                     temp_result <- table_annotation %>%
                       dplyr::filter(!is.na(id_reverse_summary)) %>%
                       dplyr::rename(id = id_reverse_summary)
                   }

                   progress <- mapProgress(n = nrow(temp_result))
                   purrr::walk(seq_along(temp_result$name), function(i){
                     mapProgressPrint(progress = progress)
                     # cat('i', i, ' ')

                     temp_feature_name <- temp_result$name[i]
                     dir.create(file.path(path_output,
                                          '02_experimental_ms2_spec_plot',
                                          paste0(temp_feature_name, '_', x)),
                                showWarnings = FALSE, recursive = TRUE)

                     temp_id <- temp_result %>%
                       dplyr::filter(name == temp_feature_name) %>%
                       dplyr::select(name:rt, id) %>%
                       tidyr::separate_rows(id, sep = ';') %>%
                       dplyr::pull(id)

                     cpd_id <- temp_id %>%
                       stringr::str_extract(pattern = 'labid\\{[A-Za-z0-9]+\\}') %>%
                       gsub(pattern = 'labid\\{', replacement = '') %>%
                       gsub(pattern = '\\}', replacement = '')

                     cpd_name <- temp_id %>%
                       stringr::str_extract(pattern = "name\\{[^\\{]+\\}") %>%
                       gsub(pattern = 'name\\{', replacement = '') %>%
                       gsub(pattern = '\\}', replacement = '')

                     cpd_score <- temp_id %>%
                       stringr::str_extract(pattern = 'score\\{0\\.[0-9]+\\}|score\\{1\\}') %>%
                       gsub(pattern = 'score\\{', replacement = '') %>%
                       gsub(pattern = '\\}', replacement = '')

                     # temp_id %>%
                     #   stringr::str_extract(pattern = 'score\\{(([1-9]\\d*\\.?\\d*)|(0\\.\\d*[1-9]))\\}') %>%
                     #   gsub(pattern = 'score\\{', replacement = '') %>%
                     #   gsub(pattern = '\\}', replacement = '')

                     idx <- which(names(ms2_result) == temp_feature_name)
                     temp_ms2_obj <- ms2_result[[idx]]

                     idx_id <- match(cpd_id, temp_ms2_obj@info$name)

                     purrr::walk(seq_along(idx_id), function(j){
                       # cat('j', j, ' ')
                       temp_idx <- idx_id[j]
                       suppressMessages(
                         temp_plot <- plotIdMs2(obj_spec = temp_ms2_obj@matchedFragments[[temp_idx]]) +
                           ggplot2::scale_colour_manual(
                             name = 'attribute',
                             labels= c(paste0('Feature ', '(', temp_feature_name, ')'),
                                       'Unmatched fragments',
                                       paste0('Library ', '(', cpd_id, ')')),
                             values = c(
                               'experiment' = 'dodgerblue',
                               'library' = 'tomato',
                               'frag_unmatch' = 'gray'
                             )
                           ) +
                           ggplot2::ggtitle(label = paste0(cpd_id[j],
                                                           ' --- ', cpd_name[j],
                                                           ' (DP: ', cpd_score[j], ')'))
                       )

                       ggplot2::ggsave(temp_plot,
                                       filename = file.path(path_output,
                                                            '02_experimental_ms2_spec_plot',
                                                            paste0(temp_feature_name, '_', x),
                                                            paste0(cpd_id[j], '.pdf')),
                                       width = 10, height = 6)

                     })
                   })


                   cat('\n\n')
                 })
               } else {
                 cat('Plot shifted spec match...\n')

                 temp_result <- table_annotation %>%
                   dplyr::filter(!is.na(id_reverse_summary)) %>%
                   dplyr::rename(id = id_reverse_summary)

                 progress <- mapProgress(n = nrow(temp_result))
                 purrr::walk(seq_along(temp_result$name), function(i){
                   # cat(i, ' ')
                   mapProgressPrint(progress = progress)

                   temp_feature_name <- temp_result$name[i]
                   dir.create(file.path(path_output,
                                        '02_experimental_ms2_spec_plot',
                                        paste0(temp_feature_name)),
                              showWarnings = FALSE, recursive = TRUE)

                   temp_id <- temp_result %>%
                     dplyr::filter(name == temp_feature_name) %>%
                     dplyr::select(name:rt, id) %>%
                     tidyr::separate_rows(id, sep = ';') %>%
                     dplyr::pull(id)

                   cpd_id <- temp_id %>%
                     stringr::str_extract(pattern = 'labid\\{[A-Za-z0-9]+\\}') %>%
                     gsub(pattern = 'labid\\{', replacement = '') %>%
                     gsub(pattern = '\\}', replacement = '')

                   cpd_name <- temp_id %>%
                     stringr::str_extract(pattern = "name\\{[^\\{]+\\}") %>%
                     gsub(pattern = 'name\\{', replacement = '') %>%
                     gsub(pattern = '\\}', replacement = '')

                   cpd_score <- temp_id %>%
                     stringr::str_extract(pattern = 'score\\{0\\.[0-9]+\\}|score\\{1\\}') %>%
                     gsub(pattern = 'score\\{', replacement = '') %>%
                     gsub(pattern = '\\}', replacement = '')

                   # temp_id %>%
                   #   stringr::str_extract(pattern = 'score\\{(([1-9]\\d*\\.?\\d*)|(0\\.\\d*[1-9]))\\}') %>%
                   #   gsub(pattern = 'score\\{', replacement = '') %>%
                   #   gsub(pattern = '\\}', replacement = '')

                   idx <- which(names(ms2_result) == temp_feature_name)
                   temp_ms2_obj <- ms2_result[[idx]]

                   idx_id <- match(cpd_id, temp_ms2_obj@info$name)

                   purrr::walk(seq_along(idx_id), function(j){
                     temp_idx <- idx_id[j]
                     suppressMessages(
                       temp_plot <- plotIdShiftMs2(obj_spec_frag_match = temp_ms2_obj@matchedFragments[[temp_idx]],
                                                   obj_spec_frag_nl = temp_ms2_obj@nlFragments[[temp_idx]]) +
                         ggplot2::scale_colour_manual(
                           name = 'attribute',
                           labels= c(paste0('Feature ', '(', temp_feature_name, ')'),
                                     'Unmatched fragments',
                                     paste0('Library ', '(', temp_ms2_obj@info$name[temp_idx], ')'),
                                     paste0('Library shift', '(', temp_ms2_obj@info$name[temp_idx], ')')),
                           values = c(
                             'experiment' = 'dodgerblue',
                             'library' = 'tomato',
                             'frag_unmatch' = 'gray',
                             'library_shift' = 'orange'
                           )
                         ) +
                         ggplot2::ggtitle(label = paste0(cpd_id[j],
                                                         ' --- ', cpd_name[j],
                                                         ' (', scoring_approach, ':', cpd_score[j], ')'))
                     )

                     ggplot2::ggsave(temp_plot,
                                     filename = file.path(path_output,
                                                          '02_experimental_ms2_spec_plot',
                                                          paste0(temp_feature_name),
                                                          paste0(cpd_id[j], '.pdf')),
                                     width = 10, height = 6)

                   })
                 })


                 cat('\n\n')

               }

             }

           })


################################################################################
# readMs1 ----------------------------------------------------------------------

#' @title readMs1
#' @author Zhiwei Zhou
#' @param filename The name of csv file.
#' @param instrument "SciexTripleTOF", "AgilentQTOF", "BrukerQTOF", "ThermoOrbitrap", 'ThermoExploris', "AgilentDTIMMS", "BrukerTIMS", 'WatersQTOF','WatersTWIMMS'
#' @param path '.'
#' @return A list. info: a table of feature information, column 1_feature name, column 2_mz, column3_rt in seconds, column4_ccs. subject: a table of peak area.
#' @importFrom magrittr '%>%'
#' @export

# filename <- 'I:/00_projects/co_worker/CaiYuping/HILICneg_modified.csv'
# instrument <- 'lcms'
# test <- readMs1(filename = 'data.csv',
#                 instrument = "SciexTripleTOF",
#                 path = 'H:/00_projects/03_MetDNA2/00_data/20200922_metdna2_development/pos_201007_metdna2')

setGeneric('readMs1',
           def = function(filename,
                          instrument = c("SciexTripleTOF", "AgilentQTOF", "BrukerQTOF", "ThermoOrbitrap", "ThermoExploris",
                                         'WatersQTOF', 'WatersTWIMMS', "AgilentDTIMMS", "BrukerTIMS"),
                          path = '.') {

             instrument <- match.arg(instrument)

             options(readr.num_columns = 0)
             temp <- readr::read_csv(file.path(path, filename),
                                     col_types = readr::cols())

             if (instrument %in% c("AgilentQTOF", "SciexTripleTOF", 'BrukerQTOF', 'WatersQTOF',
                                  "ThermoOrbitrap", "ThermoExploris")) {
               if (sum(tolower(colnames(temp)[1:3])==c('name', "mz", "rt"))!=3){
                 stop("Please check the format of 'MS1_table.csv'. [The name of column 1-3 must be named as 'name', 'mz', 'rt']")
               }

               colnames(temp)[1:3] <- c("name", "mz", "rt")
               info <- data.frame(name = temp$name,
                                  mz = as.numeric(temp$mz),
                                  rt = as.numeric(temp$rt),
                                  ccs = NA,
                                  stringsAsFactors = F)

               if (max(info$rt) < 60) {
                 info$rt <- info$rt*60
               }

               subject <- temp[, -c(1:3), drop = F]
               rownames(subject) <- info$feature

             }

             if (instrument %in% c("AgilentDTIMMS", "BrukerTIMS", 'WatersTWIMMS')) {

               if (sum(tolower(colnames(temp)[1:3])==c("mz", "rt", "ccs"))!=3){
                 stop("Please check the format of 'MS1_table.csv'. [The name of column 1-3 must be named as 'mz', 'rt', 'ccs']")
               }

               colnames(temp)[1:3] <- c("mz", "rt", "ccs")
               info <- data.frame(mz = round(as.numeric(temp$mz), digits = 4),
                                  rt = round(as.numeric(temp$rt)*60, digits = 0),
                                  ccs = round(as.numeric(temp$ccs), digits = 1),
                                  # intmed = temp %>%
                                  #   dplyr::select(-(mz:ccs)) %>%
                                  #   apply(., 1, function(x){median(x)}) %>%
                                  #   round(digits = 0),
                                  stringsAsFactors = F)
               subject <- temp[, -c(1:3), drop = F]

               fea_name <- paste("M", round(info$mz, digits = 0), "T", round(info$rt, digits = 0), "C", round(info$ccs, digits = 0), sep = "")
               fea_name <- featureReName(name = fea_name)
               info <- data.frame(name=fea_name, info, stringsAsFactors = F)
               rownames(subject) <- fea_name

             }

             result <- list(info=info, subject=subject)

             return(result)

           })


#   featureReName --------------------------------------------------------------
#' @title featureReName
#' @author Zhiwei Zhou
#' @description define feature name of replicates.
#' @param name a vector. character.

setGeneric('featureReName',
           def = function(name){
             sapply(seq(length(name)), function(i){
               temp.name <- name[i]
               temp.idx <- which(name==temp.name)
               if (length(temp.idx) > 1) {
                 temp.name <- paste(temp.name, which(temp.idx==i), sep = "_")
               }
               temp.name
             })
           })



################################################################################
# loadSpecDB -------------------------------------------------------------------

#' @title loadSpecDB
#' @author Zhiwei Zhou
#' @param lib database name, 'zhuMetLib', 'zhuMetLib_orbitrap', 'fiehnHilicLib'. Default: 'zhuMetLib'
#' @param instrument "SciexTripleTOF", "AgilentQTOF", "BrukerQTOF", "ThermoOrbitrap", "ThermoExploris", "AgilentDTIMMS", "BrukerTIMS". Default: "SciexTripleTOF"
#' @param column 'hilic', 'rp'. Default: 'hilic'
#' @param method_lc 'Amide12min', Amide23min', 'Other', 'MetlinRP'
#' @param ce "10", "20", "30", "35,15", "40", "50"; Default: '30'
#' @param polarity 'positive' or 'negative'. Default: 'positive'
#' @param adduct_list NULL
# #' @param is_rt_score whether only reserve compounds with RT
# #' @param is_ccs_score whether only reserve compounds with CCS
#' @param is_rt_calibration TRUE
#' @param path '.'
#' @export

# lib = c('zhuMetLib'),
# instrument = c("SciexTripleTOF", "AgilentQTOF", "BrukerQTOF", "ThermoOrbitrap", "AgilentDTIMMS", "BrukerTIMS"),
# column = c('hilic', 'rp'),
# method_lc = c('Amide12min', 'Amide23min'),
# ce = c("30", "10", "20", "35,15", "40", "50"),
# polarity = c('positive', 'negative'),
# adduct_list = NULL,
# # is_rt_score = FALSE,
# # is_ccs_score = FALSE,
# is_rt_calibration = TRUE,
# path = '.'

# lib = c('zhuMetLib')
# instrument = "SciexTripleTOF"
# column = 'hilic'
# method_lc = 'Amide12min'
# ce = "30"
# polarity = 'positive'
# adduct_list = NULL
# is_rt_score = FALSE
# is_ccs_score = FALSE

# test <- loadSpecDB(lib = 'zhuMetLib',
#                    instrument = 'SciexTripleTOF',
#                    column = 'hilic',
#                    method_lc = 'Amide12min',
#                    ce = '30',
#                    adduct_list = adduct_list)
#

# adduct_list <- lib_adduct_nl$positive %>% dplyr::filter(annotation == 'Yes') %>% dplyr::pull(adduct)
# test <- loadSpecDB(lib = 'zhuMetLib_orbitrap',
#                    instrument = 'ThermoOrbitrap',
#                    polarity = 'positive',
#                    column = 'hilic',
#                    method_lc = 'Amide23min',
#                    ce = 'SCE20_30_40%',
#                    adduct_list = adduct_list,
#                    path = path,
#                    is_rt_calibration = FALSE)

setGeneric(name = 'loadSpecDB',
           def = function(
             lib = c('zhuMetLib', 'zhuMetLib_orbitrap', 'fiehnHilicLib', 'zhuRPLib'),
             instrument = c("SciexTripleTOF", "AgilentQTOF", "BrukerQTOF", "ThermoOrbitrap", "ThermoExploris",
                            'WatersQTOF','WatersTWIMMS', "AgilentDTIMMS", "BrukerTIMS"),
             column = c('hilic', 'rp'),
             method_lc = c('Amide12min', 'Amide23min', 'MetlinRP', 'Other', 'zhulabRP'),
             ce = c("30", "10", "20", "35,15", "40", "50",
                    'NCE10', 'NCE20', 'NCE30', 'NCE40', 'NCE50',
                    'SCE20_30_40%', "SNCE20_30_40%"),
             polarity = c('positive', 'negative'),
             adduct_list = NULL,
             # is_rt_score = FALSE,
             # is_ccs_score = FALSE,
             file_rt_ref = 'RT_recalibration_table.csv',
             is_rt_calibration = TRUE,
             path = '.'
           ){
             lib <- match.arg(lib)
             instrument <- match.arg(instrument)
             column <- match.arg(column)
             ce <- match.arg(ce)
             method_lc <- match.arg(method_lc)
             polarity <- match.arg(polarity)

             # browser()
             if (method_lc == 'Other') {
               cat('Warning: Other LC system cannot support RT calibration.\n')
               is_rt_calibration <- FALSE
             }

             # check adduct type
             if (length(adduct_list) == 0) {
               stop('Please input adduct_list\n')
             }

             if (polarity == 'positive') {
               temp <- all(adduct_list %in% lib_adduct_nl$positive$adduct)
               if(!temp) stop('Please check the selected adducts!\n')
             } else {
               temp <- all(adduct_list %in% lib_adduct_nl$negative$adduct)
               if(!temp) stop('Please check the selected adducts!\n')
             }

             switch(instrument,
                    "AgilentQTOF" = {
                      data("zhuMetlib", envir = environment())
                      lib_data <- zhuMetlib
                      if (ce == '35,15') {
                        ce <- '35,15'
                      } else {
                        ce <- as.character(as.numeric(ce) + 10)
                      }
                    },
                    "SciexTripleTOF" = {
                      data("zhuMetlib", envir = environment())
                      lib_data <- zhuMetlib
                    },
                    'BrukerQTOF' = {
                      data("zhuMetlib", envir = environment())
                      lib_data <- zhuMetlib
                    },
                    'WatersQTOF' = {
                      data("zhuMetlib", envir = environment())
                      lib_data <- zhuMetlib
                      ce <- as.character(as.numeric(ce) + 10)
                    },
                    "OtherQTOF" = {
                      data("zhuMetlib", envir = environment())
                      lib_data <- zhuMetlib
                      ce <- as.character(as.numeric(ce) + 10)
                    },
                    "ThermoOrbitrap" = {
                      # data("orbitrapMetlib", envir = environment())
                      # lib_data = orbitrapMetlib
                      #
                      # if (lib == 'zhuMetLib_orbitrap') {
                      #   data('zhuMetlib_orbitrap', envir = environment())
                      #   lib_data <- zhuMetlib_orbitrap
                      # } else {
                      #   data('fiehnHilicLib', envir = environment())
                      #   lib_data <- fiehnHilicLib
                      # }

                      switch (lib,
                              'zhuMetLib_orbitrap' = {
                                data('zhuMetlib_orbitrap', envir = environment())
                                lib_data <- zhuMetlib_orbitrap
                              },
                              'fiehnHilicLib' = {
                                data('fiehnHilicLib', envir = environment())
                                lib_data <- fiehnHilicLib
                              },
                              'zhuRPLib' = {
                                data('zhuRPlib', envir = environment())
                                lib_data <- zhuRPlib
                              }
                      )
                    },
                    "ThermoExploris" = {
                      # data('zhuMetlib_orbitrap', envir = environment())
                      # lib_data <- zhuMetlib_orbitrap
                      switch (lib,
                              'zhuMetLib_orbitrap' = {
                                data('zhuMetlib_orbitrap', envir = environment())
                                lib_data <- zhuMetlib_orbitrap
                              },
                              'fiehnHilicLib' = {
                                data('fiehnHilicLib', envir = environment())
                                lib_data <- fiehnHilicLib
                              },
                              'zhuRPLib' = {
                                data('zhuRPlib', envir = environment())
                                lib_data <- zhuRPlib
                              }
                      )
                    },
                    "AgilentDTIMMS" = {
                      data("zhuMetlib", envir = environment())
                      lib_data = zhuMetlib
                      ce <- as.character(as.numeric(ce) + 10)
                    },
                    'WatersTWIMMS' = {
                      data("zhuMetlib", envir = environment())
                      lib_data = zhuMetlib
                      ce <- as.character(as.numeric(ce) + 10)
                    },
                    "BrukerTIMS" = {
                      data("zhuMetlib", envir = environment())
                      lib_data = zhuMetlib
                    })


             temp_polarity <- ifelse(polarity == 'positive', 'pos', 'neg')

             # lib_meta <- zhuMetlib[['meta']][[temp_polarity]][[ce]]
             # lib_spec <- zhuMetlib[['compound']][[temp_polarity]]

             if (lib %in% c('zhuMetLib_orbitrap', 'zhuRPLib')) {
               if (ce %in% c("30", "10", "20", "40", "50")) {
                 ce <- paste('CE', ce, sep = '')
               }
               if (ce %in% c('SCE20_30_40%', "SNCE20_30_40%")) {
                 ce <- gsub('_', '-', ce)
                 ce <- gsub('%', '', ce)
               }

             }

             lib_meta <- lib_data[['meta']][[temp_polarity]][[ce]]
             lib_spec <- lib_data[['compound']][[temp_polarity]]

             switch(polarity,
                    'positive' = {
                      temp_adduct_table <- lib_adduct_nl$positive %>% dplyr::filter(adduct %in% adduct_list)
                    },
                    'negative' = {
                      temp_adduct_table <- lib_adduct_nl$negative %>% dplyr::filter(adduct %in% adduct_list)
                    }
             )

             # calculate mz
             lib_mz <- lapply(seq_along(lib_meta$labid), function(i){
               temp_mz <- lib_meta$mz[i]

               result <- sapply(seq_along(temp_adduct_table$adduct), function(j){
                 calculateMz(exact_mass = temp_mz,
                             adduct = temp_adduct_table$adduct[j],
                             delta_mz = temp_adduct_table$delta_mz[j])
               })

               result
             })

             lib_mz <- lib_mz %>% do.call(rbind, .) %>% tibble::as_tibble()
             colnames(lib_mz) <- temp_adduct_table$adduct

             # lib_meta
             lib_meta <- match(lib_meta$labid, lib_data$meta$compound$id) %>%
               lib_data$meta$compound[.,]

             lib_meta <- lib_meta %>%
               dplyr::bind_cols(lib_mz) %>%
               tidyr::pivot_longer(cols = temp_adduct_table$adduct,
                                   names_to = 'adduct',
                                   values_to = 'mz') %>%
               dplyr::select(id:formula, adduct:mz, dplyr::everything())


             # add rt and ccs
             switch(method_lc,
                    'Amide12min' = {
                      lib_rt_exp <- lib_rt[[1]] %>%
                        dplyr::filter(method == method_lc)
                    },
                    'Amide23min' = {
                      lib_rt_exp <- lib_rt[[1]] %>%
                        dplyr::filter(method == method_lc)
                    },
                    'MetlinRP' = {
                      lib_rt_exp <- lib_rt[[2]]
                    },
                    'Other' = {
                      lib_rt_exp <- lib_rt[[1]] %>%
                        dplyr::filter(method == method_lc)
                    },
                    'zhulabRP' = {
                      lib_rt_exp <- lib_rt[[3]]
                    })


             if (is_rt_calibration) {
               cat('Calibrate RTs in library...\n')
               temp_result <- calibrateRT(file_rt_ref = file_rt_ref,
                                          lib_rt = lib_rt_exp,
                                          is_rt_calibration = TRUE,
                                          is_plot = TRUE,
                                          method_lc = method_lc,
                                          column = 'hilic',
                                          path = path)

               lib_rt_exp <- temp_result$lib_rt
             }

             if (method_lc == 'zhulabRP') {
               temp_rt <- match(lib_meta$id, lib_rt_exp$id) %>%
                 lib_rt_exp$rt[.] %>%
                 round(., digits = 2)
             } else {
               temp_rt <- match(lib_meta$inchikey1, lib_rt_exp$inchikey1) %>%
                 lib_rt_exp$rt[.] %>%
                 round(., digits = 2)
             }

             lib_ccs_exp <- lib_ccs[[1]]
             temp_ccs <- match(paste(lib_meta$inchikey, lib_meta$adduct, sep = '_'),
                               paste(lib_ccs_exp$inchikey, lib_ccs_exp$adduct, sep = '_')) %>%
               lib_ccs_exp$ccs[.] %>%
               round(., digits = 2)

             # generate final lib_meta
             lib_meta <- lib_meta %>%
               dplyr::mutate(rt = temp_rt,
                             ccs = temp_ccs) %>%
               dplyr::select(id:mz, rt, ccs, dplyr::everything())

             # if (is_rt_score) {
             #   lib_meta <- lib_meta %>%
             #     dplyr::filter(!is.na(rt))
             # }
             #
             # if (is_ccs_score) {
             #   lib_meta <- lib_meta %>%
             #     dplyr::filter(!is.na(ccs))
             # }

             lib_spec <- match(unique(lib_meta$id), names(lib_spec)) %>%
               lib_spec[.] %>%
               lapply(., function(x){
                 x[[ce]]
               })

             result <- list(lib_meta = lib_meta,
                            lib_spec = lib_spec)


             return(result)
           }
)




#   calibrateRT ----------------------------------------------------------------

#' @title calibrateRT
#' @description RT calibration according to the RTQC
#' @author Zhiwei Zhou, Mingdu Luo
#' \email{zhouzw@@sioc.ac.cn}
#' @param file_rt_ref Default: 'RT_recalibration_table.csv'
#' @param lib_rt
#' @param is_rt_calibration
#' @param is_plot Default: TRUE
#' @param method_lc 'Amide12min', 'Amide23min', 'MetlinRP'. Default:  'Amide12min'
#' @param column Default: 'hilic'
#' @param path '.'
#' @export

# data("lib_rt", envir = environment())
# lib_rt <- lib_rt[[1]] %>% dplyr::filter(method == 'Amide23min')
# is_rt_calibration <- TRUE
# method_lc <- 'Amide23min'
# column <- 'hilic'
# path <- '/home/zhouzw/Data_processing/20210522_debug/POS/'
# file_rt_ref = 'RT_recalibration_table.csv'
# temp_result <- calibrateRT(file_rt_ref = 'RT_recalibration_table.csv',
#                            lib_rt = lib_rt,
#                            is_rt_calibration = TRUE,
#                            is_plot = TRUE,
#                            method_lc = 'Amide23min',
#                            column = 'hilic',
#                            path = path)
#

setGeneric(name = 'calibrateRT',
           def = function(file_rt_ref = 'RT_recalibration_table.csv',
                          lib_rt,
                          is_rt_calibration = TRUE,
                          is_plot = TRUE,
                          method_lc = c('Amide12min', 'Amide23min', 'MetlinRP', 'zhulabRP'),
                          column = 'hilic',
                          path = '.'){

             # Load in house RT library according to LC method
             switch(method_lc,
                    "Amide12min" = {
                      data('rt_ref', envir = environment())
                      ref_rtqc_table <- rt_ref[[1]]
                      lc_start <- 0
                      lc_end <- 720},

                    "Amide23min" = {
                      data('rt_ref', envir = environment())
                      ref_rtqc_table <- rt_ref[[2]]
                      lc_start <- 0
                      lc_end <- 1380},

                    "MetlinRP" = {
                      data('rt_ref', envir = environment())
                      ref_rtqc_table <- rt_ref[[3]]
                      lc_start <- 0
                      lc_end <- 1500},

                    'zhulabRP' = {
                      data('rt_ref', envir = environment())
                      ref_rtqc_table <- rt_ref[[4]]
                      lc_start <- 0
                      lc_end <- 720
                    }
             )

             # Check whether to do RT recalibration, only hilic method were allowed ###
             if (any(!is_rt_calibration, column != 'hilic')){
               cat('RT recalibration was turned off.\n')
             } else {

               ### Check existence of RT calibration table ###
               if (!("RT_recalibration_table.csv" %in% dir(path))){
                 stop("There was no 'RT_recalibration_table.csv' table in the file folder.\n")
               } else {

                 exp_rtqc_rt <- read.csv(file.path(path, file_rt_ref), stringsAsFactors = F)

                 ### Check the format of exp_rtqc_rt ###
                 if (!identical(colnames(exp_rtqc_rt)[1:6], c("compound.name", "id.zhulab", "id.pubchem", "ref.mz", "rt", "polarity"))){
                   stop("Please check the format of 'RT_recalibration_table.csv' table and correct it according to our tutorial.\n")
                 } else {

                   ### Begin RT calibration ###
                   if (max(exp_rtqc_rt$rt) < 60) {
                     exp_rtqc_rt$rt <- round(exp_rtqc_rt$rt*60, digits = 4)
                   } else {
                     exp_rtqc_rt$rt <- round(exp_rtqc_rt$rt, digits = 4)
                   }

                   # the index of match rt compound and the warning
                   idx <- match(toupper(ref_rtqc_table$name),
                                toupper(exp_rtqc_rt$id.zhulab))

                   if (sum(!is.na(idx)) < 7){
                     warning(paste0("The number of used compouds for RT recalibration with LOESS was ",
                                    sum(!is.na(idx)),
                                    " and might be insufficient for a good performance.\n\n"))
                   } else {
                     cat("The number of used compouds for RT recalibration with LOESS was ",
                         sum(!is.na(idx)),
                         ".\n\n",
                         sep ='')
                   }

                   training.data <- data.frame(ref.rt = c(lc_start, ref_rtqc_table$rt[!is.na(idx)], lc_end),
                                               exp.rt = c(lc_start, exp_rtqc_rt$rt[idx[!is.na(idx)]], lc_end))

                   rownames(training.data) <- c('Start', exp_rtqc_rt$id.zhulab[idx[!is.na(idx)]], 'End')

                   ### rt.recalibration.model <- lm(exp.rt~ref.rt, data = training.data)
                   rt.recalibration.model <- loess(exp.rt~ref.rt,
                                                   data = training.data,
                                                   span = 0.75, degree = 2)

                   new.data <- data.frame(ref.rt = lib_rt$rt, stringsAsFactors = FALSE)
                   lib.calibrated.rt <- round(predict(object =  rt.recalibration.model,
                                                      newdata = new.data),
                                              digits = 2)

                   lib_rt$rt <- lib.calibrated.rt

                   result <- list(lib_rt = lib_rt,
                                  training.data = training.data,
                                  rt.recalibration.model = rt.recalibration.model)

                   dir.create(file.path(path, "01_result_initial_seed_annotation", '00_intermediate_data'), showWarnings = FALSE, recursive = TRUE)
                   save(result,
                        file = file.path(path, "01_result_initial_seed_annotation", '00_intermediate_data', 'rt_calibration_result'),
                        version = 2)

                   cat("RT recalibration was done.\n")

                   if (is_plot) {
                     # browser()
                     dir.create(file.path(path, "01_result_initial_seed_annotation", '01_rt_calibration_plot'),
                                showWarnings = FALSE, recursive = TRUE)
                     plotRtCalibration(result_rt_calibration = result[[2]],
                                       rt_recalibration_model = result[[3]],
                                       path = file.path(path, "01_result_initial_seed_annotation", '01_rt_calibration_plot'))
                   }

                   return(result)
                   rm(c(idx, lc_start, lc_end))

                 }
               }
             }
           })





################################################################################
# matchMs1WithSpecLib ----------------------------------------------------------

#' @title matchMs1WithSpecLib
#' @author Zhiwei Zhou
#' @param ms1_data ms1 data, including info and subject.
#' @param mz_tol m/z tolerance. Default: 10ppm
#' @param pf_ccs_range penalty-free CCS range. Default: 0
#' @param tolerance_ccs_range maxmium CCS tolerance. Default: 4%
#' @param is_ccs_score logical vector. Default: TRUE
#' @param lib_meta the meta information of library
#' @export

# temp <- matchMs1WithSpecLib(ms1_data = ms1_data,
#                             mz_tol = 25,
#                             pf_rt_range = 0,
#                             tolerance_rt_range = 30,
#                             pf_ccs_range = 0,
#                             tolerance_ccs_range = 4,
#                             is_rt_score = FALSE,
#                             is_ccs_score = FALSE,
#                             lib_meta = lib_meta)
#
# temp <- matchMs1WithSpecLib(ms1_data = ms1_data,
#                             mz_tol = 25,
#                             pf_rt_range = 0,
#                             tolerance_rt_range = 30,
#                             pf_ccs_range = 0,
#                             tolerance_ccs_range = 4,
#                             is_rt_score = TRUE,
#                             is_ccs_score = FALSE,
#                             is_filter = TRUE,
#                             lib_meta = lib_meta)

setGeneric(name = 'matchMs1WithSpecLib',
           def = function(ms1_data,
                          lib_meta,
                          mz_tol=25,
                          mz_ppm_thr = 400,
                          pf_rt_range = 0,
                          tolerance_rt_range = 30,
                          pf_ccs_range = 0, # penalty-free CCS range-percentage
                          tolerance_ccs_range = 4, # percentage
                          is_filter = TRUE,
                          is_rt_score = FALSE,
                          is_ccs_score = TRUE,
                          is_check_cation = TRUE
           ){
             cat("m/z & RT & CCS match...\n\n")
             ms1_info <- ms1_data$info
             result_annotation <- convertMs1Data2SpecAnnotationClass(ms1_info = ms1_info)

             # browser()
             exp_mz_range <- getMzRange(mz = ms1_info$mz,
                                        ppm = mz_tol,
                                        mz_ppm_thr = mz_ppm_thr)

             idx_match <- lapply(seq_along(result_annotation), function(i){
               temp_idx <- which(lib_meta$mz >= exp_mz_range[i,1] & lib_meta$mz <= exp_mz_range[i,2])
               temp_idx
             })

             if (is_filter) {
               if (is_rt_score) {
                 # cat("rt match\n")
                 exp_rt_range <- getRtRange(data = ms1_info$rt,
                                            abs_rt_dev = tolerance_rt_range,
                                            rel_rt_dev = 6,
                                            rt_threshold = 500,
                                            is_combined = TRUE)

                 temp_idx <- lapply(seq_along(result_annotation), function(i){
                   temp_idx <- which(lib_meta$rt >= exp_rt_range[i,1] & lib_meta$rt <= exp_rt_range[i,2])
                   # include no experimental RT candidates
                   temp_idx2 <- which(is.na(lib_meta$rt))
                   result <- unique(c(temp_idx, temp_idx2))
                 })
                 idx_match <- mapply(function(x, y){
                   intersect(x,y)
                 },
                 x=idx_match,
                 y=temp_idx)
               }

               if (is_ccs_score) {
                 # cat("CCS match\n")
                 exp_ccs_range <- getCcsRange(data = ms1_info$ccs, rel_dev = tolerance_ccs_range)
                 temp_idx <- lapply(seq(nrow(ms1_info)), function(i){
                   temp_idx <- which(lib_meta$ccs >= exp_ccs_range[i,1] & lib_meta$ccs <= exp_ccs_range[i,2])
                   # include no experimental CCS candidates
                   temp_idx2 <- which(is.na(lib_meta$ccs))
                   result <- unique(c(temp_idx, temp_idx2))
                 })
                 idx_match <- mapply(function(x, y){
                   intersect(x,y)
                 },
                 x=idx_match,
                 y=temp_idx)
               }
             }


             # cat('RT & CCS scoring...\n')
             pbapply::pboptions(type = 'timer', char='+')
             result <- pbapply::pblapply(seq_along(result_annotation), function(i){
               # cat(i);cat(' ')
               if (length(idx_match[[i]]) > 0) {
                 temp_idx <- idx_match[[i]]
                 lib_mz <- lib_meta$mz[temp_idx]
                 mz_error <- round(abs(ms1_info$mz[i]-lib_mz)/lib_mz*10^6,
                                   digits = 1)


                 if (is_rt_score) {
                   temp_rt <- ms1_data$info$rt[i]
                   temp_lib_rt <- lib_meta$rt[temp_idx]
                   delta_rt <- abs(temp_rt-temp_lib_rt)
                   rt_score <- getTrapezoidalScore(delta = delta_rt,
                                                   pf_range = pf_rt_range,
                                                   range = tolerance_rt_range)
                   rt_score <- round(rt_score, digits = 2)
                   rt_error <- round(delta_rt, digits = 1)
                 } else {
                   rt_score <- -1
                   rt_error <- -1
                 }

                 if(is_ccs_score) {
                   temp_ccs <- ms1_data$info$ccs[i]
                   temp_lib_ccs <- lib_meta$ccs[temp_idx]
                   delta_ccs <- abs(temp_ccs-temp_lib_ccs)/temp_lib_ccs*100
                   ccs_score <- getTrapezoidalScore(delta = delta_ccs,
                                                    pf_range = pf_ccs_range,
                                                    range = tolerance_ccs_range)
                   ccs_score <- round(ccs_score, digits = 2)
                   ccs_error <- round(delta_ccs, digits = 1)
                 } else {
                   ccs_score <- -1
                   ccs_error <- -1
                 }


                 if (is_check_cation) {
                   idx_check_cation <- mapply(function(x, y){
                     if (x %in% c('[M]+', '[M-2H]-')) {
                       !is.na(y)
                     } else {
                       return(TRUE)
                     }
                   },
                   x = lib_meta$adduct[temp_idx],
                   y = lib_meta$note[temp_idx]) %>%
                     which()

                   if (length(idx_check_cation) == 0) {
                     return(NULL)
                   } else {
                     temp_idx <- temp_idx[idx_check_cation]
                   }

                   result <- tibble::tibble(idx = temp_idx,
                                            id = lib_meta$id[temp_idx],
                                            name = lib_meta$name[temp_idx],
                                            formula = lib_meta$formula[temp_idx],
                                            smiles = lib_meta$smiles[temp_idx],
                                            inchikey = lib_meta$inchikey[temp_idx],
                                            inchikey1 = lib_meta$inchikey1[temp_idx],
                                            adduct = lib_meta$adduct[temp_idx],
                                            mz_lib = lib_meta$mz[temp_idx],
                                            rt_lib = lib_meta$rt[temp_idx],
                                            ccs_lib = lib_meta$ccs[temp_idx],
                                            mz_error = mz_error[idx_check_cation],
                                            # mz_score=mz_score,
                                            rt_error = rt_error[idx_check_cation],
                                            rt_score = rt_score[idx_check_cation],
                                            ccs_error = ccs_error[idx_check_cation],
                                            ccs_score = ccs_score[idx_check_cation],
                                            msms_score_forward = -1,
                                            msms_score_reverse = -1)

                   return(result)
                 }

                 result <- tibble::tibble(idx = temp_idx,
                                          id = lib_meta$id[temp_idx],
                                          name = lib_meta$name[temp_idx],
                                          formula = lib_meta$formula[temp_idx],
                                          smiles = lib_meta$smiles[temp_idx],
                                          inchikey = lib_meta$inchikey[temp_idx],
                                          inchikey1 = lib_meta$inchikey1[temp_idx],
                                          adduct = lib_meta$adduct[temp_idx],
                                          mz_lib = lib_meta$mz[temp_idx],
                                          rt_lib = lib_meta$rt[temp_idx],
                                          ccs_lib = lib_meta$ccs[temp_idx],
                                          mz_error = mz_error,
                                          # mz_score=mz_score,
                                          rt_error = rt_error,
                                          rt_score = rt_score,
                                          ccs_error = ccs_error,
                                          ccs_score = ccs_score,
                                          msms_score_forward = -1,
                                          msms_score_reverse = -1)

                 return(result)
               }

             })


             result_annotation <- mapply(function(x, y){
               if (length(y)>0) {
                 x@annotation_result <- y
                 return(x)
               } else {
                 return(x)
               }
             },
             x = result_annotation,
             y = result,
             SIMPLIFY = FALSE)

             return(result_annotation)

           }
)



#' @title getRtRange
#' @author Zhiwei Zhou
#' @param data Numeric.
#' @param abs_rt_dev Numeric. Absolute rt deviation tolerance in seconds.
#' @param rel_rt_dev Numeric. Relative rt deviation tolerance in percent %.
#' @return a table of range. column 1: minimum value; column 2: maxmium value.

# getRtRange(data = 300, abs_rt_dev = 30)
# getRtRange(data = 600, rel_rt_dev = 6)
# getRtRange(data = 600, abs_rt_dev = 30, rel_rt_dev = 6, is_combined = TRUE)
# getRtRange(data = 300, abs_rt_dev = 30, rel_rt_dev = 6, is_combined = TRUE)

setGeneric(name = 'getRtRange',
           def = function(data,
                          abs_rt_dev,
                          rel_rt_dev,
                          rt_threshold = 500,
                          is_combined = FALSE){

             if (is_combined) {
               if (!any(missing(abs_rt_dev), missing(rel_rt_dev))) {
                 rt_range <- t(sapply(data, function(x){
                   if (x <= rt_threshold) {
                     result <- x+c(-1, 1)*abs_rt_dev
                   } else {
                     result <- x*(1+c(-1, 1)*rel_rt_dev/100)
                   }

                   result

                 }))
               }

             } else {
               if (!missing(abs_rt_dev)) {
                 rt_range <- t(sapply(data, function(x){
                   x+c(-1, 1)*abs_rt_dev
                 }))
               }

               if (!missing(rel_rt_dev)) {
                 rt_range <- t(sapply(data, function(x){
                   x*(1+c(-1, 1)*rel_rt_dev/100)
                 }))
               }
             }



             return(rt_range)

           })



#' @title getCcsRange
#' @author Zhiwei Zhou
#' @param data Numeric.
#' @param rel_dev Numeric. Relative CCS deviation tolerance in percent %.
#' @return a table of range. column 1: minimum value; column 2: maxmium value.
#' @export

setGeneric(name = 'getCcsRange',
           def = function(data, rel_dev){
             temp <- t(sapply(data, function(x){
               x*(1+c(-1, 1)*rel_dev/100)
             }))
             return(temp)
           })



#' @title getTrapezoidalScore
#' @author Zhiwei Zhou
#' @description Trapezodial scoring function
#' @param delta Numeric. The difference between experiment value and library value.
#' @param pf_range Numeric. Penalty free range.
#' @param range Numeric. The maxmium tolerance.
#' @export

setGeneric('getTrapezoidalScore',
           def = function(delta,
                          pf_range,
                          range){
             delta <- abs(delta)
             delta[delta <= pf_range] <- pf_range
             score <- 1-((delta-pf_range)/(range-pf_range))
             score[score<0] <- 0
             return(score)
           })


#' @title getLinerScore
#' @author Zhiwei Zhou
#' @description calculate liner score
#' @param delta Numeric. The difference between experiment value and library value.
#' @param tolerance The maxmium tolerance
#' @export

setGeneric(name = 'getLinerScore',
           def = function(delta,
                          tolerance){
             score <- 1 - abs(delta) / tolerance
             score
           })



################################################################################
# inteMs2 ----------------------------------------------------------------------

#' @title inteMs2
#' @description Purify and integrate multiple MS/MS files
#' @author Zhiwei Zhou
#' \email{zhouzw@@sioc.ac.cn}
#' @param ms2_data ms2 data
#' @param is_include_precursor remove precursor peaks in raw spectra or not. Default: FALSE
#' @param is_deisotope remove isotope peaks in raw spectra or not. Default: TRUE
#' @param int_ms2_min_abs the minmium intensity of ms2 fragment. Default: 10
#' @param int_ms2_min_relative the minmium relative of ms2 fragment. Default: 0.01 ((1%))
#' @param ppm_precursor_filter the m/z tolerance to rexmove precursor ion in spectra. Default: 20
#' @param extract_ms2_file extract ms2 file. Default: NULL
#' @param mz_range_ms2 the range of ms2 data.
#' @export

# path <- 'H:/00_projects/03_MetDNA2/00_data/20200922_metdna2_development/pos_201007_metdna2/'
# temp_files <- dir(path) %>% stringr::str_detect('\\.mgf') %>% dir(path)[.]
# ms2_data <- read(ms2_file = file.path(path, temp_files),
#                     ms2_type = 'mgf')
# test <- inteMs2(ms2_data = ms2_data,
#                 ms2_type = 'mgf',
#                 is_include_precursor = TRUE,
#                 is_deisotope = FALSE,
#                 int_ms2_min_abs = 0,
#                 int_ms2_min_relative = 0.01,
#                 ppm_precursor_filter = 25,
#                 mz_range_ms2 = NULL)

setGeneric(name = 'inteMs2',
           def = function(
             ms2_data,
             ms2_type = c("mgf", "msp", "cef", 'mzXML'),
             instrument = c("SciexTripleTOF", "AgilentQTOF", "BrukerQTOF", "ThermoOrbitrap", "ThermoExploris",
                            'WatersQTOF','WatersTWIMMS', "AgilentDTIMMS", "BrukerTIMS"),
             metdna_version = c('version1', 'version2'),
             is_include_precursor = TRUE,
             is_deisotope = FALSE,
             int_ms2_min_abs = 0,
             int_ms2_min_relative = 0.01,
             ppm_precursor_filter = 20,
             mz_range_ms2=NULL,
             ...
           ){
             match.arg(ms2_type)

             cat('Purify and integrate MS/MS spectra\n')
             tune_MS2_spec <- pbapply::pblapply(seq_along(ms2_data), function(i){
               # cat(i, ' ')

               # to keep same with MetDNA1

               if (metdna_version == 'version1') {
                 result <- purifyMs2(spec = ms2_data[[i]]$spec,
                                     # mz_precursor = ms2_data[[i]]$info['mz'],
                                     is_include_precursor = TRUE,
                                     is_deisotope = is_deisotope,
                                     int_ms2_min_abs = int_ms2_min_abs,
                                     int_ms2_min_relative = int_ms2_min_relative,
                                     ppm_precursor_filter = ppm_precursor_filter,
                                     mz_range_ms2 = mz_range_ms2)
               } else {
                 result <- purifyMs2(spec = ms2_data[[i]]$spec,
                                     mz_precursor = ms2_data[[i]]$info['mz'],
                                     is_include_precursor = is_include_precursor,
                                     is_deisotope = is_deisotope,
                                     int_ms2_min_abs = int_ms2_min_abs,
                                     int_ms2_min_relative = int_ms2_min_relative,
                                     ppm_precursor_filter = ppm_precursor_filter,
                                     mz_range_ms2 = mz_range_ms2)
               }


               return(result)
             })

             info <- lapply(seq(length(ms2_data)), function(i){
               ms2_data[[i]]$info
             })

             info <- do.call(rbind, info)

             if (ms2_type == "mgf") {

               if (!(instrument %in% c("AgilentDTIMMS", "BrukerTIMS", 'WatersTWIMS'))) {
                 info <- data.frame(mz=as.numeric(info[,1]),
                                    rt=as.numeric(info[,2]),
                                    stringsAsFactors = F)
               } else {
                 info <- data.frame(mz=as.numeric(info[,1]),
                                    rt=as.numeric(info[,2]),
                                    ccs=as.numeric(info[,3]),
                                    stringsAsFactors = F)
               }

             }

             if (ms2_type == "cef") {
               info <- data.frame(mz=as.numeric(info[,1]),
                                  rt=as.numeric(info[,2]),
                                  ccs=as.numeric(info[,3]),
                                  stringsAsFactors = F)
             }

             if (ms2_type == "msp") {
               temp <- rownames(info)
               info <- data.frame(name=as.character(info[,1]),
                                  mz=as.numeric(info[,2]),
                                  stringsAsFactors = F)
               rownames(info) <- temp
             }

             if (ms2_type == 'mzXML') {
               temp <- rownames(info)
               info <- data.frame(mz=as.numeric(info[,1]),
                                  rt=as.numeric(info[,2]),
                                  stringsAsFactors = F)
             }

             idx_null <- which(sapply(tune_MS2_spec, is.null))

             if (length(idx_null) > 0) {
               info <- info[-idx_null,,drop=F]
               tune_MS2_spec <- tune_MS2_spec[-idx_null]
             }


             result <- list(info=info, spec=tune_MS2_spec)
             return(result)

           }
)


#   purifyMs2 -------------------------------------------------------------------


#' @title purifyMs2
#' @author Zhiwei Zhou
#' @param spec the spectra matrix. column1: mz; column2: intensity
#' @param mz_precursor the precursor m/z. Numeric.
#' @param is_include_precursor remove precursor peak in raw spectra or not_ Default: False
#' @param mz_range_ms2 the range of m/z. Default: NULL
#' @param int_ms2_min_abs the minmium intensity of ms2 fragment_ Default: 10
#' @param int_ms2_min_relative the minmium relative of ms2 fragment_ Default: 0.03 ((3%))
#' @param ppm_precursor_filter the m/z tolerance to remove precursor ion in spectra_ Default: 20
#' @export

setGeneric(name = 'purifyMs2',
           def = function(
             spec,
             mz_precursor,
             is_include_precursor = FALSE,
             is_deisotope = FALSE,
             is_remove_ring_effect = TRUE,
             mz_range_ms2 = NULL,
             int_ms2_min_abs = 10,
             int_ms2_min_relative = 0.03,
             ppm_precursor_filter = 20
           ){
             # browser()
             spec <- spec[order(spec[, 'mz']), , drop = FALSE]

             # considering precursor ion
             if (missing(mz_precursor)) {
               mz_precursor <- max(spec[, 'mz'])
               mz_precursor <- as.numeric(mz_precursor)
             }

             mz_precursor_range <- getMzRange(mz_precursor, ppm_precursor_filter, mz_ppm_thr = 400)
             idx_mz_precursor_range <- ifelse(is_include_precursor, 2, 1)

             #change mz range depend precusor include or not
             mz_cutoff <- mz_precursor_range[idx_mz_precursor_range]
             spec <- spec[spec[,'mz'] < mz_cutoff, , drop = FALSE]

             if (!is.null(mz_range_ms2)) {
               nr_keep <- which(spec[, 'mz'] >= mz_range_ms2[1] &
                                  spec[, 'mz'] <= mz_range_ms2[2])
               if (length(nr_keep) > 0) {
                 spec <- spec[nr_keep, , drop = FALSE]
               }
               else {
                 return()
               }
             }

             # discarding low intensity spec (1% highest int and int.ms2.min.abs)
             int_cutoff <- max(max(spec[, 'intensity']) *
                                 int_ms2_min_relative,
                               int_ms2_min_abs)
             spec <- spec[spec[, 'intensity'] >= int_cutoff, , drop = FALSE]
             if (nrow(spec) == 0) {
               return()
             }

             if (is_deisotope) {
               spec <- removeIsotopes(spec)
             }

             # discarding ring effects
             if (is_remove_ring_effect) {
               spec <- removeRingEffect(spec)
             }

             return(spec)

           }
)



#' @title removeIsotopes
#' @author Zhiwei Zhou
#' \email{zhouzw@@sioc.ac.cn}
#' @param spec the matrix spectra. 1st column: "mz", 2nd column: "intensity"

setGeneric('removeIsotopes',
           def = function(
             spec
           ){
             # transform spec to data.frame
             temp_spec <- data.frame(spec, annotation='M', stringsAsFactors = F)

             is_filter_isotope <- TRUE

             while (is_filter_isotope) {
               idx_max <- which.max(temp_spec[,'intensity'])
               # cat(idx.max); cat(' ')

               mono_mz <- temp_spec[idx_max, 'mz']
               mono_int <- temp_spec[idx_max, 'intensity']

               temp_spec$intensity[idx_max] <- 0

               # generate potential isotope list
               isotopes_list <- mono_mz + c(1.0003, 2.0044)

               isotopes_list <- lapply(isotopes_list, function(x){
                 c(x, x + c(-1, 1) * 0.003)
               })

               isotopes_list <- do.call(rbind, isotopes_list)
               colnames(isotopes_list) <- c('mz', 'min', 'max')

               # m/z match
               # if multiple spectra were

               temp_idx_isotope <- lapply(seq(nrow(isotopes_list)), function(i){
                 idx_isotope <- which(temp_spec[,'mz'] >= isotopes_list[i,'min'] & temp_spec[,'mz'] <= isotopes_list[i,'max'])

                 if (length(idx_isotope) > 1) {
                   temp_mz_error <- abs(temp_spec[idx_isotope, 'mz'] - isotopes_list[i,'mz'])
                   idx_isotope <- idx_isotope[which.min(temp_mz_error)]
                 }

                 idx_isotope
               })


               # judge the intensity relationship
               # if matched M+1 and M+2, intensity relationship: M > M+1 > M+2
               # if matched M+1, intensity relationship: M > M+1

               if (length(temp_idx_isotope[[1]]) > 0 & length(temp_idx_isotope[[2]]) > 0) {
                 temp_m_1_int <- temp_spec$intensity[temp_idx_isotope[[1]]]
                 temp_m_2_int <- temp_spec$intensity[temp_idx_isotope[[2]]]

                 # intensity decrease
                 if (mono_int > temp_m_1_int & temp_m_1_int > temp_m_2_int) {
                   idx_isotope <- c(temp_idx_isotope[[1]], temp_idx_isotope[[2]])
                   temp_spec$annotation[idx_isotope] <- c('M+1', 'M+2')
                   temp_spec$intensity[idx_isotope] <- 0

                   # names(idx_isotope) <- c('M+1', 'M+2')
                 }

                 if (mono_int > temp_m_1_int & temp_m_1_int <= temp_m_2_int) {
                   idx_isotope <- c(temp_idx_isotope[[1]])
                   temp_spec$annotation[idx_isotope] <- c('M+1')
                   temp_spec$intensity[idx_isotope] <- 0
                   # names(idx_isotope) <- c('M+1', 'M+2')
                 }

               }

               if (length(temp_idx_isotope[[1]]) > 0 & length(temp_idx_isotope[[2]]) == 0) {
                 temp_m_1_int <- temp_spec$intensity[temp_idx_isotope[[1]]]

                 if (mono_int > temp_m_1_int) {
                   idx_isotope <- c(temp_idx_isotope[[1]])
                   temp_spec$annotation[idx_isotope] <- c('M+1')
                   temp_spec$intensity[idx_isotope] <- 0
                 }

               }


               if (all(temp_spec$intensity==0)) {
                 is_filter_isotope <- FALSE
               }

             }

             idx_remove <- which(temp_spec$annotation != 'M')
             spec <- spec[-idx_remove, , drop=F]
             spec
           }
)




#' @title removeRingEffect
#' @author Yandong Yin; Zhiwei Zhou
#' \email{zhouzw@@sioc.ac.cn}
#' @param spec 2-column table. column 1: mz; column 2: intensity
#' @param mz_diff_thr the threold of mz difference. Default: 0.3
#' @param int_rel_thr the threold of relative intensity. Default: 0.2
#' @export

setGeneric(name = 'removeRingEffect',
           def = function(spec,
                          mz_diff_thr = 0.3,
                          int_rel_thr = 0.2){
             nr_ring <- nrow(spec) + 1
             mz <- spec[, 'mz']

             mz_diff <- diff(mz)
             idx_mzdiff <- which(mz_diff <= mz_diff_thr)
             if (length(idx_mzdiff) == 0) {
               return(spec)
             }

             nr_ring_possible <- unique(c(idx_mzdiff, idx_mzdiff + 1))

             # remove ringeffect loop
             while (TRUE) {

               idx_int_max <- which.max(spec[nr_ring_possible, 2])
               nr_int_max <- nr_ring_possible[idx_int_max] # the index of possible Ringeffect ions with maxium intensity
               int_thr <- spec[nr_int_max, 2] * int_rel_thr # the threshold = 0_2*max_int (possible ring)

               mz_diff <- abs(mz[nr_ring_possible[-idx_int_max]] - mz[nr_int_max])
               int <- spec[nr_ring_possible[-idx_int_max], 2]
               nr_ring <- append(nr_ring, nr_ring_possible[-idx_int_max][which(mz_diff <= mz_diff_thr & int <= int_thr)])
               nr_ring_possible <- nr_ring_possible[!nr_ring_possible %in% c(nr_ring, nr_int_max)]
               if (length(nr_ring_possible) == 0) {
                 break # break loop untill satisfy the nr_ring_possible==0
               }
             }

             return(spec[-nr_ring, , drop = FALSE])
           }
)




# combineMS1MS2 ----------------------------------------------------------------
#' @title combineMS1MS2
#' @description Combine MS1 and MS2 data.
#' @author Xiaotao Shen, Zhiwei Zhou
#' \email{zhouzw@@sioc.ac.cn}
#' @param ms1.file The name of ms1 peak table. Column 1 is "name", Column 2 is "mz" and column 3 is "rt".
#' @param ms2.file The vector of names of ms2 files. MS2 file must be mzXML.
#' @param mz.tol mz tol for ms1 and ms2 data matching.
#' @param rt.tol RT tol for ms1 and ms2 data matching.
#' @param path Directory.
#' @param ms2.type The type of MS2 file, default is mzXML.
#' @return Return ms1 and ms2 data.

# combineMS1MS2(ms1.file = ms1.file, ms2.file = ms2.file)

setGeneric(name = "combineMs1Ms2",
           function(ms1_data,
                    ms2_data,
                    mz_tol_combine_ms1_ms2 = 25,
                    rt_tol_combine_ms1_ms2 = 10,
                    ccs_tol_combine_ms1_ms2 = NULL,
                    path = ".",
                    ms2_type = c("mgf", "mzXML", "msp", 'cef')){

             # browser()
             ms2_info <- ms2_data$info
             ms2_spec <- ms2_data$spec

             if (ms2_type %in% c("mzXML", 'mgf', 'cef')) {

               cat("\n")
               cat("Match MS1 and MS2 data.(mz tolerance ",
                   mz_tol_combine_ms1_ms2," ppm, RT tolerance ", rt_tol_combine_ms1_ms2, " second).\n", sep = "")

               if (length(ccs_tol_combine_ms1_ms2) == 0) {
                 temp_ms1_info <- ms1_data$info %>% dplyr::select(mz, rt)
               } else {
                 temp_ms1_info <- ms1_data$info %>% dplyr::select(mz, rt, ccs)
               }

               match_result <- SXTMTmatch(data1 = ms2_info,
                                          data2 = temp_ms1_info,
                                          mz.tol = mz_tol_combine_ms1_ms2,
                                          rt.tol = rt_tol_combine_ms1_ms2,
                                          ccs.tol = ccs_tol_combine_ms1_ms2,
                                          rt.error.type = "abs")

               # 20210322 debug: one ms2 matched multiple ms1, the ms2 was assigned the ms1 name of last index
               # ms2_name <- rep(NA, nrow(ms2_info))
               # ms2_name[match_result[,"Index1"]] <- ms1_data$info$name[match_result[,"Index2"]]
               # ms2_name <- as.character(ms2_name)

               match_result <- match_result %>%
                 tibble::as_tibble() %>% dplyr::mutate(peak_name = ms1_data$info$name[Index2])

               # # remove no matched MS2
               # remove_idx <- (!(ms2_name %in% match_result$Index1)) %>% which()
               # # remove_idx <- which(is.na(ms2_name))
               # if(length(remove_idx) > 0){
               #   ms2_name <- ms2_name[-remove_idx]
               #   ms2_info <- ms2_info[-remove_idx,]
               #   ms2_spec <- ms2_spec[-remove_idx]
               # }

               # remove duplicated MS2 spectrum, if one peak has more than 1 ms2 spectrum,
               #    select the biggest of the sum intensity of top 10 fragments.
               # unique_name <- unique(ms2_name)
               unique_name <- match_result$peak_name %>% unique()
               cat("\n")
               cat("Select the most intense MS2 spectrum for one peak.\n")
               remain_idx <- pbapply::pblapply(seq_along(unique_name), function(i){
                 name <- unique_name[i]
                 temp_idx <- match_result %>% dplyr::filter(peak_name == name) %>% dplyr::pull(Index1)
                 # temp_idx <- which(ms2_name == name)
                 if(length(temp_idx) == 1) return(temp_idx)
                 temp_ms2 <- ms2_spec[temp_idx]
                 # temp_ms2 <- lapply(temp_ms2, function(x) x[[2]])
                 temp_int <- lapply(temp_ms2, function(x) {
                   x <- as.data.frame(x)
                   x <- x[order(x[,2], decreasing = TRUE),]
                   if(nrow(x) >= 10) {sum(x[1:10,2])}else{sum(x[,2])}
                 })
                 temp_int <- unlist(temp_int)
                 temp_idx <- temp_idx[which.max(temp_int)]
                 temp_idx
               })

               names(remain_idx) <- unique_name
               remain_idx <- unlist(remain_idx)

               ms2_name <- names(remain_idx)
               ms2_info <- ms2_info[remain_idx,]
               ms2_spec <- ms2_spec[remain_idx]

               ms2 <- pbapply::pblapply(seq_along(ms2_name), function(i){
                 NAME <- ms2_name[i]
                 PRECURSORMZ <- ms2_info$mz[i]
                 PRECURSORRT <- ms2_info$rt[i]
                 result_info <- data.frame(NAME, PRECURSORMZ,
                                           PRECURSORRT, stringsAsFactors = FALSE)
                 result_info <- t(result_info)
                 rownames(result_info) <- c('NAME', 'PRECURSORMZ', 'PRECURSORRT')

                 if (length(ccs_tol_combine_ms1_ms2) > 0) {
                   PRECURSORCCS <- ms2_info$ccs[i]
                   result_info <- data.frame(NAME,
                                             PRECURSORMZ,
                                             PRECURSORRT,
                                             PRECURSORCCS,
                                             stringsAsFactors = FALSE)
                   result_info <- t(result_info)
                   rownames(result_info) <- c('NAME', 'PRECURSORMZ', 'PRECURSORRT', 'PRECURSORCCS')
                 }

                 result <- list(info = result_info,
                                spec = ms2_spec[[i]])

                 return(result)
               })

               names(ms2) <- ms2_name
               save(ms2, file = file.path(path, "ms2"), compress = "xz", version = 2)

               return(ms2)

             }


             if(ms2_type == "msp"){
               # ms2 <- readMSP(file = file.path(path, ms2_data), mode = 'all')
               # ms2_name <- sapply(ms2, function(x){x[[1]]$NAME})

               cat("\n")
               cat("Match MS1 and MS2 data according to the name.\n")

               ms1_info <- ms1_data$info

               # remove spectra not in peak table
               idx <- match(ms2_info$name, ms1_info$name) %>% is.na() %>% which()

               if (length(idx) > 0) {
                 ms2_info <- ms2_info[-idx,,drop = FALSE]
                 ms2_spec <- ms2_spec[-idx]
               }

               idx <- match(ms2_info$name, ms1_info$name)
               ms2_name <- ms2_info$name
               ms2_mz <- ms1_info$mz[idx]
               ms2_rt <- ms1_info$rt[idx]

               # change ms2 info to MetDNA1 default format
               ms2 <- lapply(seq_along(ms2_name), function(i){
                 temp_info <- data.frame('NAME' = ms2_name[i],
                                         'PRECURSORMZ' = ms2_mz[i],
                                         'PRECURSORRT' = ms2_rt[i],
                                         stringsAsFactors = FALSE)
                 temp_info <- t(temp_info)
                 rownames(temp_info) <- c('NAME', 'PRECURSORMZ', 'PRECURSORRT')
                 result <- list(info = temp_info,
                                spec = ms2_spec[[i]])

                 return(result)
               })

               names(ms2) <- ms2_name
               save(ms2, file = file.path(path, "ms2"), compress = "xz", version = 2)

               return(ms2)
             }


           })



################################################################################
# matchMs2WithSpecLib ----------------------------------------------------------
#' @title matchMs2WithSpecLib
#' @author Zhiwei Zhou
#' @param ms1_result ms1 result
#' @param exp_spec correspond msms data
#' @param lib_meta library of meta information
#' @param lib_spec library of spectrum
#' @param mz_tol_ms2 Default: 35 ppm
#' @param dp_cutoff Default: 0.8
#' @param matched_frag_cutoff Default: 1
#' @param direction 'forward', 'reverse'
#' @param scoring_approach 'dp', 'bonanza', 'hybrid', 'gnps'. Default: 'dp'
#' @export

# test <- matchMs2WithSpecLib(ms1_result = ms1_result,
#                             exp_spec = exp_spec,
#                             lib_meta = lib_meta,
#                             lib_spec = lib_spec,
#                             mz_tol_ms2 = 35,
#                             dp_cutoff = 0.8,
#                             direction = 'reverse',
#                             path = path)

setGeneric(name = 'matchMs2WithSpecLib',
           def = function(
             ms1_result,
             exp_spec,
             lib_meta,
             lib_spec,
             metdna_version = c('version2', 'version1'),
             mz_tol_ms2 = 35,
             dp_cutoff = 0.8,
             matched_frag_cutoff = 1,
             path = '.',
             is_include_precursor = TRUE,
             direction = c('reverse', 'forward'),
             scoring_approach = c('dp', 'bonanza', 'hybrid', 'gnps'),
             ...
           ){
             cat('Match ms2 with Spec library...\n')

             # browser()
             metdna_version <- match.arg(metdna_version)
             direction <- match.arg(direction)
             scoring_approach <- match.arg(scoring_approach)


             pbapply::pboptions(type='timer', char='+')
             ms2_result <- pbapply::pblapply(seq_along(ms1_result), function(i){
               # browser()
               # cat(i);cat(' ')
               # i <- 10209
               # i <- 41
               # i <- 277

               # modify the experimental ms2
               if (length(exp_spec[[i]]) > 0) {
                 if (nrow(ms1_result[[i]]@annotation_result) > 0) {
                   # modify the msms format for SpectraTools


                   if (metdna_version == 'version1') {
                     temp_peak_name <- ms1_result[[i]]@peak_info$name
                     temp_peak_mz <- ms1_result[[i]]@peak_info$mz
                     temp_peak_spec <- exp_spec[[i]]

                     # add absolute filter to keep same with metdna1
                     temp_peak_spec <- purifyMs2(spec = temp_peak_spec,
                                                 mz_precursor = temp_peak_mz,
                                                 is_include_precursor = is_include_precursor,
                                                 is_deisotope = FALSE,
                                                 mz_range_ms2 = NULL,
                                                 int_ms2_min_abs = 30,
                                                 int_ms2_min_relative = 0.01,
                                                 ppm_precursor_filter = 20)

                     if (length(temp_peak_spec) == 0) {return(NULL)}

                     rownames(temp_peak_spec) <- NULL

                     temp_peak_data <- list(info = tibble::tibble(NAME = temp_peak_name,
                                                                  PRECURSORMZ = temp_peak_mz),
                                            spec = temp_peak_spec)

                     exp_ms2 <- convertSpectraData(ms2_data = temp_peak_data)
                   } else {
                     temp_peak_name <- ms1_result[[i]]@peak_info$name
                     temp_peak_mz <- ms1_result[[i]]@peak_info$mz
                     temp_peak_spec <- exp_spec[[i]]

                     temp_peak_data <- list(info = tibble::tibble(NAME = temp_peak_name,
                                                                  PRECURSORMZ = temp_peak_mz),
                                            spec = temp_peak_spec)
                     exp_ms2 <- convertSpectraData(ms2_data = temp_peak_data)
                   }



                   # modify the library ms2
                   if (metdna_version == 'version1') {

                     # modify the candidate spec for SpectraTools
                     temp_lib_info <- ms1_result[[i]]@annotation_result %>%
                       dplyr::select(id, mz_lib) %>%
                       dplyr::rename(name = id,
                                     mz = mz_lib)

                     # metdna1 use the monoisotope mass as the precursor mz
                     temp_lib_info$mz <- match(temp_lib_info$name, lib_meta$id) %>% lib_meta$monoisotopic_mass[.]
                     temp_lib_spec <- match(temp_lib_info$name, names(lib_spec)) %>% lib_spec[.]


                     temp_lib_spec <- lapply(seq_along(temp_lib_spec), function(j){
                       temp <- purifyMs2(spec = temp_lib_spec[[j]],
                                         mz_precursor = temp_lib_info$mz[j],
                                         is_include_precursor = is_include_precursor,
                                         is_deisotope = FALSE,
                                         is_remove_ring_effect = FALSE,
                                         mz_range_ms2 = NULL,
                                         int_ms2_min_abs = 0,
                                         int_ms2_min_relative = 0,
                                         ppm_precursor_filter = 20)

                       rownames(temp) <- NULL

                       temp
                     })

                     names(temp_lib_spec) <- temp_lib_info$name

                     lib_ms2 <- new('SpectraData',
                                    info = temp_lib_info,
                                    spectra = temp_lib_spec)

                   } else {
                     # modify the candidate spec for SpectraTools
                     temp_lib_info <- ms1_result[[i]]@annotation_result %>%
                       dplyr::select(id, mz_lib) %>%
                       dplyr::rename(name = id,
                                     mz = mz_lib)

                     temp_lib_spec <- match(temp_lib_info$name, names(lib_spec)) %>% lib_spec[.]

                     lib_ms2 <- new('SpectraData',
                                    info = temp_lib_info,
                                    spectra = temp_lib_spec)

                   }


                   # # set parameters for match
                   # matchParam <- SpectraTools::MatchParam(ppm = mz_tol_ms2,
                   #                                        methodScore = 'dp',
                   #                                        methodMatch = 'direct',
                   #                                        weightIntensity = 1,
                   #                                        weightMZ = 0,
                   #                                        cutoff = 0,
                   #                                        includePrecursor = is_include_precursor,
                   #                                        intensityExpNormed = TRUE,
                   #                                        intensityLibNormed = TRUE,
                   #                                        tuneLibSpectra = FALSE)
                   #
                   # result <- try(SpectraTools::MatchSpectra(exp_ms2,
                   #                                          lib_ms2,
                   #                                          matchParam),
                   #               silent = TRUE)


                   result <- try(runSpecMatch(obj_ms2_cpd1 = exp_ms2,
                                              obj_ms2_cpd2 = lib_ms2,
                                              mz_tol_ms2 = mz_tol_ms2,
                                              scoring_approach = scoring_approach),
                                 silent = TRUE)


                   if ((class(result) == 'try-error') | (length(result) == 0)) {
                     return(NULL)
                   }

                   return(result)
                 } else {
                   return(NULL)
                 }
               } else {
                 return(NULL)
               }
             })

             names(ms2_result) <- names(ms1_result)

             dir.create(file.path(path, "00_intermediate_data"),
                        recursive = TRUE, showWarnings = FALSE)

             save(ms2_result,
                  file = file.path(path, "00_intermediate_data", 'ms2_result'),
                  compress = 'gzip',
                  version = 2)

             # for DP score, record ms2 score as its direction
             # for other score, record ms2 score as reverse
             cat('Add the ms2 score into SpecAnnotationClass...\n')
             ms2_result_annotation <- mapply(function(x, y){
               if (nrow(x@annotation_result) > 0) {
                 if (scoring_approach == 'dp') {
                   if(length(y) > 0) {
                     temp <- y@info
                     temp_forward <- temp$scoreForward
                     temp_reverse <- temp$scoreReverse
                     temp_n_frag <- temp$n_frag_total
                   } else {
                     temp_forward <- temp_reverse <- temp_n_frag <- 0
                   }
                   x@annotation_result$msms_score_forward <- temp_forward
                   x@annotation_result$msms_score_reverse <- temp_reverse
                   x@annotation_result$msms_matched_frag <- temp_n_frag
                   return(x)

                 } else {
                   if(length(y) > 0) {

                     temp <- y@info
                     temp_forward <- 0
                     temp_reverse <- temp$score
                     temp_n_frag <- temp$n_frag_total
                   } else {
                     temp_forward <- temp_reverse <- temp_n_frag <- 0
                   }
                   x@annotation_result$msms_score_forward <- temp_forward
                   x@annotation_result$msms_score_reverse <- temp_reverse
                   x@annotation_result$msms_matched_frag <- temp_n_frag
                   return(x)
                 }
               }

               return(x)

             },
             x = ms1_result,
             y = ms2_result,
             SIMPLIFY = FALSE)

             ms2_result_annotation <- pbapply::pblapply(ms2_result_annotation, function(x){
               if (nrow(x@annotation_result) > 0) {
                 if (scoring_approach == 'dp') {
                   switch (direction,
                           'reverse' = {
                             x@annotation_result <- x@annotation_result %>%
                               dplyr::filter(msms_score_reverse >= dp_cutoff & msms_matched_frag >= matched_frag_cutoff) %>%
                               dplyr::mutate(msms_score_forward = round(msms_score_forward, 4),
                                             msms_score_reverse = round(msms_score_reverse, 4))
                           },
                           'forward' = {
                             x@annotation_result <- x@annotation_result %>%
                               dplyr::filter(msms_score_forward >= dp_cutoff & msms_matched_frag >= matched_frag_cutoff) %>%
                               dplyr::mutate(msms_score_forward = round(msms_score_forward, 4),
                                             msms_score_reverse = round(msms_score_reverse, 4))
                           }
                   )
                 } else {
                   # for other score, the score was assigned as reverse score
                   x@annotation_result <- x@annotation_result %>%
                     dplyr::filter(msms_score_reverse >= dp_cutoff & msms_matched_frag >= matched_frag_cutoff) %>%
                     dplyr::mutate(msms_score_forward = round(msms_score_forward, 4),
                                   msms_score_reverse = round(msms_score_reverse, 4))
                 }

               }
               return(x)
             })


             return(ms2_result_annotation)
           }
)


#   runSpecMatch ---------------------------------------------------------------
#' @title runSpecMatch
#' @description a interphace of runing SpectraTools
#' @author Zhiwei Zhou
#' @param obj_ms2_cpd1 experimental ms2 object
#' @param obj_ms2_cpd2 library ms2 object
#' @param mz_tol_ms2 Default: 35 ppm
#' @param scoring_approach 'dp', 'bonanza', 'hybrid', 'gnps'
#' @export

# load('/home/zhouzw/Data_processing/20210224_metdna2_development_test/obj_ms2_cpd1_for_metdna2.RData')
# load('/home/zhouzw/Data_processing/20210224_metdna2_development_test/obj_ms2_cpd2_for_metdna2.RData')
# score_dp <- runSpecMatch(obj_ms2_cpd1 = obj_ms2_cpd1, obj_ms2_cpd2 = obj_ms2_cpd2, mz_tol_ms2 = 35, scoring_approach = 'dp')
# score_bonanza <- runSpecMatch(obj_ms2_cpd1 = obj_ms2_cpd1, obj_ms2_cpd2 = obj_ms2_cpd2, mz_tol_ms2 = 35, scoring_approach = 'bonanza')
# score_hybrid <- runSpecMatch(obj_ms2_cpd1 = obj_ms2_cpd1, obj_ms2_cpd2 = obj_ms2_cpd2, mz_tol_ms2 = 35, scoring_approach = 'hybrid')
# score_gnps <- runSpecMatch(obj_ms2_cpd1 = obj_ms2_cpd1, obj_ms2_cpd2 = obj_ms2_cpd2, mz_tol_ms2 = 35, scoring_approach = 'gnps')

setGeneric(name = 'runSpecMatch',
           def = function(
             obj_ms2_cpd1,
             obj_ms2_cpd2,
             mz_tol_ms2 = 35,
             scoring_approach = c('dp', 'bonanza', 'hybrid', 'gnps'),
             ...
           ){

             match.arg(scoring_approach)

             switch (scoring_approach,
                     'dp' = {
                       intensityNormedMethod <- 'maximum'
                       methodScore <- 'dp'
                     },
                     'bonanza' = {
                       intensityNormedMethod <- 'bonanza'
                       methodScore <- 'bonanza'
                     },
                     'hybrid' = {
                       intensityNormedMethod <- 'maximum'
                       methodScore <- 'hybrid'
                     },
                     'gnps' = {
                       intensityNormedMethod <- 'gnps'
                       methodScore <- 'gnps'
                     }
             )

             matchParam <- SpectraTools::MatchParam(ppm = mz_tol_ms2,
                                                    cutoff = 0,
                                                    weightIntensity = 1,
                                                    weightMZ = 0,
                                                    normIntensity = TRUE,
                                                    tuneLibSpectra = TRUE,
                                                    intensityExpNormed = TRUE,
                                                    intensityLibNormed = TRUE,
                                                    includePrecursor = TRUE,
                                                    ppmPrecursorFilter = 30,
                                                    thrIntensityAbs = 0,
                                                    thrIntensityRel = 0,
                                                    intensityNormedMethod = intensityNormedMethod,
                                                    methodMatch = 'direct',
                                                    methodScore = methodScore) %>%
               new(Class = 'MatchParam')


             result <- try(SpectraTools::MatchSpectra(dataExp = obj_ms2_cpd1,
                                                      dataRef = obj_ms2_cpd2,
                                                      matchParam),
                           silent = TRUE)




             # # if the spectra has only one fragment, and it larger than precursor, it was removed
             # if (length(result) == 0) {
             #   n_frag_cpd1 <- nrow(obj_ms2_cpd1@spectra[[1]])
             #   n_frag_cpd2 <- nrow(obj_ms2_cpd2@spectra[[1]])
             #   n_frag_match <- 0
             #   n_nl_match <- 0
             #
             #   if (scoring_approach == 'dp') {
             #
             #   }
             #   result@info <- obj_ms2_cpd2@info %>%
             #     dplyr::mutate(scoreReverse = 0,
             #                   scoreForward = 0) %>%
             #     dplyr::mutate(n_frag_cpd1 = n_frag_cpd1,
             #                   n_frag_cpd2 = n_frag_cpd2,
             #                   n_frag_match = n_frag_match,
             #                   n_frag_nl = n_nl_match) %>%
             #     dplyr::mutate(n_frag_total = n_frag_match + n_frag_nl)
             # }

             # add matched_frag and matched_nl into the result table
             stat_matched_frag <- lapply(seq_along(result@matchedFragments), function(i){
               temp_matchedFragments <- result@matchedFragments[[i]]

               if (length(temp_matchedFragments) > 0) {
                 n_frag_cpd1 <- temp_matchedFragments %>%
                   dplyr::filter(intensity > 0) %>%
                   dplyr::count() %>%
                   dplyr::pull()

                 n_frag_cpd2 <- temp_matchedFragments %>%
                   dplyr::filter(intensityExp > 0) %>%
                   dplyr::count() %>%
                   dplyr::pull()

                 n_frag_match <- temp_matchedFragments %>%
                   dplyr::filter(intensity > 0 & intensityExp > 0) %>%
                   dplyr::count() %>%
                   dplyr::pull()

               } else {
                 n_frag_cpd1 <- nrow(obj_ms2_cpd1@spectra[[1]])
                 n_frag_cpd2 <- nrow(obj_ms2_cpd2@spectra[[i]])
                 n_frag_match <- 0
               }

               if (scoring_approach == 'dp') {
                 n_nl_match <- 0
               } else {
                 n_nl_match <- result@nlFragments[[1]] %>%
                   dplyr::filter(intensity > 0 & intensityExp > 0) %>%
                   dplyr::count() %>%
                   dplyr::pull()
               }

               temp_result <- tibble::tibble(n_frag_cpd1 = n_frag_cpd1,
                                             n_frag_cpd2 = n_frag_cpd2,
                                             n_frag_match = n_frag_match,
                                             n_frag_nl = n_nl_match) %>%
                 dplyr::mutate(n_frag_total = n_frag_match + n_frag_nl)

             })

             stat_matched_frag <- stat_matched_frag %>% dplyr::bind_rows()

             result@info <- result@info %>%
               dplyr::bind_cols(stat_matched_frag)

             return(result)
           })








################################################################################
# convertSpecAnnotationClass2Table ---------------------------------------------
#' @title convertSpecAnnotationClass2Table
#' @param ms1_data
#' @param result_annotation
#' @param lib_meta
#' @param instrument
#' @param is_rt_score
#' @param is_ccs_score
#' @param is_msms_score
#' @param direction
#' @export

setGeneric(name = 'convertSpecAnnotationClass2Table',
           def = function(
             ms1_data,
             result_annotation,
             lib_meta,
             instrument = c("SciexTripleTOF", "AgilentQTOF", "BrukerQTOF", "ThermoOrbitrap", "ThermoExploris",
                            'WatersQTOF','WatersTWIMMS', "AgilentDTIMMS", "BrukerTIMS"),
             is_rt_score = TRUE,
             is_ccs_score = FALSE,
             is_msms_score = TRUE,
             # direction = c('forward', 'reverse'),
             rt_cutoff = 30, # second
             dp_cutoff = 0.8,
             matched_frag_cutoff = 1,
             ...
           ){
             instrument <- match.arg(instrument)
             if (!(instrument %in% c("AgilentDTIMMS", "BrukerTIMS", "WatersTWIMMS"))) {
               is_ccs_score <- FALSE
             }

             cat('\n');cat('Generate metabolite annotation report\n')

             pbapply::pboptions(type='timer', char='+')
             report_table <- pbapply::pblapply(seq_along(result_annotation), function(i){
               # i <- 192; i <- 28
               x <- result_annotation[[i]]

               if (nrow(x@annotation_result) > 0) {
                 result <- x@annotation_result %>%
                   dplyr::arrange(desc(msms_score_forward),
                                  desc(msms_score_reverse),
                                  desc(ccs_score),
                                  desc(rt_score)) %>%
                   dplyr::mutate(feature = x@peak_info$name) %>%
                   dplyr::select(feature, dplyr::everything())
                 return(result)
               } else {
                 return(NULL)
               }

             })

             report_table <- report_table %>% dplyr::bind_rows()

             if (length(report_table) == 0) {
               report_result <- ms1_data$info %>%
                 dplyr::bind_cols(ms1_data$subject) %>%
                 dplyr::mutate(id_reverse_summary = NA,
                               id_forward_summary = NA)

               if (!is_ccs_score) {
                 report_result <- report_result %>%
                   dplyr::select(-ccs)
               }

               return(report_result)
             }

             id_summary_forward <- report_table %>%
               dplyr::filter(msms_score_forward >= dp_cutoff) %>%
               dplyr::mutate(id = paste0('score{', msms_score_forward, '}',
                                         'frag{', msms_matched_frag, '}',
                                         'adduct{', adduct, '}',
                                         'name{', name, '}',
                                         'labid{', id, '}')) %>%
               dplyr::group_by(feature) %>%
               dplyr::summarise(id_forward_summary = paste(id, collapse = ';')) %>%
               dplyr::ungroup()

             id_summary_reverse <- report_table %>%
               # dplyr::filter(msms_score_reverse >= dp_cutoff) %>%
               dplyr::mutate(id = paste0('score{', msms_score_reverse, '}',
                                         'frag{', msms_matched_frag, '}',
                                         'adduct{', adduct, '}',
                                         'name{', name, '}',
                                         'labid{', id, '}')) %>%
               dplyr::group_by(feature) %>%
               dplyr::summarise(id_reverse_summary = paste(id, collapse = ';')) %>%
               dplyr::ungroup()


             report_result <- ms1_data$info %>%
               dplyr::bind_cols(ms1_data$subject) %>%
               dplyr::left_join(id_summary_reverse, by = c('name' = 'feature')) %>%
               dplyr::left_join(id_summary_forward, by = c('name' = 'feature'))

             if (!is_ccs_score) {
               report_result <- report_result %>%
                 dplyr::select(-ccs)
             }

             return(report_result)

           })
