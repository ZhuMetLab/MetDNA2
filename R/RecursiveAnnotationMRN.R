################################################################################
# annotateMRN ------------------------------------------------------------------
#' @title annotateMRN
#' @description metabolic reaction network based recursive annotation
#' @author Zhiwei Zhou, Xiaotao Shen
#' @param annotation_result "ms2_match_annotation_result.csv"
#' @param ms2_data 'ms2'
#' @param prefer_adduct prefer adduct for RT prediction, including "all", "M+H", "M+Na", "M-H". Default: 'M+H'
#' @param use_default_md Whether use default MD.
#' @param metdna_version MetDNA version, 'version1', 'version2'. Default: 'version2'
#' @param direction 'forward', 'reverse'. Default: 'reverse'
#' @param Use seed from reserved match or forward match. 'forward', 'reverse'. Default: 'reverse'
#' @param column 'hilic', 'rp'
#' @param method_lc 'Amide12min', 'Amide23min', 'Other', 'MetlinRP', 'RP12min'
#' @param excluded_adduct NULL
#' @param polarity 'positive', 'negative'. Default: 'positive'
#' @param extension_step The used extended MRN with n step reaction. Steps: '0', '1', '2', '3', '4', '5', '6', '7', '8'
#' @param threads Default: 4
#' @param path '.'
#' @param max_isotope Maximum considerable isotope. Default: 4
#' @param mz_tol Default: 25 ppm
#' @param rt_tol1 RT tolerance of adductAnnotation and isotopeAnnotation. Default: 3s
#' @param rt_tol2 RT tolerance of metAnnotation. Default: 30%
#' @param cor_tol Default: 0
#' @param int_tol isotope intensity tolerance. Defautl: 500
#' @param dp_tol 0.5
#' @param max_step the maximum step of reaction pair during recursive annotation
#' @param score_cutoff Total score cutoff. Default: 0
#' @param remain FALSE
#' @param remain_per 0.5
#' @param seed_neighbor_match_plot TRUE
#' @param candidate_num reserve top N candidate number. Default: 3
#' @return
#' @export


# annotation_result <- "ms2_match_annotation_result.csv"
# ms2_data <- "ms2"
# prefer_adduct <- 'all'
# ## metdna_version  <- 'version1'
# metdna_version  <- 'version2'
# scoring_approach <- 'dp'
# lib <- 'zhuMetLib'
# column <- "hilic"
# method_lc <- 'Amide12min'
# use_default_md <- TRUE
# ## column <- 'hilic'
# excluded_adduct <- NULL
# ## polarity <- 'positive'
# polarity <- 'negative'
# extension_step = '2'
# threads <- 3
# path <- '/home/zhouzw/Data_processing/20210118_nist_urine_analysis/04_20210116_pos_only_full_8_DDA_emrn_2/'
# ## path <- 'H:/00_projects/03_MetDNA2/00_data/20200922_metdna2_development/pos_201012_compare_metdna/'
# max_isotope = 4
# mz_tol = 25
# rt_tol1 = 3
# rt_tol2 = 30
# cor_tol = 0
# int_tol = 500
# dp_tol = 0.5
# max_step = 3
# score_cutoff = 0
# remain = FALSE
# remain_per = 0.5
# seed_neighbor_match_plot = FALSE
# candidate_num = 3
# # load('E:/zhiwei_create/packages/MetDNA_V1.21/data//kegg.rpair2.rda')
#
# annotation_result = "ms2_match_annotation_result.csv"
# ms2_data = "ms2"
# prefer_adduct = "all"
# use_default_md = TRUE
# metdna_version = 'version2'
# scoring_approach = 'dp'
# direction = 'reverse'
# instrument = 'SciexTripleTOF'
# lib = 'zhuMetLib'# annotation_result = "ms2_match_annotation_result.csv"
# ms2_data = "ms2"
# prefer_adduct = "all"
# use_default_md = TRUE
# metdna_version = 'version2'
# scoring_approach = 'dp'
# direction = 'reverse'
# instrument = 'IMMS'
# lib = 'zhuMetLib'
# column = "hilic"
# method_lc = 'Amide23min'
# excluded_adduct = NULL
# polarity = "positive"
# extension_step = '2'
# threads = 3
# path = '/home/zhouzw/Data_processing/20210522_debug/'
# max_isotope = 4
# mz_tol = 25
# rt_tol1 = 3 # s
# rt_tol2 = 30 # %
# ccs_tol = 4 # %
# cor_tol = 0
# int_tol = 500
# dp_tol = 0.5
# matched_frag_tol = 1
# max_step = 3
# score_cutoff = 0
# remain = FALSE
# remain_per = 0.5
# seed_neighbor_match_plot = TRUE
# candidate_num = 3

setGeneric(name = 'annotateMRN',
           def = function(
             annotation_result = "ms2_match_annotation_result.csv",
             ms2_data = "ms2",
             prefer_adduct = c("all", "M+H", "M+Na", "M-H"),
             use_default_md = TRUE,
             use_pretrained_model = FALSE,
             metdna_version = c('version2', 'version1'),
             scoring_approach = c('dp', 'bonanza', 'hybrid', 'gnps'),
             direction = c('reverse', 'forward'),
             instrument = c("SciexTripleTOF", "AgilentQTOF", "BrukerQTOF", "ThermoOrbitrap",
                            'ThermoExploris', "AgilentDTIMMS", "BrukerTIMS", 'WatersQTOF','WatersTWIMMS'),
             lib = c('zhumetlib_qtof', 'zhumetlib_orbitrap', 'fiehnHilicLib'),
             column = c("hilic", "rp"),
             method_lc = c('Amide12min', 'Amide23min', 'Other', 'MetlinRP', 'RP12min'),
             excluded_adduct = NULL,
             polarity = c("positive", "negative"),
             extension_step = c('0', '1', '2', '3', '4', '5', '6', '7', '8'),
             threads = 3,
             path = '.',
             max_isotope = 4,
             mz_tol = 25,
             mz_ppm_thr = 400,
             rt_tol1 = 3,
             rt_tol2 = 30,
             ccs_tol = 4,
             cor_tol = 0,
             int_tol = 500,
             dp_tol = 0.5,
             matched_frag_tol = 1,
             whether_link_frag = FALSE,
             max_step = 3,
             score_cutoff = 0,
             remain = FALSE,
             remain_per = 0.5,
             seed_neighbor_match_plot = TRUE,
             candidate_num = 3,
             # test_old_mrn = c('v0.3', 'v0.4'),
             test_adduct_version = c('version2', 'version1'),
             test_evaluation = c('No', '200STD', '46STD'),
             whether_use_predRT = TRUE,
             ...
           ){
             path_output  <- file.path(path, "02_result_MRN_annotation")
             dir.create(file.path(path_output, "00_intermediate_data"),
                        showWarnings = FALSE, recursive = TRUE)

             # parameters check
             column <- match.arg(column)
             polarity <- match.arg(polarity)
             prefer_adduct <- match.arg(prefer_adduct)
             extension_step <- match.arg(extension_step)
             method_lc <- match.arg(method_lc)
             scoring_approach <- match.arg(scoring_approach)
             test_adduct_version <- match.arg(test_adduct_version)
             test_evaluation <- match.arg(test_evaluation)

             switch (metdna_version,
                     'version1' = {
                       method_lc <- 'Other'
                     },
                     'version2' = {
                       method_lc <- method_lc
                     }
             )

             # data check
             file <- c(dir(file.path(path, "01_result_initial_seed_annotation")),
                       dir(file.path(path, "01_result_initial_seed_annotation", "00_intermediate_data")))
             need_file <- c(annotation_result, ms2_data)

             file_check <- which(!need_file %in% file)

             if(length(file_check) > 0) {
               stop(paste(need_file[file_check], collapse = " & "),
                    " are not in the directory.")
             }

             # browser()
             # -----------------------------------------------------------------
             # read annotation.result from ms2Annotation

             # switch (lib,
             #         'zhuMetLib' = {
             #           data("zhuMetlib", envir = environment())
             #           inHouse_compound <- zhuMetlib$meta$compound %>% as.data.frame()
             #
             #           data("md_zhumetlib", envir = environment())
             #           md_inHouse_cpd <- md_zhumetlib
             #         },
             #         'fiehnHilicLib' = {
             #           data("fiehnHilicLib", envir = environment())
             #           inHouse_compound <- fiehnHilicLib$meta$compound %>% as.data.frame()
             #
             #           data("md_fiehnHilicLib", envir = environment())
             #           md_inHouse_cpd <- md_fiehnHilicLib
             #         }
             # )

             if (lib == 'zhumetlib_qtof') {
               data("zhuMetlib", envir = environment())
               inHouse_compound <- zhuMetlib$meta$compound %>% as.data.frame()

               data("md_zhumetlib", envir = environment())
               md_inHouse_cpd <- md_zhumetlib
             }

             if (lib == 'fiehnHilicLib') {
               data("fiehnHilicLib", envir = environment())
               inHouse_compound <- fiehnHilicLib$meta$compound %>% as.data.frame()

               data("md_fiehnHilicLib", envir = environment())
               md_inHouse_cpd <- md_fiehnHilicLib
             }

             if (lib == 'zhumetlib_orbitrap') {
               data("zhuMetlib_orbitrap", envir = environment())
               inHouse_compound <- zhuMetlib_orbitrap$meta$compound %>% as.data.frame()

               data("md_zhumetlib", envir = environment())
               md_inHouse_cpd <- md_zhumetlib
             }


             # if (test_old_mrn == 'v0.3') {
             #   data("md_mrn_emrn_v03", envir = environment())
             #   data("cpd_emrn_v03", envir = environment())
             #
             #   md_mrn_emrn <- md_mrn_emrn_v03
             #   cpd_emrn <- cpd_emrn_v03
             #   rm(list = c('md_mrn_emrn_v03', 'cpd_emrn_v03'), envir = environment())
             # } else {
             #   data('md_mrn_emrn', envir = environment())
             #   data('cpd_emrn', envir = environment())
             # }


             data('md_mrn_emrn', envir = environment())
             data('cpd_emrn', envir = environment())


             cat("Read annotation result from ms2Annotation.\n")

             # read initial seed annotation data
             data <- readInitialAnnotation(data = annotation_result,
                                           direction = direction,
                                           rt_filter = FALSE,
                                           inHouse_compound = inHouse_compound,
                                           instrument = instrument,
                                           path = file.path(path, "01_result_initial_seed_annotation"))

             # RT prediction ---------------------------------------------------
             if (metdna_version == 'version1') {
               md_emrn_cpd <- md_mrn_emrn$version1
             } else {
               md_emrn_cpd <- md_mrn_emrn$version2
             }

             if (test_evaluation == '200STD') {
               rm(list = c('md_emrn_cpd', 'cpd_emrn'));gc()
               data("md_200std", envir = environment())
               data("cpd_200stdExd", envir = environment())

               cpd_emrn <- cpd_200stdExd
               md_emrn_cpd <- md_200std
               rm(list = c('cpd_200stdExd', 'md_200std'));gc()
             }

             if (test_evaluation == '46STD') {
               rm(list = c('md_emrn_cpd', 'cpd_emrn'));gc()
               data("md_46std", envir = environment())
               data("cpd_46stdExd", envir = environment())

               cpd_emrn <- cpd_46stdExd
               md_emrn_cpd <- md_46std
               rm(list = c('cpd_46stdExd', 'md_46std'));gc()
             }

             rm(list = c('zhuMetlib', "fiehnHilicLib", 'zhuRPlib'));gc()
             rm(list = c('md_mrn_emrn', 'md_zhumetlib', "md_fiehnHilicLib"));gc()

             cat("\n");cat("Construct RT prediction model.\n")

             rt_result <- predictRT(data = data,
                                    prefer_adduct = prefer_adduct,
                                    threads = threads,
                                    md_inHouse_cpd = md_inHouse_cpd,
                                    md_kegg_cpd = md_emrn_cpd,
                                    use_default_md = use_default_md,
                                    column = column,
                                    method_lc = method_lc,
                                    metdna_version = metdna_version,
                                    use_pretrained_model = use_pretrained_model)

             save(rt_result,
                  file = file.path(path_output, "00_intermediate_data", "rt_result"),
                  compress = "xz", version = 2)

             rm(list = c("data"))
             gc()

             # add prediction RT to inhouse compound and KEGG compound
             rt_inhouse <- rt_result[[1]]
             rt_emrn <- rt_result[[2]]

             idx <- match(inHouse_compound$id, row.names(rt_inhouse))

             inHouse_compound <- inHouse_compound %>%
               tibble::as_tibble() %>%
               dplyr::mutate(rt = rt_inhouse$RT[idx],
                             attribute = rt_inhouse$attribute[idx])

             # Add the predicted rt to inHouse compound and KEGG compound
             inHouse_compound$rt[is.na(inHouse_compound$rt)] <-
               median(inHouse_compound$rt[!is.na(inHouse_compound$rt)])

             inHouse_compound$attribute <- inHouse_compound$attribute %>%
               tidyr::replace_na('Median')

             # # load lib_ccs_emrn
             # data('lib_ccs', envir = environment())



             if (metdna_version == 'version1') {
               cpd_kegg <- cpd_emrn %>%
                 dplyr::filter(source == 'KEGG_MRN') %>%
                 tidyr::separate_rows(id_kegg_synonyms, sep = ';')
               temp_name <- cpd_emrn %>%
                 dplyr::filter(source == 'KEGG_MRN') %>%
                 tidyr::separate_rows(name_synonyms, sep = ';') %>%
                 dplyr::pull(name_synonyms)

               idx <- match(cpd_kegg$id, row.names(rt_emrn))
               cpd_kegg <- cpd_kegg %>%
                 dplyr::mutate(name_synonyms = temp_name,
                               rt = rt_emrn$RT[idx],
                               attribute = rt_emrn$attribute[idx])

               cpd_kegg$rt <- cpd_kegg$rt %>%
                 tidyr::replace_na(median(cpd_kegg$rt, na.rm = TRUE))

               cpd_kegg$attribute <- cpd_kegg$attribute %>%
                 tidyr::replace_na('Median')

               metabolite <- cpd_kegg %>%
                 dplyr::select(id_kegg_synonyms, name_synonyms, dplyr::everything()) %>%
                 dplyr::rename(id_raw = id,
                               name_raw = name,
                               id = id_kegg_synonyms,
                               name = name_synonyms,
                               RT = rt,
                               attribute = attribute) %>%
                 as.data.frame()

               rm(list = c('cpd_kegg'));gc()

             } else {

               idx <- match(cpd_emrn$id, row.names(rt_emrn))
               cpd_emrn <- cpd_emrn %>%
                 tibble::as_tibble() %>%
                 dplyr::mutate(rt = rt_emrn$RT[idx],
                               attribute = rt_emrn$attribute[idx])

               cpd_emrn$rt <- cpd_emrn$rt %>%
                 tidyr::replace_na(median(cpd_emrn$rt, na.rm = TRUE))

               cpd_emrn$attribute <- cpd_emrn$attribute %>%
                 tidyr::replace_na('Median')

               metabolite <- cpd_emrn %>%
                 dplyr::rename(RT = rt,
                               attribute = attribute) %>%
                 as.data.frame()
             }


             # load metabolite_ccs
             data('lib_ccs', envir = environment())
             metabolite_ccs <- lib_ccs$emrnlib_pred %>%
               dplyr::filter(status == 'Valid') %>%
               dplyr::select(name, adduct, pred_ccs) %>%
               dplyr::filter(name %in% metabolite$id) %>%
               dplyr::mutate(adduct = dplyr::case_when(adduct == '[M+H]+'~'M+H',
                                                       adduct == '[M+Na]+'~'M+Na',
                                                       adduct == '[M-H2O+H]+'~'M-H2O+H',
                                                       adduct == '[M+NH4]+'~'M+NH4',
                                                       adduct == '[M-H]-'~'M-H',
                                                       adduct == '[M+Na-2H]-'~'M+Na-2H',
                                                       adduct == '[M+HCOO]-'~'M+HCOO')) %>%
               dplyr::mutate(temp_id = paste(name, adduct, sep = '_'))

             rm(list = c("inHouse_compound.md", "cpd_kegg.md", "rt_result", 'lib_ccs'))
             gc()

             # Metabolite annotation
             load(file.path(path, "01_result_initial_seed_annotation", "00_intermediate_data", ms2_data))
             ms2_data <- ms2

             # switch (metdna_version,
             #         'version1' = {
             #           if (polarity == 'positive') {
             #             if (column == 'hilic') {
             #               adduct_list <- c("[M+H]+", "[M-H2O+H]+",  "[M-2H2O+H]+", "[M+NH4]+", "[M+Na]+",
             #                                "[M-H+2Na]+", "[M-2H+3Na]+", "[M+K]+", "[M-H+2K]+", "[M-2H+3K]+",
             #                                "[M+CH3CN+H]+", "[M+CH3CN+Na]+", "[M+CH3COO+2H]+")
             #             } else {
             #               adduct_list <- c("[M+H]+", "[M-H2O+H]+", "[M-2H2O+H]+", "[M+NH4]+", "[M+Na]+",
             #                                "[M-H+2Na]+", "[M-2H+3Na]+", "[M+K]+", "[M-H+2K]+", "[M-2H+3K]+",
             #                                "[M+CH3CN+H]+", "[M+CH3CN+Na]+", "[2M+H]+", "[M+HCOO+H+Na]+", "[M+HCOO+H+K]+",
             #                                "[M+HCOO+2H]+", "[2M+Na]+", "[2M+K]+")
             #             }
             #           } else {
             #             if (column == 'hilic') {
             #               adduct_list <- c("[M-H]-", '[M+Cl]-', "[M+Na-2H]-", "[M+K-2H]-", "[M+CH3CN-H]-",
             #                                "[M+NH3+Cl]-", "[M+NH4-2H]-", "[M-2H]2-", "[M+CH3COO]-", "[2M-H]-")
             #             } else {
             #               adduct_list <- c("[M-H]-", '[M+Cl]-', "[M+Na-2H]-", "[M+K-2H]-", "[M+CH3CN-H]-",
             #                                "[M+NH3+Cl]-", "[M+NH4-2H]-", "[M-H2O-H]-", "[M-2H]2-", "[M+CH3COO]-",
             #                                "[M+F]-", "[2M-H]-")
             #             }
             #           }
             #
             #         },
             #         'version2' = {
             #           if (polarity == 'positive') {
             #             adduct_list <- lib_adduct_nl$positive %>% dplyr::filter(annotation == 'Yes') %>% dplyr::pull(adduct)
             #           } else {
             #             adduct_list <- lib_adduct_nl$negative %>% dplyr::filter(annotation == 'Yes') %>% dplyr::pull(adduct)
             #           }
             #
             #           if (length(excluded_adduct) > 0) {
             #             adduct_list <- adduct_list[!which(adduct_list %in% excluded_adduct)]
             #           }
             #         }
             # )

             switch (test_adduct_version,
                     'version1' = {
                       if (polarity == 'positive') {
                         if (column == 'hilic') {
                           adduct_list <- c("[M+H]+", "[M-H2O+H]+",  "[M-2H2O+H]+", "[M+NH4]+", "[M+Na]+",
                                            "[M-H+2Na]+", "[M-2H+3Na]+", "[M+K]+", "[M-H+2K]+", "[M-2H+3K]+",
                                            "[M+CH3CN+H]+", "[M+CH3CN+Na]+", "[M+CH3COO+2H]+")
                         } else {
                           adduct_list <- c("[M+H]+", "[M-H2O+H]+", "[M-2H2O+H]+", "[M+NH4]+", "[M+Na]+",
                                            "[M-H+2Na]+", "[M-2H+3Na]+", "[M+K]+", "[M-H+2K]+", "[M-2H+3K]+",
                                            "[M+CH3CN+H]+", "[M+CH3CN+Na]+", "[2M+H]+", "[M+HCOO+H+Na]+", "[M+HCOO+H+K]+",
                                            "[M+HCOO+2H]+", "[2M+Na]+", "[2M+K]+")
                         }
                       } else {
                         if (column == 'hilic') {
                           adduct_list <- c("[M-H]-", '[M+Cl]-', "[M+Na-2H]-", "[M-H2O-H]-", "[M+K-2H]-", "[M+CH3CN-H]-",
                                            "[M+NH3+Cl]-", "[M+NH4-2H]-", "[M-2H]2-", "[M+CH3COO]-", "[2M-H]-")
                         } else {
                           adduct_list <- c("[M-H]-", '[M+Cl]-', "[M+Na-2H]-", "[M-H2O-H]-", "[M+K-2H]-", "[M+CH3CN-H]-",
                                            "[M+NH3+Cl]-", "[M+NH4-2H]-", "[M-2H]2-", "[M+F]-", "[2M-H]-")
                         }
                       }

                     },
                     'version2' = {
                       if (polarity == 'positive') {
                         adduct_list <- lib_adduct_nl$positive %>% dplyr::filter(annotation == 'Yes') %>% dplyr::pull(adduct)
                       } else {
                         adduct_list <- lib_adduct_nl$negative %>% dplyr::filter(annotation == 'Yes') %>% dplyr::pull(adduct)
                       }

                       if (length(excluded_adduct) > 0) {
                         adduct_list <- adduct_list[!which(adduct_list %in% excluded_adduct)]
                       }
                     }
             )

             adduct_table <- convertAdductList2AdductTable(adduct_list = adduct_list,
                                                           polarity = polarity,
                                                           column = column)

             # metabolite <- cpd_kegg %>%
             #   dplyr::rename(RT = rt,
             #                 attribute = attribute) %>%
             #   as.data.frame()

             # MRN
             # if (test_old_mrn == 'v0.3') {
             #   data("reaction_pair_network_v03", envir = environment())
             #   reaction_pair_network <- reaction_pair_network_v03
             #   rm('reaction_pair_network_v03');gc()
             # } else {
             #   data("reaction_pair_network", envir = environment())
             # }
             data("reaction_pair_network", envir = environment())

             switch (metdna_version,
                     'version1' = {
                       metabolic_network <- reaction_pair_network$version1
                     },
                     'version2' = {
                       extension_step <- paste0('step', extension_step)
                       metabolic_network <- reaction_pair_network$version2[[extension_step]]
                     }
             )

             if (test_evaluation == '200STD') {
               data("reaction_pair_network_200std", envir = environment())
               # extension_step <- paste0('step', extension_step)
               metabolic_network <- reaction_pair_network_200std[[extension_step]]
               rm(reaction_pair_network_200std);gc()
             }

             if (test_evaluation == '46STD') {
               data("reaction_pair_network_46std", envir = environment())
               # extension_step <- paste0('step', extension_step)
               metabolic_network <- reaction_pair_network_46std[[extension_step]]
               rm(reaction_pair_network_46std);gc()
             }


             rm(reaction_pair_network);gc()


             annotateNeighbor(data = annotation_result,
                              ms2 = ms2_data,
                              adduct_table = adduct_table,
                              max_isotope = max_isotope,
                              metdna_version = metdna_version,
                              direction = direction,
                              polarity = polarity,
                              mz_tol = mz_tol,
                              mz_ppm_thr = mz_ppm_thr,
                              rt_tol1 = rt_tol1,
                              rt_tol2 = rt_tol2,
                              cor_tol = cor_tol,
                              ccs_tol = ccs_tol,
                              int_tol = int_tol,
                              dp_tol = dp_tol,
                              matched_frag_tol = matched_frag_tol,
                              whether_link_frag = whether_link_frag,
                              scoring_approach = scoring_approach,
                              max_step = max_step,
                              metabolite = metabolite,
                              metabolite_ccs = metabolite_ccs,
                              metabolic_network = metabolic_network,
                              inHouse_compound = inHouse_compound,
                              threads = threads,
                              score_cutoff = score_cutoff,
                              remain = remain,
                              remain_per = remain_per,
                              path = path,
                              seed_neighbor_match_plot = seed_neighbor_match_plot,
                              candidate_num = candidate_num,
                              instrument = instrument,
                              method_lc = method_lc,
                              test_evaluation = test_evaluation,
                              whether_use_predRT = whether_use_predRT,
                              ...)

             rm(list = c("inHouse_compound", "metabolic_network"));gc()

             cat("\n");cat("annotateMRN is done.\n")

           })




################################################################################
#   readInitialAnnotation --------------------------------------------------------

#' @title readInitialAnnotation
#' @description Read the data from ms2Annotation and seperate it to tags and sample
#' @param data The file name (csv) from ms2Annotation
#' @param path The work directory
#' @param direction = c('reverse', 'forward')
#' @param rt_filter Filter annotation according or not.
#' @param rt_tol The tolerance of RT (\%).
#' @param inHouse_compound in house compound database.
#' @author Zhiwei Zhou, Xiaotao Shen
#' @export

# data <- 'ms2_match_annotation_result.csv'
# ##path <- 'H:/00_projects/03_MetDNA2/00_data/20200922_metdna2_development/pos_201012_compare_metdna/01_result_initial_seed_annotation/'
# ##path <- 'E:/zhiwei_create/02_MetDNA2/01_data_processing/20201005_metdna2_development/pos_201005_metdna2/01_result_initial_seed_annotation'
# path <- '/home/zhouzw/Data_processing/20210830_s9fraction_incubation_on_46STDs/peak_area_filterd/pos/01_result_initial_seed_annotation/'
# inHouse_compound <- zhuMetlib$meta$compound %>% tidyr::separate_rows('id_zhulab_synoyms', sep = ';')

# direction <- 'reverse'
# test <- readInitialAnnotation(data = data, path = path, metdna_version = 'version2', inHouse_compound = inHouse_compound)
#
# test <- readInitialAnnotation(data = 'ms2_match_annotation_result.csv',
#                               path = "/home/zhouzw/Data_processing/20210720_modify_MetDNA2_for_200STD_validation/POS_200std_mouse_liver/01_result_initial_seed_annotation/",
#                               metdna_version = 'version2',
#                               direction = 'reverse',
#                               instrument = 'SciexTripleTOF',
#                               inHouse_compound = inHouse_compound)

setGeneric(name = "readInitialAnnotation",
           def = function(data = "ms2_match_annotation_result.csv",
                          path = ".",
                          # metdna_version = c('version1', 'version2'),
                          rt_filter = FALSE,
                          rt_tol = 30,
                          direction = c('reverse', 'forward'),
                          instrument = c("SciexTripleTOF", "AgilentQTOF", "BrukerQTOF", "ThermoOrbitrap",
                                         'ThermoExploris', "AgilentDTIMMS", "BrukerTIMS", 'WatersQTOF','WatersTWIMMS'),
                          inHouse_compound,
                          # test_evaluation = c('No', '200STD'),
                          ...){

             # metdna_version <- match.arg(metdna_version)
             direction <- match.arg(direction)

             options(warn = -1)
             data <- read.csv(file.path(path, data), stringsAsFactors = FALSE)

             # data <- readr::read_csv(file.path(path, data),
             #                         col_types = readr::cols(),
             #                         progress = TRUE)
             # data <- as.data.frame(data)

             if (!(instrument %in% c("AgilentDTIMMS", "BrukerTIMS", 'WatersTWIMMS'))) {
               tags.idx <- match(c("name", "mz", "rt",
                                   "id_reverse_summary",
                                   "id_forward_summary"),
                                 colnames(data))

               tags <- data[,tags.idx] %>%
                 dplyr::mutate('ccs' = -1) %>%
                 dplyr::select("name", "mz", "rt", 'ccs',
                               "id_reverse_summary",
                               "id_forward_summary")
               sample <- data[,-tags.idx]

               rm(list = "data")
               gc()
             } else {
               tags.idx <- match(c("name", "mz", "rt", 'ccs',
                                   "id_reverse_summary",
                                   "id_forward_summary"),
                                 colnames(data))

               if (sum(is.na(tags.idx))) {
                 stop('Please check the format of annotation_table\n')
               }

               tags <- data[,tags.idx]
               sample <- data[,-tags.idx]

               rm(list = "data")
               gc()
             }


             if (direction == 'reverse') {
               hit.reverse <- tags$id_reverse_summary
               idx1 <- which(!is.na(hit.reverse))
             } else {
               hit.reverse <- tags$id_forward_summary
               idx1 <- which(!is.na(hit.reverse))
             }

             hit.reverse1 <- hit.reverse[idx1]
             peak.rt <- tags$rt[idx1]
             peak.mz <- tags$mz[idx1]
             peak.ccs <- tags$ccs[idx1]

             rm(list = c("hit.reverse"))
             gc()

             hit.reverse2 <- lapply(hit.reverse1,
                                    function(x) {strsplit(x, split = ";")[[1]]})

             rm(list = c("hit.reverse1"))
             gc()

             labid2 <- lapply(hit.reverse2, function(x) {
               y <- stringr::str_extract(string = x,
                                         pattern = "labid\\{[a-zA-Z0-9\\_\\,]+\\}")
               y <- gsub(pattern = "labid\\{", replacement = "", x = y)
               y <- gsub(pattern = "\\}", replacement = "", x = y)
               y
             })

             adduct2 <- lapply(hit.reverse2, function(x){
               result <- stringr::str_extract_all(string = x,
                                                  pattern = "adduct\\{.+\\}name")
               unlist(result)
             })

             adduct2 <-
               lapply(adduct2, function(x) {
                 gsub(pattern = "adduct\\{", replacement = "", x = x)})
             adduct2 <-
               lapply(adduct2, function(x) {
                 gsub(pattern = "\\}name", replacement = "", x = x)})
             adduct2 <-
               lapply(adduct2, function(x) {
                 gsub(pattern = "\\[", replacement = "", x = x)})
             adduct2 <-
               lapply(adduct2, function(x) {
                 result <- sapply(x, function(y){
                   if (y == 'M]+' | y == 'M]-') {
                     gsub(pattern = "\\]", replacement = "", x = y)
                   } else {
                     gsub(pattern = "\\][\\+|\\-]", replacement = "", x = y)
                   }
                 })

                 names(result) <- NULL
                 result
               })


             # ##remove M- adduct
             # mapply(function(x, y){
             # temp.index <- which(y != "M-")
             # if(length(temp.index) == 0) return(NA)
             # temp.index
             # },
             # x = labid2,
             # y = adduct2)

             ##remove wrong annotation, for examle n01170 (C18H22N2) can't
             ## -H2O

             right.idx <- mapply(function(x, y) {
               temp.idx <- match(x, inHouse_compound$id)
               temp.formula <- inHouse_compound$formula[temp.idx]
               temp.formula1 <- mapply(function(x, y){sumFormula(x,y)},
                                       x = temp.formula,
                                       y = y)
               which(!is.na(temp.formula1))
             },
             x = labid2,
             y = adduct2)

             remove.idx <- which(unlist(lapply(right.idx, length)) == 0)
             if(length(remove.idx) != 0){
               right.idx <- right.idx[-remove.idx]
               hit.reverse2 <- hit.reverse2[-remove.idx]
               labid2 <- labid2[-remove.idx]
               adduct2 <- adduct2[-remove.idx]
               idx1 <- idx1[-remove.idx]
             }

             peak.rt <- tags$rt[idx1]
             peak.mz <- tags$mz[idx1]
             peak.ccs <- tags$ccs[idx1]


             hit.reverse2 <- mapply(function(x, y) {
               x[y]
             },
             x = hit.reverse2,
             y = right.idx)

             labid2 <- mapply(function(x, y) {
               x[y]
             },
             x = labid2,
             y = right.idx)

             adduct2 <- mapply(function(x, y) {
               x[y]
             },
             x = adduct2,
             y = right.idx)

             rm(list = c("adduct2"))
             gc()


             if(!rt_filter){
               standard.rt <- lapply(labid2, function(x) {
                 rep(NA, length(x))
               })

               rt.error2 <- mapply(function(x,y) {
                 error <- abs(y - x)*100/y
               }, x = peak.rt, y = standard.rt)

               standard.ccs <- lapply(labid2, function(x) {
                 rep(NA, length(x))
               })

               ccs.error2 <- mapply(function(x,y) {
                 error <- abs(y - x)*100/y
               },
               x = peak.ccs,
               y = standard.ccs)

               index <- mapply(function(x,y) {
                 idx <- abs(y - x)*100/y
                 idx[is.na(idx)] <- 0
                 idx <- sapply(idx, function(x) {ifelse(x > rt_tol, FALSE, TRUE)})
                 idx
               }, x = peak.rt, y = standard.rt)
             }else{
               standard.rt <- lapply(labid2, function(x) {
                 # temp.idx <- match(x, rt.all$`in house ID`)
                 temp.idx <- match(x, inHouse_compound$id)
                 # temp.rt <- rt.all$`RT (min)`[temp.idx]
                 temp.rt <- inHouse_compound$rt[temp.idx]
                 temp.rt
               })

               rt.error2 <- mapply(function(x,y) {
                 error <- abs(y - x)*100/y
               }, x = peak.rt, y = standard.rt)

               standard.ccs <- lapply(labid2, function(x) {
                 rep(NA, length(x))
               })

               ccs.error2 <- mapply(function(x,y) {
                 error <- abs(y - x)*100/y
               },
               x = peak.ccs,
               y = standard.ccs)

               index <- mapply(function(x,y) {
                 idx <- abs(y - x)*100/y
                 idx[is.na(idx)] <- rt_tol+1
                 idx <- sapply(idx, function(x) {ifelse(x > rt_tol, FALSE, TRUE)})
                 idx
               }, x = peak.rt, y = standard.rt)

             }

             # # remove RT filters
             # # default: assign the RT error as zero, and keep all annotations
             # standard.rt <- lapply(labid2, function(x) {
             #   rep(NA, length(x))
             # })
             #
             # rt.error2 <- mapply(function(x,y) {
             #   error <- abs(y - x)*100/y
             # },
             # x = peak.rt,
             # y = standard.rt)
             #
             # standard.ccs <- lapply(labid2, function(x) {
             #   rep(NA, length(x))
             # })
             #
             # ccs.error2 <- mapply(function(x,y) {
             #   error <- abs(y - x)*100/y
             # },
             # x = peak.ccs,
             # y = standard.ccs)
             #
             # # keep all annotations for RT/CCS NA annotation
             # index <- mapply(function(x,y) {
             #   idx <- abs(y - x)*100/y
             #   idx[is.na(idx)] <- 0
             #   idx <- sapply(idx, function(x) {ifelse(x > 30, FALSE, TRUE)})
             #   idx
             # },
             # x = peak.rt,
             # y = standard.rt)

             rm(list = c("peak.rt", 'peak.ccs'))
             gc()

             hit.reverse3 <-
               mapply(function(x,y){x[y]},
                      x = hit.reverse2,
                      y = index)

             rm(list = c("hit.reverse2"))
             gc()

             hit.reverse4 <-
               lapply(hit.reverse3, function(x) {paste(x, collapse = ";")})
             hit.reverse4 <- unlist(hit.reverse4)

             rm(list = c("hit.reverse3"))
             gc()

             rt.error3 <- mapply(function(x,y){x[y]}, x = rt.error2, y = index)
             rm(list = c("rt.error2"))
             gc()
             rt.error4 <-
               lapply(rt.error3, function(x) {paste(x, collapse = ";")})
             rt.error4 <- unlist(rt.error4)
             rm(list = c("rt.error3"))
             gc()

             # remove CCS filters
             ccs.error3 <- mapply(function(x,y){x[y]}, x = ccs.error2, y = index)
             rm(list = c("ccs.error2"))
             gc()
             ccs.error4 <-
               lapply(ccs.error3, function(x) {paste(x, collapse = ";")})
             ccs.error4 <- unlist(ccs.error4)
             rm(list = c("ccs.error3"))
             gc()


             ms2.sim3 <-
               stringr::str_extract_all(string = hit.reverse4,
                                        pattern = "score\\{0\\.[0-9]+\\}|score\\{1\\}")
             ms2.sim3 <-
               lapply(ms2.sim3, function(x) {
                 gsub(pattern = "score\\{", replacement = "", x = x)})
             ms2.sim3 <-
               lapply(ms2.sim3, function(x) {
                 gsub(pattern = "\\}", replacement = "", x = x)})
             ms2.sim4 <-
               lapply(ms2.sim3, function(x) {paste(x, collapse = ";")})
             ms2.sim4 <- unlist(ms2.sim4)
             rm(list = c("ms2.sim3"))
             gc()


             nfrag3 <-
               stringr::str_extract_all(string = hit.reverse4,
                                        pattern = "frag\\{\\d+\\}")
             nfrag3 <-
               lapply(nfrag3, function(x) {
                 gsub(pattern = "frag\\{", replacement = "", x = x)})
             nfrag3 <-
               lapply(nfrag3, function(x) {
                 gsub(pattern = "\\}", replacement = "", x = x)})
             nfrag4 <-
               lapply(nfrag3, function(x) {paste(x, collapse = ";")})
             nfrag4 <- unlist(nfrag4)
             rm(list = c("nfrag3"))
             gc()

             # adduct3 <- stringr::str_extract_all(string = hit.reverse4,
             #                                     pattern = "adduct\\{[\\(0-9a-zA-Z\\+\\-]+")
             # adduct3 <-
             #   lapply(adduct3, function(x) {
             #     gsub(pattern = "adduct\\{", replacement = "", x = x)})
             # adduct3 <-
             #   lapply(adduct3, function(x) {
             #     gsub(pattern = "\\(", replacement = "", x = x)})
             # adduct4 <- lapply(adduct3, function(x) {paste(x, collapse = ";")})
             # adduct4 <- unlist(adduct4)
             #

             adduct3 <- lapply(hit.reverse4, function(x){
               result <- (x %>% stringr::str_split(pattern = ';'))[[1]]

               result <- stringr::str_extract(result, pattern = "adduct\\{.+\\}name")
               unlist(result)
             })
             adduct3 <-
               lapply(adduct3, function(x) {
                 gsub(pattern = "adduct\\{", replacement = "", x = x)})
             adduct3 <-
               lapply(adduct3, function(x) {
                 gsub(pattern = "\\}name", replacement = "", x = x)})
             adduct3 <-
               lapply(adduct3, function(x) {
                 gsub(pattern = "\\[", replacement = "", x = x)})
             adduct3 <-
               lapply(adduct3, function(x) {
                 gsub(pattern = "\\][\\+|\\-]", replacement = "", x = x)})
             # lapply(adduct3, function(x))

             adduct4 <- lapply(adduct3, function(x) {paste(x, collapse = ";")})
             adduct4 <- unlist(adduct4)

             name3 <-
               stringr::str_extract_all(string = hit.reverse4,
                                        pattern = "name\\{[^\\{]+\\}")
             name3 <- lapply(name3, function(x) {
               gsub(pattern = "name\\{", replacement = "", x = x)})
             name3 <- lapply(name3, function(x) {
               gsub(pattern = "\\}", replacement = "", x = x)})
             name4 <- lapply(name3, function(x) {paste(x, collapse = ";")})
             name4 <- unlist(name4)
             rm(list = "name3")
             gc()

             labid3 <- mapply(function(x,y){x[y]}, x = labid2, y = index)
             labid4 <- lapply(labid3, function(x) {paste(x, collapse = ";")})
             labid4 <- unlist(labid4)
             rm(list = "labid2")
             gc()

             rm(list = c("hit.reverse4"))
             gc()

             # extract kegg_ids ------------------------------------------------
             KEGG.ID3 <- lapply(labid3, function(x){
               inHouse_compound$id_kegg[match(x, inHouse_compound$id)]
             })

             ##if one zhulab compound matched more than 1 KEGG ID, remain first ID
             KEGG.ID3 <- lapply(KEGG.ID3, function(x){
               id.len <- nchar(x)
               id.len[is.na(id.len)] <- 6
               temp.idx <- which(id.len > 6)
               if(length(temp.idx) > 0){
                 x[temp.idx] <-
                   unlist(lapply(strsplit(x[temp.idx], split = "/"), function(x){
                     x[[1]]}
                   ))}
               x
             })


             KEGG.ID4 <- lapply(KEGG.ID3, function(x) {
               paste(x, collapse = ";")
             })
             KEGG.ID4 <- unlist(KEGG.ID4)

             # if (metdna_version == 'version1') {
             #   KEGG.ID4 <- lapply(KEGG.ID3, function(x) {
             #     paste(x, collapse = ";")
             #   })
             #   KEGG.ID4 <- unlist(KEGG.ID4)
             #   KEGG.ID4[grep('NA', KEGG.ID4)] <- NA
             # } else {
             #   KEGG.ID4 <- lapply(KEGG.ID3, function(x) {
             #     idx <- which(!is.na(x))
             #
             #     if (length(idx) == 0) {return(NA)}
             #     result <- paste(x[idx], collapse = ";")
             #     result
             #   })
             #   KEGG.ID4 <- unlist(KEGG.ID4)
             #   KEGG.ID4[grep('NA', KEGG.ID4)] <- NA
             # }
             #
             # rm(list = c("KEGG.ID3"))
             # gc()

             # calculate mz error ----------------------------------------------
             accurate.mass3 <- mapply(function(x, y) {
               temp.idx <- match(x, inHouse_compound$id)
               temp.formula <- inHouse_compound$formula[temp.idx]

               if (is.na(y)){
                 temp.adduct <- character()
               } else {
                 temp.adduct <- y
               }

               temp.formula1 <- mapply(function(x, y){sumFormula(x,y)},
                                       x = temp.formula,
                                       y = temp.adduct)
               a.mass <- sapply(temp.formula1, function(x) {
                 Rdisop::getMass(Rdisop::getMolecule(formula = x))})
               a.mass
             },
             x = labid3,
             y = adduct3)

             # test <- lapply(seq_along(labid3), function(i){
             #   cat(i, ' ')
             #   x <- labid3[[i]]
             #   y <- adduct3[[i]]
             #
             #   temp.idx <- match(x, inHouse_compound$id)
             #   temp.formula <- inHouse_compound$formula[temp.idx]
             #   temp.formula1 <- mapply(function(x, y){sumFormula(x,y)},
             #                           x = temp.formula,
             #                           y = y)
             #   a.mass <- sapply(temp.formula1, function(x) {
             #     Rdisop::getMass(Rdisop::getMolecule(formula = x))})
             #   a.mass
             # })

             rm(list = c("adduct3"))
             gc()
             # for(i in 1:length(labid3)){
             #   cat(i); cat(" ")
             # x <- labid3[[i]]
             # y <- adduct3[[i]]
             # temp.formula <- inHouse_compound$Formula[match(x, inHouse_compound$lab.ID)]
             # temp.formula1 <- mapply(function(x, y){sumFormula(x,y)},
             #                         x = temp.formula,
             #                         y = y)
             # a.mass <- sapply(temp.formula1, function(x) {
             #   Rdisop::getMass(Rdisop::getMolecule(formula = x))})
             # }


             mz.error3 <- mapply(function(x, y) {
               if(is.list(x)) x <- NA
               # error <- abs(x - y)*10^6/y
               error <- abs(x - y)*10^6/ifelse(y>=400,y,400)
               unname(error)
             },
             x = accurate.mass3,
             y = peak.mz)

             rm(list = c("peak.mz"))
             gc()


             # should be fixed
             mz.error3 <- lapply(mz.error3, function(x){
               x[!is.na(x)] <- 10
               x
             })

             mz.error4 <-
               lapply(mz.error3, function(x) {paste(x, collapse = ";")})
             mz.error4 <- unlist(mz.error4)
             mz.error4[mz.error4=="NA"] <- ""


             # formula
             Formula3 <- lapply(labid3, function(x){
               inHouse_compound$formula[match(x, inHouse_compound$id)]
             })

             rm(list = c("labid3"))
             gc()

             Formula4 <- lapply(Formula3,  function(x) {
               paste(x, collapse = ";")
             })
             Formula4 <- unlist(Formula4)

             tags <- as.data.frame(tags)
             tags <- tags[,c("name","mz","rt",'ccs')]
             colnames(tags) <- c("Peak.name","mz","rt",'ccs')
             ms2.sim <- rep(NA, nrow(tags))
             adduct <- rep(NA, nrow(tags))
             rt.error <- rep(NA, nrow(tags))
             mz.error <- rep(NA, nrow(tags))
             ccs.error <- rep(NA, nrow(tags))
             Metabolite.name <- rep(NA, nrow(tags))
             KEGG.ID <- rep(NA, nrow(tags))
             Formula <- rep(NA, nrow(tags))
             isotope <- rep(NA, nrow(tags))
             labid <- rep(NA, nrow(tags))
             nfrag <- rep(NA, nrow(tags))

             ms2.sim[idx1] <- ms2.sim4
             adduct[idx1] <- adduct4
             rt.error[idx1] <- rt.error4
             mz.error[idx1] <- mz.error4
             ccs.error[idx1] <- ccs.error4
             Metabolite.name[idx1] <- name4
             KEGG.ID[idx1] <- KEGG.ID4
             # if (test_evaluation == '200STD') {KEGG.ID[idx1] <- labid4}
             Formula[idx1] <- Formula4
             labid[idx1] <- labid4
             isotope[idx1] <- "[M]"
             nfrag[idx1] <- nfrag4

             rm(list = c("labid4", "adduct4", "ms2.sim4", "rt.error4",
                         "mz.error4", 'ccs.error4', "name4", "KEGG.ID4", "Formula4", 'nfrag4'))
             gc()

             tags <- data.frame(tags, adduct, isotope, ms2.sim, nfrag, rt.error, mz.error, ccs.error,
                                Formula, Metabolite.name, KEGG.ID, labid,
                                stringsAsFactors = FALSE)
             rm(list = c("adduct", "isotope", "ms2.sim", 'nfrag', "rt.error", "mz.error", 'ccs.error',
                         "Formula", "Metabolite.name", "KEGG.ID", "labid"))
             gc()

             tags[which(tags == "", arr.ind = TRUE)] <- NA

             ##change tags to tags2 style
             result <- list(tags = tibble::as_tibble(tags),
                            sample = tibble::as_tibble(sample))
             names(result) <- c("tags", "sample")
             result <- result
           })



################################################################################
#   predictRT ------------------------------------------------------------------

#' @title predictRT
#' @description Predict RTs of kegg metabolites using MS2 matched metabolites.
#' @author Zhiwei Zhou, Xiaotao Shen
#' @param data data From readAnnotation.
#' @param prefer.adduct The reliable adducts in this LC system. Default is all.
#' @param threads How many threads do you want to use? Default is the number of your PC threads - 3.
#' @param inHouse.compound.md The molecular descriptors of inhouse compounds.
#' @param kegg.compound.md The molecular descriptors of kegg compounds.
#' @return The predicted in house RT and KEGG RT.
#' @export

# data <- 'ms2_match_annotation_result.csv'
# ## path <- 'H:/00_projects/03_MetDNA2/00_data/20200922_metdna2_development/pos_200927_metdna2/01_result_initial_seed_annotation/'
# ## path <- 'E:/zhiwei_create/02_MetDNA2/01_data_processing/20201005_metdna2_development/pos_201005_metdna2/01_result_initial_seed_annotation'
# path <- '/home/zhouzw/Data_processing/20201102_metdna2_pos_development/01_result_initial_seed_annotation/'
# inHouse_compound <- zhuMetlib$meta$compound %>% tidyr::separate_rows('id_zhulab_synoyms', sep = ';')
#
# data <- readInitialAnnotation(data = 'ms2_match_annotation_result.csv',
#                               # rt_filter = FALSE,
#                               inHouse_compound = inHouse_compound,
#                               path = path)
# prefer_adduct <- 'M+H'
# threads <- 3
# data("md_inHouse_cpd", envir = environment())
# use_default_md = TRUE
# column <- 'hilic'
# method_lc <- 'Amide23min'
# md_inHouse_cpd <- md_zhumetlib
# md_kegg_cpd <- md_mrn_emrn$version2
#
# test <- predictRT(data = data,
#                   prefer_adduct = 'M+H',
#                   md_inHouse_cpd = md_inHouse_cpd,
#                   md_kegg_cpd = md_kegg_cpd,
#                   threads = 3,
#                   use_default_md = TRUE,
#                   column = 'hilic', method_lc = 'Amide23min')

setGeneric(name = "predictRT",
           def = function(data,
                          prefer_adduct = c("all", "M+Na", "M+H", "M+NH4", "M-H"),
                          threads = 3,
                          md_inHouse_cpd,
                          md_kegg_cpd,
                          use_default_md = TRUE,
                          column = c("hilic", "rp"),
                          method_lc = c('Amide12min', 'Amide23min', 'Other', 'MetlinRP', 'RP12min'),
                          metdna_version = c('version2', 'version1'),
                          use_pretrained_model = FALSE
           ){
             # browser()
             library(randomForest)
             prefer_adduct <- match.arg(prefer_adduct)
             column <- match.arg(column)
             tags <- data[[1]]
             sample <- data[[2]]
             sample.int <- apply(sample, 1, median)
             rm(list = "sample")
             gc()

             # for Amide12min and Amide23min methods, use pre-trained model
             if ((method_lc %in% c('Amide12min', 'Amide23min') & use_pretrained_model) &
                 (column == 'hilic')) {
               cat('Note: pre-trained', method_lc, 'model was used in RT prediction\n')
               data("rt_model_pretrained", envir = environment())
               switch (method_lc,
                       'Amide12min' = {
                         rf.reg <- rt_model_pretrained$Amide12min
                       },
                       'Amide23min' = {
                         rf.reg <- rt_model_pretrained$Amide23min
                       }
               )

               marker.name <- c('XLogP', "tpsaEfficiency", "WTPT.5", "khs.dsCH", "MLogP", "nAcid", "nBase", "BCUTp.1l")
             } else {
               idx <- which(!is.na(tags$labid))
               if(length(idx) == 0) stop("No metabolites are identified by MS2 library.\n")
               tags1 <- tags[idx,]
               sample.int1 <- sample.int[idx]
               rm(list = "sample.int")
               gc()

               # if one peak has multiple annotations, remain the annotation with highest MS2 score
               labid <- tags1$labid
               labid <- lapply(labid, function(x) {
                 strsplit(x, split = ";")[[1]][1]
               })
               labid <- unlist(labid)


               adduct <- tags1$adduct
               adduct <- lapply(adduct, function(x) {
                 strsplit(x, split = ";")[[1]][1]
               })
               adduct <- unlist(adduct)

               # remove multiple peaks matched one metabolite, keep the max intensity peak
               dup.id <- unique(labid[duplicated(labid)])
               if (length(dup.id) > 0){
                 for (i in 1:length(dup.id)){
                   temp.id <- dup.id[i]
                   temp.idx <- grep(temp.id, labid)
                   temp.int <- sample.int1[temp.idx]
                   # rm(list = c("sample.int1"))
                   # temp.adduct <- adduct[temp.idx]
                   # temp.rt <- tags1$rt[temp.idx]
                   labid[temp.idx[-which.max(temp.int)]]  <- NA
                 }
               }

               rt <- tags1$rt
               rm(list = "tags1")
               gc()

               data <- data.frame(labid, rt, adduct, stringsAsFactors = FALSE)
               data <- data[!is.na(data$labid),,drop = FALSE]

               # filter the metabolite who have no MD in md_inHouse_cpd
               temp.idx <- which(data$labid %in% rownames(md_inHouse_cpd))

               # if(length(temp.idx) == 0) stop("The metabolites from MS2 match have no MD in md_inHouse_cpd.\n")
               if(length(temp.idx) == 0) stop("No metabolites are identified by MS2 library.\n")
               if(length(temp.idx) < 70){
                 prefer_adduct <- 'all'
               }
               data <- data[temp.idx,]


               # filter data using adduct or not
               #  if the prefer adduct not fit the number cutoff, sequently adding the adduct according to their number
               adduct.order <- setdiff(names(sort(table(adduct), decreasing = TRUE)), prefer_adduct)
               if (prefer_adduct != 'all') {
                 temp.idx <- which(data$adduct %in% prefer_adduct)
                 count <- 1
                 while(length(temp.idx) < 50 & count < length(adduct.order)){
                   prefer_adduct <- c(prefer_adduct, adduct.order[count])
                   temp.idx <- which(data$adduct %in% prefer_adduct)
                   count <- count + 1
                 }
                 data <- data[temp.idx,,drop = FALSE]
               }

               if(nrow(data) <= 5) stop("No or less 5 metabolites are identified.")

               cat("There are ", nrow(data),
                   " metabolites are used for RT prediction.\n", sep = "")

               idx <- match(data$labid, rownames(md_inHouse_cpd))
               md <- md_inHouse_cpd[idx,]
               # rm(list = "md_inHouse_cpd")

               # remove NA which appear in more than 50% metabolites
               remove.idx1 <-
                 which(apply(md, 2, function(x) {sum(is.na(x)/nrow(md))})>0.5)
               md1 <- md[,-remove.idx1]
               rm(list = c("md"))
               gc()

               # impute NA
               md2 <- t(impute::impute.knn(data = t(md1))[[1]])
               rm(list = "md1")
               gc()

               #remove MD which are same in all metaboites
               remove.idx2 <- which(apply(md2, 2, sd) == 0)
               md3 <- md2[,-remove.idx2]
               rm(list = c("md2"))
               gc()

               ##construct RF model
               train.y <- data$rt
               train.x <- md3
               rm(list = c('data', 'md3'))
               gc()

               if (use_default_md){
                 switch(column,
                        "hilic" = {
                          marker.name <- c('XLogP', "tpsaEfficiency", "WTPT.5", "khs.dsCH", "MLogP", "nAcid", "nBase", "BCUTp.1l")
                        },
                        "rp" = {
                          marker.name <- c('XLogP', "WTPT.4", "WTPT.5", "ALogp2", "BCUTp.1l")})
               } else {
                 # not use in MetDNA2
                 #construct RF 100 times, find the MDs which are always apperar in
                 # top 5 as markers
                 temp.fun <- function(idx, x, y){
                   suppressMessages(library(randomForest))
                   temp <- idx
                   rf.class <- randomForest::randomForest(x = x, y = y,
                                                          replace = TRUE, importance = TRUE,
                                                          proximity = TRUE)
                   imp <- randomForest::importance(rf.class)
                   rm(list = c("rf.class", "temp"))
                   gc()
                   imp
                 }
                 cat("\n")
                 cat("Find the optimal moleculr descriptor\n")
                 imp <- BiocParallel::bplapply(1:100,
                                               temp.fun,
                                               BPPARAM = BiocParallel::SnowParam(workers = threads,
                                                                                 progressbar = FALSE),
                                               x = train.x,
                                               y = train.y)

                 md.name <- list()
                 for(i in 1:length(imp)){
                   md.name[[i]] <- names(sort(imp[[i]][,1], decreasing = TRUE)[1:5])
                 }

                 md.name <- unlist(md.name)
                 md.name <- table(md.name)
                 md.name <- sort(md.name)
                 marker.name <- names(md.name[which(md.name >= 50)])
                 rm(list = "imp")
                 gc()
                 #
                 save(marker.name, file = "marker.name", compress = "xz")
                 save(md.name, file = "md.name", compress = "xz")
               }

               idx <- match(marker.name, colnames(train.x))
               idx <- idx[!is.na(idx)]
               if(length(idx) == 0){
                 stop("Your markers are not in MD data.\n")
               }
               train.x <- train.x[,idx]

               # using default tuning function or tuning with caret package
               if (metdna_version == 'version1') {
                 cat('Tuning RF function using MetDNA1 approach\n')

                 para <- NULL
                 ntree1 <- seq(300,1000,by = 200)
                 mtry1 <- seq(1,length(marker.name),by = 1)
                 for(i in 1:length(ntree1)){
                   para <- rbind(para, cbind(ntree1[i],mtry1))
                 }
                 colnames(para) <- c("ntree", "mtry")
                 mse <- NULL
                 rf.reg <- list()
                 cat("\n")
                 cat("Find the optimal parameters\n")
                 for(i in 1:nrow(para)){
                   cat(i);cat(" ")
                   temp.ntree <- para[i,1]
                   temp.mtry <- para[i,2]
                   rf.reg[[i]] <- randomForest::randomForest(x = train.x, y = train.y,
                                                             ntree = temp.ntree, mtry = temp.mtry,
                                                             replace = TRUE, importance = TRUE, proximity = TRUE)
                   mse[i] <- mean(rf.reg[[i]]$mse)
                 }
                 cat("\n")
                 result <- data.frame(para, mse, stringsAsFactors = FALSE)
                 temp.idx <- which.min(result$mse)
                 ntree <- result$ntree[temp.idx]
                 mtry <- result$mtry[temp.idx]
                 ##
                 rf.reg <- randomForest::randomForest(x = train.x, y = train.y,
                                                      ntree = ntree,
                                                      mtry = mtry,
                                                      replace = TRUE,
                                                      importance = TRUE, proximity = TRUE)

                 rm(list = c("train.x"));gc()
               } else {
                 cat('Tuning RF function using caret approach\n')
                 train.data <- data.frame(RT = train.y, train.x, stringsAsFactors = FALSE)
                 rf.reg <- trainRandomForestWithCaret(x = train.data)

                 rm(list = c("train.x"));gc()
               }
             }

             ##predict RT in inhouse database
             test.x <- md_inHouse_cpd
             rm(list = "md_inHouse_cpd")
             gc()
             test.x <- test.x[,match(marker.name, colnames(test.x))]
             test.x <- as.data.frame(test.x)

             inHouse.compound.rt <- rep(NA, nrow(test.x))
             names(inHouse.compound.rt) <- rownames(test.x)

             # impute NA in test.x
             idx1 <-
               which(apply(test.x, 1, function(x) {sum(is.na(x))/ncol(test.x) < 0.5}))
             test.x1 <- test.x[idx1,]
             test.x1 <- t(impute::impute.knn(data = t(test.x1))[[1]])
             temp.rt <- predict(object = rf.reg, newdata = test.x1)
             names(temp.rt) <- rownames(test.x1)
             inHouse.compound.rt[idx1] <- temp.rt

             # NA give the median RT of all peaks, these compounds don't effective MDs (may be caused by lacking structures)
             inHouse.compound.rt[is.na(inHouse.compound.rt)] <- median(tags$rt)
             attribute <- rep(NA, length(inHouse.compound.rt))
             attribute[idx1] <- "Predicted"
             attribute[is.na(attribute)] <- "Median"
             inHouse.rt <- data.frame(inHouse.compound.rt, attribute, stringsAsFactors = FALSE)
             colnames(inHouse.rt)[1] <- "RT"

             # predict RT in KEGG database
             test.x <- md_kegg_cpd
             rm(list = "md_kegg_cpd")
             gc()
             test.x <- test.x[,match(marker.name, colnames(test.x))]
             test.x <- as.data.frame(test.x)

             kegg.compound.rt <- rep(NA, nrow(test.x))
             names(kegg.compound.rt) <- rownames(test.x)
             ##impute NA in text.x
             idx1 <-
               which(apply(test.x, 1, function(x) {sum(is.na(x))/ncol(test.x) < 0.5}))
             test.x1 <- test.x[idx1,]
             test.x1 <- t(impute::impute.knn(data = t(test.x1))[[1]])
             temp.rt <- predict(object = rf.reg, newdata = test.x1)
             names(temp.rt) <- rownames(test.x1)
             kegg.compound.rt[idx1] <- temp.rt
             ##NA give the median RT of all peaks
             kegg.compound.rt[is.na(kegg.compound.rt)] <- median(tags$rt)
             rm(list = "tags")
             gc()
             attribute <- rep(NA, length(kegg.compound.rt))
             attribute[idx1] <- "Predicted"
             attribute[is.na(attribute)] <- "Median"
             kegg.rt <- data.frame(kegg.compound.rt, attribute, stringsAsFactors = FALSE)
             colnames(kegg.rt)[1] <- "RT"

             rt <- list(inHouse.rt, kegg.rt)
             names(rt) <- c("inHouse.rt", "KEGG.rt")
             rt <- rt
           })


#' @title trainRandomForestWithCaret
#' @author Zhiwei Zhou
#' @param x the first column is RT, others are MDs
#' @return RF model

trainRandomForestWithCaret <- function(x){
  # set up train control for 10 times cross validation and random search of mtry tune parameters
  para_control <- caret::trainControl(method='cv',
                                      number=10,
                                      search = 'random',
                                      verbose=T)

  cat("Tune and train Random Forest model with caret\n")

  #Random generate mtry values with tuneLength = 10
  set.seed(100)
  model_rf <- caret::train(RT ~ .,
                           data = x,
                           method = 'rf',
                           metric = 'Rsquared',
                           tuneLength  = 10,
                           trControl = para_control,
                           importance=T,
                           allowParallel=T)



  cat("Model training is completed!\n")

  return(model_rf)
}



################################################################################
#   annotateNeighbor -------------------------------------------------------------
#' @title annotateNeighbor
#' @description Annotate peak table from seeds.
#' @author Zhiwei Zhou, Xiaotao Shen
#' @param data The data (csv) from ms2Annotation.
#' @param ms2 The ms2 data of peak table.
#' @param adduct.table Adduct table.
#' @param max.isotope The number of isotope peaks.
#' @param polarity The polarity.
#' @param mz.tol The mz tolerance.
#' @param rt.tol1 The RT tolerance for isotope and adduct annotation. (second)
#' @param rt.tol2 The RT tolerance for metabolite annotation. (\%)
#' @param cor.tol The cor tolerance.
#' @param int.tol The intensity ratio tolerance.
#' @param dp.tol The tolerance of dot product.
#' @param max.step The max number of reaction step.
#' @param remain Remain some seeds as validation or not.
#' @param remain.per The percentage of remained seeds.
#' @param metabolite The kegg compound database.
#' @param metabolic.network kegg.rpair2.
#' @param inHouse.compound The inhouse compound database with predicted RT.
#' @param threads The number of threads.
#' @param score.cutoff Score cutoff.
#' @param path The directory.
#' @param output.path The directory of outputing results.
#' @return  The tags data with annotation result.

# data <- "ms2_match_annotation_result.csv"
# ##ms2_data <- "ms2"
# # path <- 'H:/00_projects/03_MetDNA2/00_data/20200922_metdna2_development/pos_200927_metdna2'
# # path <- 'E:/zhiwei_create/02_MetDNA2/01_data_processing/20201005_metdna2_development/pos_201012_metdna2/'
# path <- '/home/zhouzw/Data_processing/20210224_metdna2_update/'
# load(file.path(path, "01_result_initial_seed_annotation", '00_intermediate_data', 'ms2'))
#
# ##adduct.table
# max_isotope = 4
# polarity = "positive"
# mz_tol = 25
# rt_tol1 = 3
# rt_tol2 = 30
# cor_tol = 0
# int_tol = 500
# dp_tol = 0.5
# max_step = 3
# remain = FALSE
# remain_per = 0.5
# threads = 3
# score_cutoff = 0
# rt_filter = FALSE # Not available
# seed_neighbor_match_plot = FALSE
# candidate_num = 3

setGeneric(name = "annotateNeighbor",
           def = function(data = "ms2_match_annotation_result.csv",
                          ms2,
                          adduct_table,
                          max_isotope = 4,
                          polarity = c("positive", "negative"),
                          metdna_version = c('version2', 'version1'),
                          instrument = c("SciexTripleTOF", "AgilentQTOF", "BrukerQTOF", "ThermoOrbitrap",
                                         'ThermoExploris', "AgilentDTIMMS", "BrukerTIMS", 'WatersQTOF','WatersTWIMMS'),
                          method_lc = c('Amide12min', 'Amide23min', 'Other', 'MetlinRP', 'RP12min'),
                          direction = c('forward', 'reverse'),
                          scoring_approach = c('dp', 'bonanza', 'hybrid', 'gnps'),
                          mz_tol = 25,
                          mz_ppm_thr = 400,
                          rt_tol1 = 3,
                          rt_tol2 = 30,
                          cor_tol = 0,
                          int_tol = 500,
                          ccs_tol = 4, # %
                          dp_tol = 0.5,
                          matched_frag_tol = 1,
                          whether_link_frag = FALSE,
                          max_step = 3,
                          remain = FALSE,
                          remain_per = 0.5,
                          metabolite,
                          metabolite_ccs,
                          metabolic_network,
                          inHouse_compound,
                          threads = 3,
                          score_cutoff = 0,
                          path = ".",
                          # output.path = ".",
                          rt_filter = TRUE, # Not available yet, need to be fix
                          seed_neighbor_match_plot = FALSE,
                          candidate_num = 3,
                          test_evaluation = c('No', '200STD', '46STD'),
                          whether_use_predRT = TRUE,
                          ...){
             # browser()


             # dir.create(output.path)
             path_output <- file.path(path, "02_result_MRN_annotation")
             dir.create(file.path(path_output, "00_intermediate_data"), showWarnings = FALSE, recursive = TRUE)
             polarity <- match.arg(polarity)

             # export MS2 match parameter
             switch (metdna_version,
                     'version1' = {
                       rt_filter <- TRUE
                       mz_tol_ms2 <- 30
                     },
                     'version2' = {
                       if (method_lc %in% c('Other', 'MetlinRP')) {
                         rt_filter <- TRUE
                       } else {
                         rt_filter <- FALSE
                       }

                       mz_tol_ms2 <- 25
                     }
             )

             # cat("\n");cat("Read annotation result from ms2Annotation.\n")
             data <- readInitialAnnotation(data = data,
                                           path = file.path(path, "01_result_initial_seed_annotation"),
                                           metdna_version = metdna_version,
                                           direction = direction,
                                           rt_filter = rt_filter,
                                           rt_tol = rt_tol2,
                                           inHouse_compound = inHouse_compound,
                                           instrument = instrument,
                                           ...)

             # data <- readInitialAnnotation(data = data,
             #                               path = file.path(path, "01_result_initial_seed_annotation"),
             #                               metdna_version = metdna_version,
             #                               direction = direction,
             #                               rt_filter = rt_filter,
             #                               rt_tol = rt_tol2,
             #                               inHouse_compound = inHouse_compound,
             #                               instrument = instrument)

             sample <- data[[2]]
             peak_int <- apply(sample, 1, median) # for matchMs2WithNeighbor (metIden)
             rm(list = "sample");gc()

             tags <- data[[1]]
             if (test_evaluation == '200STD') {tags$KEGG.ID <- tags$labid}
             if (test_evaluation == '46STD') {tags$KEGG.ID <- tags$labid}
             # tags[which(tags=="NA", arr.ind = TRUE)] <- NA
             rm(list = "data");gc()

             ms2_name <- unname(unlist(lapply(ms2, function(x) {x[[1]]["NAME",]})))

             # tags is transformed in to tags2 type data
             cat("\n")
             cat("Transform tags to tags2.\n")
             # use the tags2 data if it is in the directory
             if (any(dir(path) == "tags2")) {
               load(file.path(path, "tags2"))
               cat("=============================IMPORTANT NOTE=========================\n")
               cat("===                      Use the tags2 from local               ====\n")
               cat("====================================================================\n")
             } else {
               tags2 <- changeTags(tags = tags,
                                   ms2 = ms2,
                                   weight_mz = 0.25,
                                   weight_rt = 0.25,
                                   weight_dp = 0.5,
                                   weight_ccs = 0,
                                   mz_tol = 25,
                                   rt_tol = 30,
                                   ccs_tol = 4,
                                   dp_tol = 0.5)
             }

             rm(list = c("tags"))
             gc()
             #save raw tags2 data in local
             # raw.tags2 <- tags2
             # save(raw.tags2, file = file.path(path_output, "raw.tags2"))
             # rm(list = "raw.tags2")
             ##remain some seeds as validation dataset
             if (remain) {
               cat("\n")
               cat(remain_per*100,"% of seeds are remained as validation data.\n")
               idx <-
                 which(unlist(lapply(tags2, function(x) {length(x@annotation)})) > 0)
               save(idx, file = file.path(path_output,"idx"), compress = "xz")
               if (any(dir(path_output)=="remain.idx")) {
                 load(file.path(path_output, "remain.idx"))
               } else {
                 remain.idx <- lapply(1:30, function(x) sort(sample(idx, round(remain_per*length(idx)))))
                 remain.idx <- remain.idx[[sample(1:30, 1)]]
                 # remain.idx <- idx[sample(1:length(idx), round(remain_per@*length(idx)))]
                 save(remain.idx, file = file.path(path_output,"remain.idx"), compress = "xz")
               }

               ##remove annotation of reamin idx in tags2
               tags2[remain.idx] <- lapply(tags2[remain.idx], function(x){
                 x@annotation <- list()
                 x
               })
             }

             # get the initinal seeds for annotation
             seed_idx <- which(unlist(lapply(tags2, function(x) {length(x@annotation)})) > 0)

             # add seed information to tags2
             tags2 <- addSeedLabel(tags2 = tags2,
                                   round = 1,
                                   seed_idx = seed_idx,
                                   score_cutoff = score_cutoff)

             seed_annotation_number <- sum(unlist(showTags2(tags2, slot = "as.seed")[[1]]))
             seed_peak_number <- length(seed_idx)


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

             dir.create(file.path(path_output, "00_intermediate_data"), showWarnings = FALSE, recursive = TRUE)
             save(matchParam, file = file.path(path_output, "00_intermediate_data", 'matchParam_recursive'), version = 2)
             rm(list = c('intensityNormedMethod', 'methodScore'));gc()


             # first round is 1

             if (any(dir(file.path(path_output, "00_intermediate_data")) == "tags2_after_annotation")){
               cat('Note: load existed tags2_after_annotation, skip recursive annotation\n')
               load(file.path(path_output,"00_intermediate_data/tags2_after_annotation"))
               new_tags <- tags2_after_annotation
               rm(list = c("tags2_after_annotation"))
             } else {
               cat('Start recursive annotation using metabolic reaction network\n')
               round <- 1
               while (length(seed_idx) > 0){
                 cat("\n")
                 cat("Round", round, "annotation\n")
                 cat("Seed peak number:", seed_peak_number)
                 cat("\n")
                 cat("Seed annotation number:", seed_annotation_number)
                 cat("\n")
                 cat("Memory used: ")
                 temp <- pryr::mem_used()
                 print(temp)
                 cat("-------------------------------------------------------------------\n")

                 # Ronund 1 annotation
                 # parameter preparation for matchMs2WithNeighbor

                 new_tags <- matchMs2WithNeighbor(
                   peak_int = peak_int,
                   tags2 = tags2,
                   seed_idx = seed_idx,
                   max_isotope = max_isotope,
                   polarity = polarity,
                   rt_tol1 = rt_tol1,
                   rt_tol2 = rt_tol2,
                   mz_tol = mz_tol,
                   mz_ppm_thr = mz_ppm_thr,
                   cor_tol = cor_tol,
                   ccs_tol = ccs_tol,
                   int_tol = int_tol,
                   dp_tol = dp_tol,
                   matched_frag_tol = matched_frag_tol,
                   whether_link_frag = whether_link_frag,
                   metabolite = metabolite,
                   metabolite_ccs = metabolite_ccs,
                   metabolic_network = metabolic_network,
                   scoring_approach = scoring_approach,
                   max_step = max_step,
                   mz_tol_ms2 = mz_tol_ms2,
                   ms2 = ms2,
                   round = round,
                   remove_index = NULL,
                   iso_annotation = TRUE,
                   add_annotation = TRUE,
                   met_annotation = TRUE,
                   threads = threads,
                   adduct_table = adduct_table,
                   instrument = instrument,
                   whether_use_predRT = whether_use_predRT,
                   ...
                 )

                 # if (round == 1) {
                 #   browser()
                 # }

                 # cat('finish matchMs2WithNeighbor\n\n')
                 # remove the metAnnotation which is from the peak itself
                 annotation_type <- showTags2(new_tags, slot = "annotation.type")
                 annotation_idx <- lapply(annotation_type, function(x){
                   which(unlist(x) == "metAnnotation")
                 })

                 annotation_idx <- annotation_idx[which(unlist(lapply(annotation_idx, length))!=0)]

                 peak_index <- as.numeric(names(annotation_idx))

                 for (i in 1:length(peak_index)){
                   # cat("\n");cat("i:");cat(i);cat("\n")
                   temp_tags2 <- new_tags[[peak_index[i]]]
                   temp_annotation_idx <- annotation_idx[[i]]
                   for (j in temp_annotation_idx) {
                     # cat(j);cat(" ")
                     temp_annotation <- temp_tags2@annotation[[j]]
                     if (temp_annotation$From.peak ==  temp_tags2@name & temp_annotation$type == "metAnnotation"){
                       temp_tags2@annotation[[j]] <- list()
                       # count <- c(count, list(c(i,j)))
                     } else {
                       next()
                     }
                   }
                   temp_tags2@annotation <- temp_tags2@annotation[unlist(lapply(temp_tags2@annotation, length)) !=0]
                   if(length(temp_tags2@annotation) == 0){temp_tags2@annotation <- list()}
                   new_tags[[peak_index[i]]] <- temp_tags2
                 }

                 rm(list = c("annotation_type", "annotation_idx", "peak_index"))


                 # find seed for next round annotation
                 # round should be added 1
                 round <- round + 1

                 # find new seed index from new tags for next round annotation
                 annotation_idx <- which(showTags2(new_tags, slot = "annotation.len") > 0)

                 # temp_idx is index of seeds which have at least one annotation that is not seed before
                 # temp_idx <- which(unlist(lapply(showTags2(new_tags[annotation_idx], slot = "as.seed")[[1]], function(x) {any(!x)})))
                 temp_idx <- showTags2(new_tags[annotation_idx], slot = "as.seed")[[1]] %>%
                   lapply(., function(x){any(!x)}) %>%
                   unlist() %>%
                   which() %>%
                   unname()

                 seed_idx_new <- annotation_idx[temp_idx] # possible peak idx

                 # evaluate each new annotation with seed criteria:
                 #     1, not isotope; 2, not seed before; 3, score bigger than cutoff

                 idx_len <- unlist(lapply(seq_along(seed_idx_new), function(i) {
                   temp_idx <- seed_idx_new[i]
                   x <- new_tags[[temp_idx]]
                   annotation <- x@annotation
                   annotation_isotope <- unlist(lapply(annotation, function(x) {x$isotope}))
                   annotation_id <- unlist(lapply(annotation, function(x) {x$to}))
                   annotation_as_seed <- unlist(lapply(annotation, function(x) {x$as.seed}))
                   annotation_score <- unlist(lapply(annotation, function(x) {x$score}))
                   temp_idx <-
                     which(annotation_isotope == "[M]" &
                             !duplicated(annotation_id) &
                             !annotation_as_seed &
                             annotation_score >= score_cutoff)
                   length(temp_idx)
                 }))

                 seed_idx_new <- seed_idx_new[which(idx_len > 0)]

                 # if no new annotation obtained, breaks recurisive
                 if(length(seed_idx_new) == 0) {break()}

                 # add seed information to tags2
                 new_tags <- addSeedLabel(tags2 = new_tags,
                                          round = round,
                                          seed_idx = seed_idx_new,
                                          score_cutoff = score_cutoff)

                 tags2 <- new_tags
                 rm(list = "new_tags")
                 gc()

                 # save(tags2, file = file.path('H:/00_projects/03_MetDNA2/00_data/20200922_metdna2_development/pos_200927_metdna2/test',
                 #                              paste0('tags_round', round, '.RData')))

                 dir.create(file.path(path, 'test'),
                            showWarnings = FALSE, recursive = TRUE)
                 save(tags2,
                      file = file.path(path, 'test',
                                       paste0('tags_round_', round, '_210304.RData')))

                 # need be fixed
                 seed_annotation_number <-
                   unlist(showTags2(tags2, slot = "as.seed")[[2]])
                 seed_annotation_number <- seed_annotation_number[!is.na(seed_annotation_number)]
                 seed_annotation_number <- sum(seed_annotation_number == round)
                 seed_peak_number <- length(seed_idx_new)

                 seed_idx <- seed_idx_new

                 dir.create(file.path(path, 'test'),
                            showWarnings = FALSE, recursive = TRUE)
                 save(seed_idx, file = file.path(path, 'test',
                                                 paste0('seed_idx', round, '_210304.RData')))
               }


               # --------------------------------------------------------------------------

               # order annotation for each peak according to score
               new_tags <-
                 lapply(new_tags,
                        function(x) {filterPeak(object = x, score_thr = 0)})

               # save new_tags2 after annotation
               tags2_after_annotation <- new_tags

               save(tags2_after_annotation,
                    file = file.path(path_output, "00_intermediate_data","tags2_after_annotation"),
                    compress = "xz", version = 2)
               rm(tags2_after_annotation)
               gc()
             }

             # begin to give confidence to annotation of peaks
             result_matrix_annotation <- trans2Matrix(tags2 = new_tags, base = "peak")

             # assign confidence level of each ID -----------------------------------------
             cat('\n'); cat("Assign confidence.\n")
             unique_id <- unique(result_matrix_annotation$to)
             pbapply::pboptions(type="timer", style = 1)


             # for each unique id:
             #    1. group all annotation into peak group: RT = 3 (with same annotation)
             #    2. if seed exist, only keep the 1st annotation records
             #    3. if not seed annotation, only keep the highest score records
             # (Note: the metAnnotation records with highest score is prior to reserve after version 0.6.62)

             id_result <- lapply(seq_along(unique_id), function(i) {
               temp_id <- unique_id[i]
               temp_idx <- which(result_matrix_annotation$to == temp_id)
               temp_result <- result_matrix_annotation[temp_idx,]
               temp_rt <- temp_result$rt

               # group peaks according to RT
               switch(metdna_version,
                      'version1' = {
                        rt_class <- groupRT(rt = temp_rt, rt.tol = rt_tol1)
                      },
                      'version2' = {
                        rt_class <- groupRT2(rt_vec = temp_rt, rt_tol = rt_tol1)
                      })
               # rt_class <- groupRT(rt = temp_rt, rt.tol = rt_tol1)
               temp_result <- lapply(rt_class, function(x) {temp_result[x,]})

               temp_result <-  lapply(seq_along(temp_result), function(i) {
                 x <- temp_result[[i]]
                 x1 <- data.frame()
                 unique.peak.name <- unique(x$name)

                 # "Note": the number of annotation for one peak
                 for (j in seq_along(unique.peak.name)) {
                   temp_idx <- which(x$name == unique.peak.name[j])
                   if (length(temp_idx) == 1) {
                     Note <- 1
                     x1 <- rbind(x1, data.frame(x[temp_idx,], Note, stringsAsFactors = FALSE))
                     # colnames(x1)[ncol(x1)] <- "Note"
                   } else {
                     temp.x <- x[temp_idx,]
                     if(any(temp.x$type == "seed")){
                       temp.x <- temp.x[which(temp.x$type == "seed")[1],]
                     } else {
                       # temp.x <- temp.x[order(temp.x$score, decreasing = TRUE),]
                       # temp.x <- temp.x[1,]

                       switch(metdna_version,
                              'version1' = {
                                temp.x <- temp.x[order(temp.x$score, decreasing = TRUE),]
                                temp.x <- temp.x[1,]
                              },
                              'version2' = {
                                # modify order to keep annotation records
                                #    isotopeAnnotation > metAnnotation > adductAnnotation
                                temp.x <- temp.x %>%
                                  dplyr::arrange(match(type, c('metAnnotation', 'isotopeAnnotation', 'adductAnnotation')),
                                                 dplyr::desc(score))
                                temp.x <- temp.x[1,]
                              }
                       )
                     }
                     Note <- length(temp_idx)
                     x1 <- rbind(x1, data.frame(temp.x, Note, stringsAsFactors = FALSE))
                     # colnames(x1)[ncol(x1)] <- "Note"
                   }
                 }

                 colnames(x1)[ncol(x1)] <- "Note"
                 x1
               })

               temp_confidence <- unlist(lapply(temp_result, function(x) assignConfidence(x, polarity = polarity)))

               temp_result <- mapply(function(x,
                                              Confidence,
                                              group) {
                 result <-
                   data.frame(x,
                              Confidence,
                              paste(x$to, group, sep = "_"),
                              stringsAsFactors = FALSE)
                 colnames(result)[ncol(result)] <- "group"
                 list(result)
               },
               x = temp_result,
               Confidence = temp_confidence,
               group = c(1:length(temp_result)))

               temp_result <- do.call(rbind, temp_result)
               temp_result
             })

             id_result_with_confidence <- do.call(rbind, id_result)

             # remove redundancy ----------------------------------------------------------
             cat("\n");cat("Remove redundancy.\n")

             id_result_redun_rm <- removeRedundancy(result = id_result_with_confidence,
                                                    path = path_output,
                                                    polarity = polarity)

             load(file.path(path_output, "redun"))
             record_redun1 <- redun
             save(record_redun1, file = file.path(path_output, "00_intermediate_data","record_redun1"), compress = "xz")
             unlink(x = file.path(path_output, "redun"), recursive = TRUE)
             save(id_result_redun_rm, file = file.path(path_output, "00_intermediate_data","id_result_redun_rm"), compress = "xz")

             # filter new_tags according to annotation.result
             new_tags <- convertResultMatrix2Tags(result = id_result_redun_rm, tags2 = new_tags)
             tags2_after_redundancy_remove <- new_tags
             save(tags2_after_redundancy_remove,
                  file = file.path(path_output, "00_intermediate_data","tags2_after_redundancy_remove"), compress = "xz")

             # change tags2_after_redundancy_remove to csv like ms2Annotation
             cat("\n")
             # data("cpd_kegg", envir = environment())

             temp <- getAnnotationResult(tags2 = tags2_after_redundancy_remove,
                                         tags.result2 = id_result_redun_rm,
                                         metabolite = metabolite,
                                         candidate.num = candidate_num,
                                         instrument = instrument,
                                         score.cutoff = 0.4)

             cat("\n")
             readr::write_csv(temp, file.path(path_output, "MRN.annotation.result.csv"))

             # output MS/MS matching spectra -----------------------------------
             if(seed_neighbor_match_plot){
               cat('Plot spectra match results (Seed vs Neighbor).\n')
               # data("kegg.compound", envir = environment())
               # load(file.path(path_output, "00_intermediate_data", "ms2"))

               plotNeighbor(tags2 = tags2_after_redundancy_remove,
                            path = path_output,
                            ms2 = ms2,
                            kegg.compound = metabolite,
                            matchParam = matchParam)

               cat("\n")
             }


             rm(list = c("tags2.after.redundancy.remove", "new.tags", "temp"))
             gc()
           })


################################################################################
#     changeTags -------------------------------------------------------------------
#' @title changeTags: Chages tags (matrix or data.frame) to S4 class peakInfo
#' @description changeTages is used to change tags information (matrix or data frame) to peakInfo (S4 class).
#' @author Zhiwei Zhou, Xiaotao Shen
#' @param tags The tags information.
#' @param ms2 The ms2.
#' @param adduct_table The adduct table.
#' @param weight_mz mz weight.
#' @param weight_rt RT weight
#' @param weight_dp Dop product weight
#' @param mz_tol The mz tolerance.
#' @param rt_tol The RT tolerance.
#' @param dp_tol The dp tolerance.
#' @param ... other parameters.
#' @return Tags2 data.
#' @export


setGeneric(name = "changeTags",
           def = function(tags,
                          ms2,
                          # adduct_table,
                          weight_mz = 0.25,
                          weight_rt = 0.25,
                          weight_ccs = 0,
                          weight_dp = 0.5,
                          mz_tol = 25,
                          rt_tol = 30,
                          dp_tol = 0.5,
                          ccs_tol = 4,
                          ...) {
             #score formula
             ##calculate score for peak, using mz error, rt error
             ##and int error
             temp.x <- c(0, mz_tol)
             temp.y <- c(1, 0)
             lm.reg.mz <- lm(temp.y~temp.x)

             temp.x <- c(0, rt_tol)
             temp.y <- c(1, 0)
             lm.reg.rt <- lm(temp.y~temp.x)

             temp.x <- c(0, ccs_tol)
             temp.y <- c(1, 0)
             lm.reg.ccs <- lm(temp.y~temp.x)

             temp.x <- c(dp_tol, 1)
             temp.y <- c(0, 1)
             lm.reg.dp <- lm(temp.y~temp.x)

             # load ms2 information
             if(!is.null(ms2)){
               ms2 <- ms2
               ms2.name <-
                 unname(unlist(lapply(ms2, function(x) {x[[1]]["NAME",]})))
             }

             # change space and nan to NA
             # tags[which(tags == "", arr.ind = TRUE)] <- NA
             # tags[which(is.nan(tags), arr.ind = TRUE)] <- NA
             # change tags to a list, a element is a peak.

             tags1 <- apply(tags, 1, list)

             # change tags1 to peakInfo S4 class
             tags2 <- lapply(tags1, function(x){
               temp <- new(Class = "PeakInfo",
                           name = stringr::str_trim(as.character(x[[1]]["Peak.name"])),
                           mz = as.numeric(x[[1]]["mz"]),
                           rt = as.numeric(x[[1]]["rt"]),
                           ccs = as.numeric(x[[1]]["ccs"]),
                           ##add a empty ms2
                           ms2 = data.frame(),
                           annotation = list())

               # add seed annotation information
               if(!is.na(x[[1]][["KEGG.ID"]])){
                 kegg.id <- x[[1]]["KEGG.ID"]
                 kegg.id <- strsplit(kegg.id, split = ";")[[1]]
                 kegg.id[kegg.id == "NA"] <- NA
                 if(any(!is.na(kegg.id))){
                   adduct <- strsplit(x[[1]]["adduct"], split = ";")[[1]]
                   ms2.sim <-
                     as.numeric(strsplit(x[[1]]["ms2.sim"], split = ";")[[1]])
                   nfrag <-
                     as.numeric(strsplit(x[[1]]["nfrag"], split = ";")[[1]]) # number of matched fragments

                   rt.error <-
                     as.numeric(strsplit(x[[1]]["rt.error"], split = ";")[[1]])
                   # for the identification whose rt error is NA, make it as rt_tol
                   rt.error[is.na(rt.error)] <- rt_tol

                   ccs.error <-
                     as.numeric(strsplit(x[[1]]["ccs.error"], split = ";")[[1]])
                   # for the identification whose ccs error is NA, make it as ccs_tol
                   ccs.error[is.na(ccs.error)] <- ccs_tol

                   mz.error <- as.numeric(strsplit(x[[1]]["mz.error"], split = ";")[[1]])
                   formula <- strsplit(x[[1]]["Formula"], split = ";")[[1]]
                   temp.idx <- which(!is.na(kegg.id))

                   kegg.id <- kegg.id[temp.idx]
                   adduct <- adduct[temp.idx]
                   ms2.sim <- ms2.sim[temp.idx]
                   nfrag <- nfrag[temp.idx]
                   rt.error <- rt.error[temp.idx]
                   mz.error <- mz.error[temp.idx]
                   ccs.error <- ccs.error[temp.idx]
                   formula <- formula[temp.idx]


                   # add annotation in KEGG
                   for(i in 1:length(kegg.id)){
                     score1 <- coefficients(lm.reg.mz)[1] +
                       coefficients(lm.reg.mz)[2] * mz.error[i]
                     score2 <- coefficients(lm.reg.rt)[1] +
                       coefficients(lm.reg.rt)[2] * rt.error[i]
                     score3 <- coefficients(lm.reg.dp)[1] +
                       coefficients(lm.reg.dp)[2] * ms2.sim[i]
                     score4 <- coefficients(lm.reg.ccs)[1] +
                       coefficients(lm.reg.ccs)[2] * ccs.error[i]

                     score <-
                       score1*weight_mz + score2*weight_rt  + score4*weight_ccs + score3*weight_dp

                     temp@annotation[[i]] <-
                       list(type = "seed",
                            From = NA,
                            From.peak = NA,
                            to = stringr::str_trim(kegg.id[i]),
                            step = NA,
                            level = 1,
                            as.seed = FALSE,
                            as.seed.round = NA,
                            isotope = "[M]",
                            adduct = stringr::str_trim(adduct[i]),
                            charge = 1,
                            # as.numeric(adduct_table$charge[match(adduct[i], adduct_table$name)]),
                            formula = stringr::str_trim(formula[i]),
                            mz.error = as.numeric(mz.error[i]),
                            rt.error = as.numeric(rt.error[i]),
                            ccs.error = as.numeric(ccs.error[i]),
                            int.error = NA,
                            ms2.sim = as.numeric(ms2.sim[i]),
                            nfrag = as.numeric(nfrag[i]),
                            score = unname(score)
                       )
                   }


                 }
               }
               temp
             })

             # if has ms2, add ms2 information
             if(!is.null(ms2)){
               tags2 <- lapply(tags2, function(x){
                 x@ms2 <- data.frame(ms2[[match(x@name, ms2.name)]]$spec,
                                     stringsAsFactors = FALSE)
                 x
               })
             }


             # if score of annotation is less than 0, then change it to 0.001
             tags2 <- lapply(tags2, function(x){
               if(length(x@annotation) == 0) return(x)
               annotation <- x@annotation
               annotation <- lapply(annotation, function(y){
                 if(y$score < 0) y$score <- 0.001
                 y
               })
               x@annotation <- annotation
               x
             })

             tags2 <- tags2
           }
)


# test <- lapply(seq_along(tags2), function(i){
#   cat(i, ' ')
#   x <- tags2[[i]]
#   if (length(x@annotation) == 0) return(NULL)
#   annotation <- x@annotation
#   annotation <- lapply(seq_along(annotation), function(j){
#     cat(j, ' ')
#     y <- annotation[[j]]
#     result <- tibble::tibble(feature = x@name,
#                              id = y$to,
#                              score = y$score)
#
#     result
#   })
#
#   annotation <- annotation %>% dplyr::bind_rows()
#
#   annotation
# })


#     addSeedLabel ---------------------------------------------------------------------
#' @title addSeedLabel
#' @description Label seed in tags2 data.
#' @author Zhiwei Zhou, Xiaotao Shen
#' @param tags2 The tags2 data.
#' @param round The annotation round.
#' @param seed_idx The index of seeds.
#' @param score_cutoff The score cutoff.
#' @return The tags2 data with labeled seeds.

setGeneric(name = "addSeedLabel",
           def = function(tags2,
                          round,
                          seed_idx,
                          score_cutoff = 0,
                          matched_frag_tol = 1){

             # For each possible metabolite (annotation), select seed with 4 creteria:
             #   1. remain monoisotope, 2. unique KEGG ID (for feature), 3. never as seed, 4. score >= cutoff

             tags2[seed_idx] <- lapply(tags2[seed_idx], function(x) {
               annotation <- x@annotation
               annotation.isotope <-
                 unlist(lapply(annotation, function(x) {x$isotope}))
               annotation.id <-
                 unlist(lapply(annotation, function(x) {x$to}))
               annotation.as.seed <-
                 unlist(lapply(annotation, function(x) {x$as.seed}))
               annotation.score <-
                 unlist(lapply(annotation, function(x) {x$score}))
               annotation.nfrag <-
                 unlist(lapply(annotation, function(x) {x$nfrag}))
               temp.idx <-
                 which(annotation.isotope == "[M]" &
                         !duplicated(annotation.id) &
                         !annotation.as.seed &
                         annotation.score >= score_cutoff &
                         (annotation.nfrag >= matched_frag_tol | is.na(annotation.nfrag)))
               if (length(temp.idx) != 0) {
                 annotation[temp.idx] <-
                   lapply(annotation[temp.idx], function(x) {x$as.seed <- TRUE;x})
                 annotation[temp.idx] <-
                   lapply(annotation[temp.idx], function(x) {x$as.seed.round <- round;x})
                 x@annotation <- annotation
               } else {
                 x@annotation <- annotation
               }
               x
             })
             return(tags2)
           })


################################################################################
#     matchMs2WithNeighbor ---------------------------------------------------------

#' @title matchMs2WithNeighbor
#' @description Annotate peak table from seeds.
#' @author Zhiwei Zhou, Xiaotao Shen
#' @param peak_int The median/mean intensity of peaks.
#' @param tags2 The tags2 data of all peaks.
#' @param seed_idx The index of seeds.
#' @param max_isotope The number of isotope peaks.
#' @param polarity The polarity.
#' @param rt_tol1 The RT tolerance for isotope and adduct annotation. (second)
#' @param rt_tol2 The RT tolerance for metabolite annotation (\%).
#' @param mz_tol The mz tolerance.
#' @param cor_tol The cor tolerance.
#' @param int_tol The intensity ratio tolerance.
#' @param dp_tol The tolerance of dot product.
#' @param metabolite The kegg compound database.
#' @param metabolic_network kegg.rpair2
#' @param max_step The max number of reaction step.
#' @param ms2 The ms2 data of peak table.
#' @param round The round of annotation.
#' @param remove_index The index of peaks which you dont want to annotate.
#' @param iso_annotation Isotope annotation or not.
#' @param add_annotation Adduct annotation or not.
#' @param met_annotation Metabolite annotation or not.
#' @param threads The number of threads.
#' @param adduct_table Adduct table.
#' @return  The tags data with annotation result.

# setwd("/home/jasper/work/MetDNA/metIden")
# load("peak_int")
# load("tags2")
# load("idx")
# load("kegg.compound.rda")
# metabolite <- kegg.compound
# load("kegg.rpair2.rda")
# metabolic_network <- kegg.rpair2
# load('adduct_table')
# load("ms2")
#
# max_isotope = 4
# polarity = "positive"
# #rt_tol1 is absolute
# rt_tol1 = 3
# #rt_tol2 is relative
# rt_tol2 = 30
# mz_tol = 25
# cor_tol = 0
# #for isotope annotation
# int_tol = 500
# dp_tol = 0.5
# max_step = 3
# round <- 1
# # cor.matrix,
# remove_index = NULL
# iso_annotation = TRUE
# add_annotation = TRUE
# met_annotation = TRUE
# threads = 3
# seed_idx <- seed_idx
# adduct_table <- convertAdductList2AdductTable(polarity = 'positive',
#                                               column = 'hilic',
#                                               adduct_type = c('Common', 'Rare'))
# metabolite = cpd_kegg %>%
#   dplyr::rename('RT' = Amide23.RT,
#                 'attribute' = Amide23.RT.attribute) %>%
#   dplyr::mutate(RT = RT*60)
#
# # load('D:/01_r_package/MetDNA_V1.21/data/kegg.rpair2.rda')
# metabolic_network = kegg.rpair2

setGeneric(name = "matchMs2WithNeighbor",
           def = function(peak_int,
                          tags2,
                          seed_idx,
                          ms2,
                          round,
                          metabolite,
                          metabolite_ccs,
                          metabolic_network,
                          max_isotope = 4,
                          polarity = c("positive", "negative"),
                          instrument = c("SciexTripleTOF", "AgilentQTOF", "BrukerQTOF", "ThermoOrbitrap",
                                         'ThermoExploris', "AgilentDTIMMS", "BrukerTIMS", 'WatersQTOF','WatersTWIMMS'),
                          scoring_approach = c('dp', 'bonanza', 'hybrid', 'gnps'),
                          rt_tol1 = 3, # rt_tol1 is absolute
                          rt_tol2 = 30, # rt_tol2 is relative
                          mz_tol = 25,
                          mz_ppm_thr = 400,
                          cor_tol = 0,
                          int_tol = 500, # for isotope annotation
                          ccs_tol = 4, # for metaboliteAnnotation, %
                          dp_tol = 0.5,
                          matched_frag_tol = 1,
                          whether_link_frag = FALSE, # whehter consider propagation with matched fragment number
                          max_step = 3,
                          mz_tol_ms2 = 25, # fragment match tolerance
                          # cor.matrix,
                          remove_index = NULL,
                          iso_annotation = TRUE,
                          add_annotation = TRUE,
                          met_annotation = TRUE,
                          threads = 3,
                          whether_use_predRT = TRUE,
                          adduct_table,
                          ...){
             polarity <- match.arg(polarity)
             # browser()
             # path <- '/home/zhouzw/Data_processing/20210224_metdna2_update/'
             # scoring_approach <- match.arg(scoring_approach)

             # ms2 is the ms2 data base, ms2.name is the name of ms2 spectra
             ms2_name <- unname(unlist(lapply(ms2, function(x) {x[[1]]["NAME",]})))

             # tags2 is the S4 class peakInfo for tags
             # peak_mz and peak_rt is the mz and rt of all peaks
             peak_name <- unname(unlist(lapply(tags2, function(x) x@name)))
             peak_mz <- unname(unlist(lapply(tags2, function(x) x@mz)))
             peak_rt <- unname(unlist(lapply(tags2, function(x) x@rt)))
             peak_ccs <- unname(unlist(lapply(tags2, function(x) x@ccs)))
             names(peak_mz) <- names(peak_rt) <- names(peak_rt) <- peak_name

             # isotope annotation ----------------------------------------------
             if (iso_annotation) {

               # save(seed_idx, file = "seed_idx")
               # save(tags2, file = "tags2")
               # save(annotation_idx, file = "annotation_idx")
               cat("\n")
               cat("Isotope annotation.\n")
               annotation_idx <- 1:length(tags2)
               if(!is.null(remove_index)) {
                 annotation_idx <- setdiff(annotation_idx, remove_index)
               }
               pbapply::pboptions(type="timer", style = 1)

               # annotate isotopes for each seed annotation
               isotope_result <- lapply(seq_along(seed_idx), function(i){
                 # cat(i); cat(' ')
                 temp_idx <- seed_idx[i]
                 seed <- tags2[[temp_idx]]
                 seed.as.seed.round <- showTags2(list(seed), slot = "as.seed")[[2]][[1]]

                 # select the specific annotation as seed with specific round
                 seed.i <- which(seed.as.seed.round == round)
                 if(length(seed.i)==0){
                   return(NULL)
                 }else{
                   temp_result <- lapply(seed.i,function(j){
                     query_name <- seed@name
                     query_id <- seed@annotation[[j]]$to
                     query_charge <- seed@annotation[[j]]$charge
                     query_level <- seed@annotation[[j]]$level
                     query_formula <- stringr::str_trim(seed@annotation[[j]]$formula)
                     query_adduct <- seed@annotation[[j]]$adduct
                     query_mz <- seed@mz
                     query_rt <- seed@rt
                     query_int <- peak_int[temp_idx]
                     peak_mz_all <- peak_mz[annotation_idx]
                     peak_rt_all <- peak_rt[annotation_idx]
                     peak_int_all <- peak_int[annotation_idx]

                     # cor_all <- cor.matrix[annotation_idx, anno.i]
                     cor_all <- rep(1, length(annotation_idx))
                     isotopes_result <- annotateIsotopeMRN(id = query_id,
                                                           formula = query_formula,
                                                           adduct = query_adduct,
                                                           charge = query_charge,
                                                           mz = query_mz,
                                                           rt = query_rt,
                                                           int = query_int,
                                                           peak_mz_all = peak_mz_all,
                                                           peak_rt_all = peak_rt_all,
                                                           peak_int_all = peak_int_all,
                                                           peak_cor_all = cor_all,
                                                           rt_tol = rt_tol1,
                                                           mz_tol = mz_tol,
                                                           mz_ppm_thr = mz_ppm_thr,
                                                           cor_tol = cor_tol,
                                                           int_tol = int_tol,
                                                           max_isotope = max_isotope)
                     isotopes_result
                   })
                   temp_result
                 }
               })

               # add isotope_result to tags2, specific to round
               anno_idx <- showTags2(tags2[seed_idx], slot = "as.seed")[[2]]
               anno_idx <- lapply(anno_idx, function(x) which(x == round))
               for(j in 1:length(seed_idx)){
                 # cat(j); cat(" ")
                 temp_result <- isotope_result[[j]]
                 for(k in 1:length(temp_result)){
                   tags2 <- addIsotopeResult2PeakInfo(isotopes.result = temp_result[[k]],
                                                      mz.tol = mz_tol,
                                                      rt.tol = rt_tol1,
                                                      int.tol = int_tol,
                                                      weight.mz = 0.45,
                                                      weight.rt = 0.45,
                                                      weight.int = 0.1,
                                                      peak.info = tags2,
                                                      peak.idx = seed_idx[j],
                                                      anno.idx = anno_idx[[j]][k],
                                                      annotation.idx = annotation_idx)
                 }

               }

               rm(list = "isotope.result")
               gc()
             }

             # browser()
             # dir.create(file.path(path, 'test'), showWarnings = FALSE, recursive = TRUE)
             # save(tags2,
             #      file = file.path(path, 'test', 'temp_tags2_isotope_201110.RData'))

             # dir.create(file.path(path, 'test'), showWarnings = FALSE, recursive = TRUE)
             # save(tags2,
             #      file = file.path('D:/project/00_zhulab/01_metdna2/00_data/20220412_test_metdna2', 'test', paste0('round', round, '_temp_tags2_isotope_220420.RData')),
             #      version = 2)

             # adduct annotation -----------------------------------------------
             if (add_annotation) {
               cat("\n");cat("Adduct annotation.\n")

               # # for each seed, annotate its adducts
               # 
               # temp_fun <- function(index, # chunks from seed_idx
               #                      tags2, peak_mz, peak_rt,
               #                      polarity, mz_tol, rt_tol, mz_ppm_thr,
               #                      adduct_table,
               #                      annotation_idx,
               #                      round){
               #   annotation_idx <- 1:length(tags2)
               #   library(Rcpp, warn.conflicts = FALSE, quietly = TRUE)
               #   # suppressMessages(data(thermo, package = "CHNOSZ", envir = environment()))
               #   temp_tags2 <- tags2[index]
               #   rm(tags2)
               #   gc()
               #   pbapply::pboptions(type="timer", style = 1)
               #   # result <- pbapply::pblapply(temp_tags2, function(x){
               #   result <- lapply(seq_along(temp_tags2), function(i){
               #     seed <- temp_tags2[[i]]
               #     seed.as.seed.round <-
               #       showTags2(list(seed), slot = "as.seed")[[2]][[1]]
               #     seed.i <- which(seed.as.seed.round == round)
               #     peak_mz_all <- peak_mz[annotation_idx]
               #     peak_rt_all <- peak_rt[annotation_idx]
               #     peak_cor_all <- rep(1, length(annotation_idx))
               #     rm(list=c("peak_mz", "peak_rt"))
               #     gc()
               #     if(length(seed.i)==0){
               #       return(NULL)
               #     }else{
               #       temp_result <- lapply(seed.i,function(j){
               #         query_name <- seed@name
               #         query_id <- seed@annotation[[j]]$to
               #         query_charge <- seed@annotation[[j]]$charge
               #         query_level <- seed@annotation[[j]]$level
               #         query_formula <- stringr::str_trim(seed@annotation[[j]]$formula)
               #         query_adduct <- seed@annotation[[j]]$adduct
               #         query_mz <- seed@mz
               #         query_rt <- seed@rt
               # 
               #         adduct_result <-
               #           annotateAdductMRN(id = query_id,
               #                             formula = query_formula,
               #                             adduct = query_adduct,
               #                             polarity = polarity,
               #                             mz = query_mz,
               #                             rt = query_rt,
               #                             peak_mz_all = peak_mz_all,
               #                             peak_rt_all = peak_rt_all,
               #                             rt_tol = rt_tol,
               #                             mz_tol = mz_tol,
               #                             mz_ppm_thr = mz_ppm_thr,
               #                             adduct_table = adduct_table,
               #                             cor_tol = cor_tol)
               #         adduct_result
               #       })
               #       temp_result
               #     }
               #   })
               #   result
               # }
               # 
               # 
               # cl <- snow::makeCluster(threads, type = "SOCK")
               # 
               # snow::clusterExport(cl, list = list("showTags2",
               #                                     "annotateAdductMRN",
               #                                     "sumFormula",
               #                                     "checkElement",
               #                                     "splitFormula",
               #                                     "pasteElement"))
               # nc <- length(cl)
               # options(warn = -1)
               # if (length(seed_idx) < threads){
               #   threads <-  1
               # }
               # 
               # # split seed index into chunks according to threads
               # # ichunks pass to index
               # ichunks <- split(seed_idx, 1:threads)
               # # library(Rcpp)
               # adduct_result <-
               #   snow::clusterApply(cl = cl, ichunks,
               #                      fun = temp_fun,
               #                      tags2 = tags2,
               #                      peak_mz = peak_mz,
               #                      peak_rt = peak_rt,
               #                      polarity = polarity,
               #                      rt_tol = rt_tol1,
               #                      mz_tol = mz_tol,
               #                      mz_ppm_thr = mz_ppm_thr,
               #                      # cor_tol = cor_tol,
               #                      round = round,
               #                      adduct_table = adduct_table,
               #                      annotation_idx = annotation_idx)
               # 
               # snow::stopCluster(cl)
               # 
               # browser()
               # # combine adduct chunks to list
               # adduct_result1 <- adduct_result[[1]]
               # if(threads > 1){
               #   for(i in 2:length(adduct_result)){
               #     adduct_result1 <- c(adduct_result1, adduct_result[[i]])
               #   }
               # }
               # adduct_result1 <- adduct_result1[order(unlist(ichunks))]
               # rm(list = "adduct_result")
               # gc()

               tempFunAdduct <- function(k){
                 index <- seed_idx[k]
                 
                 annotation_idx <- 1:length(tags2)
                 # library(Rcpp, warn.conflicts = FALSE, quietly = TRUE)
                 # suppressMessages(data(thermo, package = "CHNOSZ", envir = environment()))
                 temp_tags2 <- tags2[index]
                 rm(tags2)
                 gc()
                 # pbapply::pboptions(type="timer", style = 1)
                 # result <- pbapply::pblapply(temp_tags2, function(x){
                 result <- lapply(seq_along(temp_tags2), function(i){
                   seed <- temp_tags2[[i]]
                   seed.as.seed.round <-
                     showTags2(list(seed), slot = "as.seed")[[2]][[1]]
                   seed.i <- which(seed.as.seed.round == round)
                   peak_mz_all <- peak_mz[annotation_idx]
                   peak_rt_all <- peak_rt[annotation_idx]
                   peak_cor_all <- rep(1, length(annotation_idx))
                   rm(list=c("peak_mz", "peak_rt"))
                   gc()
                   if(length(seed.i)==0){
                     return(NULL)
                   }else{
                     temp_result <- lapply(seed.i,function(j){
                       query_name <- seed@name
                       query_id <- seed@annotation[[j]]$to
                       query_charge <- seed@annotation[[j]]$charge
                       query_level <- seed@annotation[[j]]$level
                       query_formula <- stringr::str_trim(seed@annotation[[j]]$formula)
                       query_adduct <- seed@annotation[[j]]$adduct
                       query_mz <- seed@mz
                       query_rt <- seed@rt
                       
                       adduct_result <-
                         annotateAdductMRN(id = query_id,
                           formula = query_formula,
                           adduct = query_adduct,
                           polarity = polarity,
                           mz = query_mz,
                           rt = query_rt,
                           peak_mz_all = peak_mz_all,
                           peak_rt_all = peak_rt_all,
                           rt_tol = rt_tol1,
                           mz_tol = mz_tol,
                           mz_ppm_thr = mz_ppm_thr,
                           adduct_table = adduct_table,
                           cor_tol = cor_tol)
                       adduct_result
                     })
                     temp_result
                   }
                 })
                 result
                 
               }
               
               
               library(parallel)
               # cl <- makeCluster(3L)
               cl <- parallel::makeCluster(threads)
               parallel::clusterExport(cl, c('tempFunAdduct',
                 "checkElement",
                 "annotateAdductMRN",
                 "showTags2",
                 "seed_idx",
                 "tags2",
                 "peak_mz",
                 "peak_rt",
                 "polarity",
                 "rt_tol1",
                 'mz_tol',
                 'mz_ppm_thr',
                 'cor_tol',
                 'round',
                 'adduct_table',
                 'annotation_idx'),
                 envir = environment())
               
               system.time(
                 adduct_result <- parallel::parLapply(cl = cl,
                   seq_along(seed_idx),
                   function(k) {tempFunAdduct(k)}
                 )
               )
               stopCluster(cl)
               
               adduct_result1 <- lapply(adduct_result, function(x){
                 x[[1]]
               })
               rm(list = "adduct_result")
               gc()


               # # debugs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               # adduct_result <- pbapply::pblapply(seq_along(seed_idx), function(k){
               #   cat(k, ' ')
               #   index <- seed_idx[k]
               # 
               #   annotation_idx <- 1:length(tags2)
               #   # library(Rcpp, warn.conflicts = FALSE, quietly = TRUE)
               #   # suppressMessages(data(thermo, package = "CHNOSZ", envir = environment()))
               #   temp_tags2 <- tags2[index]
               #   rm(tags2)
               #   gc()
               #   pbapply::pboptions(type="timer", style = 1)
               #   # result <- pbapply::pblapply(temp_tags2, function(x){
               #   result <- lapply(seq_along(temp_tags2), function(i){
               #     seed <- temp_tags2[[i]]
               #     seed.as.seed.round <-
               #       showTags2(list(seed), slot = "as.seed")[[2]][[1]]
               #     seed.i <- which(seed.as.seed.round == round)
               #     peak_mz_all <- peak_mz[annotation_idx]
               #     peak_rt_all <- peak_rt[annotation_idx]
               #     peak_cor_all <- rep(1, length(annotation_idx))
               #     rm(list=c("peak_mz", "peak_rt"))
               #     gc()
               #     if(length(seed.i)==0){
               #       return(NULL)
               #     }else{
               #       temp_result <- lapply(seed.i,function(j){
               #         query_name <- seed@name
               #         query_id <- seed@annotation[[j]]$to
               #         query_charge <- seed@annotation[[j]]$charge
               #         query_level <- seed@annotation[[j]]$level
               #         query_formula <- stringr::str_trim(seed@annotation[[j]]$formula)
               #         query_adduct <- seed@annotation[[j]]$adduct
               #         query_mz <- seed@mz
               #         query_rt <- seed@rt
               #         
               #         adduct_result <-
               #           annotateAdductMRN(id = query_id,
               #             formula = query_formula,
               #             adduct = query_adduct,
               #             polarity = polarity,
               #             mz = query_mz,
               #             rt = query_rt,
               #             peak_mz_all = peak_mz_all,
               #             peak_rt_all = peak_rt_all,
               #             rt_tol = rt_tol1,
               #             mz_tol = mz_tol,
               #             mz_ppm_thr = mz_ppm_thr,
               #             adduct_table = adduct_table,
               #             cor_tol = cor_tol)
               #         adduct_result
               #       })
               #       temp_result
               #     }
               #   })
               #   result
               # })
               # 
               # 
               # adduct_result1 <- lapply(adduct_result, function(x){
               #   x[[1]]
               # })
               # rm(list = "adduct_result")
               # gc()

               for (j in 1:length(seed_idx)) {
                 # cat(j, ' ')
                 temp_result <- adduct_result1[[j]]
                 for (k in 1:length(temp_result)) {
                   tags2 <- addAdductResult2PeakInfo(adduct.result = temp_result[[k]],
                                                     mz.tol = mz_tol,
                                                     rt.tol = rt_tol1,
                                                     weight.mz = 0.8,
                                                     weight.rt = 0.2,
                                                     peak.info = tags2,
                                                     peak.idx = seed_idx[j],
                                                     anno.idx = anno_idx[[j]][k],
                                                     annotation.idx = annotation_idx)

                 }
               }
               rm(list = "adduct_result1")
               gc()
             }

             # dir.create(file.path(path, 'test'), showWarnings = FALSE, recursive = TRUE)
             # save(tags2,
             #      file = file.path(path, 'test', 'temp_tags2_adduct_201110.RData'))

             # dir.create(file.path(path, 'test'), showWarnings = FALSE, recursive = TRUE)
             # save(tags2,
             #      file = file.path(path, 'test', paste0('round', round, '_temp_tags2_adduct_210304.RData')),
             #      version = 2)


             # metabolite annotation -------------------------------------------
             if (met_annotation) {
               cat("\n");cat("Metabolite annotation.\n")
               # browser()

               # 20210317
               tempFun <- function(k){
                 index <- seed_idx[k]
                 annotation_idx <- 1:length(tags2)

                 library(Rcpp, quietly = TRUE,
                         logical.return = FALSE, warn.conflicts = FALSE)

                 library(dplyr, quietly = TRUE,
                         logical.return = FALSE, warn.conflicts = FALSE)

                 library(SpectraTools, quietly = TRUE,
                         logical.return = FALSE, warn.conflicts = FALSE)
                 # suppressMessages(data(thermo, package = "CHNOSZ", envir = environment()))
                 temp_tags2 <- tags2[index]
                 # rm(tags2)
                 # gc()
                 result <- lapply(seq_along(temp_tags2), function(i){
                   seed <- temp_tags2[[i]]
                   seed.as.seed.round <-
                     showTags2(list(seed), slot = "as.seed")[[2]][[1]]
                   seed.i <- which(seed.as.seed.round == round)
                   peak_mz_all <- peak_mz[annotation_idx]
                   peak_rt_all <- peak_rt[annotation_idx]
                   peak_ccs_all <- peak_ccs[annotation_idx]
                   # peak_cor_all <- rep(1, length(annotation_idx))
                   rm(list=c("peak_mz", "peak_rt"))
                   gc()
                   if (length(seed.i) == 0) {
                     return(NULL)
                   } else {
                     temp_result <- lapply(seed.i, function(j){
                       query_name <- seed@name
                       query_id <- seed@annotation[[j]]$to
                       query_charge <- seed@annotation[[j]]$charge
                       query_level <- seed@annotation[[j]]$level
                       query_formula <- stringr::str_trim(seed@annotation[[j]]$formula)
                       query_adduct <- seed@annotation[[j]]$adduct
                       query_mz <- seed@mz
                       query_rt <- seed@rt
                       query_ccs <- seed@ccs
                       metabolite_result <- NULL
                       reaction_step <- 1
                       while (is.null(metabolite_result) & reaction_step <= max_step){
                         # cat(reaction_step, ' ')
                         metabolite_result <-
                           annotateMetaboliteMRN(metabolite_name = query_name,
                                                 metabolite_id = query_id,
                                                 formula = query_formula,
                                                 adduct = query_adduct,
                                                 scoring_approach = scoring_approach,
                                                 instrument = instrument,
                                                 mz = query_mz,
                                                 rt = query_rt,
                                                 peak_mz_all = peak_mz_all,
                                                 peak_rt_all = peak_rt_all,
                                                 peak_ccs_all = peak_ccs_all,
                                                 ms2 = ms2,
                                                 mz_tol = mz_tol,
                                                 mz_ppm_thr = mz_ppm_thr,
                                                 rt_tol = rt_tol2,
                                                 ccs_tol = ccs_tol,
                                                 dp_tol = dp_tol,
                                                 matched_frag_tol = matched_frag_tol,
                                                 whether_link_frag = whether_link_frag,
                                                 step = reaction_step,
                                                 mz_tol_ms2 = mz_tol_ms2,
                                                 metabolite = metabolite,
                                                 metabolite_ccs = metabolite_ccs,
                                                 metabolic_network = metabolic_network,
                                                 adduct_table = adduct_table,
                                                 whether_use_predRT = whether_use_predRT)

                         reaction_step <- reaction_step + 1
                       }
                       metabolite_result
                     })
                     temp_result
                   }
                 })
                 result

               }

               library(parallel)
               # cl <- makeCluster(3L)
               cl <- parallel::makeCluster(threads)
               parallel::clusterExport(cl, c('tempFun',
                                             "showTags2",
                                             "annotateMetaboliteMRN",
                                             'runSpecMatch',
                                             "checkElement",
                                             "sumFormula",
                                             "splitFormula",
                                             "pasteElement",
                                             "getNeighbor",
                                             'runSpecMatch',
                                             'convertSpectraData',
                                             '%>%',
                                             'polarity',
                                             'tags2',
                                             'ms2',
                                             'seed_idx',
                                             'max_step',
                                             'peak_mz',
                                             'peak_rt',
                                             'peak_ccs',
                                             'scoring_approach',
                                             'instrument',
                                             'mz_tol',
                                             'mz_ppm_thr',
                                             'rt_tol2',
                                             'dp_tol',
                                             'ccs_tol',
                                             'metabolite',
                                             'metabolite_ccs',
                                             'metabolic_network',
                                             'round',
                                             'adduct_table',
                                             'matched_frag_tol',
                                             'whether_link_frag',
                                             'mz_tol_ms2',
                                             'scoring_approach',
                                             'whether_use_predRT'),
                                       envir = environment())

               system.time(
                 metabolite_result <- parallel::parLapply(cl = cl,
                                                          seq_along(seed_idx),
                                                          function(k) {tempFun(k)}
                 )
               )
               stopCluster(cl)



               # # debug 20220421
               # metabolite_result <- lapply(seq_along(seed_idx), function(k){
               #   # k <- 82
               #   cat(k, ' ')
               #   index <- seed_idx[k]
               #   annotation_idx <- 1:length(tags2)
               # 
               #   library(Rcpp, quietly = TRUE,
               #     logical.return = FALSE, warn.conflicts = FALSE)
               # 
               #   library(dplyr, quietly = TRUE,
               #     logical.return = FALSE, warn.conflicts = FALSE)
               # 
               #   library(SpectraTools, quietly = TRUE,
               #     logical.return = FALSE, warn.conflicts = FALSE)
               #   # suppressMessages(data(thermo, package = "CHNOSZ", envir = environment()))
               #   temp_tags2 <- tags2[index]
               #   # rm(tags2)
               #   # gc()
               #   result <- lapply(seq_along(temp_tags2), function(i){
               #     seed <- temp_tags2[[i]]
               #     seed.as.seed.round <-
               #       showTags2(list(seed), slot = "as.seed")[[2]][[1]]
               #     seed.i <- which(seed.as.seed.round == round)
               #     peak_mz_all <- peak_mz[annotation_idx]
               #     peak_rt_all <- peak_rt[annotation_idx]
               #     peak_ccs_all <- peak_ccs[annotation_idx]
               #     # peak_cor_all <- rep(1, length(annotation_idx))
               #     rm(list=c("peak_mz", "peak_rt"))
               #     gc()
               #     if (length(seed.i) == 0) {
               #       return(NULL)
               #     } else {
               #       temp_result <- lapply(seed.i, function(j){
               #         query_name <- seed@name
               #         query_id <- seed@annotation[[j]]$to
               #         query_charge <- seed@annotation[[j]]$charge
               #         query_level <- seed@annotation[[j]]$level
               #         query_formula <- stringr::str_trim(seed@annotation[[j]]$formula)
               #         query_adduct <- seed@annotation[[j]]$adduct
               #         query_mz <- seed@mz
               #         query_rt <- seed@rt
               #         query_ccs <- seed@ccs
               #         metabolite_result <- NULL
               #         reaction_step <- 1
               #         while (is.null(metabolite_result) & reaction_step <= max_step){
               #           # cat(reaction_step, ' ')
               #           metabolite_result <-
               #             annotateMetaboliteMRN(metabolite_name = query_name,
               #               metabolite_id = query_id,
               #               formula = query_formula,
               #               adduct = query_adduct,
               #               scoring_approach = scoring_approach,
               #               instrument = instrument,
               #               mz = query_mz,
               #               rt = query_rt,
               #               peak_mz_all = peak_mz_all,
               #               peak_rt_all = peak_rt_all,
               #               peak_ccs_all = peak_ccs_all,
               #               ms2 = ms2,
               #               mz_tol = mz_tol,
               #               mz_ppm_thr = mz_ppm_thr,
               #               rt_tol = rt_tol2,
               #               ccs_tol = ccs_tol,
               #               dp_tol = dp_tol,
               #               matched_frag_tol = matched_frag_tol,
               #               whether_link_frag = whether_link_frag,
               #               step = reaction_step,
               #               mz_tol_ms2 = mz_tol_ms2,
               #               metabolite = metabolite,
               #               metabolite_ccs = metabolite_ccs,
               #               metabolic_network = metabolic_network,
               #               adduct_table = adduct_table,
               #               whether_use_predRT = whether_use_predRT)
               # 
               #           reaction_step <- reaction_step + 1
               #         }
               #         metabolite_result
               #       })
               #       temp_result
               #     }
               #   })
               #   result
               # 
               # })

               metabolite_result1 <- lapply(metabolite_result, function(x){x[[1]]})

               # metabolite_result1 <- metabolite_result1[order(unlist(ichunks))]

               # metabolite_result1 <- metabolite_result
               # add metabolite_result to tags2
               for (j in 1:length(seed_idx)){
                 # cat(j, ' ')
                 if (instrument %in% c("AgilentDTIMMS", "BrukerTIMS", "WatersTWIMMS")) {
                   weight.mz <- 0.1
                   weight.rt <- weight.ccs <- 0.2
                   weight.dp <- 0.5
                 } else {
                   weight.mz <- weight.rt <- 0.25
                   weight.dp <- 0.5
                   weight.ccs <- 0
                 }

                 temp_result <- metabolite_result1[[j]]
                 for (k in 1:length(temp_result)){
                   tags2 <- addMetaboliteResult2PeakInfo(metabolite.result = temp_result[[k]],
                                                         mz.tol = mz_tol,
                                                         rt.tol = rt_tol2,
                                                         dp.tol = dp_tol,
                                                         ccs.tol = ccs_tol,
                                                         matched_frag_tol = matched_frag_tol,
                                                         whether_link_frag = whether_link_frag,
                                                         weight.mz = weight.mz,
                                                         weight.rt = weight.rt,
                                                         weight.ccs = weight.ccs,
                                                         weight.dp = weight.dp,
                                                         peak.info = tags2,
                                                         peak.idx = seed_idx[j],
                                                         anno.idx = anno_idx[[j]][k],
                                                         metabolite = metabolite,
                                                         annotation.idx = annotation_idx)
                 }
               }
               rm(list = "metabolite_result1", "metabolite_result")
               gc()
             }

             # dir.create(file.path(path, 'test'), showWarnings = FALSE, recursive = TRUE)
             # save(tags2,
             #      file = file.path(path, 'test', paste0('round', round, '_temp_tags2_met_210304.RData')),
             #      version = 2)

             # export all results after recursive annotaion --------------------
             new_tags <- tags2
             rm(list = "tags2")
             gc()
             # new_tags <- new_tags

             return(new_tags)
           }
)





################################################################################
#       annotateIsotopeMRN -----------------------------------------------------------

#' @title annotateIsotopeMRN
#' @description Find the isotopes of a known metabolite according to rt, mz and intensity. Old function in MetDNA1.
#' @author Zhiwei Zhou, Xiaotao Shen
#' @param id The KEGG ID of metabolite.
#' @param formula The formula of metabolite.
#' @param adduct The adduct of metabolite.
#' @param charge The charge of metabolite.
#' @param mz The mz of metabolite.
#' @param rt The rt of metabolite.
#' @param int The mean/median intensity of metabolite accross the samples.
#' @param peak.mz The mz of peaks.
#' @param peak.rt The rt of peaks.
#' @param peak.int The mean/median intensity of peaks accross the samples.
#' @param cor The correlations between metabolite and peaks.
#' @param rt.tol The rt tolerance (second).
#' @param mz.tol The mz tolerance.
#' @param cor.tol The cor tolerance.
#' @param int.tol The intensity ratio tolerance.
#' @param max.isotope The max number of isotopes.
#' @return  Isotope annotation information of peaks.

##test data
# setwd('/home/jasper/work/MetDNA/isotopeAnnotation')
# id = NA
# formula = "C21H27N7O14P2"
# adduct = "M+H"
# charge = 1
# mz = 664.1139
# rt = 812.143
# int = 205109.4
# data <- readr::read_csv("annotation.result.csv")
# data <- as.data.frame(data)
# sample <- data[,grep("Sample", colnames(data))]
# tags <- data[,-grep("Sample", colnames(data))]
# peak.mz <- tags$mz
# peak.rt <- tags$rt
# peak.int <- apply(sample, 1, median)
# cor <- rep(1, nrow(tags))
# rt.tol = 5
# or.tol = 0
# int.tol = 500
# max.isotope = 4
#
# system.time(a <- isotopeAnnotation(peak.mz = peak.mz, peak.rt = peak.rt, peak.int = peak.int, cor = cor))
# system.time(b <- isotopeAnnotation1(peak.mz = peak.mz, peak.rt = peak.rt, peak.int = peak.int, cor = cor))


#
# id = 'C01770'
# formula = "C4H6O2"
# adduct = "M+H"
# charge = 1
# mz = 87.04448
# rt = 685.0255
# int = 327779.1
# peak_mz_all = peak_mz.temp
# peak_rt_all = peak_rt.temp
# peak_int_all = peak_int.temp
# peak_cor_all = cor.temp
# rt_tol = 3
# mz_tol = 25
# cor_tol = 0
# int_tol = 500
# max_isotope = 4


# annotateIsotopeMRN(id = 'C01770',
#                    formula = "C4H6O2",
#                    adduct = 'M+H',
#                    charge = 1,
#                    mz = 87.0448,
#                    rt = 685.0255,
#                    int = 327779.1,
#                    peak_mz_all = peak_mz.temp,
#                    peak_rt_all = peak_rt.temp,
#                    peak_int_all = peak_int.temp,
#                    peak_cor_all = cor.temp,
#                    rt_tol = 3,
#                    mz_tol = 25,
#                    cor_tol = 0,
#                    int_tol = 500,
#                    max_isotope = 4)

setGeneric(name = "annotateIsotopeMRN",
           def = function(# metabolite information
             id = NA,
             formula = "C21H27N7O14P2",
             adduct = "M+H",
             charge = 1,
             mz = 664.1139,
             rt = 812.143,
             int = 205109.4,
             ## peak information
             peak_mz_all,
             peak_rt_all,
             peak_int_all,
             peak_cor_all,
             ## other parameters
             rt_tol = 5,
             mz_tol = 25,
             mz_ppm_thr = 400,
             cor_tol = 0,
             int_tol = 500,
             max_isotope = 4){
             # browser()


             formula1 <- sumFormula(formula = formula, adduct = adduct)

             # should be fix latter ??!!!!!
             if(is.na(formula1)) formula1 <- formula

             molecule <- Rdisop::getMolecule(formula = formula1,
                                             z = charge,
                                             maxisotopes = max_isotope + 1)
             isotopes <- t(Rdisop::getIsotope(molecule = molecule))
             rownames(isotopes) <-
               c("[M]",paste("[M","+",c(1:(nrow(isotopes)-1)),"]", sep = ""))
             isotopes <- data.frame(isotopes, rownames(isotopes),
                                    stringsAsFactors = FALSE)
             colnames(isotopes) <- c("mz", "intensity", "isotope")
             accurate_mz <- mz

             # convert to rela intensity (normalized to the peak)
             peak_int_all <- peak_int_all/int

             # if int is 0
             peak_int_all[is.na(peak_int_all)] <- 1
             peak_int_all[is.nan(peak_int_all)] <- 1
             peak_int_all[is.infinite(peak_int_all)] <- 100000

             # rt and correlation filtering
             rt_error <- abs(rt - peak_rt_all)
             index1 <- which(rt_error <= rt_tol & peak_cor_all >= cor_tol)

             if(length(index1) == 0) return(NULL)
             #all the peaks are filtered using rt and names as 1
             peak_mz_all1 <- peak_mz_all[index1]
             peak_rt_all1 <- peak_rt_all[index1]
             peak_int_all1 <- peak_int_all[index1]
             peak_cor_all1 <- peak_cor_all[index1]

             #iso_idx is the index of peak for isotopes
             iso_idx <- NULL
             mz_error <- NULL
             int_error <- NULL
             iso <- NULL
             correlation <- NULL

             for(i in 2:nrow(isotopes)){
               # cat(i, ' ')
               # browser()

               # mz and intensity infromation
               temp_mz <- as.numeric(isotopes[i,"mz"])
               temp_int <-
                 as.numeric(isotopes[i,"intensity"])/as.numeric(isotopes[1,"intensity"])

               # calculate error
               # peak_mz_all.error <- abs(temp_mz - peak_mz_all1)*10^6/temp_mz
               # peak_mz_all_error <- abs(temp_mz - peak_mz_all1)*10^6/ifelse(temp_mz>=400,temp_mz,400)

               peak_mz_all_error <- abs(temp_mz - peak_mz_all1)*10^6/ifelse(temp_mz>=mz_ppm_thr, temp_mz, mz_ppm_thr)
               peak_int_all_error <-
                 abs(temp_int - peak_int_all1)*100/temp_int

               # has peak matched the mz tolerance
               idx <- which(peak_mz_all_error <= mz_tol)
               # idx=0, no peaks matched
               if(length(idx) == 0) {
                 iso_idx[i] <- NA
                 mz_error[i] <- NA
                 int_error[i] <- NA
                 correlation[i] <- NA
                 iso[i] <- NA
               }

               # one peak matched, it is
               if(length(idx) == 1) {
                 iso_idx[i] <- idx
                 mz_error[i] <- peak_mz_all_error[idx]
                 int_error[i] <- peak_int_all_error[idx]
                 correlation[i] <- peak_cor_all1[idx]
                 iso[i] <- isotopes[i,3]
               }

               # more than one matched, see the intensity ratio error (select the intensity ratio )
               if(length(idx) > 1) {
                 idx <- idx[which.min(peak_int_all_error[idx])]
                 iso_idx[i] <- idx
                 mz_error[i] <- peak_mz_all_error[idx]
                 int_error[i] <- peak_int_all_error[idx]
                 correlation[i] <- peak_cor_all1[idx]
                 iso[i] <- isotopes[i,3]
               }
             }

             # -----------------------------------------------------------------
             if(all(is.na(iso_idx))) {
               return(NULL)
             }else{
               iso_idx <- iso_idx[!is.na(iso_idx)]
               index2 <- index1[iso_idx]
               mz_error <- mz_error[!is.na(mz_error)]
               int_error <- int_error[!is.na(int_error)]
               correlation <- correlation[!is.na(correlation)]
               rt_error <- rt_error[index2]
               iso <- iso[!is.na(iso)]
               iso_info <- data.frame(index2, rt_error, mz_error,
                                      int_error, correlation,
                                      iso, stringsAsFactors = FALSE)
             }
             colnames(iso_info) <- c("peakIndex", "rtError.s",
                                     "mzError.ppm", "IntensityRatioError.%",
                                     "correlation", "isotopes")

             # remove peaks which have large intensity ratio error. For M+1 isotope, the tolerance is set as 30%, and for other isotopes,
             # the tolerance is set as 100%.

             remove_idx1 <-
               which(iso_info$isotopes == "[M+1]" & iso_info$"IntensityRatioError.%" > int_tol)

             remove_idx2 <-
               which(iso_info$isotopes != "[M+1]" & iso_info$"IntensityRatioError.%" > int_tol)

             remove_idx <- c(remove_idx1, remove_idx2)
             if(length(remove_idx) != 0){
               iso_info <- iso_info[-c(remove_idx1, remove_idx2),]
             }
             if(nrow(iso_info) == 0) {return(NULL)}
             iso_info <- data.frame(id, iso_info, stringsAsFactors = FALSE)
             rm(list = c("peak_mz_all", "peak_rt_all", "peak_int_all", "cor"))
             gc()
             iso_info <- iso_info

             return(iso_info)
           })


#       addIsotopeResult2PeakInfo ----------------------------------------------------
#' @title addIsotopeResult2PeakInfo
#' @description Add isotope result into tags2 data.
#' @author Zhiwei Zhou, Xiaotao Shen
#' @param isotopes.result The isotope result.
#' @param mz.tol The mz tolerance.
#' @param rt.tol The RT tolerance.
#' @param int.tol The intensity ratio tolerance.
#' @param weight.mz mz weight.
#' @param weight.rt RT weight.
#' @param weight.int Intesnity ratio weight.
#' @param peak.info The tags2 data.
#' @param peak.idx The index of seeds.
#' @param anno.idx The index of annotation
#' @param annotation.idx The index of peaks which are annotated.
#' @return Tags2 data.
# export

setGeneric(name = "addIsotopeResult2PeakInfo",
           def = function(isotopes.result,
                          mz.tol = 25,
                          rt.tol = 3,
                          int.tol = 500,
                          weight.mz = 0.45,
                          weight.rt = 0.45,
                          weight.int = 0.1,
                          peak.info,
                          peak.idx,
                          anno.idx,
                          annotation.idx = c(1:length(peak.info))){
             #
             if (is.null(isotopes.result)) return(peak.info)
             index <- annotation.idx[isotopes.result[,"peakIndex"]]
             seed <- peak.info[[peak.idx]]
             seed.name <- seed@name
             #score formula
             ##calculate score for peak, using mz error, rt error
             temp.x <- c(0, mz.tol)
             temp.y <- c(1, 0)
             lm.reg.mz <- lm(temp.y~temp.x)

             temp.x <- c(0, rt.tol)
             temp.y <- c(1, 0)
             lm.reg.rt <- lm(temp.y~temp.x)

             temp.x <- c(0, int.tol)
             temp.y <- c(1, 0)
             lm.reg.int <- lm(temp.y~temp.x)


             for(i in 1:length(index)){
               if(length(peak.info[[index[i]]]@annotation) == 0){
                 k = 1
                 peak.info[[index[i]]]@annotation <- list(
                   list(type = NA,
                        From = NA,
                        From.peak = NA,
                        to = NA,
                        step = NA,
                        level = NA,
                        as.seed = FALSE,
                        as.seed.round = NA,
                        isotope = NA,
                        adduct = NA,
                        charge = NA,
                        formula = NA,
                        mz.error = NA,
                        rt.error = NA,
                        ccs.error = NA,
                        int.error = NA,
                        ms2.sim = NA,
                        nfrag = NA,
                        score = NA)
                 )
               }else{
                 k = length(peak.info[[index[i]]]@annotation) + 1
                 peak.info[[index[i]]]@annotation <-
                   c(peak.info[[index[i]]]@annotation,
                     list(peak.info[[index[i]]]@annotation[[1]]))
               }

               peak.info[[index[i]]]@annotation[[k]]$type = "isotopeAnnotation"
               peak.info[[index[i]]]@annotation[[k]]$From =
                 unname(seed@annotation[[anno.idx]]$to)
               peak.info[[index[i]]]@annotation[[k]]$From.peak =
                 unname(seed.name)
               peak.info[[index[i]]]@annotation[[k]]$to =
                 unname(seed@annotation[[anno.idx]]$to)
               peak.info[[index[i]]]@annotation[[k]]$step = NA
               peak.info[[index[i]]]@annotation[[k]]$level =
                 unname(seed@annotation[[anno.idx]]$level + 1)
               peak.info[[index[i]]]@annotation[[k]]$as.seed = FALSE
               peak.info[[index[i]]]@annotation[[k]]$as.seed.round = NA
               peak.info[[index[i]]]@annotation[[k]]$isotope =
                 unname(isotopes.result$isotopes[i])
               peak.info[[index[i]]]@annotation[[k]]$adduct =
                 unname(seed@annotation[[anno.idx]]$adduct)
               peak.info[[index[i]]]@annotation[[k]]$charge =
                 unname(seed@annotation[[anno.idx]]$charge)
               peak.info[[index[i]]]@annotation[[k]]$formula =
                 unname(seed@annotation[[anno.idx]]$formula)
               peak.info[[index[i]]]@annotation[[k]]$mz.error =
                 unname(isotopes.result$mzError.ppm[i])
               peak.info[[index[i]]]@annotation[[k]]$rt.error =
                 unname(isotopes.result$rtError.s[i])
               peak.info[[index[i]]]@annotation[[k]]$ccs.error = NA
               peak.info[[index[i]]]@annotation[[k]]$int.error =
                 unname(isotopes.result$IntensityRatioError..[i])
               peak.info[[index[i]]]@annotation[[k]]$ms2.sim = NA
               peak.info[[index[i]]]@annotation[[k]]$nfrag = NA

               score1 <- coefficients(lm.reg.mz)[1] +
                 coefficients(lm.reg.mz)[2] * isotopes.result$mzError.ppm[i]
               score2 <- coefficients(lm.reg.rt)[1] +
                 coefficients(lm.reg.rt)[2] * isotopes.result$rtError.s[i]
               score3 <- coefficients(lm.reg.int)[1] +
                 coefficients(lm.reg.int)[2] * isotopes.result$IntensityRatioError..[i]
               score <- score1*weight.mz + score2*weight.rt + score3*weight.int

               peak.info[[index[i]]]@annotation[[k]]$score = unname(score)
             }
             peak.info <- peak.info
           })










################################################################################
#       annotateAdductMRN ------------------------------------------------------------
#' @title annotateAdductMRN
#' @description adduct annotation
#' @author Zhiwei Zhou, Xiaotao Shen
#' @param id The kegg ID of metabolite.
#' @param formula The formula of metabolite.
#' @param adduct The adduct of metabolite.
#' @param polarity The mode of data.
#' @param mz The mz of metabolite.
#' @param rt The rt of metabolite.
#' @param adduct.table The adduct table.
#' @param peak.mz The mz of peaks.
#' @param peak.rt The rt of peaks.
#' @param cor The correlations between metabolite and peaks.
#' @param mz.tol The mz tol.
#' @param rt.tol The rt tol (second).
#' @param cor.tol The correlation tol.
#' @return  Isotope annotation information of peaks.
#' @export

#test data
# setwd('/home/jasper/work/MetDNA/adductAnnotation')
# id = NA
# formula = "C6H14N4O2"
# adduct = "M+H"
# polarity = "positive"
# mz = 175.118
# rt = 961.6225
# data <- readr::read_csv("annotation.result.csv")
# data <- as.data.frame(data)
# sample <- data[,grep("Sample", colnames(data))]
# tags <- data[,-grep("Sample", colnames(data))]
# peak.mz <- tags$mz
# peak.rt <- tags$rt
# cor <- rep(1, nrow(tags))
# mz.tol = 25
# rt.tol = 3
# cor.tol = 0
# load("adduct.table.hilic.rda")
# adduct.table <- adduct.table.hilic
#
# system.time(a <- adductAnnotation(adduct.table = adduct.table, peak.mz = peak.mz, peak.rt = peak.rt, cor = cor))
# system.time(b <- adductAnnotation1(adduct.table = adduct.table, peak.mz = peak.mz, peak.rt = peak.rt, cor = cor))


# system.time(lapply(1:100, function(x) adductAnnotation(adduct.table = adduct.table, peak.mz = peak.mz, peak.rt = peak.rt, cor = cor)))
# system.time(lapply(1:100, function(x) adductAnnotation1(adduct.table = adduct.table, peak.mz = peak.mz, peak.rt = peak.rt, cor = cor)))


# id = "C01770"
# formula = "C4H6O2"
# adduct = "M+H"
# polarity = 'positive'
# mz = 87.04448
# rt = 685.0255
# peak_mz_all = peak_mz.temp
# peak_rt_all = peak_rt.temp
# peal_cor_all = cor.temp
# rt_tol = 3
# mz_tol = 25
# cor_tol = 0
# adduct_table <- convertAdductList2AdductTable(polarity = 'positive',
#                                               column = 'hilic',
#                                               adduct_type = c('Common', 'Rare'))


setGeneric(name = "annotateAdductMRN",
           def = function(
             # metabolite information
             id = NA,
             formula = "C6H14N4O2",
             adduct = "M+H",
             polarity = c("positive", "negative"),
             mz = 175.118,
             mz_ppm_thr = 400,
             rt = 961.6225,
             adduct_table,
             ## peak information
             peak_mz_all,
             peak_rt_all,
             # peak_cor_all,
             ## other parameters
             mz_tol = 25,
             rt_tol = 3,
             cor_tol = 0.5){

             polarity = match.arg(polarity)

             # if (polarity == "positive") {
             #   adduct_table <- adduct_table[adduct_table[,"mode"] == "positive",]
             # }else{
             #   adduct_table <- adduct_table[adduct_table[,"mode"] == "negative",]
             # }

             # remove metabolite adduct
             idx <- which(adduct_table[,1] == adduct)
             if (length(idx) == 1) {
               adduct_table <- adduct_table[-idx,]
             }

             if (adduct != "M+") {
               adduct_table <- adduct_table[which(adduct_table[,1] != "M+"),]
             }

             if (adduct != 'M-') {
               adduct_table <- adduct_table[which(adduct_table[,1] != "M-"),]
             }

             adduct_name <- as.character(adduct_table[,1])
             adduct_charge <- as.numeric(adduct_table[,3])

             # remove invalid
             remain_index <- which(sapply(adduct_name, function(x) {
               checkElement(formula, x)}))

             if(length(remain_index) == 0) return(NULL)

             adduct_name <- adduct_name[remain_index]
             adduct_charge <- adduct_charge[remain_index]

             # formula1 is the new adduct formula
             formula1 <- sapply(adduct_name, function(x) {sumFormula(formula,x)})

             # if the metabolite don't have enough element to remove, formula is NA.
             charge <- adduct_charge[!is.na(formula1)]
             formula1 <- formula1[!is.na(formula1)]

             # get the accurate mass for different adduct format
             accurate_mass <- sapply(formula1, function(x) {Rdisop::getMolecule(x)$exactmass})
             accurate_mz <- accurate_mass/charge
             
             # process with electron mass in mz calculation
             if (polarity == 'positive') {
               accurate_mz <- accurate_mz - 0.0005
             } else {
               accurate_mz <- accurate_mz + 0.0005
             }

             # adduct_info is the adduct table of adduct
             adduct_info <- data.frame(accurate_mz, rt, charge,
                                       names(accurate_mz),
                                       formula, formula1,
                                       stringsAsFactors = FALSE)
             colnames(adduct_info) <-
               c("mz", "rt","charge", "adduct" , "Formula", "Adduct.Formula")

             # rt mz, and cor filtering
             rt_error <- lapply(adduct_info$rt, function(x){
               abs(x - peak_rt_all)
             })

             mz_error <- lapply(adduct_info$mz, function(x){
               # abs(x - peak_mz_all)*10^6/x
               # abs(x - peak_mz_all)*10^6/ifelse(x>=400,x,400)
               abs(x - peak_mz_all)*10^6/ifelse(x>=mz_ppm_thr,x,mz_ppm_thr)
             })

             index1 <- mapply(function(x, y){
               temp.idx <- which(x <= rt_tol & y <= mz_tol)
             },
             x = rt_error,
             y = mz_error)


             rm(list = c("peak_mz_all", "peak_rt_all", "peak_cor_all"))
             gc()
             if(all(unlist(lapply(index1, length)) == 0)) return(NULL)

             add_info <- mapply(function(temp_idx, temp_mz_error, temp_rt_error){
               if(length(temp_idx) == 0) return(list(rep(NA, 4)))
               if(length(temp_idx) == 1) {
                 list(c(temp_idx, temp_rt_error[temp_idx], temp_mz_error[temp_idx],
                        1))
               }else{
                 temp_idx <- temp_idx[which.min(temp_rt_error[temp_idx])][1]
                 list(c(temp_idx, temp_rt_error[temp_idx], temp_mz_error[temp_idx],
                        1))
               }

             },
             temp_idx = index1,
             temp_mz_error = mz_error,
             temp_rt_error = rt_error)

             rm(list = c("mz_error", "rt_error"))
             gc()

             add_info <- do.call(rbind, add_info)
             add_info <- data.frame(add_info, adduct_info$adduct,
                                    stringsAsFactors = FALSE)

             colnames(add_info) <-
               c("peakIndex", "rtError.s", "mzError.ppm", "correlation", "adducts")

             add_info <- add_info[which(apply(add_info, 1,
                                              function(x) all(!is.na(x)))),,drop = FALSE]
             add_info <- data.frame(id, add_info, stringsAsFactors = FALSE)

             add_info <- add_info

           })



#       addAdductResult2PeakInfo --------------------------------------------------------------

# title addAdductResult2PeakInfo
# description Add adduct result into tags2 data.
# author Xiaotao Shen
# \email{shenxt@@sioc.ac.cn}
# param adduct.result The adduct result.
# param mz.tol The mz tolerance.
# param rt.tol The RT tolerance.
# param weight.mz mz weight.
# param weight.rt RT weight.
# param peak.info The tags2 data.
# param peak.idx The index of seeds.
# param anno.idx The index of annotation
# param annotation.idx The index of peaks which are annotated.
# return Tags2 data.
# export
setGeneric(name = "addAdductResult2PeakInfo",
           def = function(adduct.result,
                          mz.tol = 25,
                          rt.tol = 3,
                          weight.mz = 0.8,
                          weight.rt = 0.2,
                          peak.info,
                          peak.idx,
                          anno.idx,
                          annotation.idx = c(1:length(peak.info))){
             if(is.null(adduct.result)) return(peak.info)
             index <- annotation.idx[adduct.result[,"peakIndex"]]
             seed <- peak.info[[peak.idx]]

             #------------------------------------------------------------------
             ##seed information
             seed.name <- seed@name
             seed.from.name <- seed@annotation[[anno.idx]]$From.peak
             seed.from <- seed@annotation[[anno.idx]]$From

             ##find the peaks which annotate the seed. For example, if seed
             #annotatte
             #if seed.from.name or seed.form is NA, this means this seed is the
             #raw seed for annotation
             if(!is.na(seed.from.name) | !is.na(seed.from)){
               #index1 is the index peak which have annotation
               index1 <-
                 index[which(showTags2(peak.info[index],
                                       slot = "annotation.len") > 0)]
               if(length(index1) == 0){
                 index <- index
               }else{
                 #remove.i is the index of annotated peaks which should not be
                 #annotated by this seed
                 remove.i <- NULL
                 for(i in 1:length(index1)){
                   name <- peak.info[[index1[i]]]@name
                   # id <- peak.info[[index1[i]]]@annotation[[1]]$to
                   id <-
                     unlist(lapply(peak.info[[210]]@annotation, function(x) x$to))
                   as.seed <-
                     unlist(lapply(peak.info[[210]]@annotation, function(x) x$as.seed))
                   id <- id[as.seed]
                   if(length(id) == 0) id <- "no"
                   if(seed.from.name == name & any(seed.from == id)){
                     remove.i[i] <- index1[i]
                   }
                 }
                 remove.i <- remove.i[!is.na(remove.i)]
                 index <- setdiff(index, remove.i)
               }
             }
             if(length(index) == 0) return(peak.info)
             #------------------------------------------------------------------

             #score formula
             ##calculate score for peak, using mz error, rt error
             ##and int error
             temp.x <- c(0, mz.tol)
             temp.y <- c(1, 0)
             lm.reg.mz <- lm(temp.y~temp.x)

             temp.x <- c(0, rt.tol)
             temp.y <- c(1, 0)
             lm.reg.rt <- lm(temp.y~temp.x)


             for(i in 1:length(index)){
               if(length(peak.info[[index[i]]]@annotation) == 0){
                 k = 1
                 peak.info[[index[i]]]@annotation <- list(
                   list(type = NA,
                        From = NA,
                        From.peak = NA,
                        to = NA,
                        step = NA,
                        level = NA,
                        as.seed = FALSE,
                        as.seed.round = NA,
                        isotope = NA,
                        adduct = NA,
                        charge = NA,
                        formula = NA,
                        mz.error = NA,
                        rt.error = NA,
                        ccs.error = NA,
                        int.error = NA,
                        ms2.sim = NA,
                        nfrag = NA,
                        score = NA)
                 )
               }else{
                 k = length(peak.info[[index[i]]]@annotation) + 1
                 peak.info[[index[i]]]@annotation <-
                   c(peak.info[[index[i]]]@annotation,
                     list(peak.info[[index[i]]]@annotation[[1]]))
               }
               peak.info[[index[i]]]@annotation[[k]]$type = "adductAnnotation"
               peak.info[[index[i]]]@annotation[[k]]$From =
                 unname(seed@annotation[[anno.idx]]$to)
               peak.info[[index[i]]]@annotation[[k]]$From.peak =
                 unname(seed.name)
               peak.info[[index[i]]]@annotation[[k]]$to =
                 unname(seed@annotation[[anno.idx]]$to)
               peak.info[[index[i]]]@annotation[[k]]$step = NA
               peak.info[[index[i]]]@annotation[[k]]$level =
                 unname(seed@annotation[[anno.idx]]$level + 1)
               peak.info[[index[i]]]@annotation[[k]]$as.seed = FALSE
               peak.info[[index[i]]]@annotation[[k]]$as.seed.round = NA
               peak.info[[index[i]]]@annotation[[k]]$isotope =
                 '[M]'
               peak.info[[index[i]]]@annotation[[k]]$adduct =
                 unname(adduct.result$adducts[i])
               peak.info[[index[i]]]@annotation[[k]]$charge =
                 unname(seed@annotation[[anno.idx]]$charge)
               peak.info[[index[i]]]@annotation[[k]]$formula =
                 unname(seed@annotation[[anno.idx]]$formula)
               peak.info[[index[i]]]@annotation[[k]]$mz.error =
                 unname(adduct.result$mzError.ppm[i])
               peak.info[[index[i]]]@annotation[[k]]$rt.error =
                 unname(adduct.result$rtError.s[i])
               peak.info[[index[i]]]@annotation[[k]]$ccs.error = NA
               peak.info[[index[i]]]@annotation[[k]]$int.error = NA
               peak.info[[index[i]]]@annotation[[k]]$ms2.sim = NA
               peak.info[[index[i]]]@annotation[[k]]$nfrag = NA

               score1 <- coefficients(lm.reg.mz)[1] +
                 coefficients(lm.reg.mz)[2] * adduct.result$mzError.ppm[i]
               score2 <- coefficients(lm.reg.rt)[1] +
                 coefficients(lm.reg.rt)[2] * adduct.result$rtError.s[i]

               score <- score1*weight.mz + score2*weight.rt

               peak.info[[index[i]]]@annotation[[k]]$score = unname(score)
             }
             peak.info <- peak.info
           })





#       convertAdductList2AdductTable ----------------------------------------------------------
#' @title convertAdductList2AdductTable
#' @description convert adduct_list to adduct_table (MetDNA format)
#' @author Zhiwei Zhou
#' @param adduct_list Default: NULL
#' @param polarity 'positive', 'negative'
#' @param column 'hilic', 'rp'

setGeneric(name = 'convertAdductList2AdductTable',
           def = function(
             adduct_list = NULL,
             polarity = c('positive', 'negative'),
             column = c('hilic', 'rp')
           ){

             switch(polarity,
                    'positive' = {
                      data('lib_adduct_nl', envir = environment())
                      if (length(adduct_list) > 0) {
                        adduct_list <- lib_adduct_nl$positive %>%
                          dplyr::filter(adduct %in% adduct_list) %>%
                          dplyr::select(adduct, delta_mz)
                      } else {
                        adduct_list <- lib_adduct_nl$positive %>%
                          dplyr::filter(annotation == 'Yes') %>%
                          dplyr::select(adduct, delta_mz)
                      }

                      adduct_list <- adduct_list %>%
                        dplyr::mutate(name = gsub(adduct, pattern = '\\[|(\\]\\+|\\]\\-)|(\\]2\\-)', replacement = '')) %>%
                        dplyr::mutate(name = dplyr::case_when(
                          name == 'M' ~ 'M+',
                          name != 'M' ~ name
                        )) %>%
                        dplyr::mutate(nmol = 1, charge = 1) %>%
                        dplyr::mutate(nmol = dplyr::case_when(
                          stringr::str_detect(name, pattern = '2M')~2,
                          stringr::str_detect(name, pattern = '3M')~3,
                          !stringr::str_detect(name, pattern = '2M|3M')~1)
                        ) %>%
                        dplyr::select(name, nmol, charge, delta_mz) %>%
                        dplyr::rename(massdiff = delta_mz) %>%
                        dplyr::mutate(charge = dplyr::case_when(name == 'M-2H' ~ 2,
                                                                name != 'M-2H' ~ 1)) %>%
                        as.data.frame()


                    },
                    'negative' = {
                      data('lib_adduct_nl', envir = environment())
                      if (length(adduct_list) > 0) {
                        adduct_list <- lib_adduct_nl$negative %>%
                          dplyr::filter(adduct %in% adduct_list) %>%
                          dplyr::select(adduct, delta_mz)
                      } else {
                        adduct_list <- lib_adduct_nl$negative %>%
                          dplyr::filter(annotation == 'Yes') %>%
                          dplyr::select(adduct, delta_mz)
                      }

                      adduct_list <- adduct_list %>%
                        dplyr::mutate(name = gsub(adduct, pattern = '\\[|(\\]\\+|\\]\\-)|(\\]2\\-)', replacement = '')) %>%
                        dplyr::mutate(name = dplyr::case_when(
                          name == 'M' ~ 'M-',
                          name != 'M' ~ name
                        )) %>%
                        dplyr::mutate(nmol = 1, charge = 1) %>%
                        dplyr::mutate(nmol = dplyr::case_when(
                          stringr::str_detect(name, pattern = '2M')~2,
                          stringr::str_detect(name, pattern = '3M')~3,
                          !stringr::str_detect(name, pattern = '2M|3M')~1)
                        ) %>%
                        dplyr::select(name, nmol, charge, delta_mz) %>%
                        dplyr::rename(massdiff = delta_mz) %>%
                        dplyr::mutate(charge = dplyr::case_when(name == 'M-2H' ~ 2,
                                                                name != 'M-2H' ~ 1)) %>%
                        as.data.frame()
                    })


             return(adduct_list)

           })


################################################################################
#       annotateMetaboliteMRN --------------------------------------------------------
#' @title annotateMetaboliteMRN
#' @description Annotate peak table from one metabolite.
#' @author Zhiwei Zhou
#' \email{zhouzw@@sioc.ac.cn}
#' @param metabolite.name The metabolite name.
#' @param metabolite.id The metabolite ID.
#' @param formula The formula of metabolite.
#' @param adduct The adduct of metabolite.
#' @param polarity The polarity.
#' @param mz The mz of metabolite.
#' @param rt The RT of metabolite.
#' @param peak.mz The mz of all peaks.
#' @param peak.rt The RT of all peaks.
#' @param cor The correlation of metabolite between all peaks.
#' @param ms2 The ms2 data of peak table.
#' @param mz.tol The mz tolerance.
#' @param rt.tol The RT tolerance for metabolite annotation. (\%)
#' @param cor.tol The cor tolerance.
#' @param dp.tol The tolerance of dot product.
#' @param step The reaction step.
#' @param metabolite The kegg compound database.
#' @param metabolic.network kegg.rpair2
#' @param adduct.table Adduct table.
#' @return  The metabolite annotation information.
#' @export

# 20210316
# load('/home/zhouzw/Data_processing/20210315_tims_data_processing_metdna/test/210316_before_metAnnotation.RData')
#
# metabolite_name = "M136T80C125" # seed name
# metabolite_id = "C03758" # seed id
# formula = "C8H11NO2" # not use
# adduct = "M-H2O+H"
# polarity = 'positive'
# scoring_approach = 'dp'
# mz = 136.0726
# rt = 80
# # ## peak information
# load('/home/zhouzw/Data_processing/20210315_tims_data_processing_metdna/test/210316_peak_mz_all.RData')
# load('/home/zhouzw/Data_processing/20210315_tims_data_processing_metdna/test/210316_peak_rt_all.RData')
# load('/home/zhouzw/Data_processing/20210315_tims_data_processing_metdna/test/210316_peak_ccs_all.RData')
# load('/home/zhouzw/Data_processing/20210315_tims_data_processing_metdna/test/210316_peak_cor_all.RData')
# load('/home/zhouzw/Data_processing/20210315_tims_data_processing_metdna/01_result_initial_seed_annotation/00_intermediate_data/ms2')
#
# # peak_mz_all
# # peak_rt_all
# # peak_cor_all # not use
# # ms2
# ## other parameters
# mz_tol = 25
# #r tol is relative(%)
# rt_tol = 30
# ccs_tol = 4
# cor_tol = 0 # not use
# dp_tol = 0.5
# matched_frag_tol = 1
# step = 1
# metabolite
# metabolic_network
# adduct_table
#
# load('/home/zhouzw/Data_processing/20210315_tims_data_processing_metdna/test/210316_peak_mz_all.RData')
# load('/home/zhouzw/Data_processing/20210315_tims_data_processing_metdna/test/210316_peak_rt_all.RData')
# load('/home/zhouzw/Data_processing/20210315_tims_data_processing_metdna/test/210316_peak_ccs_all.RData')
# load('/home/zhouzw/Data_processing/20210315_tims_data_processing_metdna/test/210316_peak_cor_all.RData')

# metabolite_ccs <- lib_ccs$emrnlib_pred %>%
#   dplyr::filter(status == 'Valid') %>%
#   dplyr::select(name, adduct, pred_ccs) %>%
#   dplyr::filter(name %in% metabolite$id) %>%
#   dplyr::mutate(adduct = dplyr::case_when(adduct == '[M+H]+'~'M+H',
#                                           adduct == '[M+Na]+'~'M+Na',
#                                           adduct == '[M-H2O+H]+'~'M-H2O+H',
#                                           adduct == '[M+NH4]+'~'M+NH4',
#                                           adduct == '[M-H]-'~'M-H',
#                                           adduct == '[M+Na-2H]-'~'M+Na-2H',
#                                           adduct == '[M+HCOO]-'~'M+HCOO')) %>%
#   dplyr::mutate(temp_id = paste(name, adduct, sep = '_'))
#
# annotateMetaboliteMRN(
#   metabolite_name = "M136T80C125", # seed name
#   metabolite_id = "C03758", # seed id
#   formula = "C8H11NO2", # not use
#   adduct = "M-H2O+H",
#   scoring_approach = 'dp',
#   # instrument = 'IMMS',
#   instrument = "SciexTripleTOF",
#   mz = 136.0726,
#   rt = 80,
#   peak_mz_all = peak_mz_all,
#   peak_rt_all = peak_rt_all,
#   peak_ccs_all = peak_ccs_all,
#   ms2 = ms2,
#   mz_tol = 25,
#   rt_tol = 30,
#   ccs_tol = 4,
#   dp_tol = 0.5,
#   matched_frag_tol = 1,
#   step = 1,
#   metabolite = metabolite,
#   metabolite_ccs = metabolite_ccs,
#   metabolic_network = metabolic_network,
#   adduct_table = adduct_table
# )

# annotateMetaboliteMRN(
#   metabolite_name = "M329T79C179", # seed name
#   metabolite_id = "C03758", # seed id
#   formula = "C8H11NO2", # not use
#   adduct = "2M+Na",
#   scoring_approach = 'dp',
#   # instrument = 'IMMS',
#   instrument = "IMMS",
#   mz = 329.1472,
#   rt = 79,
#   peak_mz_all = peak_mz_all,
#   peak_rt_all = peak_rt_all,
#   peak_ccs_all = peak_ccs_all,
#   ms2 = ms2,
#   mz_tol = 25,
#   rt_tol = 30,
#   ccs_tol = 4,
#   dp_tol = 0.5,
#   matched_frag_tol = 1,
#   step = 1,
#   metabolite = metabolite,
#   metabolite_ccs = metabolite_ccs,
#   metabolic_network = metabolic_network,
#   adduct_table = adduct_table
# )

setGeneric(name = "annotateMetaboliteMRN",
           def = function(metabolite_name = "M175T962", # seed name
                          metabolite_id = "C00062", # seed id
                          formula = "C6H14N4O2", # not use
                          adduct = "M+H",
                          scoring_approach = c('dp', 'bonanza', 'hybrid', 'gnps'),
                          instrument = c("SciexTripleTOF", "AgilentQTOF", "BrukerQTOF", "ThermoOrbitrap",
                                         'ThermoExploris', "AgilentDTIMMS", "BrukerTIMS", 'WatersQTOF','WatersTWIMMS'),
                          mz = 175.118,
                          rt = 961.6225,
                          ## peak information
                          peak_mz_all,
                          peak_rt_all,
                          peak_ccs_all,
                          ms2,
                          ## other parameters
                          mz_tol = 25,
                          mz_ppm_thr = 400,
                          rt_tol = 30, # %
                          ccs_tol = 4, # %
                          # cor_tol = 0, # not use
                          dp_tol = 0.5,
                          matched_frag_tol = 1,
                          whether_link_frag = FALSE,
                          step = 1,
                          mz_tol_ms2 = 25,
                          whether_use_predRT = TRUE,
                          metabolite,
                          metabolite_ccs,
                          metabolic_network,
                          adduct_table,
                          ...){

             # browser()
             # only the peaks with MS2 are remainded
             peak_name <- names(peak_rt_all)
             raw_peak_name <- peak_name
             ms2_name <- unname(unlist(lapply(ms2, function(x) x[[1]][1,1])))
             temp_index <- match(ms2_name, peak_name)
             if(length(temp_index) == 0) return(NULL)

             peak_mz_all <- peak_mz_all[temp_index]
             peak_rt_all <- peak_rt_all[temp_index]
             peak_ccs_all <- peak_ccs_all[temp_index]
             peak_name <- peak_name[temp_index]
             # peak_cor_all <- peak_cor_all[temp_index]


             # if the annotation is a adduct annotaion without MS2, return NA
             if (!(metabolite_name %in% ms2_name)) {return(NULL)}

             # retrieve neighbor metabolite entries
             met_result <- getNeighbor(metabolite.id = metabolite_id,
                                       step = step,
                                       graph = metabolic_network)

             if (is.null(met_result)) return(NULL)
             if (nrow(met_result) > 100){met_result <- met_result[1:100,]}

             # remove the metabolite which has no standard formula
             remove_idx <- match(met_result$Node.ID, metabolite$id) %>%
               metabolite$monoisotopic_mass[.] %>%
               is.na() %>%
               which()
             if(length(remove_idx) != 0) met_result <- met_result[-remove_idx,,drop = FALSE]
             if(nrow(met_result) == 0) return(NULL)

             # add adduct mz, RT, and CCS information
             temp_idx <- match(met_result$Node.ID, metabolite$id)
             met_result <- met_result %>%
               dplyr::mutate(monoisotopic_mass = metabolite$monoisotopic_mass[temp_idx],
                             formula = metabolite$formula[temp_idx],
                             rt = metabolite$RT[temp_idx])

             # for IM-MS data, only keep available adduct in AllCCS
             if (instrument %in% c("AgilentDTIMMS", "BrukerTIMS", 'WatersTWIMMS')) {
               temp_adduct_table <- adduct_table %>%
                 dplyr::filter(name %in% c('M+H', 'M+Na', 'M+NH4', 'M-H2O+H', 'M-H', 'M+Na-2H'))
             } else {
               temp_adduct_table <- adduct_table
             }

             met_result_mz <- lapply(seq_along(met_result$Node.ID), function(i){
               temp_mz <- met_result$monoisotopic_mass[i]
               result <- sapply(seq_along(temp_adduct_table$name), function(j){
                 calculateMz(exact_mass = temp_mz,
                             adduct = temp_adduct_table$name[j],
                             delta_mz = temp_adduct_table$massdiff[j],
                             nmol = temp_adduct_table$nmol[j],
                             ncharge = temp_adduct_table$charge[j])
               })
               result
             })

             met_result_mz <- met_result_mz %>% do.call(rbind, .) %>% tibble::as_tibble()
             colnames(met_result_mz) <- temp_adduct_table$name
             met_result_meta <- met_result %>%
               dplyr::bind_cols(met_result_mz) %>%
               tidyr::pivot_longer(cols = temp_adduct_table$name,
                                   names_to = 'adduct',
                                   values_to = 'mz') %>%
               dplyr::select(Metabolite.ID, Node.ID, Step, formula, monoisotopic_mass, adduct, mz, rt)

             # add ccs information for met_result_meta
             if (instrument %in% c("AgilentDTIMMS", "BrukerTIMS", 'WatersTWIMMS')) {
               temp_idx <- match(paste(met_result_meta$Node.ID, met_result_meta$adduct, sep = '_'),
                                 paste(metabolite_ccs$temp_id))
               met_result_meta <- met_result_meta %>% dplyr::mutate(ccs = metabolite_ccs$pred_ccs[temp_idx])
             } else {
               met_result_meta <- met_result_meta %>% dplyr::mutate(ccs = -1)
             }

             # check element for adducts: formula without O, remove adduct
             check_element <- mapply(function(x,y){
               checkElement(formula = x, adduct = y)
             },
             x = met_result_meta$formula,
             y = met_result_meta$adduct)

             met_result_meta <- met_result_meta %>%
               dplyr::mutate(check_element = check_element) %>%
               dplyr::filter(check_element) %>%
               dplyr::select(-check_element)

             if(nrow(met_result_meta) == 0) return(NULL)

             # m/z, RT, CCS match
             neighbor_mz_range <- getMzRange(mz = met_result_meta$mz,
                                             ppm = mz_tol,
                                             mz_ppm_thr = mz_ppm_thr)
             idx_match <- lapply(seq_along(met_result_meta$Node.ID), function(i){
               temp_idx <- which(peak_mz_all >= neighbor_mz_range[i,1] & peak_mz_all <= neighbor_mz_range[i,2])
               temp_idx
             })

             if (whether_use_predRT) {
               neighbor_rt_range <- getRtRange(data = met_result_meta$rt,
                                               rel_rt_dev = rt_tol,
                                               rt_threshold = 0,
                                               is_combined = FALSE)

               rt_idx <- lapply(seq_along(met_result_meta$Node.ID), function(i){
                 temp_idx <- which(peak_rt_all >= neighbor_rt_range[i,1] & peak_rt_all <= neighbor_rt_range[i,2])
                 result <- temp_idx
               })

               idx_match <- mapply(function(x, y){
                 intersect(x,y)
               },
               x=idx_match,
               y=rt_idx)
             }

             if (instrument %in% c("AgilentDTIMMS", "BrukerTIMS", 'WatersTWIMMS')) {
               neighbor_ccs_range <- getCcsRange(data = met_result_meta$ccs, rel_dev = ccs_tol)
               ccs_idx <- lapply(seq(nrow(met_result_meta)), function(i){
                 temp_idx <- which(peak_ccs_all >= neighbor_ccs_range[i,1] & peak_ccs_all <= neighbor_ccs_range[i,2])
                 # include no experimental CCS candidates
                 temp_idx2 <- which(is.na(peak_ccs_all))
                 result <- unique(c(temp_idx, temp_idx2))
               })
               idx_match <- mapply(function(x, y){
                 intersect(x,y)
               },
               x=idx_match,
               y=ccs_idx)
             }

             # keep the matched peaks and neighbor metabolites
             neighbor_result <- lapply(seq_along(idx_match), function(i){
               temp_idx <- idx_match[[i]]
               if (length(temp_idx) == 0) {return(NULL)}

               temp_met_result <- met_result_meta[i,]

               temp_peak_name <- peak_name[temp_idx]
               temp_peak_mz <- peak_mz_all[temp_idx]
               temp_peak_rt <- peak_rt_all[temp_idx]
               temp_peak_ccs <- peak_ccs_all[temp_idx]

               temp_met_result <- tibble::tibble(peak_name = temp_peak_name,
                                                 peak_mz = temp_peak_mz,
                                                 peak_rt = temp_peak_rt,
                                                 peak_ccs = temp_peak_ccs) %>%
                 dplyr::bind_cols(temp_met_result)

               return(temp_met_result)
             })

             # if no neighbor matched, return NULL
             temp <- sapply(neighbor_result, length) %>% sum()
             if (sum(temp > 0) == 0) {return(NULL)}
             neighbor_result <- neighbor_result %>% dplyr::bind_rows()


             # run spectra match for each neighbor
             # retrieve ms2 for seeds and neighbors
             seed_ms2 <- match(metabolite_name, ms2_name) %>% ms2[[.]]
             neighbor_ms2 <- match(neighbor_result$peak_name, ms2_name) %>% ms2[.]

             # convert seed ms2 to SpectraData
             obj_seed_ms2 <- list(info = tibble::tibble(NAME = seed_ms2$info[1],
                                                        # PRECURSORMZ = as.numeric(seed_ms2$info[2])
                                                        PRECURSORMZ = as.numeric(mz)),
                                  spec = seed_ms2$spec)
             obj_seed_ms2 <- convertSpectraData(ms2_data = obj_seed_ms2)

             ms2_result <- lapply(seq_along(neighbor_ms2), function(i){
               temp_ms2 <- neighbor_ms2[[i]]
               obj_neighbor_ms2 <- list(info = tibble::tibble(NAME = temp_ms2$info[1],
                                                              # PRECURSORMZ = as.numeric(temp_ms2$info[2])
                                                              PRECURSORMZ = as.numeric(neighbor_result$peak_mz[i])),
                                        spec = temp_ms2$spec)
               obj_neighbor_ms2 <- convertSpectraData(ms2_data = obj_neighbor_ms2)

               # the lib spec always use spec from the peak with smaller mz
               if (obj_neighbor_ms2@info$mz >= obj_seed_ms2@info$mz) {
                 ms2_result <- try(runSpecMatch(obj_ms2_cpd1 = obj_neighbor_ms2,
                                                obj_ms2_cpd2 = obj_seed_ms2,
                                                mz_tol_ms2 = mz_tol_ms2,
                                                scoring_approach = scoring_approach),
                                   silent = TRUE)
               } else {
                 ms2_result <- try(runSpecMatch(obj_ms2_cpd1 = obj_seed_ms2,
                                                obj_ms2_cpd2 = obj_neighbor_ms2,
                                                mz_tol_ms2 = mz_tol_ms2,
                                                scoring_approach = scoring_approach),
                                   silent = TRUE)
               }

               if (scoring_approach != 'dp') {ms2_result@info <- ms2_result@info %>% dplyr::rename(scoreReverse = score)}

               if (class(ms2_result) == 'try-error') {return(c(0, 0))}
               if (length(ms2_result) == 0) {return(c(0, 0))}

               ms2_result_score <- as.numeric(ms2_result@info$scoreReverse)
               ms2_result_frag <- as.numeric(ms2_result@info$n_frag_total)
               ms2_result <- c(ms2_result_score, ms2_result_frag)
               return(ms2_result)
             })

             ms2_result_score <- sapply(ms2_result, function(x){x[1]})
             ms2_result_n_frag <- sapply(ms2_result, function(x){x[2]})

             # merge spectra match result
             neighbor_result <- neighbor_result %>%
               dplyr::mutate(mz_error = abs(peak_mz - mz)*10^6/ifelse(mz >= mz_ppm_thr, mz, mz_ppm_thr),
                             rt_error = abs(peak_rt - rt)*100/rt,
                             ccs_error = abs(peak_ccs - ccs)*100/ccs,
                             ms2_score = ms2_result_score,
                             n_frag = ms2_result_n_frag)

             if (whether_link_frag) {
               neighbor_result <- neighbor_result %>%
                 dplyr::filter((ms2_score >= dp_tol) |
                               (n_frag >= matched_frag_tol))
             } else {
               neighbor_result <- neighbor_result %>%
                 dplyr::filter(ms2_score >= dp_tol,
                               n_frag >= matched_frag_tol)
             }

             if (nrow(neighbor_result) == 0) {
               return(NULL)
             }

             neighbor_result <- neighbor_result %>%
               # dplyr::rename(seed_id = Metabolite.ID,
               #               neighbor_id = Node.ID,
               #               step = Step)
               dplyr::rename(peakName = peak_name,
                             peakMz = peak_mz,
                             peakRT = peak_rt,
                             peakCcs = peak_ccs,
                             mzError.ppm = mz_error,
                             rtError = rt_error,
                             ccsError = ccs_error,
                             dotProduct = ms2_score,
                             nodeID = Metabolite.ID,
                             peakID = Node.ID,
                             theoreticalMZ = mz,
                             RT = rt,
                             CCS = ccs,
                             step = Step,
                             matchedFragNum = n_frag) %>%
               dplyr::mutate(peakIndex = NA,
                             correlation = 1,
                             charge = 1,
                             attribute = 'Predicted') %>%
               dplyr::select(peakName, peakIndex, peakMz, peakRT, peakCcs, mzError.ppm, rtError, ccsError, correlation,
                             dotProduct, matchedFragNum, step, nodeID, peakID, adduct, charge, theoreticalMZ, RT, CCS, attribute)

             peakIndex <- match(neighbor_result$peakName, raw_peak_name)
             neighbor_result$peakIndex <- peakIndex

             if (!(instrument %in% c("AgilentDTIMMS", "BrukerTIMS", "WatersTWIMMS"))) {
               neighbor_result$ccsError <- -1
             }

             return(neighbor_result)

           })



#       addMetaboliteResult2PeakInfo -------------------------------------------------
#' @title addMetaboliteResult2PeakInfo
#' @description Add annotation information form metAnnotation into peakInfo.
#' @author Zhiwei Zhou, Xiaotao Shen
#' @param metabolite.result metabolite.result from metAnnotation
#' @param mz.tol m/z tolerrance
#' @param rt.tol RT tolerrance
#' @param dp.tol Dot product tolerrance
#' @param weight.mz m/z weight
#' @param weight.rt RT weight
#' @param weight.dp Dot product weight
#' @param peak.info Peak info data
#' @param peak.idx The index of seed peaks.
#' @param anno.idx The index of seed annotations.
#' @param metabolite KEGG compound database
#' @param annotation.idx The index of peaks which are used for annotation
#' @return The tags2 data.
# #' @export

setGeneric(name = "addMetaboliteResult2PeakInfo",
           def = function(metabolite.result,
                          mz.tol = 25,
                          rt.tol = 30,
                          ccs.tol = 4,
                          dp.tol = 0.5,
                          matched_frag_tol = 1,
                          whether_link_frag = FALSE,
                          weight.mz = 0.25,
                          weight.rt = 0.25,
                          weight.ccs = 0,
                          weight.dp = 0.5,
                          peak.info,
                          peak.idx = 333,
                          anno.idx = 1,
                          metabolite = metabolite,
                          annotation.idx = c(1:length(peak.info))){
             # check parameters
             if(!is.list(peak.info)) stop("peak.info provided is invalid.")

             if(is.null(metabolite.result)) return(peak.info)

             # index is the index of peaks which are annotated by seed
             index <- annotation.idx[metabolite.result$peakIndex]
             # seed is the seed which annotate peaks
             seed <- peak.info[[peak.idx]]

             #------------------------------------------------------------------------------
             ##seed information
             seed.name <- seed@name
             seed.from.name <- seed@annotation[[anno.idx]]$From.peak
             seed.from <- seed@annotation[[anno.idx]]$From

             # find the peaks which annotate the seed. For example, if seed annotate
             #if seed.from.name or seed.form is NA, this means this seed is the
             #raw seed for annotation
             if(!is.na(seed.from.name) & !is.na(seed.from)){
               #index1 is the index peak which have annotation
               index1 <-
                 index[which(showTags2(peak.info[index], slot = "annotation.len") > 0)]
               if (length(index1) == 0){
                 index <- index
               } else {
                 #remove.i is the index of annotated peaks which should not be
                 #annotated by this seed
                 remove.i <- NULL
                 for(i in 1:length(index1)){
                   name <- peak.info[[index1[i]]]@name
                   # id <- peak.info[[index1[i]]]@annotation[[1]]$to
                   id <-
                     unlist(lapply(peak.info[[index1[i]]]@annotation, function(x) x$to))
                   as.seed <-
                     unlist(lapply(peak.info[[index1[i]]]@annotation, function(x) x$as.seed))
                   id <- id[as.seed]
                   if(length(id) == 0) id <- "no"
                   if(seed.from.name == name & any(seed.from == id)){
                     remove.i[i] <- index1[i]
                   }
                 }
                 remove.i <- remove.i[!is.na(remove.i)]
                 remove.i <- unique(remove.i)
                 if (length(remove.i) != 0) {
                   metabolite.result <- metabolite.result[-which(remove.i == index),]
                   if (nrow(metabolite.result) != 0) {
                     index <- metabolite.result$peakIndex}else{
                       index <- NULL
                     }
                 }
               }
             }
             if(length(index) == 0) return(peak.info)
             #------------------------------------------------------------------
             #
             #score formula
             ##calculate score for peak, using mz error, rt error
             ##and int error
             temp.x <- c(0, mz.tol)
             temp.y <- c(1, 0)
             lm.reg.mz <- lm(temp.y~temp.x)

             temp.x <- c(0, rt.tol)
             temp.y <- c(1, 0)
             lm.reg.rt <- lm(temp.y~temp.x)

             temp.x <- c(0, ccs.tol)
             temp.y <- c(1, 0)
             lm.reg.ccs <- lm(temp.y~temp.x)

             temp.x <- c(dp.tol, 1)
             temp.y <- c(0, 1)
             lm.reg.dp <- lm(temp.y~temp.x)

             if (whether_link_frag) {
               temp.x <- c(matched_frag_tol, 10)
               temp.y <- c(0.3, 1)
               lm.reg.frag <- lm(temp.y~temp.x)
             }

             for(i in 1:length(index)){
               if(length(peak.info[[index[i]]]@annotation) == 0){
                 k = 1
                 peak.info[[index[i]]]@annotation <- list(
                   list(type = NA,
                        From = NA,
                        From.peak = NA,
                        to = NA,
                        step = NA,
                        level = NA,
                        as.seed = FALSE,
                        as.seed.round = NA,
                        isotope = NA,
                        adduct = NA,
                        charge = NA,
                        formula = NA,
                        mz.error = NA,
                        rt.error = NA,
                        ccs.error = NA,
                        int.error = NA,
                        ms2.sim = NA,
                        nfrag = NA,
                        score = NA))
               }else{
                 k = length(peak.info[[index[i]]]@annotation) + 1
                 peak.info[[index[i]]]@annotation <-
                   c(peak.info[[index[i]]]@annotation,
                     list(peak.info[[index[i]]]@annotation[[1]]))
               }

               ##remove the annotation which is from the peak that it annotate
               ##for example, if peak A is metabolte 1, and the it annotate peak B as
               ##metabolite 2, and then Peak B is used as seed and annotate peak A, so to
               #peak A, all the annotation from Peak B (metabolite 1) should be removed

               peak.info[[index[i]]]@annotation[[k]]$type = "metAnnotation"
               peak.info[[index[i]]]@annotation[[k]]$From =
                 unname(seed@annotation[[anno.idx]]$to)
               peak.info[[index[i]]]@annotation[[k]]$From.peak =
                 unname(seed.name)
               peak.info[[index[i]]]@annotation[[k]]$to =
                 unname(metabolite.result$peakID[i])
               peak.info[[index[i]]]@annotation[[k]]$step =
                 unname(metabolite.result$step[i])
               peak.info[[index[i]]]@annotation[[k]]$level =
                 unname(seed@annotation[[anno.idx]]$level + 1)
               peak.info[[index[i]]]@annotation[[k]]$as.seed = FALSE
               peak.info[[index[i]]]@annotation[[k]]$as.seed.round = NA
               peak.info[[index[i]]]@annotation[[k]]$isotope = '[M]'
               peak.info[[index[i]]]@annotation[[k]]$adduct =
                 unname(metabolite.result$adduct[i])
               peak.info[[index[i]]]@annotation[[k]]$charge =
                 unname(metabolite.result$charge[i])
               peak.info[[index[i]]]@annotation[[k]]$formula =
                 unname(metabolite$formula[match(metabolite.result$peakID[i], metabolite$id)])
               peak.info[[index[i]]]@annotation[[k]]$mz.error =
                 unname(metabolite.result$mzError.ppm[i])
               peak.info[[index[i]]]@annotation[[k]]$rt.error =
                 unname(metabolite.result$rtError[i])
               peak.info[[index[i]]]@annotation[[k]]$ccs.error =
                 unname(metabolite.result$ccsError[i])
               peak.info[[index[i]]]@annotation[[k]]$int.error = NA
               peak.info[[index[i]]]@annotation[[k]]$ms2.sim =
                 unname(metabolite.result$dotProduct[i])
               peak.info[[index[i]]]@annotation[[k]]$nfrag =
                 unname(metabolite.result$matchedFragNum[i])

               score1 <- coefficients(lm.reg.mz)[1] +
                 coefficients(lm.reg.mz)[2] * metabolite.result$mzError.ppm[i]
               score2 <- coefficients(lm.reg.rt)[1] +
                 coefficients(lm.reg.rt)[2] * metabolite.result$rtError[i]
               if (score2 < 0) {score2 <- 0} # when not use predRT, the score2 may be less than 0. This scoe assign 0
               score3 <- coefficients(lm.reg.ccs)[1] +
                 coefficients(lm.reg.ccs)[2] * metabolite.result$ccsError[i]

               if (whether_link_frag) {
                 temp_score1 <- coefficients(lm.reg.dp)[1] +
                   coefficients(lm.reg.dp)[2] * metabolite.result$dotProduct[i]
                 temp_score2 <- coefficients(lm.reg.frag)[1] +
                   coefficients(lm.reg.frag)[2] * metabolite.result$matchedFragNum[i]
                 score4 <- max(c(temp_score1, temp_score2))
               } else {
                 score4 <- coefficients(lm.reg.dp)[1] +
                   coefficients(lm.reg.dp)[2] * metabolite.result$dotProduct[i]
               }

               score <- score1*weight.mz + score2*weight.rt + score3 * weight.ccs + score4 * weight.dp
               peak.info[[index[i]]]@annotation[[k]]$score = unname(score)
             }
             peak.info <- peak.info
           })
#       getNeighbor ------------------------------------------------------------------
#' @title getNeighbor
#' @description Get the neighbors of one metabolite in graph.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param metabolite.id The metabolite ID.
#' @param step The reaction step.
#' @param graph The kegg.rpair2.
#' @return  The neighbor information.
#' @export

setGeneric(name = "getNeighbor",
           def = function(metabolite.id,
                          step = 1,
                          graph){
             # load("kegg.rpair2")
             # if(!igraph::is.igraph(kegg.rpair1)) {
             #   kegg.rpair2 <- igraph::graph_from_graphnel(kegg.rpair1)
             # }
             if(!metabolite.id %in% igraph::V(graph)$name) return(NULL)

             temp.step <- 1
             temp.id <- metabolite.id
             cum.temp.id <- NULL
             while (temp.step <= step) {
               result <-
                 unique(unlist(lapply(temp.id, function(x) names(igraph::neighbors(graph, x)))))
               cum.temp.id <- c(cum.temp.id, temp.id)
               temp.step <- temp.step + 1
               temp.id <- result
             }

             result <- setdiff(result, cum.temp.id)
             result <- sort(result)

             if(length(result) == 0) return(NULL)
             result <-
               data.frame(metabolite.id, result, step,
                          stringsAsFactors = FALSE)
             colnames(result) <- c('Metabolite.ID', 'Node.ID', "Step")
             result <- result
           })

################################################################################
#     trans2Matrix -----------------------------------------------------------

#' @title trans2Matrix
#' @description Transform tags2 data to matrix.
#' @author Zhiwei Zhou, Xiaotao Shen
#' @param tags2 The tags2 data.
#' @param base Sort according to annotation or peak.
#' @return The matrix of tags2 data.
#' @export

setGeneric(name = "trans2Matrix",
           def = function(tags2,
                          base = c("annotation", "peak")) {

             base <- match.arg(base)
             result <- lapply(seq_along(tags2), function(i) {
               # cat(i, ' ')
               x <- tags2[[i]]
               name <- x@name
               mz <- x@mz
               rt <- x@rt
               ccs <- x@ccs
               annotation <- x@annotation
               if(length(annotation) == 0){
                 NULL
               }else{
                 # annotation <- lapply(annotation, function(x) {x <- x[names(x) != "addcut"]})
                 annotation <- do.call(rbind, lapply(annotation, unlist))
                 if(any(colnames(annotation) == "int.ratio.error") & any(colnames(annotation) == "int.error")){
                   annotation[,"int.error"] <- annotation[,"int.ratio.error"]
                   annotation <- annotation[,c(1:19)]
                   colnames(annotation)[16] <- "int.error"
                 }
                 if(any(colnames(annotation) == "int.ratio.error") & all(colnames(annotation) != "int.error")){
                   colnames(annotation)[16] <- "int.error"
                 }

                 data.frame(name, mz, rt, annotation,
                            stringsAsFactors = FALSE)
               }

             })
             result <- do.call(what = "rbind", args = result)
             if(base == "annotation"){
               result <- result[order(result$to),]
             }else{
               result <- result[order(result$name),]
             }
             return(result)
           }
)

################################################################################
#     groupRT ------------------------------------------------------------------

#' @title groupRT
#' @description Group peaks according to RT.
#' @author Zhiwei Zhou, Xiaotao Shen
#' @param rt The RTs of peaks.
#' @param rt.tol The tolerance of RT.
#' @return  The index of peaks in each group.
#' @export

# rt = temp_rt
# rt.tol = 3

setGeneric(name = "groupRT",
           def = function(rt,
                          rt.tol = 3) {
             rt <- round(rt, 2)
             breaks <- seq(min(rt), max(rt) + rt.tol, by = rt.tol)
             w <- hist(rt, breaks = breaks, plot = FALSE)
             counts <- w$counts
             idx <- which(counts != 0)
             interval.range <- data.frame(breaks[idx],breaks[idx+1], stringsAsFactors = FALSE)
             interval.range <-
               data.frame(interval.range, c(1:nrow(interval.range)), stringsAsFactors = FALSE)
             colnames(interval.range) <- c("From", "To", "Module")

             rt.class <-
               sapply(rt, function(x) {
                 a <- x - round(interval.range$From, 2)
                 b <- round(interval.range$To, 2) - x
                 idx <- which(a >= 0 & b >= 0)
                 if(length(idx) > 0) {idx <- idx[1]}
                 idx
               })

             rt.class <-
               lapply(sort(unique(rt.class)), function(x) {
                 unname(which(rt.class == x))
               })
             return(rt.class)
           }
)


#' @title groupRT2
#' @description Group peaks according to RT with recursive function.
#' @author Zhiwei Zhou
#' @param rt The RTs of peaks.
#' @param rt.tol The tolerance of RT.
#' @return  The index of peaks in each group.
#' @export
groupRT2 <- function(rt_vec, rt_tol = 3) {
  # browser()
  group_label <- rep(-1, length(rt_vec))

  i <- 1
  temp_idx <- which(group_label < 0)
  while(length(temp_idx) > 0) {
    # cat(i, ' ')
    # browser()
    temp_idx <- temp_idx[1]
    rt_centric <- rt_vec[temp_idx]
    idx_label <- which(rt_vec < rt_centric + rt_tol & rt_vec >= rt_centric - rt_tol)
    group_label[idx_label] <- i

    temp_idx <- which(group_label < 0)
    i <- i + 1
  }

  result <- lapply(unique(group_label), function(x){
    which(group_label == x)
  })

  return(result)
}

#     assignConfidence ---------------------------------------------------------
#' @title assignConfidence
#' @description Assign confidence to metabolite ID group.
#' @author Zhiwei Zhou, Xiaotao Shen
#' @param x The metabolite ID group.
#' @return  The confidence of metabolite ID group.

setGeneric(name = "assignConfidence",
           def = function(x,
                          polarity = c("positive", "negative")){
             polarity <- match.arg(polarity)

             if(polarity == "positive"){
               prefer.adduct <- c("M+H", "M+Na", "M+NH4")
             }else{
               prefer.adduct <- c("M-H", "M+CH3COO", "M+Cl")
             }

             adduct <- x$adduct
             isotope <- x$isotope
             type <- x$type
             if(any(type == "seed")) return("grade1")
             temp.idx1 <- which(isotope %in% c("[M+1]","[M+2]","[M+3]","[M+4]"))
             if(length(temp.idx1) > 0) return("grade2")
             if(length(temp.idx1) == 0){
               temp.idx2 <- which(adduct %in% prefer.adduct)
               if(length(temp.idx2) > 0) {return("grade3")}
               if(length(temp.idx2) == 0) {return("grade4")}
             }
           })


#     removeRedundancy ---------------------------------------------------------
#' @title removeRedundancy
#' @description Remove annotation and peak redundancy.
#' @author Zhiwei Zhou, Xiaotao Shen
#' @param result The result from metA.
#' @param path The directory.
#' @return  A result after removing redundancy.
#' @export

# result = id_result_with_confidence
# path = path_output
# polarity = polarity
# test <- removeRedundancy(result = result,
#                          path = path_output,
#                          polarity = 'positive')

setGeneric(name = "removeRedundancy",
           def = function(result,
                          polarity = c("positive", "negative"),
                          path = "."){
             polarity <- match.arg(polarity)

             # redundancy
             redun <- list()
             redun[[1]] <- calculateRedundancy(result = result)
             count <- 1
             delta.redun <- 1
             cat("Round: ")
             while(delta.redun > 0.01 & count <= 5){
               cat(count); cat(" ")
               unique.id <- unique(result$to)
               result <- lapply(seq_along(unique.id), function(i){
                 temp.id <- unique.id[i]
                 temp.idx <- which(result$to == temp.id)
                 temp.result <- result[temp.idx,]
                 temp.group <- unique(temp.result$group)
                 temp.idx <- lapply(temp.group, function(x) which(temp.result$group == x))
                 temp.result <- lapply(temp.idx, function(x) temp.result[x,])

                 # confidence
                 temp.confidence <- unlist(lapply(temp.result, function(x) assignConfidence(x, polarity = polarity)))

                 # if the metabolite has more than two groups, and one group confidence is grade4,
                 #     then the group with the grade 4 is removed.

                 if(length(unique(temp.confidence)) > 1 & any(temp.confidence == "grade4")){
                   temp.result <- temp.result[-which(temp.confidence == "grade4")]
                   temp.confidence <- temp.confidence[-which(temp.confidence == "grade4")]
                 }

                 temp.result <- mapply(function(x, Confidence) {
                   # result <- data.frame(x, Confidence, stringsAsFactors = FALSE)
                   result <- x
                   result$Confidence <- Confidence
                   result
                   list(result)
                 },
                 x = temp.result,
                 Confidence = temp.confidence)

                 temp.result <- do.call(rbind, temp.result)
                 temp.result
               })

               result <- do.call(rbind, result)

               # remove the one peak vs many annotation:
               #    if one peak match more than two annotations,
               #    only the annotation which has the biggest grade is remained

               unique.peak <- unique(result$name)

               remain.idx <- lapply(seq_along(unique.peak), function(i){
                 temp.peak <- unique.peak[i]
                 temp.idx <- which(temp.peak == result$name)
                 if(length(temp.idx) == 1) return(temp.idx)
                 temp.result <- result[temp.idx,]
                 temp.confidence <- temp.result$Confidence
                 if(length(unique(temp.confidence)) > 1){
                   temp.grade <- as.numeric(substr(temp.confidence, 6, 6))
                   remain.idx <- temp.idx[which(temp.grade == temp.grade[which.min(temp.grade)])]
                   return(remain.idx)
                 }
                 return(temp.idx)
               })

               remain.idx <- unlist(remain.idx)

               result <- result[remain.idx,]

               # assign confidence to result again
               group <- unique(result$group)

               result <- lapply(group, function(temp.group){
                 temp.result <- result[result$group == temp.group,]
                 new.confidence <- assignConfidence(temp.result, polarity = polarity)
                 temp.result$Confidence <- new.confidence
                 temp.result
               })

               result <- do.call(rbind, result)
               redun <- c(redun, list(calculateRedundancy(result = result)))
               delta.redun <- abs(mean(redun[[length(redun)-1]]) - mean(redun[[length(redun)]]))
               count <- count + 1
             }
             cat("\n")
             save(redun, file = file.path(path, "redun"), compress = "xz")
             result
           })


#        calculateRedundancy --------------------------------------------------
#' @title calculateRedundancy
#' @description Calculate redundancy of result from metA.
#' @author Xiaotao Shen
#' @param result The result from metA.
#' @return  peak and metabolite redundancy.
# #' @export

setGeneric(name = "calculateRedundancy",
           def = function(result){
             unique.name <- unique(result$name)
             unique.group <- unique(result$group)
             unique.id <- unique(result$to)
             peak.redun <- nrow(result)/length(unique.name)
             id.redun <- length(unique.group)/length(unique.id)
             redundancy <- c(peak.redun, id.redun)
             names(redundancy) <- c("Peak.redundancy", "Metabolite.redundancy")
             redundancy
           })

#     convertResultMatrix2Tags ------------------------------------------------


#' @title convertResultMatrix2Tags
#' @description Filter tags2 data according to result.
#' @author Zhiwei Zhou, Xiaotao Shen
#' @param result The result.
#' @param tags2 The tags2 data.
#' @return The tags2 data.
#' @export

# test <- result2Tags(result = id_result_redun_rm, tags2 = new_tags)

setGeneric(name = "convertResultMatrix2Tags",
           def = function(result,
                          tags2 = tags2){

             peak.name <- showTags2(tags2, slot = "name")
             unique.name <- unique(result$name)

             # remove the annotation of peaks which are not in result
             idx <- which(!peak.name %in% unique.name)
             tags2[idx] <- lapply(tags2[idx], function(x) {
               x@annotation <- list()
               x
             })

             for(i in 1:length(unique.name)){
               temp.name <- unique.name[i]
               temp.idx1 <- which(temp.name == peak.name)
               temp.idx2 <- which(result$name == temp.name)
               temp.result <- result[temp.idx2,]
               remain.annotation <- apply(temp.result, 1, list)
               remain.annotation <- lapply(remain.annotation, function(x) {
                 paste(x[[1]][c("type", "From", "From.peak", "to", "step", "level", "as.seed",
                                "as.seed.round","isotope", "adduct", "charge", "Formula",
                                "mz.error", "rt.error","int.error","ms2.sim", "score")],
                       collapse = ";")
               })
               remain.annotation <- unlist(remain.annotation)

               annotation <- tags2[[temp.idx1]]@annotation
               annotation <- lapply(annotation, function(x) {
                 x <- unlist(x)
                 paste(x[c("type", "From", "From.peak", "to", "step", "level", "as.seed",
                           "as.seed.round","isotope", "adduct", "charge", "Formula",
                           "mz.error", "rt.error","int.error","ms2.sim", "score")],
                       collapse = ";")
               })
               annotation <- unlist(annotation)

               remain.idx <- which(annotation %in% remain.annotation)
               tags2[[temp.idx1]]@annotation <- tags2[[temp.idx1]]@annotation[remain.idx]
             }

             tags2
           })



#     getAnnotationResult ------------------------------------------------------

#' @title getAnnotationResult
#' @author Zhiwei Zhou, Xiaotao Shen
#' @param tags2
#' @param tags.result2
#' @param metabolite
#' @param instrument
#' @param candidate.num
#' @param score.cutoff
# tags2 = tags2_after_redundancy_remove
# tags.result2 = id_result_redun_rm
# kegg.compound = cpd_kegg
# candidate.num = candidate_num
#
# temp <- getAnnotationResult(tags2 = tags2_after_redundancy_remove,
#                             tags.result2 = id_result_redun_rm,
#                             kegg.compound = cpd_kegg,
#                             candidate.num = candidate_num)

setGeneric(name = "getAnnotationResult",
           def = function(tags2,
                          tags.result2,
                          metabolite,
                          instrument = c("SciexTripleTOF", "AgilentQTOF", "BrukerQTOF", "ThermoOrbitrap",
                                         'ThermoExploris', "AgilentDTIMMS", "BrukerTIMS", 'WatersQTOF','WatersTWIMMS'),
                          candidate.num = 3000,
                          score.cutoff = 0.4,
                          ...){

             temp <- convertTags2ResultMatrix(tags2 = tags2,
                                              score.cutoff = score.cutoff,
                                              candidate.num = candidate.num)

             colnames(temp)[8] <- "ID"

             confidence <- mapply(function(x, y){
               result <- tags.result2 %>%
                 dplyr::filter(name == x & score >= score.cutoff)

               if(length(result) == 0) return(NA)

               temp_cpd <- stringr::str_split(y, pattern = ';')[[1]]
               idx <- match(temp_cpd, result$to)
               result <- result$Confidence[idx] %>% stringr::str_c(collapse = ';')
               result
             },
             x = temp$name,
             y = temp$ID)

             peak.group <- mapply(function(x, y){
               result <- tags.result2 %>%
                 dplyr::filter(name == x & score >= score.cutoff)

               if(length(result) == 0) return(NA)

               temp_cpd <- stringr::str_split(y, pattern = ';')[[1]]
               idx <- match(temp_cpd, result$to)
               result <- result$group[idx] %>% stringr::str_c(collapse = ';')
               result
             },
             x = temp$name,
             y = temp$ID)

             confidence[is.na(temp$ID)] <- NA
             peak.group[is.na(temp$ID)] <- NA

             id <- temp$ID
             compound.name <- unlist(lapply(id, function(x){
               # cat(x, ' ')
               if (is.na(x)){
                 return(NA)
               } else {
                 x <- strsplit(x = x, split = ";")[[1]]
                 temp.name <- metabolite$name[match(x, metabolite$id)]
                 # if (length(temp.name) == 0) {return(NA)}
                 temp.name <- unlist(lapply(strsplit(temp.name, split = ";"), function(x) x[1]))
                 paste(temp.name, collapse = ";")
               }
             }))


             exported_report <- temp %>%
               tibble::as_tibble() %>%
               dplyr::mutate(compound.name = compound.name,
                             peak.group = peak.group,
                             confidence = confidence)

             if (instrument %in% c("AgilentDTIMMS", "BrukerTIMS", "WatersTWIMMS")) {
               exported_report <- exported_report %>%
                 dplyr::rename('Annotation.type' = type,
                               'annotated.from.ID' = From,
                               'annotated.from.peak' = From.peak) %>%
                 dplyr::select(name, mz, rt, ccs, Annotation.type, annotated.from.ID, annotated.from.peak, ID,
                               score, isotope, adduct, Formula, matched.frag, compound.name, peak.group, confidence)


             } else {
               exported_report <- exported_report %>%
                 dplyr::rename('Annotation.type' = type,
                               'annotated.from.ID' = From,
                               'annotated.from.peak' = From.peak) %>%
                 dplyr::select(name, mz, rt, Annotation.type, annotated.from.ID, annotated.from.peak, ID,
                               score, isotope, adduct, Formula, matched.frag, compound.name, peak.group, confidence)

             }

             return(exported_report)

           })





#        convertTags2ResultMatrix ----------------------------------------------

#' @title convertTags2ResultMatrix
#' @description Transform tags2 data to matrix like ms2Annotation.
#' @author Zhiwei Zhou, Xiaotao Shen
#' @param tags2 The tags2 data.
#' @param score.cutoff The score cutoff.
#' @param candidate.num The number of candidates.
#' @return The matrix of tags2 data.
#' @export

setGeneric(name = "convertTags2ResultMatrix",
           def = function(tags2,
                          score.cutoff = 0,
                          candidate.num = 3000){
             cat("Transform tags2 to matrix:\n")
             pbapply::pboptions(type = "timer", style = 1)
             result <- pbapply::pblapply(seq_along(tags2), function(i) {
               # cat(i, ' ')
               x <- tags2[[i]]
               name <- x@name
               mz <- x@mz
               rt <- x@rt
               ccs <- x@ccs
               annotation <- x@annotation
               if (length(annotation) == 0){
                 # c(name, mz, rt, rep(NA, 17))
                 result <- data.frame(name, mz, rt, ccs,
                                      NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, ##18 NAs
                                      stringsAsFactors = FALSE)
                 colnames(result) <- c("name", "mz", "rt", 'ccs', "type", "From", "From.peak", "to", "step",
                                       "level", "as.seed", "as.seed.round", "isotope", "adduct", "charge", "Formula",
                                       "mz.error", "rt.error", 'ccs.error', "int.error", "ms2.sim", 'matched.frag', "score")
                 result
               } else {
                 type <- unlist(lapply(annotation, function(x) x$type))
                 From <- unlist(lapply(annotation, function(x) x$From))
                 From.peak <- unlist(lapply(annotation, function(x) x$From.peak))
                 to <- unlist(lapply(annotation, function(x) x$to))
                 step <- unlist(lapply(annotation, function(x) x$step))
                 level <- unlist(lapply(annotation, function(x) x$level))
                 as.seed <- unlist(lapply(annotation, function(x) x$as.seed))
                 as.seed.round <- unlist(lapply(annotation, function(x) x$as.seed.round))
                 isotope <- unlist(lapply(annotation, function(x) x$isotope))
                 adduct <- unlist(lapply(annotation, function(x) x$adduct))
                 charge <- unlist(lapply(annotation, function(x) x$charge))
                 Formula <- unlist(lapply(annotation, function(x) x$formula))
                 mz.error <- unlist(lapply(annotation, function(x) x$mz.error))
                 rt.error <- unlist(lapply(annotation, function(x) x$rt.error))
                 ccs.error <- unlist(lapply(annotation, function(x) x$ccs.error))
                 int.error <- unlist(lapply(annotation, function(x) x$int.error))
                 if (is.null(int.error)) {int.error <- unlist(lapply(annotation, function(x) x$int.ratio.error))}
                 ms2.sim <- unlist(lapply(annotation, function(x) x$ms2.sim))
                 matched.frag <- unlist(lapply(annotation, function(x) x$nfrag))
                 score <- unlist(lapply(annotation, function(x) x$score))
                 result <- data.frame(name, mz, rt, ccs, type, From, From.peak, to, step,
                                      level, as.seed,
                                      as.seed.round, isotope, adduct, charge, Formula,
                                      mz.error, rt.error, ccs.error, int.error, ms2.sim, matched.frag,
                                      score, stringsAsFactors = FALSE)
                 result <- result[!duplicated(result),]
                 result <- result[order(result$score, decreasing = TRUE),]##order result according to score
                 result <- result[result$score >= score.cutoff,]
                 if(nrow(result) > candidate.num) result <- result[1:candidate.num,]
                 if(nrow(result) == 0){return(c(name, mz, rt, ccs, rep(NA, 19)))}
                 type <- paste(result$type, collapse = ";")
                 From <- paste(result$From, collapse = ";")
                 From.peak <- paste(result$From.peak, collapse = ";")
                 to <- paste(result$to, collapse = ";")
                 step <- paste(result$step, collapse = ";")
                 level <- paste(result$level, collapse = ";")
                 as.seed <- paste(result$ as.seed, collapse = ";")
                 as.seed.round <- paste(result$as.seed.round, collapse = ";")
                 as.seed.round <- paste(result$as.seed.round, collapse = ";")
                 isotope <- paste(result$isotope, collapse = ";")
                 adduct <- paste(result$adduct, collapse = ";")
                 charge <- paste(result$charge, collapse = ";")
                 Formula <- paste(result$Formula, collapse = ";")
                 mz.error <- paste(result$mz.error, collapse = ";")
                 rt.error <- paste(result$rt.error, collapse = ";")
                 ccs.error <- paste(result$ccs.error, collapse = ";")
                 int.error <- paste(result$int.error, collapse = ";")
                 ms2.sim <- paste(result$ms2.sim, collapse = ";")
                 matched.frag <- paste(result$matched.frag, collapse = ';')
                 score <- paste(result$score, collapse = ";")
                 result <- data.frame(name, mz, rt, ccs, type, From, From.peak, to, step,
                                      level, as.seed,
                                      as.seed.round, isotope, adduct, charge, Formula,
                                      mz.error, rt.error, ccs.error, int.error, ms2.sim, matched.frag,
                                      score, stringsAsFactors = FALSE)
                 rm(list = c("type", "From", "From.peak", "to", "step", "level",
                             "as.seed", "as.seed.round", "isotope", "adduct",
                             "charge", "Formula", "mz.error", "rt.error", 'ccs.error',
                             "int.error", "ms2.sim", 'matched.frag', "score"))
                 return(result)
               }
             })

             rm(list = "tags2")
             result <- do.call(rbind, result)
             result <- as.data.frame(result)
             result <- result
           })
