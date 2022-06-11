################################################################################
# MetDNA2 ----------------------------------------------------------------------

#' @title MetDNA2
#' @description Metabolite annotation and dysregulated network analysis.
#' @author Zhiwei Zhou
#' \email{zhouzw@@sioc.ac.cn}
#' @param ms1_file The name of ms1 data
#' @param ms2_file The names of ms2 data
#' @param sample_info_file 'sample.info.csv'
#' @param metdna_version MetDNA parameter 'version1' and 'version2'
#' @param ms2_type 'mgf', 'mzXML', 'msp', 'cef' are supported. Default: 'mgf'
#' @param path Default: '.'
#' @param thread The number of threads. Default: 4
#' @param is_check_data Whether check data quality. Default: TRUE
#' @param is_anno_initial_seed Whether annotate initial seeds with spectral DB. Default: TRUE
#' @param is_anno_mrn Whether perform recursive annotation. Default: TRUE
#' @param is_credential Whether perform annotation credential. Default: TRUE
#' @param is_cred_pg_filer Whether filter error annotations by peak group annotaion in final table. Default: FALSE
#' @param is_cred_formula_filter Whether filter error annotations by formula prediction in final table. Default: FALSE
#' @param lib Default: 'zhumetlib_qtof'
#' @param polarity 'positive', 'negative'. Default: 'positive'
#' @param instrument The instrument you used to acquire data. "SciexTripleTOF", "AgilentQTOF", "BrukerQTOF", "ThermoOrbitrap",
#' 'ThermoExploris', "AgilentDTIMMS", "BrukerTIMS", 'WatersQTOF','WatersTWIMMS'. Default: "SciexTripleTOF"
#' @param column "hilic", "rp". Default: 'hilic'
#' @param ce "30", "10", "20", "35,15", "40", "50". Default: "30"
#' @param method_lc 'Amide12min', 'Amide23min', 'Other', 'MetlinRP', 'RP12min'. Default: "Amide12min"
#' @param is_rt_calibration Only valid with 'Amide12min' and 'Amide23min' methods. Default: TRUE.
#' @param mz_tol The tolerance of ms1 match. Default: 25 ppm
#' @param tolerance_rt_range The maxmium rt tolerance. Default: 30 s
#' @param dp_cutoff MS2 score cutoff. Default: 0.8
#' @param direction Reserve which direct dot product. 'reverse', 'forward'. Default: 'reverse'
#' @param is_plot_ms2 Whether plot initial-seed ms2 match plot. Default: TRUE
#' @param extension_step The used extended MRN with n step reaction. Steps: '0', '1', '2', '3', '4', '5', '6', '7', '8'
#' @param max_isotope Maximum considerable isotope. Default: 4
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
#' @param candidate_num reserve top N candidate number in recursive annotation. Default: 3
#' @param dir_GenForm
#' @param is_pred_formula_all whether predict formula for all peak group (TRUE) or concised peak groups (FALSE). Default: FALSE
#' @param formula_mz_tol_ms1 Default: 25 ppm
#' @param formula_mz_tol_ms2 Default: 25 ppm
#' @param formula_elements Default: "CHNOPS"
#' @param formula_candidate_num Reserve top 3 formula. Default: 3
#' @param platform 'windows', 'linux'. Default: 'linux'
#' @export


# ms1_file = "data.csv"
# ms2_file = NULL
# sample_info_file = "sample.info.csv"
# metdna_version = 'version2'
# ms2_type = 'mgf'
# ### path = "H:/00_projects/03_MetDNA2/00_data/20200922_metdna2_development/pos_201009_compare_metdna2_raw"
# path = "/home/zhouzw/Data_processing/20201214_metdna2_development/"

# # para of loadDB
# polarity = 'positive'
# instrument = "SciexTripleTOF"
# lib = 'zhuMetLib'
# column = 'hilic'
# ce = '30'
# method_lc = 'Amide23min'
# excluded_adduct = NULL
# is_rt_calibration = TRUE

# # para of matchMs1WithSpecLib
# mz_tol = 25
# pf_rt_range = 0
# tolerance_rt_range = 30
# pf_ccs_range = 0
# tolerance_ccs_range = 2
# is_filter = FALSE
# is_rt_score = TRUE
# is_ccs_score = FALSE
# is_msms_score = TRUE

# # parameters of matchMs2WithSpecLib
# is_include_precursor = TRUE
# is_deisotope = FALSE
# int_ms2_min_abs = 0
# int_ms2_min_relative = 0.01
# ppm_precursor_filter = 20
# mz_range_ms2 = c(0, 1700)
#
# mz_tol_combine_ms1_ms2 = 25 # ppm
# rt_tol_combine_ms1_ms2 = 10 # s
# ccs_tol_combine_ms1_ms2 = NULL # %
# mz_tol_ms2 = 35
# dp_cutoff = 0.8
# direction = 'reverse'
#
# is_plot_ms2 <- TRUE

##
# polarity = c("positive", "negative")
# extension_step = '0'
# threads = 3
# max_isotope = 4
# rt_tol1 = 3
# rt_tol2 = 30
# cor_tol = 0
# int_tol = 500
# dp_tol = 0.5
# max_step = 3
# score_cutoff = 0
# remain = FALSE
# remain_per = 0.5
# seed_neighbor_match_plot = TRUE
# candidate_num = 3
# scoring_approach_recursive = 'bonanza'
# matched_frag_cutoff = 1


# is_cred_pg_filter = FALSE
# is_cred_formula_filter = FALSE


# 210304
# ms1_file = "data.csv"
# ms2_file = NULL
# sample_info_file = "sample.info.csv"
# metdna_version = 'version2'
# ms2_type = "mgf"
# path = "/home/zhouzw/Data_processing/20210224_metdna2_update/test_dp_bonanza/"
# thread = 3
# is_check_data = TRUE
# is_anno_initial_seed = TRUE
# is_anno_mrn = TRUE
# is_credential = TRUE
# is_cred_pg_filter = FALSE
# is_cred_formula_filter = FALSE
#
# # parameters of loadDB
# lib = 'zhuMetLib'
# polarity = "positive"
# instrument = "SciexTripleTOF"
# column = "hilic"
# ce = "30"
# method_lc = 'Amide23min'
# is_rt_calibration = TRUE
#
# # parameters of matchMs1WithSpecLib
# mz_tol = 25
# tolerance_rt_range = 30
# dp_cutoff = 0.8
# direction = 'reverse'
# is_plot_ms2 = TRUE
#
# ## para of RecursiveAnnotationMRN
# extension_step = '2'
# max_isotope = 4
# rt_tol1 = 3
# rt_tol2 = 30
# cor_tol = 0
# int_tol = 500
# dp_tol = 0.5
# max_step = 3
# score_cutoff = 0
# remain = FALSE
# remain_per = 0.5
# seed_neighbor_match_plot = TRUE
# candidate_num = 3
# scoring_approach_recursive = 'bonanza'
# matched_frag_cutoff = 1
#
# # annotation credential
# dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release'
# is_pred_formula_all = FALSE
# formula_mz_tol_ms1 = 25 # ms1 tolerance
# formula_mz_tol_ms2 = 25 # ms2 tolerance
# formula_elements = "CHNOPS"
# formula_candidate_num = 3 # return top 3 candidate
# platform = 'linux'


# MetDNA2(ms1_file = 'data.csv',
#         sample_info_file = "sample.info.csv",
#         metdna_version = 'version2',
#         ms2_type = 'mgf',
#         path = '/home/zhouzw/Data_processing/20210224_metdna2_update/',
#         thread = 3,
#         polarity = 'positive',
#         instrument = "SciexTripleTOF",
#         column = "hilic",
#         lib = 'zhuMetLib',
#         ce = "30",
#         method_lc = 'Amide23min',
#         extension_step = '2',
#         scoring_approach_recursive = 'bonanza',
#         is_check_data = FALSE,
#         is_anno_initial_seed = FALSE,
#         is_anno_mrn = FALSE,
#         is_credential = TRUE)
#
# MetDNA2(ms1_file = 'data.csv',
#         sample_info_file = "sample.info.csv",
#         metdna_version = 'version2',
#         ms2_type = 'mgf',
#         path = '/home/zhouzw/Data_processing/20210402_aging_fly_report_demonstration/',
#         thread = 3,
#         polarity = 'positive',
#         instrument = "SciexTripleTOF",
#         column = "hilic",
#         lib = 'zhuMetLib',
#         ce = "30",
#         method_lc = 'Amide23min',
#         extension_step = '0',
#         scoring_approach_recursive = 'dp',
#         is_check_data = FALSE,
#         is_anno_initial_seed = FALSE,
#         is_anno_mrn = FALSE,
#         is_credential = TRUE,
#         is_bio_interpret = TRUE)


setGeneric(name = 'MetDNA2',
           def = function(
             ms1_file = "data.csv",
             ms2_file = NULL,
             sample_info_file = "sample.info.csv",
             metdna_version = c('version2', 'version1'),
             ms2_type = c("mgf", "mzXML", "msp", "cef") ,
             path = ".",
             thread = 4,
             is_check_data = TRUE,
             is_anno_initial_seed = TRUE,
             is_anno_mrn = TRUE,
             is_credential = TRUE,
             is_bio_interpret = TRUE,
             is_exported_report = TRUE,
             is_cred_pg_filter = TRUE,
             is_cred_formula_filter = FALSE,
             is_rm_intermediate_data = FALSE,

             # parameters of loadDB
             lib = c('zhumetlib_qtof', 'zhumetlib_orbitrap', 'fiehnHilicLib'),
             polarity = c("positive", "negative"),
             instrument = c("SciexTripleTOF", "AgilentQTOF", "BrukerQTOF", 'WatersQTOF', "ThermoOrbitrap", 'ThermoExploris',
                            "AgilentDTIMMS", "BrukerTIMS", 'WatersTWIMMS'),
             column = c("hilic", "rp"),
             ce = c("30", "10", "20", "35,15", "40", "50",
                    'NCE10', 'NCE20', 'NCE30', 'NCE40', 'NCE50',
                    'SCE20_30_40%', "SNCE20_30_40%"),
             method_lc = c('Amide12min', 'Amide23min', 'Other', 'MetlinRP', 'RP12min'),
             is_rt_calibration = TRUE,

             # parameters of matchMs1WithSpecLib
             mz_tol = 25, # ppm
             mz_ppm_thr = 400,
             tolerance_rt_range = 30, # s
             tolerance_ccs_range = 4, # %
             int_ms2_min_abs = 50,
             dp_cutoff = 0.8,
             direction = c('reverse', 'forward'),
             mz_tol_combine_ms1_ms2 = 25,
             rt_tol_combine_ms1_ms2 = 10,
             is_plot_ms2 = TRUE,
             adduct_limit_specLib = 'all', # 220115
             test_force_filtering_rt = NULL,

             ## para of RecursiveAnnotationMRN
             extension_step = c('0', '1', '2', '3', '4', '5', '6', '7', '8'),
             max_isotope = 4,
             rt_tol1 = 3, # s
             rt_tol2 = 30, # %
             ccs_tol = 4, # %
             cor_tol = 0,
             int_tol = 500,
             dp_tol = 0.5,
             max_step = 3,
             score_cutoff = 0,
             remain = FALSE,
             remain_per = 0.5,
             seed_neighbor_match_plot = TRUE,
             candidate_num = 5,
             scoring_approach_recursive = c('dp', 'bonanza', 'hybrid', 'gnps'),
             matched_frag_cutoff = 1,
             whether_link_frag = FALSE,
             use_pretrained_model = FALSE,
             # whether_use_predRT = TRUE,
             # test_old_mrn = c('v0.5', 'v0.3'),
             # test_adduct_version = c('version2', 'version1'),

             # annotation credential
             use_redun_rm_result = TRUE,
             dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
             is_pred_formula_all = FALSE,
             formula_mz_tol_ms1 = 25, # ms1 tolerance
             formula_mz_tol_ms2 = 25, # ms2 tolerance
             formula_elements = "CHNOPS",
             formula_candidate_num = 3, # return top 3 candidate
             platform = 'linux',
             test_rm_extra_anno_from_ini_seed = FALSE, # whether assign confidence level2 for adduct and isotopes from initial seed
             is_plot_pseudo_MS1 = TRUE,

             # result table export
             remove_conflict_seed_from_final_table = FALSE,

             # biology interpretation
             comp_group = c("W03", "W30"),
             uni_test = c("t","wilcox","anova"),
             correct_p = FALSE,
             p_cutoff = 0.05,
             fc_cutoff = 1,
             species = c("hsa", "mmu", "rno", "bta", "gga",
                         "dre", "dme", "cel", "sce", "osa",
                         "ath", "smm", "pfa", "tbr", "eco",
                         "bsu", "ppu", "sau", "tma", "syf", "mln"),
             quanti_pathway_method = c('mean', 'sum', 'median'),

             test_evaluation = c('No', '200STD', '46STD')

           ){
             # parameter check and decision
             metdna_version <- match.arg(metdna_version)
             ms2_type <- match.arg(ms2_type)
             instrument <- match.arg(instrument)
             polarity <- match.arg(polarity)
             column <- match.arg(column)
             ce <- match.arg(ce)
             lib <- match.arg(lib)
             method_lc <- match.arg(method_lc)
             direction <- match.arg(direction)
             scoring_approach_recursive <- match.arg(scoring_approach_recursive)
             uni_test <- match.arg(uni_test)
             species <- match.arg(species)
             quanti_pathway_method <- match.arg(quanti_pathway_method)
             test_evaluation <- match.arg(test_evaluation)
             # test_old_mrn <- match.arg(test_old_mrn)
             # test_adduct_version <- match.arg(test_adduct_version)

             if (test_evaluation == '200STD') {
               is_bio_interpret <- FALSE
             }
             if (test_evaluation == '46STD') {
               is_bio_interpret <- FALSE
             }
             
             # create run.log file
             cat("MetDNA2 run log\n", file = file.path(path, "run.log.txt"))
             cat('MetDNA2 verison:', as.character(packageVersion('MetDNA2')), '\n', append = TRUE)
             cat(as.character(Sys.time()), file = file.path(path, "run.log.txt"), append = TRUE)
             cat("\n", file = file.path(path, "run.log.txt"), append = TRUE)

             # check parameters data
             cat("1. Check parameters and data.\n")
             cat("1. Check parameters and data.\n", file = file.path(path, "run.log.txt"), append = TRUE)
             cat("------------------------------------------------------------\n")

             # check parameters: required
             check_result_para <- try({
               checkPara(lib = lib,
                         polarity = polarity,
                         instrument = instrument,
                         column = column,
                         ce = ce,
                         method_lc = method_lc,
                         is_rt_calibration = is_rt_calibration,
                         path = path)
             }, silent = TRUE)

             if(class(check_result_para) == "try-error"){
               cat(check_result_para[[1]], file = file.path(path, "run.log.txt"), append = TRUE)
               stop(check_result_para[[1]])
             }

             # save the parameters
             list_para <- list(
               'package_verison' = packageVersion('MetDNA2'),
               'ms1_file' = ms1_file,
               'ms2_file' = paste(ms2_file, collapse = ';'),
               'sample_info_file' = sample_info_file,
               'metdena_version' = metdna_version,
               'path' = path,
               'thread' = thread,
               'is_check_data' = is_check_data,
               'is_anno_initial_seed' = is_anno_initial_seed,
               'is_anno_mrn' = is_anno_mrn,
               'is_credential' = is_credential,
               'is_bio_interpret' = is_bio_interpret,
               'is_exported_report' = is_exported_report,
               'is_cred_pg_filter' = is_cred_pg_filter,
               'is_cred_formula_filter' = is_cred_formula_filter,

               'lib' = lib,
               'polarity' = polarity,
               'instrument' = instrument,
               'column' = column,
               'ce' = ce,
               'method_lc' = method_lc,
               'is_rt_calibration' = is_rt_calibration,

               # parameters of matchMs1WithSpecLib
               'mz_tol' = mz_tol,
               'mz_ppm_thr' = mz_ppm_thr,
               'tolerance_rt_range' = tolerance_rt_range,
               'tolerance_ccs_range' = tolerance_ccs_range,
               'dp_cutoff' = dp_cutoff,
               'direction' = direction,
               'is_plot_ms2' = is_plot_ms2,
               'mz_tol_combine_ms1_ms2' = mz_tol_combine_ms1_ms2,
               'rt_tol_combine_ms1_ms2' = rt_tol_combine_ms1_ms2,
               'int_ms2_min_abs' = int_ms2_min_abs,
               'adduct_limit_specLib' = adduct_limit_specLib,

               ## para of RecursiveAnnotationMRN
               'extension_step' = extension_step,
               'max_isotope' = max_isotope,
               'rt_tol1' = rt_tol1,
               'rt_tol2' = rt_tol2,
               'cor_tol' = cor_tol,
               'int_tol' = int_tol,
               'dp_tol' = dp_tol,
               'max_step' = max_step,
               'score_cutoff' = score_cutoff,
               'remain' = remain,
               'remain_per' = remain_per,
               'seed_neighbor_match_plot' = seed_neighbor_match_plot,
               'candidate_num' = candidate_num,
               'scoring_approach_recursive' = scoring_approach_recursive,
               'matched_frag_cutoff' = matched_frag_cutoff,
               'whether_link_frag' = whether_link_frag,
               'use_pretrained_model' = use_pretrained_model,


               # annotation credential
               'use_redun_rm_result' = use_redun_rm_result,
               'dir_GenForm' = dir_GenForm,
               'is_pred_formula_all' = is_pred_formula_all,
               'formula_mz_tol_ms1' = formula_mz_tol_ms1, # ms1 tolerance
               'formula_mz_tol_ms2' = formula_mz_tol_ms2, # ms2 tolerance
               'formula_elements' = formula_elements,
               'formula_candidate_num' = formula_candidate_num, # return top 3 candidate
               'platform' = platform,
               'is_plot_pseudo_MS1' = is_plot_pseudo_MS1,

               # result table export
               'remove_conflict_seed_from_final_table' = remove_conflict_seed_from_final_table,

               # biology interpretation
               'comp_group' = paste(comp_group, collapse = ';'),
               'uni_test' = uni_test,
               'correct_p' = correct_p,
               'p_cutoff' = p_cutoff,
               'fc_cutoff' = fc_cutoff,
               'species' = species,
               'quanti_pathway_method' = quanti_pathway_method,

               # test parameters
               # 'test_old_mrn' = test_old_mrn,
               # 'test_adduct_version' = test_adduct_version,
               'test_rm_extra_anno_from_ini_seed' = test_rm_extra_anno_from_ini_seed,
               'test_evaluation' = test_evaluation
               ) %>%
               tibble::as_tibble() %>%
               dplyr::mutate_all(as.character) %>%
               tidyr::pivot_longer(cols = dplyr::everything(), names_to = 'para')

             dir.create(file.path(path),
                        showWarnings = FALSE,
                        recursive = TRUE)

             readr::write_tsv(list_para, path = file.path(path, 'para_list.txt'))




             if (is_check_data){
               # change the ms2_type to mgf if there are mgf files, and change ms2_type to msp if there are msp files.

               # switch(ms2_type,
               #        "mgf" = {ms2.file <- grep("mgf", list.files(path), value = TRUE)},
               #        "msp" = {ms2.file <- grep("msp", list.files(path), value = TRUE)},
               #        "mzXML" = {ms2.file <- grep("mzXML", list.files(path), value = TRUE)},
               #        "cef" = {ms2.file <- grep("cef", list.files(path), value = TRUE)})

               switch(ms2_type,
                      'mgf' = {
                        ms2_file <- list.files(path = path,
                                               pattern = '.+\\.mgf$',
                                               recursive = TRUE,
                                               ignore.case = TRUE)
                        # ms2_file <- file.path(path, ms2_file)
                      },
                      'msp' = {
                        ms2_file <- list.files(path = path,
                                               pattern = '.+\\.msp$',
                                               recursive = TRUE,
                                               ignore.case = TRUE)
                        # ms2_file <- file.path(path, ms2_file)
                      },
                      'mzXML' = {
                        ms2_file <- list.files(path = path,
                                               pattern = '.+\\.mzXML$',
                                               recursive = TRUE,
                                               ignore.case = TRUE)
                        # ms2_file <- file.path(path, ms2_file)
                      },
                      'cef' = {
                        ms2_file <- list.files(path = path,
                                               pattern = '.+\\.cef$',
                                               recursive = TRUE,
                                               ignore.case = TRUE)
                        # ms2_file <- file.path(path, ms2_file)
                      })


               # check group
               if (is_bio_interpret) {
                 group <- checkQualitySampleGroup(ms1_file = ms1_file, sample_info_file = sample_info_file, group = comp_group, path = path)
               }


               check_result <- try({checkQuality(ms1_file = ms1_file,
                                                 sample_info_file = sample_info_file,
                                                 ms2_type = ms2_type,
                                                 instrument = instrument,
                                                 path = path)}, silent = TRUE)

               if(class(check_result) == "try-error"){
                 cat(check_result[[1]], file = file.path(path, "run.log.txt"), append = TRUE)
                 stop(check_result[[1]])
               }

               if(any(as.numeric(check_result[,4]) > 0)){
                 cat('Error: Please check your data to make sure that they are valid.\n', file = file.path(path, "run.log.txt"), append = TRUE)
                 stop('Please check your data to make sure that they are valid.\n')
               }
             } else {
               cat('Skip data check.\n', file = file.path(path, "run.log.txt"), append = TRUE)
               cat("Skip data check.\n")
             }


             # annotateInitialSeed
             cat("\n")
             cat('2. Initial seed annotation.\n', file = file.path(path, "run.log.txt"), append = TRUE)
             cat("2. Initial seed annotation.\n")
             cat("------------------------------------------------------------\n")
             if (is_anno_initial_seed) {

               # ms2Annotation(ms1.file = ms1.data,
               #               sample.info = old.sample.info.name,
               #               mz.tol = mz.tol,
               #               rt.tol = rt.tol.for.ms1.ms2.match,
               #               path = path,
               #               output.path = file.path(path, "MS2_match_result"),
               #               instrument = instrument,
               #               ms2.type = ms2.type,
               #               polarity = polarity,
               #               column = column,
               #               ce = ce,
               #               ms2.match.plot = ms2.match.plot)

               # annotateInitialSeed(ms1_file = ms1_file,
               #                     ms2_file = ms2_file,
               #                     metdna_version = metdna_version,
               #                     ms2_type = ms2_type,
               #                     path = path,
               #                     polarity = polarity,
               #                     instrument = instrument,
               #                     lib = lib,
               #                     column = column,
               #                     ce = ce,
               #                     method_lc = method_lc,
               #                     is_rt_calibration = is_rt_calibration,
               #                     mz_tol = mz_tol,
               #                     tolerance_rt_range = tolerance_rt_range,
               #                     is_filter = FALSE,
               #                     is_rt_score = TRUE,
               #                     dp_cutoff = 0.8,
               #                     is_plot_ms2 = TRUE)
               if (instrument %in% c("AgilentDTIMMS", "BrukerTIMS", "WatersTWIMMS")) {
                 ccs_tol_combine_ms1_ms2 <- 0.5
                 is_ccs_score <- TRUE
               } else {
                 ccs_tol_combine_ms1_ms2 <- NULL
                 is_ccs_score <- FALSE
               }

               if (adduct_limit_specLib == 'all') {
                 excluded_adduct <- NULL
               } else {
                 if (polarity == 'positive') {
                   excluded_adduct <- lib_adduct_nl$positive %>% dplyr::filter(annotation == 'Yes' | adduct %in% c('[M]+')) %>% dplyr::pull(adduct)
                   excluded_adduct <- excluded_adduct[!(excluded_adduct %in% adduct_limit_specLib)]
                 } else {
                   excluded_adduct <- lib_adduct_nl$negative %>% dplyr::filter(annotation == 'Yes' | adduct %in% c('[M-2H]-')) %>% dplyr::pull(adduct)
                   excluded_adduct <- excluded_adduct[!(excluded_adduct %in% adduct_limit_specLib)]
                 }
               }

               temp_error <- try({annotateInitialSeed(ms1_file = ms1_file,
                                                      ms2_file = ms2_file,
                                                      metdna_version = metdna_version,
                                                      ms2_type = ms2_type,
                                                      path = path,
                                                      polarity = polarity,
                                                      instrument = instrument,
                                                      lib = lib,
                                                      column = column,
                                                      ce = ce,
                                                      method_lc = method_lc,
                                                      is_rt_calibration = is_rt_calibration,
                                                      mz_tol = mz_tol,
                                                      mz_ppm_thr = mz_ppm_thr,
                                                      tolerance_rt_range = tolerance_rt_range,
                                                      tolerance_ccs_range = tolerance_ccs_range,
                                                      mz_tol_combine_ms1_ms2 = mz_tol_combine_ms1_ms2,
                                                      rt_tol_combine_ms1_ms2 = rt_tol_combine_ms1_ms2,
                                                      ccs_tol_combine_ms1_ms2 = ccs_tol_combine_ms1_ms2,
                                                      is_filter = TRUE,
                                                      is_rt_score = TRUE,
                                                      is_ccs_score = is_ccs_score,
                                                      int_ms2_min_abs = int_ms2_min_abs,
                                                      dp_cutoff = dp_cutoff,
                                                      direction = direction,
                                                      scoring_approach = 'dp',
                                                      is_plot_ms2 = is_plot_ms2,
                                                      test_adduct_version = ifelse(metdna_version == 'version2', 'version2', 'version1'),
                                                      test_evaluation = test_evaluation,
                                                      excluded_adduct = excluded_adduct,
                                                      test_force_filtering_rt = test_force_filtering_rt)},
                                 silent = TRUE)

               if(class(temp_error) == "try-error") {
                 cat(temp_error[[1]], file = file.path(path, "run.log.txt"), append = TRUE)
                 stop(temp_error[[1]])
               }

             }else{
               cat('Skip MS/MS match annotation.\n', file = file.path(path, "run.log.txt"), append = TRUE)
               cat("Skip MS/MS match annotation.\n")
             }

             # metabolic reaction network based recursive annotation
             if(is_anno_initial_seed) cat("\n")
             cat("\n")
             cat('3. Metabolic reaction network based metabolite annotation.\n', file = file.path(path, "run.log.txt"), append = TRUE)
             cat("3. Metabolic reaction network based metabolite annotation.\n")
             cat("------------------------------------------------------------\n")
             if(is_anno_mrn){

               # annotateMRN(annotation_result = "ms2_match_annotation_result.csv",
               #             ms2_data = 'ms2',
               #             prefer_adduct = "M+H",
               #             metdna_version = 'version2',
               #             lib = 'zhuMetLib_obitrap',
               #             direction = 'reverse',
               #             column = 'hilic',
               #             instrument = "ThermoOrbitrap",
               #             method_lc = 'Amide12min',
               #             polarity = 'positive',
               #             extension_step = '2',
               #             threads = 3,
               #             max_isotope = 4,
               #             path = '/home/zhouzw/Data_processing/20210319_different_instrument_test/QE/',
               #             mz_tol = 25,
               #             rt_tol1 = 3,
               #             rt_tol2 = 30,
               #             ccs_tol = 4,
               #             cor_tol = 0,
               #             int_tol = 500,
               #             dp_tol = 0.5,
               #             max_step = 3,
               #             score_cutoff = 0,
               #             remain = FALSE,
               #             remain_per = 0.5,
               #             scoring_approach = 'dp',
               #             matched_frag_tol = 1,
               #             seed_neighbor_match_plot = TRUE)


               temp_error <- try({annotateMRN(annotation_result = "ms2_match_annotation_result.csv",
                                              ms2_data = 'ms2',
                                              prefer_adduct = ifelse(polarity == "positive", "M+H", "M-H"),
                                              metdna_version = metdna_version,
                                              lib = lib,
                                              direction = direction,
                                              column = column,
                                              instrument = instrument,
                                              method_lc = method_lc,
                                              polarity = polarity,
                                              extension_step = extension_step,
                                              threads = thread,
                                              max_isotope = max_isotope,
                                              path = path,
                                              mz_tol = mz_tol,
                                              mz_ppm_thr = mz_ppm_thr,
                                              rt_tol1 = rt_tol1,
                                              rt_tol2 = rt_tol2,
                                              ccs_tol = ccs_tol,
                                              cor_tol = cor_tol,
                                              int_tol = int_tol,
                                              dp_tol = dp_tol,
                                              max_step = 3,
                                              score_cutoff = 0,
                                              remain = FALSE,
                                              remain_per = 0.5,
                                              scoring_approach = scoring_approach_recursive,
                                              matched_frag_tol = matched_frag_cutoff,
                                              whether_link_frag = whether_link_frag,
                                              seed_neighbor_match_plot = seed_neighbor_match_plot,
                                              candidate_num = candidate_num,
                                              use_pretrained_model = use_pretrained_model,
                                              # test_old_mrn = test_old_mrn,
                                              test_adduct_version = ifelse(metdna_version == 'version2', 'version2', 'version1'),
                                              test_evaluation = test_evaluation,
                                              whether_use_predRT = TRUE)},
                                 silent = TRUE)

               if (class(temp_error) == "try-error") {
                 cat(temp_error[[1]], file = file.path(path, "run.log.txt"), append = TRUE)
                 stop(temp_error[[1]])
               }
             } else {
               cat('Skip metabolic reaction network based annotation.\n', file = file.path(path, "run.log.txt"), append = TRUE)
               cat("Skip metabolic reaction network based annotation.\n")
               ## return blank results
               # read initial annotation results and transfer tags2 results (S4 MetDNA2::PeakInfo)
               initialAnnotationResult = read.csv(file.path(file.path(path, "01_result_initial_seed_annotation"),
                                                            "ms2_match_annotation_result.csv"), stringsAsFactors = FALSE)
               tags2_after_annotation = lapply(1:nrow(initialAnnotationResult), function(idx) {
                 data_PeakInfo = new(Class = "PeakInfo",
                                     name = initialAnnotationResult$name[idx] %>% as.character(),
                                     mz = initialAnnotationResult$mz[idx] %>% as.numeric() %>% round(x = ., digits = 4),
                                     rt = initialAnnotationResult$rt[idx] %>% as.numeric() %>% round(x = ., digits = 2),
                                     ccs = -1,
                                     ms2 = data.frame(),
                                     annotation = list()
                 )
               })
               tags2_after_redundancy_remove = tags2_after_annotation
               # create dir
               path_output <- file.path(path, "02_result_MRN_annotation")
               dir.create(file.path(path_output, "00_intermediate_data"), showWarnings = FALSE, recursive = TRUE)
               # save blank mrn data
               save(tags2_after_annotation,
                    file = file.path(path_output, "00_intermediate_data","tags2_after_annotation"),
                    compress = "xz", version = 2)
               save(tags2_after_redundancy_remove,
                    file = file.path(path_output, "00_intermediate_data","tags2_after_redundancy_remove"), 
                    compress = "xz", version = 2)
               rm(list = c('initialAnnotationResult', 
                           'tags2_after_annotation', 
                           'tags2_after_redundancy_remove'));gc()
             }


             if(is_anno_mrn) cat("\n")
             cat("\n")
             cat('4. Annotaion Credential.\n', file = file.path(path, "run.log.txt"), append = TRUE)
             cat("4. Annotaion Credential.\n")
             cat("------------------------------------------------------------\n")
             if(is_credential){
               # merge initial seeds and recursive results for credential
               temp <- mergeRecursiveAnnotation(path = path,
                                                thread = thread,
                                                lib = lib,
                                                direction = direction,
                                                tolerance_rt_range = tolerance_rt_range,
                                                use_redun_rm_result = use_redun_rm_result,
                                                test_evaluation = test_evaluation,
                                                # for skip functions
                                                is_anno_mrn = is_anno_mrn,
                                                is_credential = is_credential)

               rm('temp');gc()

               cat('Change formats for annotation credential.\n')
               convertAnnotationTable2InitialId(peak_table_file = ms1_file,
                                                sample_info_file =  sample_info_file,
                                                path_dir = path,
                                                polarity = polarity,
                                                # test_assi_confidence_initial_seed = test_assi_confidence_initial_seed,
                                                tool = 'MetDNA2')

               # cat('Start annotation credential.\n')
               # authenticateAnnotationCredential(peak_table_file = "peak_table.csv",
               #                                  annotation_initial_file = "annotation_initial.csv",
               #                                  ms2_data_file = "ms2_data.RData",
               #                                  path_dir = path,
               #                                  polarity = polarity,
               #                                  thread = thread,
               #                                  isotope_int_ratio_check = TRUE, # para annotateIsotope
               #                                  isotope_int_ratio_cutoff = 500,
               #                                  is_ms2_check = TRUE, # compare ms2 similarity of neutral loss
               #                                  ms2_score_cutoff = -1, # -1 represent not filter
               #                                  is_plot_pseudo_MS1 = TRUE,
               #                                  # formula prediction
               #                                  dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
               #                                  ppm = 25, # ms1 tolerance
               #                                  acc = 25, # ms2 tolerance
               #                                  elements = "CHNOPS",
               #                                  num_formula_candidate = 3, # return top 3 candidate
               #                                  platform = 'linux')

               temp_error <- try({authenticateAnnotationCredential(peak_table_file = "peak_table.csv",
                                                                   annotation_initial_file = "annotation_initial.csv",
                                                                   ms2_data_file = "ms2_data.RData",
                                                                   path_dir = path,
                                                                   polarity = polarity,
                                                                   thread = thread,
                                                                   isotope_int_ratio_check = TRUE, # para annotateIsotope
                                                                   isotope_int_ratio_cutoff = 500,
                                                                   is_ms2_check = TRUE, # compare ms2 similarity of neutral loss
                                                                   ms2_score_cutoff = -1, # -1 represent not filter
                                                                   is_plot_pseudo_MS1 = is_plot_pseudo_MS1,
                                                                   # formula prediction
                                                                   dir_GenForm = dir_GenForm,
                                                                   is_pred_formula_all = is_pred_formula_all,
                                                                   ppm = formula_mz_tol_ms1, # ms1 tolerance
                                                                   acc = formula_mz_tol_ms2, # ms2 tolerance
                                                                   elements = formula_elements,
                                                                   num_formula_candidate = formula_candidate_num, # return top 3 candidate
                                                                   platform = platform)},
                                 silent = TRUE)

               if (class(temp_error) == "try-error") {
                 cat(temp_error[[1]], file = file.path(path, "run.log.txt"), append = TRUE)
                 stop(temp_error[[1]])
               }
             } else {
               cat('Skip annotaion credential.\n', file = file.path(path, "run.log.txt"), append = TRUE)
               cat("Skip annotaion credential.\n")
               # merge initial seeds and recursive results for credential (and skip credential)
               temp <- mergeRecursiveAnnotation(path = path,
                                                thread = thread,
                                                lib = lib,
                                                direction = direction,
                                                tolerance_rt_range = tolerance_rt_range,
                                                use_redun_rm_result = use_redun_rm_result,
                                                test_evaluation = test_evaluation,
                                                is_anno_mrn = is_anno_mrn,
                                                is_credential = is_credential)
               rm('temp');gc()
               # save annot_all data
               load(file.path(path, "03_annotation_credential/00_intermediate_data/annot_all_before_credential.RData"))
               path_output <- file.path(path, '00_annotation_table')
               dir.create(file.path(path_output, '00_intermediate_data'), showWarnings = FALSE, recursive = TRUE)
               save(annot_all, file = file.path(path_output, '00_intermediate_data', 'annot_all'), version = 2)
               rm('annot_all');gc()
             }


             cat("\n")
             cat('5. Merge and export result tables.\n', file = file.path(path, "run.log.txt"), append = TRUE)
             cat("5. Merge and export result tables.\n")
             cat("------------------------------------------------------------\n")
             if(TRUE){
               temp_error <- try({
                 generateMetDNA2AnnotationResult(path = path,
                                                 thread = thread,
                                                 is_cred_pg_filter = is_cred_pg_filter,
                                                 is_cred_formula_filter = is_cred_formula_filter,
                                                 is_pred_formula_all = is_pred_formula_all,
                                                 tolerance_rt_range = tolerance_rt_range,
                                                 direction = direction,
                                                 instrument = instrument,
                                                 mz_tol = mz_tol,
                                                 rt_tol = rt_tol2,
                                                 ccs_tol = ccs_tol,
                                                 candidate_num = candidate_num,
                                                 remove_conflict_seed_from_final_table = remove_conflict_seed_from_final_table,
                                                 test_evaluation = test_evaluation,
                                                 # for skip functions
                                                 is_anno_mrn = is_anno_mrn,
                                                 is_credential = is_credential)},
                 silent = TRUE)


               # temp_error <- try({generateMetDNA2AnnotationResult(path = path,
               #                                                    thread = thread,
               #                                                    is_cred_pg_filter = is_cred_pg_filter,
               #                                                    is_cred_formula_filter = is_cred_formula_filter,
               #                                                    is_pred_formula_all = is_pred_formula_all,
               #                                                    direction = direction,
               #                                                    tolerance_rt_range = tolerance_rt_range,
               #                                                    test_rm_extra_anno_from_ini_seed = test_rm_extra_anno_from_ini_seed)},
               #                   silent = TRUE)

               if (class(temp_error) == "try-error") {
                 cat(temp_error[[1]], file = file.path(path, "run.log.txt"), append = TRUE)
                 stop(temp_error[[1]])
               }

               # temp_error <- try({generateMetDNA2AnnotationResultWithCredential(path = path,
               #                                                                  type_order = type_order)},
               #                   silent = TRUE)
               #
               # if (class(temp_error) == "try-error") {
               #   cat(temp_error[[1]], file = file.path(path, "run.log.txt"), append = TRUE)
               #   stop(temp_error[[1]])
               # }

             } else {
               cat('Skip merge and export annotation credential.\n', file = file.path(path, "run.log.txt"), append = TRUE)
               cat('Skip merge and export annotation credential.\n')
             }


             cat("\n")
             cat('6. Biology interpretation.\n', file = file.path(path, "run.log.txt"), append = TRUE)
             cat("6. Biology interpretation.\n")
             cat("------------------------------------------------------------\n")
             if (is_bio_interpret){
               # browser()

               temp_error <- try({
                 if (polarity == 'positive') {
                   interpretBiology(sample_file_pos = ms1_file,
                                    sample_info_file_pos = sample_info_file,
                                    table_identification_pair_file_pos = 'table3_identification_pair.csv',
                                    polarity = polarity,
                                    path = path,
                                    metdna_version = metdna_version,
                                    comp_group = comp_group,
                                    uni_test = uni_test,
                                    correct_p = correct_p,
                                    p_cutoff = p_cutoff,
                                    fc_cutoff = fc_cutoff,
                                    species = species,
                                    quanti_pathway_method = quanti_pathway_method,
                                    organze_basis = 'kegg_id',
                                    extension_step = extension_step)
                 } else {
                   interpretBiology(sample_file_neg = ms1_file,
                                    sample_info_file_neg = sample_info_file,
                                    table_identification_pair_file_neg = 'table3_identification_pair.csv',
                                    polarity = polarity,
                                    path = path,
                                    metdna_version = metdna_version,
                                    comp_group = comp_group,
                                    uni_test = uni_test,
                                    correct_p = correct_p,
                                    p_cutoff = p_cutoff,
                                    fc_cutoff = fc_cutoff,
                                    species = species,
                                    quanti_pathway_method = quanti_pathway_method,
                                    organze_basis = 'kegg_id',
                                    extension_step = extension_step)
                 }


               }, silent = TRUE)

               if (class(temp_error) == "try-error") {
                 cat(temp_error[[1]], file = file.path(path, "run.log.txt"), append = TRUE)
                 # is_bio_interpret <- FALSE
                 # stop(temp_error[[1]])
               }

             } else {
               cat('Skip bioloy interpretation.\n', file = file.path(path, "run.log.txt"), append = TRUE)
               cat("Skip bioloy interpretation.\n")
             }


             cat("\n")
             cat('7. MetDNA2 analysis report generation.\n', file = file.path(path, "run.log.txt"), append = TRUE)
             cat("7. MetDNA2 analysis report generation.\n")
             cat("------------------------------------------------------------\n")
             if (is_exported_report){
               temp_error <- try({
                 generateMetDNA2Report(sample_info =  sample_info_file,
                                       path = path,
                                       polarity = polarity,
                                       export_type = 'html',
                                       extension_step = extension_step,
                                       is_rt_calibration = is_rt_calibration,
                                       # for skip functions
                                       is_anno_mrn = is_anno_mrn,
                                       is_credential = is_credential,
                                       is_bio_interpret = FALSE)
               }, silent = TRUE)

               if (class(temp_error) == "try-error") {
                 cat(temp_error[[1]], file = file.path(path, "run.log.txt"), append = TRUE)
                 # stop(temp_error[[1]])
               }

             } else {
               cat('Skip MetDNA2 analysis report generation.\n', file = file.path(path, "run.log.txt"), append = TRUE)
               cat("Skip MetDNA2 analysis report generation.\n")
             }

             # remove intermediate data
             if (is_rm_intermediate_data) {
               temp_file_list <- c('00_annotation_table/00_intermediate_data',
                                   '01_result_initial_seed_annotation',
                                   '02_result_MRN_annotation/00_intermediate_data',
                                   '02_result_MRN_annotation/01_surrogate_ms2_spec_plot',
                                   '03_annotation_credential',
                                   '04_biology_intepretation/00_intermediate_data',
                                   '04_biology_intepretation/data.csv',
                                   '04_biology_intepretation/sample.info.csv',
                                   '04_biology_intepretation/stat_pathway_quantitative_result.csv',
                                   '04_biology_intepretation/table3_identification_pair.csv',
                                   '05_analysis_report/01_template',
                                   'test')

               removeFiles(files = temp_file_list, path = path)

               if (metdna_version == 'version1') {
                 removeFiles(files = '00_annotation_table', path = path)
               }
             }



           })



################################################################################
# runMetDNA2 -------------------------------------------------------------------

#' @title runMetDNA2
#' @description a wrapper to run MetDNA2
#' @author Zhiwei Zhou
#' @param ms1_file_pos
#' @param ms2_file_neg
#' @param ms2_file_pos
#' @param ms2_file_neg
#' @param sample_info_file_pos
#' @param sample_info_file_neg
#' @param path_pos directory of positive data
#' @param path_neg directory of negative data
#' @param metdna_version MetDNA parameter 'version1' and 'version2'
#' @param ms2_type 'mgf', 'mzXML', 'msp', 'cef' are supported. Default: 'mgf'
#' @param thread The number of threads. Default: 4
#' @param is_check_data Whether check data quality. Default: TRUE
#' @param is_anno_initial_seed Whether annotate initial seeds with spectral DB. Default: TRUE
#' @param is_anno_mrn Whether perform recursive annotation. Default: TRUE
#' @param is_credential Whether perform annotation credential. Default: TRUE
#' @param is_cred_pg_filer Whether filter error annotations by peak group annotaion in final table. Default: FALSE
#' @param is_cred_formula_filter Whether filter error annotations by formula prediction in final table. Default: FALSE
#' @param lib Default: 'zhumetlib_qtof'
#' @param polarity 'positive', 'negative'. Default: 'positive'
#' @param instrument The instrument you used to acquire data. "SciexTripleTOF", "AgilentQTOF", "BrukerQTOF", "ThermoOrbitrap",
#' ThermoExploris', "AgilentDTIMMS", "BrukerTIMS", 'WatersQTOF','WatersTWIMMS'. Default: "SciexTripleTOF"
#' @param column "hilic", "rp". Default: 'hilic'
#' @param ce "30", "10", "20", "35,15", "40", "50". Default: "30"
#' @param method_lc 'Amide12min', 'Amide23min', 'Other', 'MetlinRP', 'RP12min'. Default: "Amide12min"
#' @param is_rt_calibration Only valid with 'Amide12min' and 'Amide23min' methods. Default: TRUE.
#' @param mz_tol The tolerance of ms1 match. Default: 25 ppm
#' @param tolerance_rt_range The maxmium rt tolerance. Default: 30 s
#' @param dp_cutoff MS2 score cutoff. Default: 0.8
#' @param direction Reserve which direct dot product. 'reverse', 'forward'. Default: 'reverse'
#' @param is_plot_ms2 Whether plot initial-seed ms2 match plot. Default: TRUE
#' @param extension_step The used extended MRN with n step reaction. Steps: '0', '1', '2', '3', '4', '5', '6', '7', '8'
#' @param max_isotope Maximum considerable isotope. Default: 4
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
#' @param candidate_num reserve top N candidate number in recursive annotation. Default: 3
#' @param dir_GenForm
#' @param is_pred_formula_all whether predict formula for all peak group (TRUE) or concised peak groups (FALSE). Default: FALSE
#' @param formula_mz_tol_ms1 Default: 25 ppm
#' @param formula_mz_tol_ms2 Default: 25 ppm
#' @param formula_elements Default: "CHNOPS"
#' @param formula_candidate_num Reserve top 3 formula. Default: 3
#' @param platform 'windows', 'linux'. Default: 'linux'
#' @export

setGeneric(name = 'runMetDNA2',
           function(
             ms1_file_pos = "data.csv",
             ms1_file_neg = "data.csv",
             ms2_file_pos = NULL,
             ms2_file_neg = NULL,
             sample_info_file_pos = "sample.info.csv",
             sample_info_file_neg = "sample.info.csv",
             path_pos = ".",
             path_neg = ".",
             metdna_version = c('version2', 'version1'),
             ms2_type = NULL,
             thread = 3,
             is_check_data = TRUE,
             is_anno_initial_seed = TRUE,
             is_anno_mrn = TRUE,
             is_credential = TRUE,
             is_bio_interpret = TRUE,
             is_exported_report = TRUE,
             is_cred_pg_filter = FALSE,
             is_cred_formula_filter = FALSE,
             is_rm_intermediate_data = TRUE,
             is_webserver = TRUE,

             # parameters of loadDB
             lib = c('zhumetlib_qtof', 'zhumetlib_orbitrap', 'fiehnHilicLib'),
             polarity = c("positive", "negative", 'both'),
             instrument = c("SciexTripleTOF", "AgilentQTOF", "BrukerQTOF", "ThermoOrbitrap",
                            'ThermoExploris', "AgilentDTIMMS", "BrukerTIMS", 'WatersQTOF','WatersTWIMMS'),
             column = c("hilic", "rp"),
             ce = c("30", "10", "20", "35,15", "40", "50",
                    'NCE10', 'NCE20', 'NCE30', 'NCE40', 'NCE50',
                    'SCE20_30_40%', "SNCE20_30_40%"),
             method_lc = c('Amide12min', 'Amide23min', 'Other', 'MetlinRP', 'RP12min'),
             is_rt_calibration = TRUE,

             # parameters of matchMs1WithSpecLib
             mz_tol = 15, # ppm
             mz_ppm_thr = 400,
             tolerance_rt_range = 25, # s
             tolerance_ccs_range = 4, # %
             int_ms2_min_abs = 50,
             dp_cutoff = 0.8,
             direction = c('reverse', 'forward'),
             mz_tol_combine_ms1_ms2 = 25,
             rt_tol_combine_ms1_ms2 = 10,
             is_plot_ms2 = TRUE,
             adduct_limit_specLib = 'limited', # 220115
             test_force_filtering_rt = NULL,

             ## para of RecursiveAnnotationMRN
             extension_step = c('0', '1', '2', '3', '4', '5', '6', '7', '8'),
             max_isotope = 4,
             rt_tol1 = 3, # s
             rt_tol2 = 30, # %
             ccs_tol = 4, # %
             cor_tol = 0,
             int_tol = 500,
             dp_tol = 0.5,
             max_step = 3,
             score_cutoff = 0,
             remain = FALSE,
             remain_per = 0.5,
             seed_neighbor_match_plot = TRUE,
             candidate_num = 10,
             scoring_approach_recursive = c('dp', 'bonanza', 'hybrid', 'gnps'),
             matched_frag_cutoff = 4,
             whether_link_frag = TRUE,
             use_pretrained_model = FALSE,
             # test_old_mrn = c('v0.5', 'v0.3'),
             # test_adduct_version = c('version2', 'version1'),

             # annotation credential
             use_redun_rm_result = TRUE,
             # dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
             dir_GenForm = '/GenForm',
             is_pred_formula_all = FALSE,
             formula_mz_tol_ms1 = 25, # ms1 tolerance
             formula_mz_tol_ms2 = 25, # ms2 tolerance
             formula_elements = "CHNOPS",
             formula_candidate_num = 3, # return top 3 candidate
             platform = 'linux',
             test_rm_extra_anno_from_ini_seed = FALSE, # whether assign confidence level2 for adduct and isotopes from initial seed
             is_plot_pseudo_MS1 = TRUE,

             # result table export
             remove_conflict_seed_from_final_table = TRUE,

             # biology interpretation
             comp_group = c("W03", "W30"),
             uni_test = c("t","wilcox","anova"),
             correct_p = FALSE,
             p_cutoff = 0.05,
             fc_cutoff = 1,
             species = c("hsa", "mmu", "rno", "bta", "gga",
                         "dre", "dme", "cel", "sce", "osa",
                         "ath", "smm", "pfa", "tbr", "eco",
                         "bsu", "ppu", "sau", "tma", "syf", "mln"),
             quanti_pathway_method = c('mean', 'sum', 'median')
           ){
             # path decision
             if (polarity == "both"){
               path <- retieveUpperPath(path_pos)

               if (!all(c('POS', 'NEG') %in% list.files(path))) {
                 stop('Positive data need to be in POS directory, Negeative data need to be in NEG directory\n')
               }

               path <- file.path(path, "BOTH")
               dir.create(path, showWarnings = FALSE, recursive = TRUE)
             } else {
               path <- ifelse(polarity == 'positive', path_pos, path_neg)
             }

             if (is_webserver) {
               # extension_step <- '0'
               is_plot_ms2 <- FALSE
               seed_neighbor_match_plot <- FALSE
               is_plot_pseudo_MS1 <- FALSE
             }
             
             # browser()
             # decision adduct
             # temp_adduct_pos <- lib_adduct_nl$positive %>% dplyr::filter(annotation == 'Yes') %>% dplyr::pull(adduct)
             # temp_adduct_neg <- lib_adduct_nl$negative %>% dplyr::filter(annotation == 'Yes') %>% dplyr::pull(adduct)
             
             temp_adduct_pos <- lib_adduct_nl$positive %>% dplyr::pull(adduct)
             temp_adduct_neg <- lib_adduct_nl$negative %>% dplyr::pull(adduct)
             temp_adduct <- c(temp_adduct_pos, temp_adduct_neg)
             
             if (adduct_limit_specLib %in% c('limited', 'all')) {
               switch(adduct_limit_specLib,
                      'limited' = {
                        list_adduct_limit_specLib <- list('pos' = c('[M+H]+', '[M]+', '[M-H2O+H]+'),
                                                         'neg' = c('[M-H]-', '[M-2H]-', '[M-H2O-H]-')) # 20220415 add -H2O adducts (by ZHS)
                      },
                      'all' = {
                        list_adduct_limit_specLib <- list('pos' = 'all',
                                                         'neg' = 'all')
                      })
             } else if (all(adduct_limit_specLib %in% temp_adduct)) {
               list_adduct_limit_specLib <- list('pos' = adduct_limit_specLib[which(adduct_limit_specLib %in% temp_adduct_pos)],
                                                 'neg' = adduct_limit_specLib[which(adduct_limit_specLib %in% temp_adduct_neg)])
             } else {
               stop('Please check the parameter [adduct_limit_specLib], it should be support ')
             }
             
             rm(list = c('temp_adduct', 'temp_adduct_pos', 'temp_adduct_neg'));gc()

             # ms2_type decision
             if (length(ms2_type) == 0) {
               # change the ms2_type to mgf if there are mgf files, and change ms2_type to msp if there are msp files.
               if (polarity %in% c('positive', 'both')) {
                 postfix <- list.files(path_pos) %>% getPostfix
                 postfix <- unique(tolower(postfix[!is.na(postfix)]))

                 if(any(postfix == "mgf")){
                   ms2_type <- "mgf"
                 }

                 if(any(postfix == "msp")){
                   ms2_type <- "msp"
                 }

                 if(any(postfix == "mzxml")){
                   ms2_type <- "mzXML"
                 }

                 if(any(postfix == "cef")){
                   ms2_type <- "cef"
                 }

               } else {
                 postfix <- list.files(path_neg) %>% getPostfix
                 postfix <- unique(postfix[!is.na(postfix)])

                 if(any(postfix == "mgf")){
                   ms2_type <- "mgf"
                 }

                 if(any(postfix == "msp")){
                   ms2_type <- "msp"
                 }

                 if(any(postfix == "mzxml")){
                   ms2_type <- "mzXML"
                 }

                 if(any(postfix == "cef")){
                   ms2_type <- "cef"
                 }
               }
             }

             # lib decision
             if (length(lib) > 1) {
               if (instrument %in% c("SciexTripleTOF", "AgilentQTOF", "BrukerQTOF",
                                     "AgilentDTIMMS", "BrukerTIMS", 'WatersQTOF','WatersTWIMMS')) {
                 lib <- 'zhumetlib_qtof'
               } else {
                 lib <- 'zhumetlib_orbitrap'
               }
             }
             
             # adjust mz tolerance and thresold according to orbitrap instrument
             if (instrument %in% c("ThermoOrbitrap", 'ThermoExploris')) {
               mz_tol <- 15 # ppm
               mz_ppm_thr <- 0
             }
             

             # parameter check and decision
             metdna_version <- match.arg(metdna_version)
             # ms2_type <- match.arg(ms2_type)
             instrument <- match.arg(instrument)
             polarity <- match.arg(polarity)
             column <- match.arg(column)
             ce <- match.arg(ce)
             lib <- match.arg(lib)
             method_lc <- match.arg(method_lc)
             direction <- match.arg(direction)
             scoring_approach_recursive <- match.arg(scoring_approach_recursive)
             uni_test <- match.arg(uni_test)
             species <- match.arg(species)
             quanti_pathway_method <- match.arg(quanti_pathway_method)
             # test_adduct_version <- match.arg(test_adduct_version)

             # browser()
             if (polarity == 'both') {
               list_para <- list(
                 'package_verison' = packageVersion('MetDNA2'),
                 'ms1_file' = paste(ms1_file_pos, ms1_file_neg, sep = ';'),
                 'ms2_file' = paste(paste(ms2_file_pos, collapse = ';'), paste(ms2_file_neg, collapse = ';'), sep = ';'),
                 'sample_info_file' = paste(sample_info_file_pos, sample_info_file_neg, sep = ';'),
                 'metdena_version' = metdna_version,
                 'path' = path,
                 'thread' = thread,
                 'is_check_data' = is_check_data,
                 'is_anno_initial_seed' = is_anno_initial_seed,
                 'is_anno_mrn' = is_anno_mrn,
                 'is_credential' = is_credential,
                 'is_bio_interpret' = is_bio_interpret,
                 'is_exported_report' = is_exported_report,
                 'is_cred_pg_filter' = is_cred_pg_filter,
                 'is_cred_formula_filter' = is_cred_formula_filter,
                 'is_rm_intermediate_data' = is_rm_intermediate_data,

                 'lib' = lib,
                 'polarity' = polarity,
                 'instrument' = instrument,
                 'column' = column,
                 'ce' = ce,
                 'method_lc' = method_lc,
                 'is_rt_calibration' = is_rt_calibration,

                 # parameters of matchMs1WithSpecLib
                 'mz_tol' = mz_tol,
                 'mz_ppm_thr' = mz_ppm_thr,
                 'tolerance_rt_range' = tolerance_rt_range,
                 'tolerance_ccs_range' = tolerance_ccs_range,
                 'dp_cutoff' = dp_cutoff,
                 'direction' = direction,
                 'is_plot_ms2' = is_plot_ms2,
                 'mz_tol_combine_ms1_ms2' = mz_tol_combine_ms1_ms2,
                 'rt_tol_combine_ms1_ms2' = rt_tol_combine_ms1_ms2,
                 'int_ms2_min_abs' = int_ms2_min_abs,
                 'adduct_limit_specLib' = paste(adduct_limit_specLib, collapse = ';'), # 220115

                 ## para of RecursiveAnnotationMRN
                 'extension_step' = extension_step,
                 'max_isotope' = max_isotope,
                 'rt_tol1' = rt_tol1,
                 'rt_tol2' = rt_tol2,
                 'cor_tol' = cor_tol,
                 'int_tol' = int_tol,
                 'dp_tol' = dp_tol,
                 'max_step' = max_step,
                 'score_cutoff' = score_cutoff,
                 'remain' = remain,
                 'remain_per' = remain_per,
                 'seed_neighbor_match_plot' = seed_neighbor_match_plot,
                 'candidate_num' = candidate_num,
                 'scoring_approach_recursive' = scoring_approach_recursive,
                 'matched_frag_cutoff' = matched_frag_cutoff,
                 'whether_link_frag' = whether_link_frag,
                 'use_pretrained_model' = use_pretrained_model,


                 # annotation credential
                 'use_redun_rm_result' = use_redun_rm_result,
                 'dir_GenForm' = dir_GenForm,
                 'is_pred_formula_all' = is_pred_formula_all,
                 'formula_mz_tol_ms1' = formula_mz_tol_ms1, # ms1 tolerance
                 'formula_mz_tol_ms2' = formula_mz_tol_ms2, # ms2 tolerance
                 'formula_elements' = formula_elements,
                 'formula_candidate_num' = formula_candidate_num, # return top 3 candidate
                 'platform' = platform,
                 'is_plot_pseudo_MS1' = is_plot_pseudo_MS1,

                 # result table export
                 'remove_conflict_seed_from_final_table' = remove_conflict_seed_from_final_table,

                 # biology interpretation
                 'comp_group' = paste(comp_group, collapse = ';'),
                 'uni_test' = uni_test,
                 'correct_p' = correct_p,
                 'p_cutoff' = p_cutoff,
                 'fc_cutoff' = fc_cutoff,
                 'species' = species,
                 'quanti_pathway_method' = quanti_pathway_method,

                 # test parameters
                 # 'test_adduct_version' = test_adduct_version,
                 'test_rm_extra_anno_from_ini_seed' = test_rm_extra_anno_from_ini_seed
               ) %>%
                 tibble::as_tibble() %>%
                 dplyr::mutate_all(as.character) %>%
                 tidyr::pivot_longer(cols = dplyr::everything(), names_to = 'para')

               dir.create(file.path(path),
                          showWarnings = FALSE,
                          recursive = TRUE)

               readr::write_tsv(list_para, path = file.path(path, 'para_list.txt'))
             }

             # create run.log file
             cat("MetDNA2 run log\n", file = file.path(path, "run.log.txt"))
             cat('MetDNA2 verison:', as.character(packageVersion('MetDNA2')), '\n', append = TRUE)
             cat(as.character(Sys.time()), file = file.path(path, "run.log.txt"), append = TRUE)
             cat("\n", file = file.path(path, "run.log.txt"), append = TRUE)

             # check data
             if (polarity == "positive" | polarity == "both") {
               cat('\n\n');cat("========================================================\n")
               cat("MetDNA2 analysis for positive mode data.\n")
               cat("========================================================\n",
                   file = file.path(path, "run.log.txt"), append = TRUE)
               cat("MetDNA2 analysis for positive mode data.\n",
                   file = file.path(path, "run.log.txt"), append = TRUE)

               # default set for MetDNA2 web server 20220215
               if (is_webserver) {
                 temp_adduct_limit_specLib <- c('[M+H]+', '[M]+')
               } else {
                 temp_adduct_limit_specLib <- list_adduct_limit_specLib$pos
               }

               MetDNA2(ms1_file = ms1_file_pos,
                       ms2_file = ms2_file_pos,
                       sample_info_file = sample_info_file_pos,
                       metdna_version = metdna_version,
                       ms2_type = ms2_type,
                       path = path_pos,
                       thread = thread,
                       is_check_data = is_check_data,
                       is_anno_initial_seed = is_anno_initial_seed,
                       is_anno_mrn = is_anno_mrn,
                       is_credential = is_credential,
                       is_bio_interpret =  ifelse(polarity == "both", FALSE, TRUE),
                       is_exported_report = ifelse(polarity == "both", FALSE, TRUE),
                       is_cred_pg_filter = is_cred_pg_filter,
                       is_cred_formula_filter = is_cred_formula_filter,
                       is_rm_intermediate_data = ifelse(polarity == "both", FALSE, is_rm_intermediate_data),
                       lib = lib,
                       polarity = 'positive',
                       instrument = instrument,
                       column = column,
                       ce = ce,
                       method_lc = method_lc,
                       is_rt_calibration = is_rt_calibration,
                       mz_tol = mz_tol,
                       mz_ppm_thr = mz_ppm_thr,
                       tolerance_rt_range = tolerance_rt_range,
                       tolerance_ccs_range = tolerance_ccs_range,
                       int_ms2_min_abs = int_ms2_min_abs,
                       dp_cutoff = dp_cutoff,
                       direction = direction,
                       mz_tol_combine_ms1_ms2 = mz_tol_combine_ms1_ms2,
                       rt_tol_combine_ms1_ms2 = rt_tol_combine_ms1_ms2,
                       adduct_limit_specLib = temp_adduct_limit_specLib,
                       test_force_filtering_rt = test_force_filtering_rt,
                       is_plot_ms2 = is_plot_ms2,
                       extension_step = extension_step,
                       max_isotope = max_isotope,
                       rt_tol1 = rt_tol1,
                       rt_tol2 = rt_tol2,
                       ccs_tol = ccs_tol,
                       cor_tol = cor_tol,
                       int_tol = int_tol,
                       dp_tol = dp_tol,
                       max_step = max_step,
                       score_cutoff = score_cutoff,
                       remain = remain,
                       remain_per = remain_per,
                       seed_neighbor_match_plot = seed_neighbor_match_plot,
                       candidate_num = candidate_num,
                       scoring_approach_recursive = scoring_approach_recursive,
                       matched_frag_cutoff = matched_frag_cutoff,
                       whether_link_frag = whether_link_frag,
                       use_pretrained_model = use_pretrained_model,
                       # test_adduct_version = test_adduct_version,
                       use_redun_rm_result = use_redun_rm_result,
                       dir_GenForm = dir_GenForm,
                       is_pred_formula_all = is_pred_formula_all,
                       formula_mz_tol_ms1 = formula_mz_tol_ms1,
                       formula_mz_tol_ms2 = formula_mz_tol_ms2,
                       formula_elements = formula_elements,
                       formula_candidate_num = formula_candidate_num,
                       platform = platform,
                       remove_conflict_seed_from_final_table = remove_conflict_seed_from_final_table,
                       is_plot_pseudo_MS1 = is_plot_pseudo_MS1,
                       test_rm_extra_anno_from_ini_seed = test_rm_extra_anno_from_ini_seed,
                       comp_group = comp_group,
                       uni_test = uni_test,
                       correct_p = correct_p,
                       p_cutoff = p_cutoff,
                       species = species,
                       quanti_pathway_method = quanti_pathway_method)

             }


             if (polarity == "negative" | polarity == "both") {
               cat('\n\n');cat("========================================================\n")
               cat("========================================================\n", file = file.path(path, "run.log.txt"), append = TRUE)
               cat("MetDNA2 analysis for negative mode data.\n", file = file.path(path, "run.log.txt"), append = TRUE)
               cat("MetDNA2 analysis for negative mode data.\n")

               # default set for MetDNA2 web server 20220215
               if (is_webserver) {
                 temp_adduct_limit_specLib <- c('[M-H]-', '[M-2H]-')
               } else {
                 temp_adduct_limit_specLib <- list_adduct_limit_specLib$neg
               }

               MetDNA2(ms1_file = ms1_file_neg,
                       ms2_file = ms2_file_neg,
                       sample_info_file = sample_info_file_neg,
                       metdna_version = metdna_version,
                       ms2_type = ms2_type,
                       path = path_neg,
                       thread = thread,
                       is_check_data = is_check_data,
                       is_anno_initial_seed = is_anno_initial_seed,
                       is_anno_mrn = is_anno_mrn,
                       is_credential = is_credential,
                       is_bio_interpret =  ifelse(polarity == "both", FALSE, TRUE),
                       is_exported_report = ifelse(polarity == "both", FALSE, TRUE),
                       is_cred_pg_filter = is_cred_pg_filter,
                       is_cred_formula_filter = is_cred_formula_filter,
                       is_rm_intermediate_data = ifelse(polarity == "both", FALSE, is_rm_intermediate_data),
                       lib = lib,
                       polarity = 'negative',
                       instrument = instrument,
                       column = column,
                       ce = ce,
                       method_lc = method_lc,
                       is_rt_calibration = is_rt_calibration,
                       mz_tol = mz_tol,
                       mz_ppm_thr = mz_ppm_thr,
                       tolerance_rt_range = tolerance_rt_range,
                       tolerance_ccs_range = tolerance_ccs_range,
                       int_ms2_min_abs = int_ms2_min_abs,
                       dp_cutoff = dp_cutoff,
                       direction = direction,
                       mz_tol_combine_ms1_ms2 = mz_tol_combine_ms1_ms2,
                       rt_tol_combine_ms1_ms2 = rt_tol_combine_ms1_ms2,
                       adduct_limit_specLib = temp_adduct_limit_specLib,
                       test_force_filtering_rt = test_force_filtering_rt,
                       is_plot_ms2 = is_plot_ms2,
                       extension_step = extension_step,
                       max_isotope = max_isotope,
                       rt_tol1 = rt_tol1,
                       rt_tol2 = rt_tol2,
                       ccs_tol = ccs_tol,
                       cor_tol = cor_tol,
                       int_tol = int_tol,
                       dp_tol = dp_tol,
                       max_step = max_step,
                       score_cutoff = score_cutoff,
                       remain = remain,
                       remain_per = remain_per,
                       seed_neighbor_match_plot = seed_neighbor_match_plot,
                       candidate_num = candidate_num,
                       scoring_approach_recursive = scoring_approach_recursive,
                       matched_frag_cutoff = matched_frag_cutoff,
                       whether_link_frag = whether_link_frag,
                       use_pretrained_model = use_pretrained_model,
                       # test_adduct_version = test_adduct_version,
                       use_redun_rm_result = use_redun_rm_result,
                       dir_GenForm = dir_GenForm,
                       is_pred_formula_all = is_pred_formula_all,
                       formula_mz_tol_ms1 = formula_mz_tol_ms1,
                       formula_mz_tol_ms2 = formula_mz_tol_ms2,
                       formula_elements = formula_elements,
                       formula_candidate_num = formula_candidate_num,
                       platform = platform,
                       remove_conflict_seed_from_final_table = remove_conflict_seed_from_final_table,
                       is_plot_pseudo_MS1 = is_plot_pseudo_MS1,
                       test_rm_extra_anno_from_ini_seed = test_rm_extra_anno_from_ini_seed,
                       comp_group = comp_group,
                       uni_test = uni_test,
                       correct_p = correct_p,
                       p_cutoff = p_cutoff,
                       species = species,
                       quanti_pathway_method = quanti_pathway_method)

             }


             if (polarity == "both") {
               cat("\n")
               cat("========================================================\n")
               cat("========================================================\n",
                   file = file.path(path, "run.log.txt"), append = TRUE)
               cat("Biology interpretation for positive and negative mode data.\n",
                   file = file.path(path, "run.log.txt"), append = TRUE)
               cat("Biology interpretation for positive and negative mode data.\n")

               temp_path <- retieveUpperPath(path)

               if (is_bio_interpret){
                 temp_error <- try({
                   interpretBiology(sample_file_pos = ms1_file_pos,
                                    sample_info_file_pos = sample_info_file_pos,
                                    table_identification_pair_file_pos = 'table3_identification_pair.csv',
                                    sample_file_neg = ms1_file_neg,
                                    sample_info_file_neg = sample_info_file_neg,
                                    table_identification_pair_file_neg = 'table3_identification_pair.csv',
                                    polarity = polarity,
                                    path = temp_path,
                                    metdna_version = metdna_version,
                                    comp_group = comp_group,
                                    uni_test = uni_test,
                                    correct_p = correct_p,
                                    p_cutoff = p_cutoff,
                                    fc_cutoff = fc_cutoff,
                                    species = species,
                                    quanti_pathway_method = quanti_pathway_method,
                                    organze_basis = 'kegg_id',
                                    extension_step = extension_step)
                 }, silent = TRUE)

                 if (class(temp_error) == "try-error") {
                   cat(temp_error[[1]], file = file.path(path, "run.log.txt"), append = TRUE)
                   # is_bio_interpret <- FALSE
                   # stop(temp_error[[1]])
                 }

               } else {
                 cat('Skip bioloy interpretation.\n', file = file.path(path, "run.log.txt"), append = TRUE)
                 cat("Skip bioloy interpretation.\n")
               }




               cat("\n")
               cat("========================================================\n")
               cat("========================================================\n",
                   file = file.path(path, "run.log.txt"), append = TRUE)
               cat("MetDNA2 analysis report generation.\n",
                   file = file.path(path, "run.log.txt"), append = TRUE)
               cat("MetDNA2 analysis report generation.\n")

               if (is_exported_report){
                 temp_error <- try({
                   temp_path <- retieveUpperPath(path)
                   generateMetDNA2Report(sample_info =  sample_info_file_pos,
                                         path = temp_path,
                                         polarity = 'both',
                                         export_type = 'html',
                                         extension_step = extension_step,
                                         is_rt_calibration = is_rt_calibration,
                                         is_bio_interpret = FALSE)
                 }, silent = TRUE)

                 if (class(temp_error) == "try-error") {
                   cat(temp_error[[1]], file = file.path(path, "run.log.txt"), append = TRUE)
                   # stop(temp_error[[1]])
                 }

               } else {
                 cat('Skip MetDNA2 analysis report generation.\n', file = file.path(path, "run.log.txt"), append = TRUE)
                 cat("Skip MetDNA2 analysis report generation.\n")
               }

               # remove intermediate data
               if (is_rm_intermediate_data) {
                 temp_file_list <- c('00_annotation_table/00_intermediate_data',
                                     '01_result_initial_seed_annotation',
                                     '02_result_MRN_annotation/00_intermediate_data',
                                     '02_result_MRN_annotation/01_surrogate_ms2_spec_plot',
                                     # '03_annotation_credential',
                                     '04_biology_intepretation/00_intermediate_data',
                                     '04_biology_intepretation/data.csv',
                                     '04_biology_intepretation/data_pos.csv',
                                     '04_biology_intepretation/data_neg.csv',
                                     '04_biology_intepretation/sample.info.csv',
                                     '04_biology_intepretation/sample.info_pos.csv',
                                     '04_biology_intepretation/sample.info_neg.csv',
                                     '04_biology_intepretation/stat_pathway_quantitative_result.csv',
                                     '04_biology_intepretation/table3_identification_pair.csv',
                                     '04_biology_intepretation/table3_identification_pair_pos.csv',
                                     '04_biology_intepretation/table3_identification_pair_neg.csv',
                                     '04_biology_intepretation/table3_identification_pair_pos_and_neg.csv',
                                     '05_analysis_report/01_template',
                                     'test')

                 temp_path <- retieveUpperPath(path)
                 removeFiles(files = temp_file_list, path = file.path(temp_path, 'POS'))
                 removeFiles(files = temp_file_list, path = file.path(temp_path, 'NEG'))
                 removeFiles(files = temp_file_list, path = file.path(temp_path, 'BOTH'))

                 if (metdna_version == 'version1') {
                   removeFiles(files = '00_annotation_table', path = file.path(temp_path, 'POS'))
                   removeFiles(files = '00_annotation_table', path = file.path(temp_path, 'NEG'))
                 }
               }


               cat("\n")
               cat("========================================================\n")
               cat("========================================================\n",
                   file = file.path(path, "run.log.txt"), append = TRUE)
               cat("MetDNA2 analysis done.\n",
                   file = file.path(path, "run.log.txt"), append = TRUE)
               cat("MetDNA2 analysis done.\n")

             }


           })



################################################################################
# startup massage --------------------------------------------------------------
.onAttach <- function(libname, pkgname){
  packageStartupMessage("
If you have any questions, please send email to zhouzw@sioc.ac.cn or jiangzhu@sioc.ac.cn.
Authors: Zhiwei Zhou, Mingdu Luo and Dr. Zhengjiang Zhu (jiangzhu@sioc.ac.cn).
Maintainer: Zhiwei Zhou, Mingdu Luo

Version 1.0.00 (20220215)
-------------
o modify default parameters to consistent with submitted manuscript
o open the parameter extended_step on web server
o open global peak correlation network function on web server

Version 1.0.01 (20220221)
-------------
o modify method_lc 'zhulabRP' to 'RP12min'
o update rt_ref for RP12min

Version 1.0.02 (20220224)
-------------
o update some typo in report

Version 1.1.00 (20220226)
-------------
o use MetLib functions to replace loadDB
o add skip credential and MRN sections

Version 1.1.01 (20220302)
-------------
o add adduct_limit_specLib decision for local runing, and default use 'limited' (POS:  c('[M+H]+', '[M]+'), NEG: c('[M-H]-', '[M-2H]-'))
o adjust parameter mz_tol and mz_ppm_thr according to instrument. (Default: QTOF: 15ppm & 400; QE: 15ppm & 0)
o fix some bugs
o modify logic of analysis report generation (Only all function enabled)

Version 1.1.11 (20220303)
-------------
o modify default tolerance_rt_range as 25s
o remove test_force_filtering_rt parameter

Version 1.1.12 (20220322)
-------------
o modify bugs: manually given adduct list 

Version 1.1.13 (20220415)
-------------
o modify bugs: addRecursiveResult2MetDNA2AnnoClass: add structure information into recursive annotation (id_kegg redirection)
o set adduct_limit_specLib default 'limited' (POS:  c('[M+H]+', '[M]+', '[M-H2O+H]+), NEG: c('[M-H]-', '[M-2H]-', '[M-H2O-H]-))

Version 1.1.15 (20220421)
-------------
o fix bugs: mz_ppm_thr parameter was not passed in matchMs2WithNeighbor
o change annotateAdductMRN as single thread

Version 1.1.16 (20220421)
-------------
o recode the function to call annotateAdductMRN with multiple threads 

Version 1.1.17 (20220423)
-------------
o modify bugs: addRecursiveResult2MetDNA2AnnoClass: id_kegg redirection

Version 1.2 / 1.1.18 (20220427)
-------------
o modify bugs: addRecursiveResult2MetDNA2AnnoClass: id_kegg redirection from lib_kegg (181214)
o cancel 1.1.13[1] and 1.1.17 modification

Version 1.2.10 (20220506)
-------------
o fix bugs: ignore electron mass in adduct m/z calculation during recursive annotation 
")
}

