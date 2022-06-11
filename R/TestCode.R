
################################################################################
# test code for feature clustering ---------------------------------------------

# 01 convert annotation table to PeakGroup -------------------------------------
# root_path <- './inst/tempdata/'
#
# peak_table <- readr::read_csv('./inst/extdata/peak_table_200STD_neg_200805.csv')
# peak_table_annotated <- readr::read_csv('./inst/extdata/peak_table_annotated_200STD_neg_200805.csv')
#
# progress <- MetDNA2::mapProgress(n = length(peak_table_annotated$name))
#
# list_peak_group <- map(seq_along(peak_table_annotated$name), function(i){
#
#   MetDNA2::mapProgressPrint(progress)
#   temp_name <- peak_table_annotated$name[i]
#   temp_rt <- peak_table_annotated$rt[i]
#   temp_adduct <- peak_table_annotated$adduct[i]
#
#   result <- generatePeakGroup(peak_table = peak_table,
#                               base_peak_name = temp_name,
#                               base_peak_adduct = temp_adduct,
#                               tol_rt = 3,
#                               cutoff_ssc = 0)
#   return(result)
#
# })
#
# names(list_peak_group) <- paste(peak_table_annotated$name,
#                                 peak_table_annotated$adduct,
#                                 sep = '_')
#
# save(list_peak_group,
#      file = file.path('./inst/tempdata', 'list_peak_group_200805.RData'),
#      version = 2)

# raw_msms <- ImmsTools::readMSP('./inst/extdata/spectra_200STD_neg_200805.msp', mode = 'all')
# temp <- sapply(raw_msms, function(x){
#   x[[1]]$NAME
# })
# names(raw_msms) <- temp
#
# dir.create('./inst/tempdata', showWarnings = FALSE, recursive = TRUE)
# save(raw_msms,
#      file = file.path('./inst/tempdata', 'raw_msms_200805.RData'),
#      version = 2)

# 02 annotate peak groups ------------------------------------------------------
# root_path <- './inst/tempdata/'
# load(file.path(root_path, 'list_peak_group_200805.RData'))
# load(file.path(root_path, 'raw_msms_200805.RData'))
# progress <- MetDNA2::mapProgress(length(list_peak_group))
# list_peak_group_annotation <- map(seq_along(list_peak_group), function(i){
#   MetDNA2::mapProgressPrint(progress)
#   cat(i, ' ')
#   result <- annotatePeakGroup(peak_group = list_peak_group[[i]],
#                               ms2_data = raw_msms,
#                               polarity = 'negative',
#                               tol_mz = 10,
#                               is_ms2_check = TRUE,
#                               ms2_score_cutoff = -1)
#
#   return(result)
#
# })
#
# names(list_peak_group_annotation) <- names(list_peak_group)
#
# save(list_peak_group_annotation,
#      file = file.path(root_path, 'list_peak_group_annotation_200805.RData'),
#      version = 2)
#

# 03 peak group concise --------------------------------------------------------
# # remove conflict peaks
# root_path <- './inst/tempdata/'
# load('./inst/tempdata/list_peak_group_annotation_200805.RData')
# list_peak_group_annotation_concised <- concisePeak2PeakGroup(list_peak_group_annotation = list_peak_group_annotation)
#
# save(list_peak_group_annotation_concised,
#      file = file.path(root_path, 'list_peak_group_annotation_concised_p2pg_200805.RData'),
#      version = 2)

# # remove overlap peakgroups
# load('./inst/tempdata/list_peak_group_annotation_concised_p2pg_200805.RData')
# list_peak_group_annotation_concised <- concisePeakGroup2PeakGroup(list_peak_group_annotation = list_peak_group_annotation_concised)
# save(list_peak_group_annotation_concised,
#      file = file.path(root_path, 'list_peak_group_annotation_concised_pg2pg_200805.RData'),
#      version = 2)


# # convert peak groups to semi-targeted annotation table ----------------------


################################################################################
# path_dir <- "I:/00_projects/03_MetDNA2/00_data/20200813_metdna_validation_set_data/20200815_liver_200std/pos"
# authenticateAnnotationCredential(peak_table_file = 'peak_table.csv',
#                                  annotation_initial_file = 'annotation_initial.csv',
#                                  ms2_data_file = 'ms2_data.RData',
#                                  path_dir = path_dir,
#                                  polarity = 'positive',
#                                  tol_mz = 10,
#                                  isotope_int_ratio_check = TRUE,
#                                  isotope_int_ratio_cutoff = 500,
#                                  is_ms2_check = TRUE,
#                                  ms2_score_cutoff = -1,
#                                  is_plot_pseudo_MS1 = TRUE,
#                                  dir_GenForm = "H:/software/GeneForm",
#                                  ppm = 20, acc = 20, elements = "CHNOPS",
#                                  num_formula_candidate = 3)

# peak_table_file = 'peak_table.csv'
# annotation_initial_file = 'annotation_initial.csv'
# ms2_data_file = 'ms2_data.RData'
# path_dir = path_dir
# polarity = 'positive'
# tol_mz = 10
# isotope_int_ratio_check = TRUE
# isotope_int_ratio_cutoff = 500
# is_ms2_check = TRUE
# ms2_score_cutoff = -1
# is_plot_pseudo_MS1 = TRUE
# dir_GenForm = "H:/software/GeneForm"
# ppm = 20
# acc = 20
# elements = "CHNOPS"
# num_formula_candidate = 3


################################################################################
# setwd('/home/zhouzw/Data_processing/20201019_annotation_credential_ubuntu/')


# peak_table_file = 'idresults.csv',
# ms2_file = 'spectra.msp',
# sample_info_file = 'sample.info.csv',
# # sample_name,
# path_dir = '.',
# lib = c('zhulib', 'zhumetlib', 'kegg'),
# match_method = c('forward', 'reverse'),
# polarity = c('positive', 'negative'),
# tool = c('MetAnalyzer1', 'MetAnalyzer2', 'MetDNA', 'MetDNA2'),

# convertAnnotationTable2InitialId(peak_table_file = 'idresults.csv',
#                                  ms2_file = 'spectra.msp',
#                                  sample_info_file =  'sample.info.csv',
#                                  path_dir = '/home/zhouzw/Data_processing/20201019_annotation_credential_ubuntu/',
#                                  lib = 'zhumetlib',
#                                  polarity = 'positive',
#                                  tool = 'MetAnalyzer1')

# authenticateAnnotationCredential(peak_table_file = "peak_table.csv",
#                                  annotation_initial_file = "annotation_initial.csv",
#                                  ms2_data_file = "ms2_data.RData",
#                                  path_dir = '/home/zhouzw/Data_processing/20201019_annotation_credential_ubuntu/',
#                                  polarity = 'positive',
#                                  thread = 4,
#                                  platform = 'linux')

# load('./annot_credential/00_intermediate_data/list_peak_group.RData')
# load('./ms2_data.RData')
#
# load('/home/zhouzw/Data_processing/20201019_fly_aging_pos/annot_credential/00_intermediate_data/list_peak_group_annotation_concised.RData')
# load('/home/zhouzw/Data_processing/20201019_fly_aging_pos/annot_credential/00_intermediate_data/list_peak_group.RData')
# load('/home/zhouzw/Data_processing/20201019_fly_aging_pos/ms2_data.RData')


# load('/home/zhouzw/Data_processing/20201020_fly_aging_pos/MetDNA2_pos/MRN_annotation_result/intermediate_data/tags2_after_annotation')
# load('/home/zhouzw/Data_processing/20201020_fly_aging_pos/MetDNA2_pos/MRN_annotation_result/intermediate_data/')


# MetDNA2(ms1_file = 'data.csv',
#         sample_info_file = "sample.info.csv",
#         metdna_version = 'version2',
#         ms2_type = 'mgf',
#         path = '/home/zhouzw/Data_processing/20201020_fly_aging_pos/MetDNA2_pos_step2/',
#         thread = 3,
#         polarity = 'positive',
#         instrument = "SciexTripleTOF",
#         column = "hilic",
#         lib = 'zhuMetLib',
#         ce = "30",
#         method_lc = 'Amide23min',
#         extension_step = '2')


# MetDNA2(ms1_file = 'data.csv',
#         sample_info_file = "sample.info.csv",
#         metdna_version = 'version2',
#         ms2_type = 'mgf',
#         path = '/home/zhouzw/Data_processing/20201022_fly_aging_neg/20201022_fly_neg/MetDNA2_step0/',
#         thread = 3,
#         polarity = 'negative',
#         instrument = "SciexTripleTOF",
#         column = "hilic",
#         lib = 'zhuMetLib',
#         ce = "30",
#         method_lc = 'Amide23min',
#         extension_step = '0')


# MetDNA2(ms1_file = 'data.csv',
#         sample_info_file = "sample.info.csv",
#         metdna_version = 'version2',
#         ms2_type = 'mgf',
#         path = '/home/zhouzw/Data_processing/20201022_fly_aging_neg/20201022_fly_neg/MetDNA2_step2/',
#         thread = 3,
#         polarity = 'negative',
#         instrument = "SciexTripleTOF",
#         column = "hilic",
#         lib = 'zhuMetLib',
#         ce = "30",
#         method_lc = 'Amide23min',
#         extension_step = '2')

# load('/home/zhouzw/Data_processing/20201022_fly_aging_neg/20201022_fly_neg/MetDNA2_step2/annot_credential/00_intermediate_data/list_peak_group.RData')
# load('/home/zhouzw/Data_processing/20201022_fly_aging_neg/20201022_fly_neg/MetDNA2_step2/annot_credential/00_intermediate_data/list_peak_group_annotation.RData')
# path_dir <- '/home/zhouzw/Data_processing/20201022_fly_aging_neg/20201022_fly_neg/MetDNA2_step2'
# type_order = c('seed', 'none_seed')
# which(unlist(list_concised_peak_group) == "M210T543_[M+F]-")
#
# sapply(list_concised_peak_group, function(x){
#   "M210T543_[M+F]-" %in% x
# }) %>% which()



# MetDNA2(ms1_file = 'data.csv',
#         sample_info_file = "sample.info.csv",
#         metdna_version = 'version1',
#         ms2_type = 'mgf',
#         path = '/home/zhouzw/Data_processing/20201022_fly_aging_neg/20201022_fly_neg/MetDNA1/',
#         thread = 3,
#         polarity = 'negative',
#         instrument = "SciexTripleTOF",
#         column = "hilic",
#         lib = 'zhuMetLib',
#         ce = "30",
#         method_lc = 'Amide23min',
#         extension_step = '2')


# MetDNA2(ms1_file = 'data.csv',
#         sample_info_file = "sample.info.csv",
#         metdna_version = 'version1',
#         ms2_type = 'mgf',
#         path = '/home/zhouzw/Data_processing/20201020_fly_aging_pos/MetDNA1_pos/',
#         thread = 3,
#         polarity = 'positive',
#         instrument = "SciexTripleTOF",
#         column = "hilic",
#         lib = 'zhuMetLib',
#         ce = "30",
#         method_lc = 'Amide23min',
#         extension_step = '0')


# path <- '/home/zhouzw/Data_processing/20201020_fly_aging_pos/MetDNA2_pos_step2/'
# convertAnnotationTable2InitialId(peak_table_file = "peak_table.csv",
#                                  sample_info_file =  "sample.info.csv",
#                                  path_dir = path,
#                                  polarity = 'positive',
#                                  tool = 'MetDNA2')

# authenticateAnnotationCredential(peak_table_file = "peak_table.csv",
#                                  annotation_initial_file = "annotation_initial.csv",
#                                  ms2_data_file = "ms2_data.RData",
#                                  path_dir = path,
#                                  polarity = 'positive',
#                                  thread = 3,
#                                  isotope_int_ratio_check = TRUE, # para annotateIsotope
#                                  isotope_int_ratio_cutoff = 500,
#                                  is_ms2_check = TRUE, # compare ms2 similarity of neutral loss
#                                  ms2_score_cutoff = -1, # -1 represent not filter
#                                  is_plot_pseudo_MS1 = TRUE,
#                                  # formula prediction
#                                  dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/bin/Debug',
#                                  ppm = 25, # ms1 tolerance
#                                  acc = 25, # ms2 tolerance
#                                  elements = "CHNOPS",
#                                  num_formula_candidate = 3, # return top 3 candidate
#                                  platform = 'linux')
#
#
# MetDNA2(ms1_file = 'data.csv',
#         sample_info_file = "sample.info.csv",
#         metdna_version = 'version2',
#         ms2_type = 'mgf',
#         path = '/home/zhouzw/Data_processing/20201102_metdna2_pos_development',
#         thread = 3,
#         polarity = 'positive',
#         instrument = "SciexTripleTOF",
#         column = "hilic",
#         lib = 'zhuMetLib',
#         ce = "30",
#         method_lc = 'Amide23min',
#         extension_step = '0',
#         dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
#         platform = 'linux',
#         is_cred_pg_filter = FALSE,
#         is_cred_formula_filter = FALSE)


# path <- '/home/zhouzw/Data_processing/20201113_predict_for_all_peakgroup'
# convertAnnotationTable2InitialId(peak_table_file = "peak_table.csv",
#                                  sample_info_file =  "sample.info.csv",
#                                  path_dir = path,
#                                  polarity = 'positive',
#                                  tool = 'MetDNA2')
#
# authenticateAnnotationCredential(peak_table_file = "peak_table.csv",
#                                  annotation_initial_file = "annotation_initial.csv",
#                                  ms2_data_file = "ms2_data.RData",
#                                  path_dir = path,
#                                  polarity = 'positive',
#                                  thread = 3,
#                                  isotope_int_ratio_check = TRUE, # para annotateIsotope
#                                  isotope_int_ratio_cutoff = 500,
#                                  is_ms2_check = TRUE, # compare ms2 similarity of neutral loss
#                                  ms2_score_cutoff = -1, # -1 represent not filter
#                                  is_plot_pseudo_MS1 = TRUE,
#                                  # formula prediction
#                                  dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
#                                  is_use_annotated_peak_group = FALSE,
#                                  ppm = 25, # ms1 tolerance
#                                  acc = 25, # ms2 tolerance
#                                  elements = "CHNOPS",
#                                  num_formula_candidate = 3, # return top 3 candidate
#                                  platform = 'linux')

# dir(file.path(path, 'annot_credential', "00_intermediate_data"))
#
# MetDNA2(ms1_file = 'data.csv',
#         sample_info_file = "sample.info.csv",
#         metdna_version = 'version2',
#         ms2_type = 'mgf',
#         path = '/home/zhouzw/Data_processing/20201102_metdna2_pos_development',
#         thread = 3,
#         polarity = 'positive',
#         instrument = "SciexTripleTOF",
#         column = "hilic",
#         lib = 'zhuMetLib',
#         ce = "30",
#         method_lc = 'Amide23min',
#         extension_step = '0',
#         dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
#         platform = 'linux',
#         is_check_data = FALSE,
#         is_anno_initial_seed = FALSE,
#         is_anno_mrn = FALSE,
#         is_credential = FALSE,
#         is_cred_pg_filter = FALSE,
#         is_cred_formula_filter = FALSE)


# # generate merged table
# MetDNA2(ms1_file = 'data.csv',
#         sample_info_file = "sample.info.csv",
#         metdna_version = 'version2',
#         ms2_type = 'mgf',
#         path = '/home/zhouzw/Data_processing/20201113_rerun_200std_pos/modified_package_result',
#         thread = 3,
#         polarity = 'positive',
#         instrument = "SciexTripleTOF",
#         column = "hilic",
#         lib = 'zhuMetLib',
#         ce = "30",
#         method_lc = 'Amide23min',
#         extension_step = '0',
#         dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
#         platform = 'linux',
#         is_check_data = FALSE,
#         is_anno_initial_seed = FALSE,
#         is_anno_mrn = FALSE,
#         is_credential = TRUE,
#         is_rt_calibration = FALSE,
#         is_cred_pg_filter = FALSE,
#         is_cred_formula_filter = FALSE)

# path_dir = '/home/zhouzw/Data_processing/20201113_rerun_200std_pos/modified_package_result'


# # generate merged table
# MetDNA2(ms1_file = 'data.csv',
#         sample_info_file = "sample.info.csv",
#         metdna_version = 'version2',
#         ms2_type = 'mgf',
#         path = '/home/zhouzw/Data_processing/20201113_rerun_200std_pos/modified_package_result_all_formula',
#         thread = 3,
#         polarity = 'positive',
#         instrument = "SciexTripleTOF",
#         column = "hilic",
#         lib = 'zhuMetLib',
#         ce = "30",
#         method_lc = 'Amide23min',
#         extension_step = '0',
#         dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
#         platform = 'linux',
#         is_check_data = FALSE,
#         is_anno_initial_seed = FALSE,
#         is_anno_mrn = FALSE,
#         is_credential = TRUE,
#         is_pred_formula_all = TRUE,
#         is_rt_calibration = FALSE,
#         is_cred_pg_filter = FALSE,
#         is_cred_formula_filter = FALSE)


# # generate merged table
# MetDNA2(ms1_file = 'data.csv',
#         sample_info_file = "sample.info.csv",
#         metdna_version = 'version2',
#         ms2_type = 'mgf',
#         path = '/home/zhouzw/Data_processing/20201113_rerun_200std_pos/modified_package_result_all_formula',
#         thread = 3,
#         polarity = 'positive',
#         instrument = "SciexTripleTOF",
#         column = "hilic",
#         lib = 'zhuMetLib',
#         ce = "30",
#         method_lc = 'Amide23min',
#         extension_step = '0',
#         dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
#         platform = 'linux',
#         is_check_data = FALSE,
#         is_anno_initial_seed = FALSE,
#         is_anno_mrn = FALSE,
#         is_credential = FALSE,
#         is_pred_formula_all = TRUE,
#         is_rt_calibration = FALSE,
#         is_cred_pg_filter = FALSE,
#         is_cred_formula_filter = FALSE)


# MetDNA2(ms1_file = 'data.csv',
#         sample_info_file = "sample.info.csv",
#         metdna_version = 'version2',
#         ms2_type = 'msp',
#         path = '/home/zhouzw/Data_processing/20201130_debug_for_XZD_problem/POS/',
#         thread = 3,
#         polarity = 'positive',
#         instrument = "SciexTripleTOF",
#         column = "hilic",
#         lib = 'zhuMetLib',
#         ce = "30",
#         method_lc = 'Amide12min',
#         extension_step = '0',
#         direction = 'reverse',
#         dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
#         platform = 'linux',
#         is_check_data = FALSE,
#         is_anno_initial_seed = TRUE,
#         is_anno_mrn = TRUE,
#         is_credential = TRUE,
#         is_rt_calibration = FALSE,
#         is_cred_pg_filter = FALSE,
#         is_cred_formula_filter = FALSE)



# MetDNA2(ms1_file = 'data.csv',
#         sample_info_file = "sample.info.csv",
#         metdna_version = 'version2',
#         ms2_type = 'msp',
#         path = '/home/zhouzw/Data_processing/20201102_metdna2_pos_development/',
#         thread = 3,
#         polarity = 'positive',
#         instrument = "SciexTripleTOF",
#         column = "hilic",
#         lib = 'zhuMetLib',
#         ce = "30",
#         method_lc = 'Amide12min',
#         extension_step = '2',
#         direction = 'reverse',
#         dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
#         platform = 'linux',
#         is_check_data = FALSE,
#         is_anno_initial_seed = FALSE,
#         is_anno_mrn = FALSE,
#         is_credential = TRUE,
#         is_rt_calibration = FALSE,
#         is_cred_pg_filter = FALSE,
#         is_cred_formula_filter = FALSE)


# MetDNA2(ms1_file = 'data.csv',
#         sample_info_file = "sample.info.csv",
#         metdna_version = 'version2',
#         ms2_type = 'mgf',
#         path = '/home/zhouzw/Data_processing/20201126_fly_step0_for_known_improvement/MetDNA2_update/',
#         thread = 3,
#         polarity = 'positive',
#         instrument = "SciexTripleTOF",
#         column = "hilic",
#         lib = 'zhuMetLib',
#         ce = "30",
#         method_lc = 'Amide23min',
#         extension_step = '0',
#         dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
#         platform = 'linux',
#         is_check_data = FALSE,
#         is_anno_initial_seed = FALSE,
#         is_anno_mrn = FALSE,
#         is_credential = TRUE,
#         is_rt_calibration = TRUE,
#         is_cred_pg_filter = FALSE,
#         is_cred_formula_filter = FALSE)

# load('/home/zhouzw/Data_processing/20201126_fly_step0_for_known_improvement/MetDNA2_update/annot_credential/00_intermediate_data/list_peak_group_annotation_concised.RData')
# load('/home/zhouzw/Data_processing/20201126_fly_step0_for_known_improvement/MetDNA2_update/annot_credential/00_intermediate_data/list_peak_group_formula.RData')


# load(file.path(path_dir, 'MRN_annotation_result/intermediate_data/'))



# MetDNA2(ms1_file = 'data.csv',
#         sample_info_file = "sample.info.csv",
#         metdna_version = 'version2',
#         ms2_type = 'mgf',
#         path = '/home/zhouzw/Data_processing/20201201_debug_mouse_liver_200std/02_POS/modified_package_result/',
#         thread = 3,
#         polarity = 'positive',
#         instrument = "SciexTripleTOF",
#         column = "hilic",
#         lib = 'zhuMetLib',
#         ce = "30",
#         method_lc = 'Amide23min',
#         extension_step = '0',
#         dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
#         platform = 'linux',
#         is_check_data = FALSE,
#         is_anno_initial_seed = FALSE,
#         is_anno_mrn = FALSE,
#         is_credential = FALSE,
#         is_rt_calibration = TRUE,
#         is_cred_pg_filter = FALSE,
#         is_cred_formula_filter = FALSE)



# path <- '/home/zhouzw/Data_processing/20201214_metdna2_development/'
# direction <- 'reverse'
# thread <- 3
# tolerance_rt_range <- 30
# ms1_file <- "data.csv"
# sample_info_file <- "sample.info.csv"
# polarity <- 'positive'

# load(file.path(path_dir,
#                'annot_credential',
#                "00_intermediate_data",
#                'list_peak_group_annotation.RData'))
# type_order = c('level1', 'level2', 'level3')

# is_cred_pg_filter <- FALSE
# is_cred_formula_filter <- FALSE


################################################################################
# 20201216 test ----------------------------------------------------------------
#
# path <- '/home/zhouzw/Data_processing/20201214_metdna2_development/'
# load(file.path(path, "01_result_initial_seed_annotation", '00_intermediate_data', 'result_annotation'))
#
# test <- lapply(result_annotation, function(x){
#   x@annotation_result
# })
#
# test <- test %>% dplyr::bind_rows()

# load('/home/zhouzw/Data_processing/20201214_metdna2_development/01_result_initial_seed_annotation/00_intermediate_data/')

# type_order <- c('level1', 'level2', 'level3')



################################################################################
# 20201229 debug for LWB data --------------------------------------------------
# path <- '/home/zhouzw/Data_processing/20201223_debug_msp/pos/'
#
# load('/home/zhouzw/Data_processing/20201223_debug_msp/pos/03_annotation_credential/00_intermediate_data/list_peak_group_annotation_concised.RData')
# load('/home/zhouzw/Data_processing/20201223_debug_msp/pos/03_annotation_credential/00_intermediate_data/list_peak_group_formula.RData')
# load('/home/zhouzw/Data_processing/20201223_debug_msp/pos/ms2_data.RData')


################################################################################
# 20201229

# path <- '/home/zhouzw/Data_processing/20201223_debug/20201217_fly_pos_main_server'
# load('/home/zhouzw/Data_processing/20201223_debug/20201217_fly_pos_main_server/00_annotation_table/00_intermediate_data/id_merge')
# load('/home/zhouzw/Data_processing/20201223_debug/20201217_fly_pos_main_server/00_annotation_table/00_intermediate_data/annot_all')


################################################################################
# 20210106
#
# path <- '/home/zhouzw/Data_processing/20210105_obitrap480_demo_data'
# MetDNA2(ms1_file = 'data.csv',
#         sample_info_file = "sample.info.csv",
#         metdna_version = 'version2',
#         ms2_type = 'msp',
#         path = '/home/zhouzw/Data_processing/20210105_obitrap480_demo_data',
#         thread = 3,
#         polarity = 'positive',
#         instrument = "ThermoOrbitrap",
#         column = "hilic",
#         lib = 'fiehnHilicLib',
#         ce = "SNCE20_30_40%",
#         method_lc = 'Amide12min',
#         extension_step = '0',
#         dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
#         platform = 'linux',
#         is_check_data = FALSE,
#         is_anno_initial_seed = FALSE,
#         is_anno_mrn = FALSE,
#         is_credential = FALSE,
#         is_rt_calibration = FALSE,
#         is_cred_pg_filter = FALSE,
#         is_cred_formula_filter = FALSE)



################################################################################
# 2021011 ----------------------------------------------------------------------
#
# path <- '/home/zhouzw/Data_processing/20201223_debug/20201217_fly_pos_main_server/'
# load('/home/zhouzw/Data_processing/20201223_debug/20201217_fly_pos_main_server/00_annotation_table/00_intermediate_data/id_merge')
# load('/home/zhouzw/Data_processing/20201223_debug/20201217_fly_pos_main_server/00_annotation_table/00_intermediate_data/annot_all')
#
# test <- generateResultTable4AllAnnotation(annot_all = annot_all,
#                                           path = path,
#                                           lib = 'zhuMetLib',
#                                           direction = 'reverse',
#                                           is_cred_pg_filter = FALSE,
#                                           is_cred_formula_filter = FALSE)
#
# id_merge <- test


################################################################################
# 20210111 debug LMD neg -------------------------------------------------------
# # problem: NL is same as ISF, but assign FALSE in is_credential_pg
# path <- '/home/zhouzw/Data_processing/20201223_debug/20201217_fly_neg_main_server'
# load('/home/zhouzw/Data_processing/20201223_debug/20201217_fly_neg_main_server/00_annotation_table/00_intermediate_data/annot_all')
# load('/home/zhouzw/Data_processing/20201223_debug/20201217_fly_neg_main_server/00_annotation_table/00_intermediate_data/id_merge')

# # problem: when M+_[M+H]+和M-_[M-H]- pair existed, [M+H]+ may be falsely replaced as [M+1] of M+
# path <- '/home/zhouzw/Data_processing/20201223_debug/20201217_fly_pos_main_server/'
# load('/home/zhouzw/Data_processing/20201223_debug/20201217_fly_pos_main_server/00_annotation_table/00_intermediate_data/id_merge')
# load('/home/zhouzw/Data_processing/20201223_debug/20201217_fly_pos_main_server/00_annotation_table/00_intermediate_data/annot_all')
# load('/home/zhouzw/Data_processing/20201223_debug/20201217_fly_pos_main_server/03_annotation_credential/00_intermediate_data/list_peak_group_annotation_concised.RData')
# load('/home/zhouzw/Data_processing/20201223_debug/20201217_fly_pos_main_server/03_annotation_credential/00_intermediate_data/list_peak_group.RData')
# load('/home/zhouzw/Data_processing/20201223_debug/20201217_fly_pos_main_server/ms2_data.RData')

# path <- '/home/zhouzw/Data_processing/20201223_debug/05_20201218_NEG_main_server/'
# load('/home/zhouzw/Data_processing/20201223_debug/05_20201218_NEG_main_server/00_annotation_table/00_intermediate_data/annot_all')
# load('/home/zhouzw/Data_processing/20201223_debug/05_20201218_NEG_main_server/00_annotation_table/00_intermediate_data/id_merge')



# path <- '/home/zhouzw/Data_processing/20201223_debug/20201217_fly_neg_main_server'
# load('/home/zhouzw/Data_processing/20201223_debug/20201217_fly_neg_main_server/00_annotation_table/00_intermediate_data/annot_all')
# load('/home/zhouzw/Data_processing/20201223_debug/20201217_fly_neg_main_server/00_annotation_table/00_intermediate_data/id_merge')
#
# MetDNA2(ms1_file = 'data.csv',
#         sample_info_file = "sample.info.csv",
#         metdna_version = 'version2',
#         ms2_type = 'mgf',
#         path = '/home/zhouzw/Data_processing/20201223_debug/20201217_fly_pos_main_server/',
#         thread = 3,
#         polarity = 'positive',
#         instrument = "SciexTripleTOF",
#         column = "hilic",
#         lib = 'zhuMetLib',
#         ce = "30",
#         method_lc = 'Amide23min',
#         extension_step = '0',
#         dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
#         platform = 'linux',
#         is_check_data = FALSE,
#         is_anno_initial_seed = FALSE,
#         is_anno_mrn = FALSE,
#         is_credential = FALSE,
#         is_rt_calibration = TRUE,
#         is_cred_pg_filter = FALSE,
#         is_cred_formula_filter = FALSE)



################################################################################
# path <- '/home/zhouzw/Data_processing/20210118_nist_urine_analysis/04_20210116_pos_only_full_8_DDA_emrn_2/'
# load('/home/zhouzw/Data_processing/20210118_nist_urine_analysis/04_20210116_pos_only_full_8_DDA_emrn_2/02_result_MRN_annotation/00_intermediate_data/tags2_after_annotation')
# load('/home/zhouzw/Data_processing/20210118_nist_urine_analysis/04_20210116_pos_only_full_8_DDA_emrn_2/02_result_MRN_annotation/00_intermediate_data/tags2_after_redundancy_remove')


################################################################################
# load('/home/zhouzw/Data_processing/20210224_metdna2_update/01_result_initial_seed_annotation/00_intermediate_data/ms1_result')
# load('/home/zhouzw/Data_processing/20210224_metdna2_update/01_result_initial_seed_annotation/00_intermediate_data/ms2')
# load('/home/zhouzw/Data_processing/20210224_metdna2_update/01_result_initial_seed_annotation/00_intermediate_data/ms1_data')
# load('/home/zhouzw/Data_processing/20210224_metdna2_update/01_result_initial_seed_annotation/00_intermediate_data/lib_meta')
# load('/home/zhouzw/Data_processing/20210224_metdna2_update/01_result_initial_seed_annotation/00_intermediate_data/lib_spec')

# load('/home/zhouzw/Data_processing/20210224_metdna2_update/01_result_initial_seed_annotation/00_intermediate_data/ms2_result')

# save(exp_ms2, file = '/home/zhouzw/Data_processing/20210224_metdna2_development_test/exp_ms2_210225.RData', version = 2)
# save(lib_ms2, file = '/home/zhouzw/Data_processing/20210224_metdna2_development_test/lib_ms2_210225.RData', version = 2)

# load('/home/zhouzw/Data_processing/20210224_metdna2_development_test/exp_ms2_210225.RData')
# load('/home/zhouzw/Data_processing/20210224_metdna2_development_test/lib_ms2_210225.RData')
#
#
# annotateInitialSeed(ms1_file = "data.csv",
#                     ms2_file = NULL,
#                     metdna_version = 'version2',
#                     ms2_type = 'mgf',
#                     path = "/home/zhouzw/Data_processing/20210224_metdna2_update/",
#
#                     polarity = 'positive',
#                     instrument = "SciexTripleTOF",
#                     lib = 'zhuMetLib',
#                     column = 'hilic',
#                     ce = '30',
#                     method_lc = 'Amide23min',
#                     excluded_adduct = NULL,
#                     is_rt_calibration = TRUE,
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
#                     scoring_approach = 'bonanza',
#                     # direction = 'reverse',
#
#                     is_plot_ms2 <- TRUE
#                     )

# load('/home/zhouzw/Data_processing/20210224_metdna2_update/01_result_initial_seed_annotation/00_intermediate_data/ms1_data')
# load('/home/zhouzw/Data_processing/20210224_metdna2_update/01_result_initial_seed_annotation/00_intermediate_data/ms2_result')
# load('/home/zhouzw/Data_processing/20210224_metdna2_update/01_result_initial_seed_annotation/00_intermediate_data/result_annotation')
# path_output <- '/home/zhouzw/Data_processing/20210224_metdna2_update/01_result_initial_seed_annotation/'
# table_annotation <- readr::read_csv('/home/zhouzw/Data_processing/20210224_metdna2_update/01_result_initial_seed_annotation/ms2_match_annotation_result.csv')

# load('/home/zhouzw/Data_processing/20210224_metdna2_update/01_result_initial_seed_annotation/00_intermediate_data/ms1_result')
# load('/home/zhouzw/Data_processing/20210224_metdna2_update/01_result_initial_seed_annotation/00_intermediate_data/result_annotation')
# load('/home/zhouzw/Data_processing/20210224_metdna2_update/01_result_initial_seed_annotation/00_intermediate_data/lib_meta')
# instrument = "SciexTripleTOF"
# is_rt_score = TRUE
# is_ccs_score = FALSE
# is_ms2_score = TRUE
# dp_cutoff <- 0.8

# table_annotation <- convertSpecAnnotationClass2Table(ms1_data = ms1_data,
#                                                      result_annotation = result_annotation,
#                                                      lib_meta = lib_meta,
#                                                      instrument = instrument,
#                                                      is_rt_score = is_rt_score,
#                                                      is_ccs_score = is_ccs_score,
#                                                      is_msms_score = is_msms_score,
#                                                      dp_cutoff = dp_cutoff)
#
# test_idx <- which(unlist(lapply(tags2, function(x) {length(x@annotation)})) > 0)
# tags2[[38]]@annotation
# tags2[[28]]@annotation
# tags2[[1205]]
# tags2[[786]]
#
#
# test <- tags2[test_idx] %>% convertTags2ResultMatrix()
# tags2[[test_idx[206]]]

#
# test_idx <- which(unlist(lapply(new_tags, function(x) {length(x@annotation)})) > 0)
# test <- new_tags[test_idx] %>% convertTags2ResultMatrix()

# test_idx <- which(unlist(lapply(tags2, function(x) {length(x@annotation)})) > 0)
# tags2[[38]]@annotation

# test_idx2 <- which(unlist(lapply(tags2, function(x) {length(x@annotation)})) > 0)
# tags2[[220]]@annotation

# test_idx3 <- which(unlist(lapply(tags2, function(x) {length(x@annotation)})) > 0)
# tags2[[195]]@annotation

# test_idx <- which(unlist(lapply(new_tags, function(x) {length(x@annotation)})) > 0)

# 20210308
# load('/home/zhouzw/Data_processing/20210224_metdna2_update/test_dp_bonanza/03_annotation_credential/00_intermediate_data/list_peak_group_annotation_concised.RData')
# load('/home/zhouzw/Data_processing/20210224_metdna2_update/test_dp_bonanza/ms2_data.RData')



# # 20210309
# MetDNA2(ms1_file = 'data.csv',
#         sample_info_file = "sample.info.csv",
#         metdna_version = 'version2',
#         ms2_type = 'mgf',
#         path = '/home/zhouzw/Data_processing/20210309_test_metdna_v0.40/aging_fly_POS',
#         thread = 3,
#         polarity = 'positive',
#         instrument = "SciexTripleTOF",
#         column = "hilic",
#         lib = 'zhuMetLib',
#         ce = "30",
#         method_lc = 'Amide23min',
#         is_rt_calibration = FALSE,
#         extension_step = '0',
#         scoring_approach_recursive = 'bonanza',
#         is_check_data = TRUE,
#         is_anno_initial_seed = TRUE,
#         is_anno_mrn = TRUE,
#         is_credential = TRUE,
#         is_bio_interpret = TRUE,
#         comp_group = c("W03", "W30"),
#         uni_test = "t",
#         correct_p = FALSE,
#         p_cutoff = 0.05,
#         species = 'dme',
#         quanti_pathway_method = 'mean')

################################################################################
# # 20210312 export cpd_emrn for LMD to predict CCS

#
# dir.create('/home/zhouzw/Data_processing/20210312_pred_ccs_for_metdna2', showWarnings = FALSE, recursive = TRUE)
#
# temp_data <- cpd_emrn %>% select(id, smiles, inchikey)
#
# readr::write_csv(temp_data,
#                  file = '/home/zhouzw/Data_processing/20210312_pred_ccs_for_metdna2/cpd_emrn_ccs_pred_210312.csv')



################################################################################
# # 20210315 MetDNA2 test for IM-MS data ---------------------------------------

# # # initial seed
# annotateInitialSeed(ms1_file = "data.csv",
#                     ms2_file = NULL,
#                     metdna_version = 'version2',
#                     ms2_type = 'mgf',
#                     path = "/home/zhouzw/Data_processing/20210315_tims_data_processing_metdna/",
#
#                     polarity = 'positive',
#                     instrument = "IMMS",
#                     lib = 'zhuMetLib',
#                     column = 'hilic',
#                     ce = '30',
#                     method_lc = 'Amide23min',
#                     excluded_adduct = NULL,
#                     is_rt_calibration = TRUE,
#
#                     mz_tol = 25,
#                     pf_rt_range = 0,
#                     tolerance_rt_range = 30,
#                     pf_ccs_range = 0,
#                     tolerance_ccs_range = 2,
#                     is_filter = TRUE,
#                     is_rt_score = TRUE,
#                     is_ccs_score = TRUE,
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

# recursive annotation
#
# path = "/home/zhouzw/Data_processing/20210315_tims_data_processing_metdna/01_result_initial_seed_annotation/"
# test_idx <- which(unlist(lapply(tags2, function(x) {length(x@annotation)})) > 0)
# tags2[[449]]@annotation
# tags2[[502]]@annotation
#
# test_idx2 <- which(unlist(lapply(tags2, function(x) {length(x@annotation)})) > 0)
# tags2[[1205]]
# tags2[[786]]


# test_idx3 <- which(unlist(lapply(tags2, function(x) {length(x@annotation)})) > 0)
# tags2[[93]]
# tags2[[356]]

# test_idx4 <- which(unlist(lapply(tags2, function(x) {length(x@annotation)})) > 0)
# tags2[[207]]
# tags2[[356]]


# test_idx_round2 <- which(unlist(lapply(tags2, function(x) {length(x@annotation)})) > 0)
# tags2[[7952]]
# tags2[[356]]
#
# tags2[test_idx_round2] %>% convertTags2ResultMatrix()

#
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
# showTags2(tags2, slot = "annotation.type")
# test <- sapply(result, ncol) %>% sapply(., length)
# which(test > 0)

#
#
# test <- convertTags2ResultMatrix(tags2 = tags2_after_redundancy_remove)
#
# which(sapply(annot_all, function(x)nrow(x@initial_seed_annotation)) > 0)
#
# test <- lapply(annot_all, function(x){
#   x@initial_seed_annotation
# })
#
# test <- test %>% dplyr::bind_rows()


################################################################################
# 20210326 test MetDNA MRN v0.3 ------------------------------------------------
#
# MetDNA2(ms1_file = 'data.csv',
#         sample_info_file = "sample.info.csv",
#         metdna_version = 'version2',
#         ms2_type = 'mgf',
#         path = '/home/zhouzw/Data_processing/20210324_evaluate_metdna1_2_mrn_change_effect/',
#         thread = 3,
#         polarity = 'positive',
#         int_ms2_min_abs = 50,
#         instrument = "SciexTripleTOF",
#         column = "hilic",
#         lib = 'zhuMetLib',
#         ce = "30",
#         method_lc = 'Amide23min',
#         is_rt_calibration = TRUE,
#         extension_step = '0',
#         scoring_approach_recursive = 'dp',
#         dp_tol = 0.5,
#         matched_frag_cutoff = 1,
#         use_redun_rm_result = TRUE,
#         is_check_data = TRUE,
#         is_anno_initial_seed = FALSE,
#         is_anno_mrn = TRUE,
#         is_credential = TRUE,
#         is_bio_interpret = FALSE,
#         # dir_GenForm = '/mnt/data4/program/share/bin',
#         comp_group = c("A", "B"),
#         uni_test = "t",
#         correct_p = FALSE,
#         p_cutoff = 0.05,
#         species = 'hsa',
#         quanti_pathway_method = 'mean',
#         test_old_mrn = 'v0.3')

################################################################################
# 20210329 ---------------------------------------------------------------------
#
# MetDNA2(ms1_file = 'data.csv',
#         sample_info_file = "sample.info.csv",
#         metdna_version = 'version2',
#         ms2_type = 'mgf',
#         # path = '/mnt/data5/memberdata/ZHOUzhiwei/01_data/20210324_MRN_t# rim_evaluation/04_metdna2_mrn03',
#         path = '/home/zhouzw/Data_processing/20210324_evaluate_metdna1_2_mrn_change_effect/',
#         thread = 3,
#         polarity = 'positive',
#         int_ms2_min_abs = 30,
#         instrument = "SciexTripleTOF",
#         column = "hilic",
#         lib = 'zhuMetLib',
#         ce = "30",
#         method_lc = 'Amide23min',
#         is_rt_calibration = TRUE,
#         extension_step = '0',
#         scoring_approach_recursive = 'dp',
#         dp_tol = 0.5,
#         matched_frag_cutoff = 1,
#         use_redun_rm_result = TRUE,
#         is_check_data = TRUE,
#         is_anno_initial_seed = TRUE,
#         is_anno_mrn = FALSE,
#         is_credential = FALSE,
#         is_bio_interpret = FALSE,
#         # dir_GenForm = '/mnt/data4/program/share/bin',
#         comp_group = c("A", "B"),
#         uni_test = "t",
#         correct_p = FALSE,
#         p_cutoff = 0.05,
#         species = 'hsa',
#         quanti_pathway_method = 'mean')


################################################################################
# 20210407 debug ---------------------------------------------------------------

# MetDNA2(ms1_file = 'data.csv',
#         sample_info_file = "sample.info.csv",
#         metdna_version = 'version2',
#         ms2_type = 'mgf',
#         path = '/home/zhouzw/Data_processing/20210330_debug/',
#         thread = 20,
#         polarity = 'positive',
#         int_ms2_min_abs = 50,
#         instrument = "SciexTripleTOF",
#         column = "hilic",
#         lib = 'zhuMetLib',
#         ce = "30",
#         method_lc = 'Amide23min',
#         is_rt_calibration = FALSE,
#         extension_step = '0',
#         scoring_approach_recursive = 'dp',
#         dp_tol = 0.5,
#         matched_frag_cutoff = 1,
#         use_redun_rm_result = TRUE,
#         is_check_data = FALSE,
#         is_anno_initial_seed = TRUE,
#         is_anno_mrn = TRUE,
#         is_credential = TRUE,
#         is_bio_interpret = TRUE,
#         is_exported_report = TRUE,
#         # dir_GenForm = '/share/home/zhuzhengjiang/share/bin',
#         comp_group = c("24W", "78W"),
#         uni_test = "t",
#         correct_p = FALSE,
#         p_cutoff = 0.05,
#         species = 'mmu',
#         quanti_pathway_method = 'mean')


# ms1_file = "data.csv"
# ms2_file = NULL
# metdna_version = 'version2'
# ms2_type = 'mgf'
# path = "/home/zhouzw/Data_processing/20210330_debug/"
#
# polarity = 'positive'
# instrument = "SciexTripleTOF"
# lib = 'zhuMetLib'
# column = 'hilic'
# ce = '30'
# method_lc = 'Amide23min'
# excluded_adduct = NULL
# is_rt_calibration = FALSE
#
# mz_tol = 25
# pf_rt_range = 0
# tolerance_rt_range = 30
# pf_ccs_range = 0
# tolerance_ccs_range = 2
# is_filter = TRUE
# is_rt_score = TRUE
# is_ccs_score = FALSE
# is_ms2_score = TRUE
#
# is_include_precursor = TRUE
# int_ms2_min_abs = 50
# int_ms2_min_relative = 0.01
# mz_tol_combine_ms1_ms2 = 25 # ppm
# rt_tol_combine_ms1_ms2 = 10 # s
# ccs_tol_combine_ms1_ms2 = NULL # %
# mz_tol_ms2 = 35
# dp_cutoff = 0.8
# matched_frag_cutoff = 1
# direction = 'reverse'
# scoring_approach = 'dp'
# # direction = 'reverse'
#
# is_plot_ms2 <- TRUE

# checkQuality(ms1_file = "data.csv",
#              sample_info_file = "sample.info.csv",
#              ms2_type = "mgf",
#              path = "/home/zhouzw/Data_processing/20210330_debug/")


# MetDNA2(ms1_file = 'data.csv',
#         sample_info_file = "sample.info.csv",
#         metdna_version = 'version2',
#         ms2_type = 'mgf',
#         path = '/home/zhouzw/Data_processing/20210330_debug/',
#         thread = 20,
#         polarity = 'positive',
#         int_ms2_min_abs = 50,
#         instrument = "IMMS",
#         column = "hilic",
#         lib = 'zhuMetLib',
#         ce = "30",
#         method_lc = 'Amide23min',
#         is_rt_calibration = FALSE,
#         extension_step = '0',
#         scoring_approach_recursive = 'dp',
#         dp_tol = 0.5,
#         matched_frag_cutoff = 1,
#         use_redun_rm_result = TRUE,
#         is_check_data = TRUE,
#         is_anno_initial_seed = TRUE,
#         is_anno_mrn = TRUE,
#         is_credential = TRUE,
#         is_bio_interpret = FALSE,
#         is_exported_report = FALSE,
#         # dir_GenForm = '/share/home/zhuzhengjiang/share/bin',
#         comp_group = c("24W", "78W"),
#         uni_test = "t",
#         correct_p = FALSE,
#         p_cutoff = 0.05,
#         species = 'mmu',
#         quanti_pathway_method = 'mean')


################################################################################
# # 20210413 test different instruments for v0.5.12
# MetDNA2(ms1_file = 'data.csv',
#         sample_info_file = "sample.info.csv",
#         metdna_version = 'version2',
#         ms2_type = 'msp',
#         path = '/home/zhouzw/Data_processing/20210319_different_instrument_test/QE/',
#         thread = 3,
#         polarity = 'positive',
#         instrument = "ThermoExploris",
#         column = "hilic",
#         lib = 'zhuMetLib_orbitrap',
#         ce = "SNCE20_30_40%",
#         method_lc = 'Amide12min',
        # extension_step = '0',
        # scoring_approach_recursive = 'dp',
        # use_redun_rm_result = FALSE,
        # is_check_data = TRUE,
        # is_anno_initial_seed = TRUE,
        # is_anno_mrn = TRUE,
        # is_credential = TRUE,
        # is_bio_interpret = FALSE,
        # is_exported_report = TRUE)

# MetDNA2(ms1_file = 'data.csv',
#         sample_info_file = "sample.info.csv",
#         metdna_version = 'version2',
#         ms2_type = 'mgf',
#         path = '/home/zhouzw/Data_processing/20210319_different_instrument_test/TIMS/',
#         thread = 3,
#         polarity = 'positive',
#         int_ms2_min_abs = 50,
#         instrument = "BrukerTIMS",
#         column = "hilic",
#         lib = 'zhuMetLib',
#         ce = "30",
#         extension_step = '0',
#         scoring_approach_recursive = 'dp',
#         use_redun_rm_result = FALSE,
#         is_check_data = TRUE,
#         is_anno_initial_seed = TRUE,
#         is_anno_mrn = TRUE,
#         is_credential = TRUE,
#         is_bio_interpret = FALSE,
#         is_exported_report = TRUE)

################################################################################

# MetDNA2(ms1_file = 'data.csv',
#         sample_info_file = "sample.info.csv",
#         metdna_version = 'version2',
#         ms2_type = 'mgf',
#         path = '/home/zhouzw/Data_processing/20210421_urine_xcms_pos_for_gnps/MetDNA2/',
#         thread = 3,
#         polarity = 'positive',
#         int_ms2_min_abs = 50,
#         instrument = "SciexTripleTOF",
#         column = "hilic",
#         lib = 'zhuMetLib',
#         ce = "30",
#         extension_step = '2',
#         scoring_approach_recursive = 'dp',
#         use_redun_rm_result = TRUE,
#         is_check_data = FALSE,
#         is_anno_initial_seed = TRUE,
#         is_anno_mrn = TRUE,
#         is_credential = TRUE,
#         is_bio_interpret = FALSE,
#         is_exported_report = FALSE)

################################################################################
# 20210423 test change test_assi_confidence_initial_seed -----------------------
#
#
# convertAnnotationTable2InitialId(peak_table_file = 'data.csv',
#                                  sample_info_file =  'sample.info.csv',
#                                  path_dir = '/home/zhouzw/Data_processing/20210330_debug/int30_met2_addct1/',
#                                  polarity = 'positive',
#                                  test_assi_confidence_initial_seed = TRUE,
#                                  tool = 'MetDNA2')
#
# path <- '/home/zhouzw/Data_processing/20210330_debug/int30_met2_addct1/'

################################################################################
# data('cpd_emrn', envir = environment())
# metabolite <- cpd_emrn
#
# # export MS2 match parameter
# switch (scoring_approach,
#         'dp' = {
#           intensityNormedMethod <- 'maximum'
#           methodScore <- 'dp'
#         },
#         'bonanza' = {
#           intensityNormedMethod <- 'bonanza'
#           methodScore <- 'bonanza'
#         },
#         'hybrid' = {
#           intensityNormedMethod <- 'maximum'
#           methodScore <- 'hybrid'
#         },
#         'gnps' = {
#           intensityNormedMethod <- 'gnps'
#           methodScore <- 'gnps'
#         }
# )
#
# matchParam <- SpectraTools::MatchParam(ppm = 25,
#                                        cutoff = 0,
#                                        weightIntensity = 1,
#                                        weightMZ = 0,
#                                        normIntensity = TRUE,
#                                        tuneLibSpectra = TRUE,
#                                        intensityExpNormed = TRUE,
#                                        intensityLibNormed = TRUE,
#                                        includePrecursor = TRUE,
#                                        ppmPrecursorFilter = 30,
#                                        thrIntensityAbs = 0,
#                                        thrIntensityRel = 0,
#                                        intensityNormedMethod = intensityNormedMethod,
#                                        methodMatch = 'direct',
#                                        methodScore = methodScore) %>%
#   new(Class = 'MatchParam')
#
# # data("kegg.compound", envir = environment())
# load(file.path('/home/zhouzw/Data_processing/20210330_debug/unlabeling_ecoli_neg_bonanza/01_result_initial_seed_annotation/', "00_intermediate_data", "ms2"))
#
# path_output <- file.path('/home/zhouzw/Data_processing/20210330_debug/unlabeling_ecoli_neg_bonanza/01_result_initial_seed_annotation',
#                          "02_result_MRN_annotation")
#
# load('/home/zhouzw/Data_processing/20210330_debug/unlabeling_ecoli_neg_bonanza/02_result_MRN_annotation/00_intermediate_data/tags2_after_redundancy_remove')
#
# plotNeighbor(tags2 = tags2_after_redundancy_remove,
#              path = path_output,
#              ms2 = ms2,
#              kegg.compound = metabolite,
#              matchParam = matchParam)


################################################################################
# 20210429
# load('/home/zhouzw/Data_processing/20210330_debug/aging_fly_pos_version_msp/01_result_initial_seed_annotation/00_intermediate_data/ms1_data')
# save(exp_ms2, file = '/home/zhouzw/Data_processing/20210330_debug/exp_ms2_210429.RData', version = 2)
# load('/home/zhouzw/Data_processing/20210330_debug/aging_fly_pos_version_msp/01_result_initial_seed_annotation/00_intermediate_data/lib_meta')
# load('/home/zhouzw/Data_processing/20210330_debug/aging_fly_pos_version_msp/01_result_initial_seed_annotation/00_intermediate_data/lib_spec')
# save(lib_ms2, file = '/home/zhouzw/Data_processing/20210330_debug/lib_ms2_210429.RData', version = 2)


################################################################################
# 20210508 debug MS2 match
# load('/home/zhouzw/Data_processing/20210330_debug/20210429_aging_fly_neg_initial_seed_test/01_result_initial_seed_annotation/00_intermediate_data/ms2')
# load('/home/zhouzw/Data_processing/20210330_debug/20210429_aging_fly_neg_initial_seed_test/01_result_initial_seed_annotation/00_intermediate_data/ms2_result')
# load('/home/zhouzw/Data_processing/20210330_debug/20210429_aging_fly_neg_initial_seed_test/01_result_initial_seed_annotation/00_intermediate_data/ms1_result')
# load('/home/zhouzw/Data_processing/20210330_debug/20210429_aging_fly_neg_initial_seed_test/01_result_initial_seed_annotation/00_intermediate_data/ms2_result')
# load('/home/zhouzw/Data_processing/20210330_debug/20210429_aging_fly_neg_initial_seed_test/01_result_initial_seed_annotation/00_intermediate_data/lib_meta')
# load('/home/zhouzw/Data_processing/20210330_debug/20210429_aging_fly_neg_initial_seed_test/01_result_initial_seed_annotation/00_intermediate_data/lib_spec')
# load('/home/zhouzw/Data_processing/20210330_debug/20210429_aging_fly_neg_initial_seed_test/01_result_initial_seed_annotation/00_intermediate_data/ms1_data')


################################################################################
# 20210515 modify MetDNA1 and MetDNA2 to keep same in 0 round ------------------

# ms1_file = "data.csv"
# ms2_file = NULL
# sample_info_file = "sample.info.csv"
# metdna_version = 'version1'
# ms2_type = 'mgf'
# # # ms2_type = 'msp'
# # path = "/home/zhouzw/Data_processing/20201102_metdna2_pos_development/"
# path = 'H:/00_projects/03_MetDNA2/00_data/20210508_debug/'
#
# # para of loadDB
# polarity = 'positive'
# instrument = "SciexTripleTOF"
# lib = 'zhuMetLib'
# column = 'hilic'
# ce = '30'
# # method_lc = 'Amide23min'
# method_lc = 'Other'
# excluded_adduct = NULL
# is_rt_calibration = FALSE
#
# # para of matchMs1WithSpecLib
# mz_tol = 25
# pf_rt_range = 0
# tolerance_rt_range = 30
# pf_ccs_range = 0
# tolerance_ccs_range = 4
# is_filter = TRUE
# is_rt_score = FALSE
# is_ccs_score = FALSE
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
# # ccs_tol_combine_ms1_ms2 = NULL # %
# mz_tol_ms2 = 35
# dp_cutoff = 0.8
# matched_frag_cutoff = 1
# direction = 'reverse'
# scoring_approach = 'dp'
#
# is_plot_ms2 <- TRUE
# test_adduct_version <- 'version1'


# annotateInitialSeed(ms1_file = ms1_file,
#                     ms2_file = ms2_file,
#                     metdna_version = metdna_version,
#                     ms2_type = ms2_type,
#                     path = 'H:/00_projects/03_MetDNA2/00_data/20210508_debug/metdna2/',
#                     polarity = polarity,
#                     instrument = instrument,
#                     lib = lib,
#                     column = column,
#                     ce = ce,
#                     method_lc = method_lc,
#                     is_rt_calibration = is_rt_calibration,
#                     mz_tol = mz_tol,
#                     tolerance_rt_range = tolerance_rt_range,
#                     tolerance_ccs_range = tolerance_ccs_range,
#                     mz_tol_combine_ms1_ms2 = mz_tol_combine_ms1_ms2,
#                     rt_tol_combine_ms1_ms2 = rt_tol_combine_ms1_ms2,
#                     ccs_tol_combine_ms1_ms2 = NULL,
#                     is_filter = FALSE,
#                     is_rt_score = FALSE,
#                     is_ccs_score = is_ccs_score,
#                     int_ms2_min_abs = int_ms2_min_abs,
#                     dp_cutoff = dp_cutoff,
#                     direction = direction,
#                     scoring_approach = 'dp',
#                     is_plot_ms2 = is_plot_ms2,
#                     test_adduct_version = test_adduct_version)
#
#
# load('H:/00_projects/03_MetDNA2/00_data/20210508_debug/metdna2/01_result_initial_seed_annotation/00_intermediate_data/result_annotation')
# test <- lapply(result_annotation, function(x){x@annotation_result})
# test <- test %>% dplyr::bind_rows()


# modify readInitialAnnotation -------------------------------------------------
# if (lib == 'zhuMetLib') {
#   data("zhuMetlib", envir = environment())
#   inHouse_compound <- zhuMetlib$meta$compound %>% as.data.frame()
#
#   data("md_zhumetlib", envir = environment())
#   md_inHouse_cpd <- md_zhumetlib
# }
# rm(zhuMetlib, md_zhumetlib)
#
# data <- readInitialAnnotation(data = "ms2_match_annotation_result.csv",
#                               metdna_version = 'version1',
#                               direction = 'reverse',
#                               inHouse_compound = inHouse_compound,
#                               instrument = "SciexTripleTOF",
#                               path = 'H:/00_projects/03_MetDNA2/00_data/20210508_debug/metdna2/01_result_initial_seed_annotation/')
#
# data$tags %>% dplyr::filter(!is.na(labid)) %>% nrow()
#
#
# if (metdna_version == 'version1') {
#   md_emrn_cpd <- md_mrn_emrn$version1
# } else {
#   md_emrn_cpd <- md_mrn_emrn$version2
# }
#
# data("md_zhumetlib", envir = environment())
# md_inHouse_cpd <- md_zhumetlib
# rm(md_zhumetlib)
#
# md_emrn_cpd <- md_mrn_emrn$version2
#
# rm(list = c('zhuMetlib', "fiehnHilicLib"));gc()
# rm(list = c('md_mrn_emrn', 'md_zhumetlib', "md_fiehnHilicLib"));gc()
#
# cat("\n");cat("Construct RT prediction model.\n")
#
# rt_result <- predictRT(data = data,
#                        prefer_adduct = "all",
#                        threads = 3,
#                        md_inHouse_cpd = md_inHouse_cpd,
#                        md_kegg_cpd = md_emrn_cpd,
#                        use_default_md = TRUE,
#                        column = 'hilic',
#                        method_lc = 'Other')
#
# path <- 'H:/00_projects/03_MetDNA2/00_data/20210508_debug/metdna2'
# metabolite
# metabolite_ccs
# metabolic_network
# inHouse_compound
# dir.create('H:/00_projects/03_MetDNA2/00_data/20210508_debug/metdna2/temp', showWarnings = FALSE, recursive = TRUE)
# save(metabolite, file = file.path(path, 'temp', 'metabolite'))
# save(metabolite_ccs, file = file.path(path, 'temp', 'metabolite_ccs'))
# save(metabolic_network, file = file.path(path, 'temp', 'metabolic_network'))
# save(inHouse_compound, file = file.path(path, 'temp', 'inHouse_compound'))
# save(ms2_data, file = file.path(path, 'temp', 'ms2_data'))
# save(adduct_table, file = file.path(path, 'temp', 'adduct_table'))
#

# path <- 'H:/00_projects/03_MetDNA2/00_data/20210508_debug/metdna2'
# load(file.path(path, 'temp', 'metabolite'))
# load(file.path(path, 'temp', 'metabolic_network'))
# load(file.path(path, 'temp', 'metabolite_ccs'))
# load(file.path(path, 'temp', 'inHouse_compound'))
# load(file.path(path, 'temp', 'ms2_data'))
# load(file.path(path, 'temp', 'adduct_table'))

# path <- 'H:/00_projects/03_MetDNA2/00_data/20210508_debug/metdna2'
# data = "ms2_match_annotation_result.csv"
# # ms2
# # adduct_table
# max_isotope = 4
# polarity = 'positive'
# metdna_version = 'version1'
# instrument = "SciexTripleTOF"
# direction = 'reverse'
# scoring_approach = 'dp'
# mz_tol = 25
# rt_tol1 = 3
# rt_tol2 = 30
# cor_tol = 0
# int_tol = 500
# ccs_tol = 4 # %
# dp_tol = 0.5
# matched_frag_tol = 1
# max_step = 3
# remain = FALSE
# remain_per = 0.5
# # metabolite
# # metabolite_ccs
# # metabolic_network
# # inHouse_compound
# threads = 3
# score_cutoff = 0
# # path = "."
# # output.path = ".",
# rt_filter = TRUE # Not available yet, need to be fix
# seed_neighbor_match_plot = FALSE
# candidate_num = 3
#
#
#
#
#
# data("zhuMetlib", envir = environment())
# inHouse_compound <- zhuMetlib$meta$compound %>% as.data.frame()
#
# data("md_zhumetlib", envir = environment())
# md_inHouse_cpd <- md_zhumetlib
#
# rm(zhuMetlib, md_zhumetlib)
# data <- readInitialAnnotation(data = data,
#                               path = file.path(path, "01_result_initial_seed_annotation"),
#                               metdna_version = metdna_version,
#                               direction = direction,
#                               rt_filter = rt_filter,
#                               rt_tol = rt_tol2,
#                               inHouse_compound = inHouse_compound,
#                               instrument = instrument)
#
# data$tags %>% dplyr::filter(!is.na(labid)) %>% nrow()
#
# save(data, file = 'H:/00_projects/03_MetDNA2/00_data/20210508_debug/metdna2/temp/data_rt_rm')
#
# test <- tags2 %>% trans2Matrix()

################################################################################
#
# MetDNA2(ms1_file = 'data.csv',
#         sample_info_file = "sample.info.csv",
#         metdna_version = 'version1',
#         test_adduct_version = 'version1',
#         ms2_type = 'mgf',
#         path = '/home/zhouzw/Data_processing/20210330_debug/metdna2/',
#         thread = 3,
#         polarity = 'positive',
#         instrument = "SciexTripleTOF",
#         column = "hilic",
#         lib = 'zhuMetLib',
#         ce = "30",
#         method_lc = 'Other',
#         extension_step = '0',
#         scoring_approach_recursive = 'dp',
#         is_rt_calibration = FALSE,
#         is_check_data = TRUE,
#         is_anno_initial_seed = TRUE,
#         is_anno_mrn = TRUE,
#         is_credential = FALSE,
#         is_bio_interpret = FALSE)


################################################################################
# load('/home/zhouzw/Data_processing/20210330_debug/aging_fly_pos_version1_msp/02_result_MRN_annotation/00_intermediate_data/tags2_after_redundancy_remove')
# load('/home/zhouzw/Data_processing/20210330_debug/aging_fly_pos_version1_msp/02_result_MRN_annotation/00_intermediate_data/id_result_redun_rm')
#
#
# test <- getAnnotationResult(tags2 = tags2_after_redundancy_remove,
#                             tags.result2 = id_result_redun_rm,
#                             metabolite = cpd_emrn,
#                             score.cutoff = 0.4,
#                             candidate.num = 5)

# tags2_after_redundancy_remove


################################################################################
# 20210520 debug ---------------------------------------------------------------

# path = '.'
# thread = 4
# is_cred_pg_filter = TRUE
# is_cred_formula_filter = TRUE
# is_pred_formula_all = FALSE
# direction = 'reverse'
# tolerance_rt_range = 30
# instrument = "SciexTripleTOF"
# test_rm_extra_anno_from_ini_seed = FALSE

# path = '/home/zhouzw/Data_processing/20210330_debug/200STD_mouse_liver_pos/'
# export_type = 'html'
# extension_step = '2'
# is_rt_calibration = TRUE
# is_bio_interpret = TRUE

#
#
# setGeneric(name = "trans2Matrix",
#            def = function(tags2,
#                           base = c("annotation", "peak")) {
#
#              base <- match.arg(base)
#              result <- lapply(seq_along(tags2), function(i) {
#                # cat(i, ' ')
#                x <- tags2[[i]]
#                name <- x@name
#                mz <- x@mz
#                rt <- x@rt
#                ccs <- x@ccs
#                annotation <- x@annotation
#                if(length(annotation) == 0){
#                  NULL
#                }else{
#                  # annotation <- lapply(annotation, function(x) {x <- x[names(x) != "addcut"]})
#                  annotation <- do.call(rbind, lapply(annotation, unlist))
#                  if(any(colnames(annotation) == "int.ratio.error") & any(colnames(annotation) == "int.error")){
#                    annotation[,"int.error"] <- annotation[,"int.ratio.error"]
#                    annotation <- annotation[,c(1:19)]
#                    colnames(annotation)[16] <- "int.error"
#                  }
#                  if(any(colnames(annotation) == "int.ratio.error") & all(colnames(annotation) != "int.error")){
#                    colnames(annotation)[16] <- "int.error"
#                  }
#
#                  data.frame(name, mz, rt, annotation,
#                             stringsAsFactors = FALSE)
#                }
#
#              })
#              result <- do.call(what = "rbind", args = result)
#              if(base == "annotation"){
#                result <- result[order(result$to),]
#              }else{
#                result <- result[order(result$name),]
#              }
#              return(result)
#            }
# )

#
#
# setGeneric(name = 'statTypeListPeakGroup',
#            function(
#              list_peak_group
#            ){
#              table_annotated <- purrr::map(list_peak_group, function(x){
#                x@peak_list_annotated %>%
#                  dplyr::mutate(peak_group = paste0(x@base_peak_name,
#                                                    '_',
#                                                    x@base_peak_adduct))
#              }) %>% dplyr::bind_rows()
#
#              if (nrow(table_annotated) == 0) {
#                result <- tibble::tibble('No. of unique peak' = 0,
#                                         'No. of peak' = 0,
#                                         'No. of isotope' = 0,
#                                         'No. of adduct' = 0,
#                                         'No. of NL' = 0,
#                                         'No. of ISF' = 0)
#
#                return(result)
#              }
#
#              table_annotated <- table_annotated %>%
#                tidyr::pivot_longer(cols = c('adductAnnotation',
#                                             'neutralLossAnnotation',
#                                             'isfAnnotation'),
#                                    names_to = 'adduct_nl_isf',
#                                    values_to = 'label') %>%
#                dplyr::filter(!is.na(label) | isotopeAnnotation %in% c('[M+1]',
#                                                                       '[M+2]',
#                                                                       '[M+3]',
#                                                                       '[M+4]'))
#
#              temp_isotope <- table_annotated %>%
#                dplyr::filter(isotopeAnnotation %in% c('[M+1]',
#                                                       '[M+2]',
#                                                       '[M+3]',
#                                                       '[M+4]')) %>%
#                dplyr::distinct(peak_name,
#                                isotopeAnnotation,
#                                peak_group,
#                                .keep_all = TRUE)
#
#              temp_adduct_nl_isf <- table_annotated %>%
#                dplyr::filter(!is.na(label))
#
#
#              num_unique_peak <- table_annotated %>%
#                dplyr::distinct(peak_name) %>%
#                dplyr::count() %>%
#                dplyr::pull(n)
#
#              num_peak <- nrow(temp_isotope) + nrow(temp_adduct_nl_isf)
#              num_isotope <- nrow(temp_isotope)
#              num_adduct <- temp_adduct_nl_isf %>%
#                dplyr::filter(adduct_nl_isf == 'adductAnnotation') %>% nrow()
#              num_nl <- temp_adduct_nl_isf %>%
#                dplyr::filter(adduct_nl_isf == 'neutralLossAnnotation') %>% nrow()
#              num_isf <- temp_adduct_nl_isf %>%
#                dplyr::filter(adduct_nl_isf == 'isfAnnotation') %>% nrow()
#
#              result <- tibble::tibble('No. of unique peak' = num_unique_peak,
#                                       'No. of peak' = num_peak,
#                                       'No. of isotope' = num_isotope,
#                                       'No. of adduct' = num_adduct,
#                                       'No. of NL' = num_nl,
#                                       'No. of ISF' = num_isf)
#
#              return(result)
#            })
#
#

################################################################################
# 20210522 ---------------------------------------------------------------------
# MetDNA2(ms1_file = 'data.csv',
#         sample_info_file = "sample.info.csv",
#         metdna_version = 'version2',
#         ms2_type = 'mgf',
#         path = '/home/zhouzw/Data_processing/20210522_debug/',
#         polarity = 'negative',
#         instrument = "ThermoExploris",
#         column = "hilic",
#         lib = 'zhuMetLib_orbitrap',
#         ce = "SNCE20_30_40%", # '20' for ce 20
#         method_lc = 'Amide12min',
#         extension_step = '2',
#         # dir_GenForm = '/share/home/zhuzhengjiang/share/bin',
#         thread = 2,
#         is_rt_calibration = TRUE,
#         is_check_data = TRUE,
#         is_anno_initial_seed = TRUE,
#         is_anno_mrn = TRUE,
#         is_credential = TRUE,
#         is_bio_interpret = FALSE,
#         is_exported_report = TRUE,
#         comp_group = c("QC", "Sample"))


################################################################################
# 20210526 debug ---------------------------------------------------------------
#
# load('/home/zhouzw/Data_processing/20210522_debug/02_result_MRN_annotation/00_intermediate_data/tags2_after_annotation')
#
# tags2_after_annotation

# tags2_after_annotation[[1]]@annotation
# names(tags2_after_annotation) <- sapply(tags2_after_annotation, function(x){
#   x@name
# })
#
# showTag2AnnotationMatrix <- function(tags2_after_annotation,
#                                      round) {
#
#   round_idx <- showTags2(tags2_after_annotation, "annotation.level") %>% lapply(unlist)
#
# }


################################################################################
# 20210527 add checkElement ----------------------------------------------------


# temp_error <- try({annotateMRN(annotation_result = "ms2_match_annotation_result.csv",
#                                ms2_data = 'ms2',
#                                prefer_adduct = 'M+H',
#                                metdna_version = 'version1',
#                                lib = 'zhuMetLib',
#                                direction = 'reverse',
#                                column = 'hilic',
#                                instrument = 'SciexTripleTOF',
#                                method_lc = 'Amide23min',
#                                polarity = 'positive',
#                                extension_step = '0',
#                                threads = 3,
#                                max_isotope = 4,
#                                path = '/home/zhouzw/Data_processing/20210522_debug/',
#                                mz_tol = 25,
#                                rt_tol1 = 3,
#                                rt_tol2 = 30,
#                                ccs_tol = 4,
#                                cor_tol = 0,
#                                int_tol = 500,
#                                dp_tol = 0.5,
#                                max_step = 3,
#                                score_cutoff = 0,
#                                remain = FALSE,
#                                remain_per = 0.5,
#                                scoring_approach = 'dp',
#                                matched_frag_tol = 1,
#                                seed_neighbor_match_plot = TRUE,
#                                candidate_num = 3,
#                                # test_old_mrn = test_old_mrn,
#                                test_adduct_version = 'version1')},
#                   silent = TRUE)

################################################################################
# # 20210601
# temp_adduct_table <- adduct_table
#
# met_result_mz <- lapply(seq_along(met_result$Node.ID), function(i){
#   temp_mz <- met_result$monoisotopic_mass[i]
#   result <- sapply(seq_along(temp_adduct_table$name), function(j){
#     calculateMz(exact_mass = temp_mz,
#                 adduct = temp_adduct_table$name[j],
#                 delta_mz = temp_adduct_table$massdiff[j])
#   })
#   result
# })
#
#
# calculateMz(exact_mass = temp_mz,
#             adduct = temp_adduct_table$name[j],
#             delta_mz = temp_adduct_table$massdiff[j],
#             nmol = temp_adduct_table$)
#
  # result <- sapply(seq_along(temp_adduct_table$name), function(j){
  #   calculateMz(exact_mass = temp_mz,
  #               adduct = temp_adduct_table$name[j],
  #               delta_mz = temp_adduct_table$massdiff[j],
  #               nmol = temp_adduct_table$nmol[j],
  #               ncharge = temp_adduct_table$charge[j])
  # })
#
#   calculateMz(exact_mass = 100,
#               adduct = "M+NH3+Cl",
#               delta_mz = 51.99600)

# calculateMz(exact_mass = 893.2197,
#             adduct = "M+NH3+Cl",
#             delta_mz = 51.99600)


################################################################################
# 20210604 for lwb -------------------------------------------------------------
#
# library(MetDNA2)
# MetDNA2(ms1_file = 'data.csv',
#         sample_info_file = "sample.info.csv",
#         metdna_version = 'version2',
#         ms2_type = 'msp',
#         path = '/home/zhouzw/Data_processing/20210522_debug/',
#         polarity = 'positive',
#         instrument = "SciexTripleTOF",
#         column = "hilic",
#         lib = 'zhuMetLib',
#         ce = "30",
#         method_lc = 'Amide12min',
#         extension_step = '0',
#         # dir_GenForm = '/share/home/zhuzhengjiang/share/bin',
#         thread = 20,
#         is_check_data = FALSE,
#         is_anno_initial_seed = FALSE,
#         is_anno_mrn = FALSE,
#         is_credential = TRUE,
#         is_rt_calibration = TRUE)

# load('/home/zhouzw/Data_processing/20210522_debug/03_annotation_credential/00_intermediate_data/annot_all_before_credential.RData')
#
# temp <- mergeRecursiveAnnotation(path = '/home/zhouzw/Data_processing/20210522_debug/',
#                                  thread = 3,
#                                  lib = 'zhuMetLib',
#                                  direction = 'reverse',
#                                  tolerance_rt_range = 30,
#                                  use_redun_rm_result = TRUE)

#
# cat('Change formats for annotation credential.\n')
# convertAnnotationTable2InitialId(peak_table_file = "data.csv",
#                                  sample_info_file =  "sample.info.csv",
#                                  path_dir = '/home/zhouzw/Data_processing/20210522_debug/',
#                                  polarity = 'positive',
#                                  # test_assi_confidence_initial_seed = test_assi_confidence_initial_seed,
#                                  tool = 'MetDNA2')

# path = '/home/zhouzw/Data_processing/20210522_debug/'
# lib = 'zhuMetLib'
# thread = 4
# direction = 'reverse'
# is_cred_pg_filter = TRUE
# is_cred_formula_filter = TRUE
# tolerance_rt_range = 30

# authenticateAnnotationCredential(peak_table_file = "peak_table.csv",
#                                  annotation_initial_file = "annotation_initial.csv",
#                                  ms2_data_file = "ms2_data.RData",
#                                  path_dir = '/home/zhouzw/Data_processing/20210522_debug/',
#                                  polarity = 'positive',
#                                  thread = 3,
#                                  isotope_int_ratio_check = TRUE, # para annotateIsotope
#                                  isotope_int_ratio_cutoff = 500,
#                                  is_ms2_check = TRUE, # compare ms2 similarity of neutral loss
#                                  ms2_score_cutoff = -1, # -1 represent not filter
#                                  is_plot_pseudo_MS1 = TRUE,
#                                  # formula prediction
#                                  dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
#                                  is_pred_formula_all = FALSE,
#                                  ppm = 25, # ms1 tolerance
#                                  acc = 25, # ms2 tolerance
#                                  elements = 'CHONPS',
#                                  num_formula_candidate = 3, # return top 3 candidate
#                                  platform = 'linux')

################################################################################
#
# path_dir <- '/home/zhouzw/Data_processing/20210522_debug/01_aging_fly_pos_version2/'
# load('/home/zhouzw/Data_processing/20210522_debug/01_aging_fly_pos_version2/03_annotation_credential/00_intermediate_data/list_peak_group_annotation_concised.RData')

# load('/home/zhouzw/Data_processing/20210522_debug/01_aging_fly_pos_version2/03_annotation_credential/00_intermediate_data/list_peak_group_annotation.RData')
#
# path <- '/home/zhouzw/Data_processing/20210522_debug/01_aging_fly_pos_version2/'
#
# load('/home/zhouzw/Data_processing/20210522_debug/01_aging_fly_pos_version2/00_annotation_table/00_intermediate_data/annot_all')
# load('/home/zhouzw/Data_processing/20210522_debug/01_aging_fly_pos_version2/00_annotation_table/00_intermediate_data/')



# path <- '/home/zhouzw/Data_processing/20210522_debug/test/'
# load('/home/zhouzw/Data_processing/20210522_debug/test/00_annotation_table/00_intermediate_data/annot_all')
#
################################################################################
# path = '/home/zhouzw/Data_processing/20210522_debug/01_aging_fly_pos_version2/'
# thread = 4
# is_cred_pg_filter = TRUE
# is_cred_formula_filter = FALSE
# is_pred_formula_all = FALSE
# tolerance_rt_range = 30
# direction = 'reverse'
# instrument = 'SciexTripleTOF'
# mz_tol = 25
# rt_tol2 = 30
# ccs_tol = 4
# candidate_num = 5
# lib <- 'zhuMetLib'

################################################################################
#
# library(MetDNA2)
# MetDNA2(ms1_file = 'data.csv',
#         sample_info_file = "sample.info.csv",
#         metdna_version = 'version2',
#         ms2_type = 'mgf',
#         path = '/home/zhouzw/Data_processing/20210522_debug/test/',
#         polarity = 'positive',
#         instrument = "SciexTripleTOF",
#         column = "hilic",
#         lib = 'zhuMetLib',
#         ce = "30",
#         method_lc = 'Amide23min',
#         extension_step = '0',
#         # dir_GenForm = '/share/home/zhuzhengjiang/share/bin',
#         thread = 3,
#         is_check_data = FALSE,
#         is_anno_initial_seed = FALSE,
#         is_anno_mrn = FALSE,
#         is_credential = FALSE,
#         is_bio_interpret = FALSE,
#         is_exported_report = TRUE,
#         is_rt_calibration = FALSE)

# load('/home/zhouzw/Data_processing/20210522_debug/test/01_result_initial_seed_annotation/00_intermediate_data/ms1_result')
# load('/home/zhouzw/Data_processing/20210522_debug/test/01_result_initial_seed_annotation/00_intermediate_data/lib_meta')
# load('/home/zhouzw/Data_processing/20210522_debug/test/01_result_initial_seed_annotation/00_intermediate_data/ms1_data')

# /home/zhouzw/Data_processing/20210522_debug/extended_step2/

#
# MetDNA2(ms1_file = 'data.csv',
#         sample_info_file = "sample.info.csv",
#         metdna_version = 'version2',
#         ms2_type = 'mgf',
#         path = '/home/zhouzw/Data_processing/20210522_debug/extended_step2/',
#         polarity = 'positive',
#         instrument = "SciexTripleTOF",
#         column = "hilic",
#         lib = 'zhuMetLib',
#         ce = "30",
#         method_lc = 'Amide23min',
#         extension_step = '2',
#         # dir_GenForm = '/share/home/zhuzhengjiang/share/bin',
#         thread = 3,
#         is_check_data = TRUE,
#         is_anno_initial_seed = TRUE,
#         is_anno_mrn = TRUE,
#         is_credential = TRUE,
#         is_bio_interpret = FALSE,
#         is_exported_report = FALSE,
#         is_rt_calibration = FALSE)

################################################################################
#
# library(MetDNA2)
# MetDNA2(ms1_file = 'data.csv',
#         sample_info_file = "sample.info.csv",
#         metdna_version = 'version2',
#         ms2_type = 'mgf',
#         path = '/home/zhouzw/Data_processing/20210522_debug/POS/',
#         thread = 3,
#         polarity = 'positive',
#         int_ms2_min_abs = 50,
#         instrument = "SciexTripleTOF",
#         column = "hilic",
#         lib = 'zhuMetLib',
#         ce = "30",
#         method_lc = 'Amide23min',
#         is_rt_calibration = TRUE,
#         extension_step = '2',
#         scoring_approach_recursive = 'dp',
#         dp_tol = 0.5,
#         matched_frag_cutoff = 1,
#         use_redun_rm_result = TRUE,
#         is_check_data = FALSE,
#         is_anno_initial_seed = FALSE,
#         is_anno_mrn = FALSE,
#         is_credential = FALSE,
#         is_bio_interpret = FALSE,
#         is_exported_report = FALSE,
#         # dir_GenForm = '/share/home/zhuzhengjiang/share/bin',
#         comp_group = c("W03", "W30"),
#         uni_test = "t",
#         correct_p = TRUE,
#         p_cutoff = 0.05,
#         species = 'dme',
#         quanti_pathway_method = 'mean',
#         test_adduct_version = 'version2',
#         candidate_num = 5)
#
# path = '/home/zhouzw/Data_processing/20210522_debug/POS/'

################################################################################

# load('/home/zhouzw/Data_processing/20210626_metdna2_merge_pos_neg/POS_and_NEG/04_biology_intepretation/00_intermediate_data/result_pathway_enrichment.RData')
# path_output <- '/home/zhouzw/Data_processing/20210626_metdna2_merge_pos_neg/POS_and_NEG/04_biology_intepretation'

################################################################################
# load('/home/zhouzw/Data_processing/20210522_debug/WHM/20210624_WHM_metdna2_CRC_capRT_neg_V2/00_annotation_table/00_intermediate_data/annot_all')
# load('/home/zhouzw/Data_processing/20210522_debug/WHM/20210624_WHM_metdna2_CRC_capRT_neg_V2/00_annotation_table/00_intermediate_data/list_identification')
# load('/home/zhouzw/Data_processing/20210522_debug/WHM/20210624_WHM_metdna2_CRC_capRT_neg_V2/00_annotation_table/00_intermediate_data/table_identification')

################################################################################
#
# load('/home/zhouzw/Data_processing/20210522_debug/neg_mode/03_annotation_credential/00_intermediate_data/list_peak_group_annotation_concised.RData')
# load('/home/zhouzw/Data_processing/20210522_debug/neg_mode/03_annotation_credential/ms2_data.RData')
# load(file.path(path_dir, '03_annotation_credential', "00_intermediate_data", 'list_peak_group_formula.RData'))
# load('/home/zhouzw/Data_processing/20210522_debug/neg_mode/02_result_MRN_annotation/00_intermediate_data/tags2_after_redundancy_remove')
#
# path_dir <- "/home/zhouzw/Data_processing/20210522_debug/neg_mode/"
# authenticateFormulaCredential(list_peak_group = list_peak_group_annotation_concised,
#                               ms2_data = raw_msms,
#                               polarity = 'negative',
#                               path_dir = path_dir,
#                               thread = 3,
#                               dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
#                               ppm = 25,
#                               acc = 25,
#                               platform = 'linux',
#                               elements = 'CHONPS',
#                               num_formula_candidate = 3)
#
# setwd(temp_path)


################################################################################
# library(MetDNA2)
# runMetDNA2(ms1_file_pos = "data.csv",
#            sample_info_file_pos = "sample.info.csv",
#            path_pos = "/home/zhouzw/Data_processing/20210522_debug/NEG",
#            ms1_file_neg = "data.csv",
#            sample_info_file_neg = "sample.info.csv",
#            path_neg = "/home/zhouzw/Data_processing/20210522_debug/POS",
#            metdna_version = 'version2',
#            ms2_type = 'mgf',
#            thread = 3,
#            lib = 'zhuMetLib',
#            polarity = 'positive',
#            instrument = "SciexTripleTOF",
#            column = "hilic",
#            ce = "30",
#            method_lc = 'Amide23min',
#            is_rt_calibration = FALSE,
#            extension_step = '0',
#            comp_group = c("W03", "W30"),
#            # dir_GenForm = '/share/home/zhuzhengjiang/share/bin',
#            species = 'dme',
#            quanti_pathway_method = 'mean')

################################################################################
# temp_path <- '/home/zhouzw/Data_processing/20210522_debug/both_polarity/'
# runMetDNA2(ms1_file_pos = "data.csv",
#            ms1_file_neg = "data.csv",
#            sample_info_file_pos = "sample.info.csv",
#            sample_info_file_neg = "sample.info.csv",
#            path_pos = file.path(temp_path, "POS"),
#            path_neg = file.path(temp_path, "NEG"),
#            metdna_version = 'version2',
#            ms2_type = 'mgf',
#            thread = 20,
#            lib = 'zhuMetLib',
#            polarity = 'both',
#            instrument = "SciexTripleTOF",
#            column = "hilic",
#            ce = "30",
#            method_lc = 'Amide23min',
#            is_rt_calibration = FALSE,
#            # dir_GenForm = '/share/home/zhuzhengjiang/share/bin',
#            extension_step = '0',
#            comp_group = c("W03", "W30"),
#            species = 'dme',
#            quanti_pathway_method = 'mean')

################################################################################
# runMetDNA2(ms1_file_pos = "data.csv",
#            sample_info_file_pos = "sample.info.csv",
#            path_pos = "/home/zhouzw/Data_processing/20210522_debug/pos_version2",
#            metdna_version = 'version2',
#            ms2_type = 'msp',
#            thread = 20,
#            lib = 'zhuMetLib',
#            polarity = 'positive',
#            instrument = "SciexTripleTOF",
#            column = "hilic",
#            ce = "30",
#            method_lc = 'Amide23min',
#            is_rt_calibration = TRUE,
#            extension_step = '0',
#            # dir_GenForm = '/share/home/zhuzhengjiang/share/bin',
#            comp_group = c("W03", "W30"),
#            species = 'dme',
#            quanti_pathway_method = 'mean')
#
# convertAnnotationTable2InitialId(peak_table_file = 'data.csv',
#                                  sample_info_file =  'sample.info.csv',
#                                  path_dir = '/home/zhouzw/Data_processing/20210522_debug/pos_version2/',
#                                  polarity = 'positive',
#                                  # test_assi_confidence_initial_seed = test_assi_confidence_initial_seed,
#                                  tool = 'MetDNA2')
#
# generateMetDNA2AnnotationResult(path = '/home/zhouzw/Data_processing/20210522_debug/pos_version2/',
#                                 thread = 3,
#                                 is_cred_pg_filter = TRUE,
#                                 is_cred_formula_filter = FALSE,
#                                 is_pred_formula_all = FALSE,
#                                 tolerance_rt_range = 30,
#                                 direction = 'reverse',
#                                 instrument = 'SciexTripleTOF',
#                                 mz_tol = 25,
#                                 rt_tol = 30,
#                                 ccs_tol = 3,
#                                 candidate_num = 5)
#
# temp <- mergeRecursiveAnnotation(path = '/home/zhouzw/Data_processing/20210522_debug/pos_version2/',
#                                  thread = 3,
#                                  lib = 'zhuMetLib',
#                                  direction = 'reverse',
#                                  tolerance_rt_range = 30,
#                                  use_redun_rm_result = TRUE)
#
# len1 <- sapply(annot_all, function(x){
#   nrow(x@initial_seed_annotation)
# })
# len2 <- sapply(annot_all, function(x){
#   nrow(x@recursive_annotation)
# })
#
# which(len1> 0 & len2>0)
#
# which(len1> 0 & len2==0)
#
# which(len1 == 0 & len2 > 0)

################################################################################
#
# runMetDNA2(ms1_file_pos = "data.csv",
#            sample_info_file_pos = "sample.info.csv",
#            path_pos = "/home/zhouzw/Data_processing/20210522_debug/metdna2/",
#            metdna_version = 'version2',
#            ms2_type = 'mgf',
#            thread = 20,
#            lib = 'zhuMetLib',
#            polarity = 'positive',
#            instrument = "SciexTripleTOF",
#            column = "hilic",
#            ce = "30",
#            method_lc = 'Amide12min',
#            is_rt_calibration = TRUE,
#            extension_step = '0',
#            is_check_data = TRUE,
#            is_anno_initial_seed = FALSE,
#            is_bio_interpret = FALSE,
#            # dir_GenForm = '/share/home/zhuzhengjiang/share/bin',
#            comp_group = c("W03", "W30"),
#            species = 'dme',
#            quanti_pathway_method = 'mean')
#
# interpretBiology(sample_file_pos =  "data.csv",
#                  sample_info_file_pos = "sample.info.csv",
#                  table_identification_pair_file_pos = 'table3_identification_pair.csv',
#                  polarity = 'positive',
#                  path = "/home/zhouzw/Data_processing/20210522_debug/metdna2/",
#                  metdna_version = 'version2',
#                  comp_group = c("W03", "W30"),
#                  uni_test = 't',
#                  correct_p = TRUE,
#                  p_cutoff = 0.05,
#                  fc_cutoff = 1,
#                  species = 'dme',
#                  quanti_pathway_method = 'mean',
#                  organze_basis = 'kegg_id',
#                  extension_step = '0')
#
# checkQualitySampleGroup(ms1_file = "data.csv",
#                         sample_info_file = "sample.info.csv",
#                         group = c("W03", "W30"),
#                         path = "/home/zhouzw/Data_processing/20210522_debug/metdna2/")

# library(MetDNA2)
# runMetDNA2(ms1_file_pos = "data.csv",
#            ms1_file_neg = "data.csv",
#            sample_info_file_pos = "sample.info.csv",
#            sample_info_file_neg = "sample.info.csv",
#            path_pos = "/home/zhouzw/Data_processing/20210522_debug/both_emrn0/POS/",
#            path_neg = "./NEG/",
#            metdna_version = 'version2',
#            ms2_type = 'msp',
#            thread = 3,
#            lib = 'zhuMetLib',
#            polarity = 'positive',
#            instrument = "SciexTripleTOF",
#            column = "hilic",
#            ce = "30",
#            method_lc = 'Amide23min',
#            # dir_GenForm = '/share/home/zhuzhengjiang/share/bin',
#            is_rt_calibration = FALSE,
#            is_anno_initial_seed = FALSE,
#            is_anno_mrn = FALSE,
#            extension_step = '0',
#            comp_group = c("W03", "W30"),
#            species = 'dme',
#            quanti_pathway_method = 'mean')
#
#
# runMetDNA2(ms1_file_pos = "data.csv",
#            ms1_file_neg = "data.csv",
#            sample_info_file_pos = "sample.info.csv",
#            sample_info_file_neg = "sample.info.csv",
#            path_pos = "/home/zhouzw/Data_processing/20210522_debug/both_emrn0/POS/",
#            path_neg = "/home/zhouzw/Data_processing/20210522_debug/both_emrn0/NEG/",
#            metdna_version = 'version2',
#            ms2_type = 'msp',
#            thread = 3,
#            lib = 'zhuMetLib',
#            polarity = 'both',
#            instrument = "SciexTripleTOF",
#            column = "hilic",
#            ce = "30",
#            method_lc = 'Amide23min',
#            # dir_GenForm = '/share/home/zhuzhengjiang/share/bin',
#            is_rt_calibration = FALSE,
#            is_anno_initial_seed = FALSE,
#            is_anno_mrn = FALSE,
#            is_credential = FALSE,
#            extension_step = '0',
#            comp_group = c("W03", "W30"),
#            species = 'dme',
#            quanti_pathway_method = 'mean')
#
#
################################################################################
#
# runMetDNA2(ms1_file_pos = "data.csv",
#            sample_info_file_pos = "sample.info.csv",
#            path_pos = "/home/zhouzw/Data_processing/20210522_debug/01_aging_fly_pos/",
#            metdna_version = 'version1',
#            ms2_type = 'mgf',
#            thread = 3,
#            lib = 'zhuMetLib',
#            polarity = 'positive',
#            instrument = "SciexTripleTOF",
#            column = "hilic",
#            ce = "30",
#            method_lc = 'Other',
#            is_rt_calibration = FALSE,
#            extension_step = '0',
#            is_check_data = TRUE,
#            is_anno_initial_seed = TRUE,
#            is_bio_interpret = TRUE,
#            # dir_GenForm = '/share/home/zhuzhengjiang/share/bin',
#            comp_group = c("W03", "W30"),
#            species = 'dme',
#            quanti_pathway_method = 'mean')

################################################################################
# path_dir <- '/home/zhouzw/Data_processing/20210522_debug/01_aging_fly_pos'
# load('/home/zhouzw/Data_processing/20210522_debug/01_aging_fly_pos/03_annotation_credential/00_intermediate_data/list_peak_group.RData')
# test2 <- sapply(list_peak_group, function(x){x@base_peak_adduct})
# which(!(test2 %in% c(lib_adduct_nl$positive$adduct, lib_adduct_nl$negative$adduct)))
# load('/home/zhouzw/Data_processing/20210522_debug/01_aging_fly_pos/03_annotation_credential/peak_table.csv')
# load('/home/zhouzw/Data_processing/20210522_debug/01_aging_fly_pos/03_annotation_credential/00_intermediate_data/annotation_initial.RData')
# load('/home/zhouzw/Data_processing/20210522_debug/01_aging_fly_pos/01_result_initial_seed_annotation/00_intermediate_data/result_annotation')
# load('/home/zhouzw/Data_processing/20210522_debug/01_aging_fly_pos/03_annotation_credential/00_intermediate_data/annot_all_before_credential.RData')
# test_name <- sapply(tags2_after_redundancy_remove, function(x){x@name})
# names(tags2_after_redundancy_remove) <- test_name
# load('/home/zhouzw/Data_processing/20210522_debug/01_aging_fly_pos/03_annotation_credential/00_intermediate_data/id_merge_before_credential.RData')

################################################################################
#
# sample_file_pos = 'data.csv'
# sample_info_file_pos = 'sample.info.csv'
# table_identification_pair_file_pos = 'table3_identification_pair.csv'
# # sample_file_neg = 'data.csv'
# # sample_info_file_neg = 'sample.info.csv'
# # table_identification_pair_file_neg = 'table3_identification_pair.csv'
# path = '/home/zhouzw/Data_processing/20210705_compare_metdna1_metdna2_pathway/20210702_LMD_pathway_comparsion/02_metdna2/01_aging_fly_pos/'
# correct_p = TRUE
# comp_group <- c('W03','W30')
#
# load('/home/zhouzw/Data_processing/20210705_compare_metdna1_metdna2_pathway/20210702_LMD_pathway_comparsion/temp/MetDNA1_sig_ids_210705.RData')
# MetBioInterpretation::enrichPathway(met_sig = id,
#                                     species = species,
#                                     lib_pathway = lib_pathway,
#                                     test_method = 'hypergeometric')
#

################################################################################
#
# sample_file_pos = 'data.csv'
# sample_info_file_pos = 'sample.info.csv'
# table_identification_pair_file_pos = 'table3_identification_pair.csv'
# sample_file_neg = 'data.csv'
# sample_info_file_neg = 'sample.info.csv'
# table_identification_pair_file_neg = 'table3_identification_pair.csv'
# path = '/home/zhouzw/Data_processing/20210705_compare_metdna1_metdna2_pathway/20210702_LMD_pathway_comparsion/02_metdna2/03_both'
# polarity = 'both'
# metdna_version = 'version1'
# comp_group = c("W30", "W03")
# uni_test = "t"
# correct_p = TRUE
# p_cutoff = 0.05
# fc_cutoff = 1
# species = 'dme'
# quanti_pathway_method = 'mean'
# organze_basis = 'kegg_id'
# extension_step = '0'
#
# runMetDNA2(ms1_file_pos = "data.csv",
#            ms1_file_neg = "data.csv",
#            sample_info_file_pos = "sample.info.csv",
#            sample_info_file_neg = "sample.info.csv",
#            path_pos = '/home/zhouzw/Data_processing/20210705_compare_metdna1_metdna2_pathway/20210702_LMD_pathway_comparsion/02_metdna2/03_both/POS',
#            path_neg = '/home/zhouzw/Data_processing/20210705_compare_metdna1_metdna2_pathway/20210702_LMD_pathway_comparsion/02_metdna2/03_both/NEG',
#            metdna_version = 'version1',
#            ms2_type = 'mgf',
#            thread = 3,
#            lib = 'zhuMetLib',
#            is_anno_initial_seed = FALSE,
#            is_anno_mrn = FALSE,
#            is_credential = FALSE,
#            polarity = 'both',
#            instrument = "SciexTripleTOF",
#            column = "hilic",
#            ce = "30",
#            method_lc = 'Other',
#            # dir_GenForm = '/share/home/zhuzhengjiang/share/bin',
#            is_rt_calibration = FALSE,
#            extension_step = '0',
#            comp_group = c("W03", "W30"),
#            species = 'dme',
#            quanti_pathway_method = 'mean')


################################################################################

# runMetDNA2(
#   path_pos = "/home/zhouzw/Data_processing/20210522_debug/01_metanalyzer_pos/",
#   # path_neg = "/home/zhouzw/Data_processing/20210706_metdna2test_dev/02_demo_data/01_mgf/03_aging_fly_both/NEG/",
#   metdna_version = 'version2',
#   polarity = 'positive',
#   instrument = "SciexTripleTOF",
#   column = "hilic",
#   ce = "30",
#   method_lc = 'Amide23min',
#   is_rt_calibration = TRUE,
#   extension_step = '0',
#   comp_group = c("W03", "W30"),
#   species = 'dme',
#   p_cutoff = 0.05,
#   fc_cutoff = 1,
#   correct_p = TRUE)

# runMetDNA2(
#   path_pos = "/home/zhouzw/Data_processing/20210522_debug/03_msdial_pos/",
#   # path_neg = "/home/zhouzw/Data_processing/20210706_metdna2test_dev/02_demo_data/01_mgf/03_aging_fly_both/NEG/",
#   metdna_version = 'version2',
#   polarity = 'positive',
#   instrument = "SciexTripleTOF",
#   column = "hilic",
#   ce = "30",
#   method_lc = 'Amide23min',
#   is_rt_calibration = TRUE,
#   extension_step = '0',
#   comp_group = c("W03", "W30"),
#   species = 'dme',
#   p_cutoff = 0.05,
#   fc_cutoff = 1,
#   correct_p = TRUE)


################################################################################
# 20210709 Bug2 ----------------------------------------------------------------
# load('/home/zhouzw/Data_processing/20210522_debug/POS_bug2/00_annotation_table/00_intermediate_data/annot_all')
# annot_all$M258T301@initial_seed_annotation
# path <- '/home/zhouzw/Data_processing/20210522_debug/POS_bug2'
# load('/home/zhouzw/Data_processing/20210522_debug/POS_bug2/00_annotation_table/00_intermediate_data/peak_group_id_table')
#
# Bug1 -------------------------------------------------------------------------
# load('/home/zhouzw/Data_processing/20210522_debug/POS_bug1/03_annotation_credential/00_intermediate_data/list_peak_group.RData')
# load('/home/zhouzw/Data_processing/20210522_debug/POS_bug1/03_annotation_credential/ms2_data.RData')
# load('/home/zhouzw/Data_processing/20210522_debug/POS_bug1/03_annotation_credential/00_intermediate_data/list_peak_group_annotation_concised.RData')
# load('/home/zhouzw/Data_processing/20210522_debug/POS_bug1/03_annotation_credential/00_intermediate_data/list_peak_group_annotation.RData')
# ms2_data <- raw_msms
# path_dir <- '/home/zhouzw/Data_processing/20210522_debug/POS_bug1'
#
# authenticatePeakGroupCredential(list_peak_group = list_peak_group,
#                                 ms2_data = raw_msms,
#                                 path_dir = '/home/zhouzw/Data_processing/20210522_debug/POS_bug1',
#                                 polarity = 'positive',
#                                 tol_mz = 25,
#                                 # para annotateIsotope
#                                 isotope_int_ratio_check = TRUE,
#                                 isotope_int_ratio_cutoff = 500,
#                                 is_ms2_check = TRUE, # compare ms2 similarity of neutral loss
#                                 ms2_score_cutoff = -1, # -1 represent not filter
#                                 cutoff_ssc = 0.3,
#                                 cutoff_ssc_int = 3000,
#                                 is_rule_limitation = TRUE,
#                                 cutoff_topN = 5,
#                                 is_plot_pseudo_MS1 = TRUE,
#                                 thread = 4,
#                                 type_order = c('level1', 'level2', 'level3'))

#
# result <- annotatePeakGroup(peak_group = list_peak_group[[i]],
#                             ms2_data = ms2_data,
#                             polarity = 'positive',
#                             tol_mz = 25,
#                             # para annotateIsotope
#                             isotope_int_ratio_check = TRUE,
#                             isotope_int_ratio_cutoff = 500,
#                             is_ms2_check = TRUE, # compare ms2 similarity of neutral loss
#                             ms2_score_cutoff = -1, # -1 represent not filter
#                             cutoff_ssc = 0.3,
#                             cutoff_ssc_int = 3000,
#                             is_rule_limitation = TRUE,
#                             cutoff_topN = 5)

################################################################################
# 20210710 debug local-server --------------------------------------------------
#
# runMetDNA2(
#   path_pos = "/home/zhouzw/Data_processing/20210522_debug/POS",
#   path_neg = "/home/zhouzw/Data_processing/20210522_debug/POS",
#   metdna_version = "version2",
#   polarity = "positive",
#   instrument = "SciexTripleTOF",
#   column = "hilic",
#   ce = "30",
#   method_lc = "Amide23min",
#   correct_p = TRUE,
#   extension_step = "0",
#   comp_group = c("W03", "W30"),
#   species = "hsa",
#   p_cutoff = 0.050000,
#   fc_cutoff = 1.000000,
#   is_rt_calibration = TRUE)

# lib <- 'zhuMetLib'
# instrument <- 'SciexTripleTOF'
# polarity <- 'positive'
# column <- 'hilic'
# ce <- '30'
# is_rt_calibration <- TRUE

# test <- loadSpecDB(lib = 'zhuMetLib_orbitrap',
#                    instrument = 'ThermoOrbitrap',
#                    polarity = 'positive',
#                    column = 'hilic',
#                    method_lc = 'Amide23min',
#                    ce = 'SCE20_30_40%',
#                    adduct_list = adduct_list,
#                    path = path,
#                    is_rt_calibration = FALSE)

################################################################################
# load('/home/zhouzw/Data_processing/20210522_debug/01_aging_fly_pos/03_annotation_credential/00_intermediate_data/annotation_initial.RData')
# load('/home/zhouzw/Data_processing/20210522_debug/01_aging_fly_pos/03_annotation_credential/00_intermediate_data/list_peak_group_formula.RData')
# load('/home/zhouzw/Data_processing/20210522_debug/01_aging_fly_pos/03_annotation_credential/00_intermediate_data/list_peak_group_annotation_concised.RData')
# path_dir <- '/home/zhouzw/Data_processing/20210522_debug/01_aging_fly_pos'
#
# refineAnnotation(annotation_initial,
#                  list_peak_group_annotation_concised,
#                  list_peak_group_formula,
#                  path_dir = path_dir)
#
# runMetDNA2(ms1_file_pos = "data.csv",
#            sample_info_file_pos = "sample.info.csv",
#            path_pos = "/home/zhouzw/Data_processing/20210522_debug/01_aging_fly_pos/",
#            metdna_version = 'version2',
#            ms2_type = 'mgf',
#            thread = 4,
#            lib = 'zhuMetLib',
#            polarity = 'positive',
#            instrument = "SciexTripleTOF",
#            column = "hilic",
#            ce = "30",
#            method_lc = 'Amide23min',
#            is_rt_calibration = FALSE,
#            extension_step = '0',
#            is_check_data = TRUE,
#            is_anno_initial_seed = FALSE,
#            is_anno_mrn = FALSE,
#            is_credential = TRUE,
#            is_bio_interpret = FALSE,
#            dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
#            comp_group = c("W03", "W30"),
#            p_cutoff = 0.01,
#            correct_p = TRUE,
#            species = 'dme',
#            quanti_pathway_method = 'mean')
#
# load('/home/zhouzw/Data_processing/20210522_debug/01_aging_fly_pos/00_annotation_table/00_intermediate_data/annot_all')
# path <- '/home/zhouzw/Data_processing/20210522_debug/01_aging_fly_pos'

################################################################################
# 20210713 ---------------------------------------------------------------------
#
# runMetDNA2(ms1_file_pos = "data.csv",
#            sample_info_file_pos = "sample.info.csv",
#            path_pos = "/home/zhouzw/Data_processing/20210522_debug/01_aging_fly_pos/",
#            metdna_version = 'version2',
#            ms2_type = 'mgf',
#            thread = 20,
#            lib = 'zhuMetLib',
#            polarity = 'positive',
#            instrument = "SciexTripleTOF",
#            column = "hilic",
#            ce = "30",
#            method_lc = 'Amide23min',
#            is_rt_calibration = FALSE,
#            extension_step = '0',
#            is_check_data = TRUE,
#            is_anno_initial_seed = TRUE,
#            is_anno_mrn = TRUE,
#            is_credential = TRUE,
#            is_bio_interpret = TRUE,
#            dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
#            comp_group = c("W03", "W30"),
#            p_cutoff = 0.01,
#            correct_p = TRUE,
#            species = 'dme',
#            quanti_pathway_method = 'mean')
#
#
# runMetDNA2(
#   path_pos = "/home/zhouzw/Data_processing/20210522_debug/7aaebfdbafb2393074745fde658a4096/POS",
#   path_neg = "/home/zhouzw/Data_processing/20210522_debug/7aaebfdbafb2393074745fde658a4096/NEG",
#   metdna_version = "version2",
#   polarity = "positive",
#   instrument = "SciexTripleTOF",
#   column = "hilic",
#   ce = "30",
#   method_lc = "Amide23min",
#   correct_p = TRUE,
#   extension_step = "0",
#   comp_group = c("W30", "W03"),
#   species = "dme",
#   p_cutoff = 0.050000,
#   fc_cutoff = 1.500000,
#   is_rt_calibration = TRUE)

################################################################################
# interpretBiology(sample_file_pos = 'data.csv',
#                  sample_info_file_pos = 'sample.info.csv',
#                  table_identification_pair_file_pos = 'table3_identification_pair.csv',
#                  sample_file_neg = 'data.csv',
#                  sample_info_file_neg = 'sample.info.csv',
#                  table_identification_pair_file_neg = 'table3_identification_pair.csv',
#                  path = '/home/zhouzw/Data_processing/20210522_debug/web_bug/907bd62c2668fe3e2a035ca9ab12a1c3/',
#                  polarity = 'both',
#                  metdna_version = 'version1',
#                  comp_group = c("W", "R"),
#                  uni_test = 't',
#                  correct_p = TRUE,
#                  p_cutoff = 0.01,
#                  fc_cutoff = 1,
#                  species = 'mmu',
#                  quanti_pathway_method = 'mean',
#                  organze_basis = 'kegg_id',
#                  extension_step = '0')
#
# interpretBiology(sample_file_pos = ms1_file_pos,
#                  sample_info_file_pos = sample_info_file_pos,
#                  table_identification_pair_file_pos = 'table3_identification_pair.csv',
#                  sample_file_neg = ms1_file_neg,
#                  sample_info_file_neg = sample_info_file_neg,
#                  table_identification_pair_file_neg = 'table3_identification_pair.csv',
#                  polarity = polarity,
#                  path = temp_path,
#                  metdna_version = metdna_version,
#                  comp_group = comp_group,
#                  uni_test = uni_test,
#                  correct_p = correct_p,
#                  p_cutoff = p_cutoff,
#                  fc_cutoff = fc_cutoff,
#                  species = species,
#                  quanti_pathway_method = quanti_pathway_method,
#                  organze_basis = 'kegg_id',
#                  extension_step = extension_step)
#
# interpretBiology(sample_file_pos = 'data.csv',
#                  sample_info_file_pos = 'sample.info.csv',
#                  table_identification_pair_file_pos = 'table3_identification_pair.csv',
#                  sample_file_neg = 'data.csv',
#                  sample_info_file_neg = 'sample.info.csv',
#                  table_identification_pair_file_neg = 'table3_identification_pair.csv',
#                  path = '/home/zhouzw/Data_processing/20210522_debug/web_bug/907bd62c2668fe3e2a035ca9ab12a1c3/',
#                  polarity = 'both',
#                  metdna_version = 'version1',
#                  comp_group = c("W03", "W30"),
#                  uni_test = 't',
#                  correct_p = TRUE,
#                  p_cutoff = 0.01,
#                  fc_cutoff = 1,
#                  species = 'mmu',
#                  quanti_pathway_method = 'mean',
#                  organze_basis = 'kegg_id',
#                  extension_step = '0')

################################################################################
#
# runMetDNA2(
#   path_pos = "/home/zhouzw/Data_processing/20210522_debug/web_bug/907bd62c2668fe3e2a035ca9ab12a1c3/POS",
#   path_neg = "/home/zhouzw/Data_processing/20210522_debug/web_bug/907bd62c2668fe3e2a035ca9ab12a1c3/NEG",
#   metdna_version = "version1",
#   polarity = "both",
#   instrument = "SciexTripleTOF",
#   column = "rp",
#   ce = "30",
#   method_lc = "Other",
#   correct_p = TRUE,
#   extension_step = "0",
#   comp_group = c("W", "R"),
#   species = "mmu",
#   p_cutoff = 0.010000,
#   fc_cutoff = 1.000000,
#   is_rt_calibration = FALSE)


# grep(pattern = "mgf", list.files(path = '/home/zhouzw/Data_processing/20210713_MetProcess_data/data_process/POS/', pattern = 'pattern=".+\\.mgf$', recursive = TRUE), value = TRUE)
#
#
# # # list.files(pattern="(?i)\\.mgf$", recursive = TRUE, full.names = TRUE)
# list.files(path = '/home/zhouzw/Data_processing/20210713_MetProcess_data/data_process/POS/',
#            pattern = '.+\\.mgf$',
#            recursive = TRUE)
#
# file.path('.', list.files(path = '/home/zhouzw/Data_processing/20210713_MetProcess_data/data_process/POS/',
#                           pattern = '.+\\.mgf$',
#                           recursive = TRUE,
#                           ignore.case = TRUE))

# list.files(path = '/home/zhouzw/Data_processing/20210713_MetProcess_data/data_process/POS/', recursive = TRUE)

################################################################################
# runMetDNA2(
#   path_pos = "/home/zhouzw/Data_processing/20210522_debug/01_metanalyzer_pos/",
#   path_neg = "/home/zhouzw/Data_processing/20210522_debug/7aaebfdbafb2393074745fde658a4096/NEG",
#   metdna_version = "version2",
#   polarity = "positive",
#   instrument = "SciexTripleTOF",
#   column = "hilic",
#   ce = "30",
#   method_lc = "Amide23min",
#   correct_p = TRUE,
#   extension_step = "0",
#   comp_group = c("W30", "W03"),
#   species = "dme",
#   p_cutoff = 0.050000,
#   fc_cutoff = 1.500000,
#   is_rt_calibration = TRUE)
#
# runMetDNA2(
#   path_pos = "/home/zhouzw/Data_processing/20210522_debug/03_msdial_pos/",
#   path_neg = "/home/zhouzw/Data_processing/20210522_debug/04_msdial_neg/",
#   metdna_version = "version2",
#   polarity = "negative",
#   instrument = "SciexTripleTOF",
#   column = "hilic",
#   ce = "30",
#   method_lc = "Amide23min",
#   correct_p = TRUE,
#   extension_step = "0",
#   comp_group = c("W30", "W03"),
#   species = "dme",
#   p_cutoff = 0.050000,
#   fc_cutoff = 1.500000,
#   is_rt_calibration = TRUE)

################################################################################
#
# annotateInitialSeed(ms1_file = "data.csv",
#                     ms2_file = NULL,
#                     metdna_version = 'version2',
#                     ms2_type = 'mgf',
#                     path = "/home/zhouzw/Data_processing/20210720_modify_MetDNA2_for_200STD_validation/POS_200std_mouse_liver/",
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
#                     test_adduct_version = 'version2',
#                     test_evaluation = '200STD',
#                     # direction = 'reverse',
#
#                     is_plot_ms2 <- TRUE
#                     )
#
#
# annotateMRN(annotation_result = "ms2_match_annotation_result.csv",
#             ms2_data = 'ms2',
#             prefer_adduct = "M+H",
#             metdna_version = 'version2',
#             lib = 'zhuMetLib',
#             direction = 'reverse',
#             column = 'hilic',
#             instrument = "SciexTripleTOF",
#             method_lc = 'Amide23min',
#             polarity = 'positive',
#             extension_step = '2',
#             threads = 3,
#             max_isotope = 4,
#             path = "/home/zhouzw/Data_processing/20210720_modify_MetDNA2_for_200STD_validation/POS_200std_mouse_liver/",
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
#             seed_neighbor_match_plot = TRUE,
#             test_adduct_version = 'version2',
#             test_evaluation = '200STD')


# annotation_result = "ms2_match_annotation_result.csv"
# ms2_data = 'ms2'
# prefer_adduct = "M+H"
# metdna_version = 'version2'
# lib = 'zhuMetLib'
# direction = 'reverse'
# column = 'hilic'
# instrument = "SciexTripleTOF"
# method_lc = 'Amide23min'
# polarity = 'positive'
# extension_step = '2'
# threads = 3
# max_isotope = 4
# path = "/home/zhouzw/Data_processing/20210720_modify_MetDNA2_for_200STD_validation/POS_200std_mouse_liver/"
# mz_tol = 25
# rt_tol1 = 3
# rt_tol2 = 30
# ccs_tol = 4
# cor_tol = 0
# int_tol = 500
# dp_tol = 0.5
# max_step = 3
# score_cutoff = 0
# remain = FALSE
# remain_per = 0.5
# scoring_approach = 'dp'
# matched_frag_tol = 1
# seed_neighbor_match_plot = TRUE
# test_adduct_version = 'version2'
# test_evaluation = '200STD'
#
# save(new_tags, file = '/home/zhouzw/Data_processing/20210720_modify_MetDNA2_for_200STD_validation/POS_200std_mouse_liver/temp_new_tags_210721.RData')

# load('/home/zhouzw/Data_processing/20210720_modify_MetDNA2_for_200STD_validation/POS_200std_mouse_liver/02_result_MRN_annotation/00_intermediate_data/tags2_after_redundancy_remove')
# load('/home/zhouzw/Data_processing/20210720_modify_MetDNA2_for_200STD_validation/POS_200std_mouse_liver/02_result_MRN_annotation/00_intermediate_data/id_result_redun_rm')
# metabolite
# temp <- getAnnotationResult(tags2 = tags2_after_redundancy_remove,
#                             tags.result2 = id_result_redun_rm,
#                             metabolite = metabolite,
#                             candidate.num = 5,
#                             instrument = 'SciexTripleTOF',
#                             score.cutoff = 0.4)


# temp <- mergeRecursiveAnnotation(path = path,
#                                  thread = thread,
#                                  lib = lib,
#                                  direction = direction,
#                                  tolerance_rt_range = tolerance_rt_range,
#                                  use_redun_rm_result = use_redun_rm_result)
#
# temp <- mergeRecursiveAnnotation(path = path,
#                                  thread = 3,
#                                  lib = 'zhuMetLib',
#                                  direction = 'reverse',
#                                  tolerance_rt_range = 30,
#                                  use_redun_rm_result = TRUE,
#                                  test_evaluation = '200STD')
#
# load('/home/zhouzw/Data_processing/20210720_modify_MetDNA2_for_200STD_validation/POS_200std_mouse_liver/03_annotation_credential/00_intermediate_data/id_merge_before_credential.RData')


# temp <- mergeRecursiveAnnotation(path = "/home/zhouzw/Data_processing/20210720_modify_MetDNA2_for_200STD_validation/POS_200std_mouse_liver/",
#                                  thread = 3,
#                                  lib = 'zhuMetLib',
#                                  direction = 'reverse',
#                                  tolerance_rt_range = 30,
#                                  use_redun_rm_result = TRUE,
#                                  test_evaluation = '200STD')
#
# rm('temp');gc()
#
# cat('Change formats for annotation credential.\n')
# convertAnnotationTable2InitialId(peak_table_file = "data.csv",
#                                  sample_info_file =  "sample.info.csv",
#                                  path_dir = "/home/zhouzw/Data_processing/20210720_modify_MetDNA2_for_200STD_validation/POS_200std_mouse_liver/",
#                                  polarity = 'positive',
#                                  # test_assi_confidence_initial_seed = test_assi_confidence_initial_seed,
#                                  tool = 'MetDNA2')

# authenticateAnnotationCredential(peak_table_file = "peak_table.csv",
#                                  annotation_initial_file = "annotation_initial.csv",
#                                  ms2_data_file = "ms2_data.RData",
#                                  path_dir = "/home/zhouzw/Data_processing/20210720_modify_MetDNA2_for_200STD_validation/POS_200std_mouse_liver/",
#                                  polarity = 'positive',
#                                  thread = 3,
#                                  isotope_int_ratio_check = TRUE, # para annotateIsotope
#                                  isotope_int_ratio_cutoff = 500,
#                                  is_ms2_check = TRUE, # compare ms2 similarity of neutral loss
#                                  ms2_score_cutoff = -1, # -1 represent not filter
#                                  is_plot_pseudo_MS1 = TRUE,
#                                  # formula prediction
#                                  dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
#                                  is_pred_formula_all = FALSE,
#                                  ppm = 25, # ms1 tolerance
#                                  acc = 25, # ms2 tolerance
#                                  elements = 'CHONPS',
#                                  num_formula_candidate = 3, # return top 3 candidate
#                                  platform = 'linux')
#
#
# generateMetDNA2AnnotationResult(path = "/home/zhouzw/Data_processing/20210720_modify_MetDNA2_for_200STD_validation/POS_200std_mouse_liver/",
#                                 thread = 3,
#                                 is_cred_pg_filter = TRUE,
#                                 is_cred_formula_filter = FALSE,
#                                 is_pred_formula_all = FALSE,
#                                 tolerance_rt_range = 30,
#                                 direction = 'reverse',
#                                 instrument = 'SciexTripleTOF',
#                                 mz_tol = 25,
#                                 rt_tol = 30,
#                                 ccs_tol = 3,
#                                 candidate_num = 5,
#                                 test_evaluation = '200STD')

#
#
# generateMetDNA2AnnotationResult(path = path,
#                                 thread = thread,
#                                 is_cred_pg_filter = is_cred_pg_filter,
#                                 is_cred_formula_filter = is_cred_formula_filter,
#                                 is_pred_formula_all = is_pred_formula_all,
#                                 tolerance_rt_range = tolerance_rt_range,
#                                 direction = direction,
#                                 instrument = instrument,
#                                 mz_tol = mz_tol,
#                                 rt_tol = rt_tol2,
#                                 ccs_tol = ccs_tol,
#                                 candidate_num = candidate_num,
#                                 test_evaluation = test_evaluation)

# MetDNA2(ms1_file = 'data.csv',
#         sample_info_file = "sample.info.csv",
#         metdna_version = 'version2',
#         ms2_type = 'mgf',
#         path = '/home/zhouzw/Data_processing/20210720_modify_MetDNA2_for_200STD_validation/POS',
#         thread = 3,
#         polarity = 'positive',
#         instrument = "SciexTripleTOF",
#         column = "hilic",
#         lib = 'zhuMetLib',
#         ce = "30",
#         method_lc = 'Amide23min',
#         extension_step = '2',
#         scoring_approach_recursive = 'dp',
#         is_check_data = TRUE,
#         is_anno_initial_seed = TRUE,
#         is_anno_mrn = TRUE,
#         is_credential = TRUE,
#         is_bio_interpret = FALSE,
#         test_evaluation = '200STD')


# test <- readMGF('/home/zhouzw/Data_processing/20210522_debug/Data_from_Youyou/s6_1.mgf')

# test <- readMSP(file = '/home/zhouzw/Data_processing/20210522_debug/20210721/spectra_ext0.msp', mode = 'all')
# test <- readMSP(file = '/home/zhouzw/Data_processing/20210522_debug/20210721/aging_fly_pos.msp', mode = 'all')

# cpd_zhulab <- zh
# idx <- match()
# lib_rt$zhumetlib_exp

# cpd_zhulab <- zhuMetlib$meta$compound
# temp_inchikey <- match(lib_rt$zhumetlib_exp$id, cpd_zhulab$id) %>% cpd_zhulab$inchikey[.]
# temp_inchikey1 <- match(lib_rt$zhumetlib_exp$id, cpd_zhulab$id) %>% cpd_zhulab$inchikey1[.]
#
# temp_lib_rt <- lib_rt$zhumetlib_exp %>%
#   mutate(new_inchikey = temp_inchikey, new_inchikey1 = temp_inchikey1) %>%
#   dplyr::mutate(same_inchikey = (inchikey == new_inchikey),
#                 same_inchikey1 = (inchikey1 == new_inchikey1))
#
# temp_lib_rt %>% dplyr::filter(!same_inchikey) %>% View()
#
#
# MetDNA2(ms1_file = 'data.csv',
#         sample_info_file = "sample.info.csv",
#         metdna_version = 'version2',
#         ms2_type = 'mgf',
#         path = '/home/zhouzw/Data_processing/20210720_modify_MetDNA2_for_200STD_validation/POS',
#         thread = 3,
#         polarity = 'positive',
#         instrument = "SciexTripleTOF",
#         column = "hilic",
#         lib = 'zhuMetLib',
#         ce = "30",
#         method_lc = 'Amide23min',
#         extension_step = '2',
#         scoring_approach_recursive = 'dp',
#         is_check_data = TRUE,
#         is_anno_initial_seed = TRUE,
#         is_anno_mrn = TRUE,
#         is_credential = TRUE,
#         is_bio_interpret = FALSE,
#         test_evaluation = '200STD')
#
# MetDNA2(ms1_file = 'data.csv',
#         sample_info_file = "sample.info.csv",
#         metdna_version = 'version2',
#         ms2_type = 'mgf',
#         path = '/home/zhouzw/Data_processing/20210720_modify_MetDNA2_for_200STD_validation/POS_with_formula/',
#         thread = 3,
#         polarity = 'positive',
#         instrument = "SciexTripleTOF",
#         column = "hilic",
#         lib = 'zhuMetLib',
#         ce = "30",
#         method_lc = 'Amide23min',
#         extension_step = '2',
#         scoring_approach_recursive = 'dp',
#         is_check_data = TRUE,
#         is_anno_initial_seed = TRUE,
#         is_anno_mrn = TRUE,
#         is_credential = TRUE,
#         is_bio_interpret = FALSE,
#         is_cred_formula_filter = TRUE,
#         test_evaluation = '200STD')


################################################################################
# 20210724 ---------------------------------------------------------------------

#   bug reported by LuoMingdu
# interpretBiology(sample_file_neg = "data.csv",
#                  sample_info_file_neg = "sample.info.csv",
#                  table_identification_pair_file_neg = 'table3_identification_pair.csv',
#                  polarity = 'negative',
#                  path = '/home/zhouzw/Data_processing/20210522_debug/06_20210723_200std_in_mice_liver_neg/',
#                  metdna_version = 'version2',
#                  comp_group = c('a', 'b'),
#                  uni_test = 't',
#                  correct_p = TRUE,
#                  p_cutoff = 0.05,
#                  fc_cutoff = 1,
#                  species = 'mmu',
#                  quanti_pathway_method = 'mean',
#                  organze_basis = 'kegg_id',
#                  extension_step = '0')


#   bug reported by web
# interpretBiology(sample_file_neg = "data.csv",
#                  sample_info_file_neg = "sample.info.csv",
#                  table_identification_pair_file_neg = 'table3_identification_pair.csv',
#                  polarity = 'positive',
#                  path = '/home/zhouzw/Data_processing/20210522_debug/45a44036a8e777bb37b9a9fb3237056e/POS/',
#                  metdna_version = 'version2',
#                  comp_group = c('1', '2'),
#                  uni_test = 't',
#                  correct_p = TRUE,
#                  p_cutoff = 0.05,
#                  fc_cutoff = 1,
#                  species = 'ath',
#                  quanti_pathway_method = 'mean',
#                  organze_basis = 'kegg_id',
#                  extension_step = '0')
#
#
# interpretBiology(sample_file_neg = "data.csv",
#                  sample_info_file_neg = "sample.info.csv",
#                  table_identification_pair_file_neg = 'table3_identification_pair.csv',
#                  polarity = 'positive',
#                  path = '/home/zhouzw/Data_processing/20210522_debug/660266eeac841ccf2fc65d62faeac965/POS/',
#                  metdna_version = 'version2',
#                  comp_group = c('A','B'),
#                  uni_test = 't',
#                  correct_p = TRUE,
#                  p_cutoff = 0.05,
#                  fc_cutoff = 1,
#                  species = 'hsa',
#                  quanti_pathway_method = 'mean',
#                  organze_basis = 'kegg_id',
#                  extension_step = '0')
#
#
# interpretBiology(sample_file_neg = "data.csv",
#                  sample_info_file_neg = "sample.info.csv",
#                  table_identification_pair_file_neg = 'table3_identification_pair.csv',
#                  polarity = 'negative',
#                  path = '/home/zhouzw/Data_processing/20210522_debug/54fc634bbf6464d141a6d8bb14603577/NEG/',
#                  metdna_version = 'version2',
#                  comp_group = c('NO2', 'DO2'),
#                  uni_test = 't',
#                  correct_p = FALSE,
#                  p_cutoff = 0.05,
#                  fc_cutoff = 1,
#                  species = 'hsa',
#                  quanti_pathway_method = 'mean',
#                  organze_basis = 'kegg_id',
#                  extension_step = '0')

# path <- '/home/zhouzw/Data_processing/20210522_debug/54fc634bbf6464d141a6d8bb14603577/NEG/04_biology_intepretation/'
#
# temp <- try(MetBioInterpretation::quantiPathway(quanti_table = 'pathway_metabolite_quantitative_result.csv',
#                                                 sample_info = 'sample.info.csv',
#                                                 comp_group = c('NO2', 'DO2'),
#                                                 is_scale = FALSE, # peak area has been scaled in selectQuantiMetabolite
#                                                 correct_p = FALSE,
#                                                 uni_test = 't',
#                                                 path = '/home/zhouzw/Data_processing/20210522_debug/54fc634bbf6464d141a6d8bb14603577/NEG/04_biology_intepretation/',
#                                                 by_what = 'mean',
#                                                 species = 'hsa',
#                                                 is_consider_stereo_isomer = FALSE),
#             silent = TRUE)


################################################################################
# load('/home/zhouzw/Data_processing/20210522_debug/01_20210723_lps_cell_pos/00_annotation_table/00_intermediate_data/annot_all')

# generateExportTable(annot_all = annot_all,
#                     direction = 'reverse',
#                     tolerance_rt_range = 30,
#                     lib = 'zhuMetLib',
#                     path = '/home/zhouzw/Data_processing/20210522_debug/01_20210723_lps_cell_pos/',
#                     instrument = 'SciexTripleTOF',
#                     is_cred_formula_filter = FALSE,
#                     candidate_num = 5,
#                     mz_tol = 25,
#                     rt_tol = 30,
#                     ccs_tol = 3,
#                     test_evaluation = 'No')
#
# interpretBiology(sample_file_pos = "data.csv",
#                  sample_info_file_pos = "sample.info.csv",
#                  table_identification_pair_file_pos = 'table3_identification_pair.csv',
#                  polarity = 'positive',
#                  path = '/home/zhouzw/Data_processing/20210522_debug/01_20210723_lps_cell_pos/',
#                  metdna_version = 'version2',
#                  comp_group = c('lps', 'control'),
#                  uni_test = 't',
#                  correct_p = TRUE,
#                  p_cutoff = 0.05,
#                  fc_cutoff = 1,
#                  species = 'mmu',
#                  quanti_pathway_method = 'mean',
#                  organze_basis = 'kegg_id',
#                  extension_step = '0')
# load('/home/zhouzw/Data_processing/20210522_debug/01_20210723_lps_cell_pos/04_biology_intepretation/00_intermediate_data/result_stat.RData')
# table_class_enrich <- MetBioInterpretation::generateClassEnrichTable(table_identification_pair = 'table3_identification_pair.csv',
#                                                                      result_stat = result_stat)
# load('/home/zhouzw/Data_processing/20210522_debug/01_20210723_lps_cell_pos/test.Rdata')
# path_output <- '/home/zhouzw/Data_processing/20210522_debug/01_20210723_lps_cell_pos/04_biology_intepretation'

#################################################################################

# interpretBiology(sample_file_neg = "data.csv",
#                  sample_info_file_neg = "sample.info.csv",
#                  table_identification_pair_file_neg = 'table3_identification_pair.csv',
#                  polarity = 'negative',
#                  path = '/home/zhouzw/Data_processing/20210522_debug/06_20210723_200std_in_mice_liver_neg/',
#                  metdna_version = 'version2',
#                  comp_group = c('a', 'b'),
#                  uni_test = 't',
#                  correct_p = TRUE,
#                  p_cutoff = 0.05,
#                  fc_cutoff = 1,
#                  species = 'mmu',
#                  quanti_pathway_method = 'mean',
#                  organze_basis = 'kegg_id',
#                  extension_step = '0')
#
# runMetDNA2(sample_file_neg = "data.csv",
#            sample_info_file_neg = "sample.info.csv",
#            table_identification_pair_file_neg = 'table3_identification_pair.csv',
#            polarity = 'negative',
#            path = '/home/zhouzw/Data_processing/20210522_debug/06_20210723_200std_in_mice_liver_neg/',
#            metdna_version = 'version2',
#            comp_group = c('a', 'b'), )


# runMetDNA2(
#   path_neg = '/home/zhouzw/Data_processing/20210522_debug/06_20210723_200std_in_mice_liver_neg/',
#   metdna_version = "version2",
#   polarity = "negative",
#   instrument = "SciexTripleTOF",
#   column = "hilic",
#   ce = "30",
#   dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
#   method_lc = "Amide12min",
#   correct_p = TRUE,
#   extension_step = "0",
#   comp_group = c('a', 'b'),
#   species = "mmu",
  # is_check_data = FALSE,
  # is_anno_initial_seed = FALSE,
  # is_anno_mrn = FALSE,
  # is_credential = FALSE,
  # is_webserver = FALSE,
  # is_rm_intermediate_data = FALSE,
#   p_cutoff = 0.05,
#   fc_cutoff = 1.000000,
#   is_rt_calibration = TRUE)


################################################################################
# 20210726 ---------------------------------------------------------------------
#
# runMetDNA2(ms1_file_pos = "data.csv",
#            sample_info_file_pos = "sample.info.csv",
#            path_pos = '/home/zhouzw/Data_processing/20210522_debug/POS/',
#            metdna_version = 'version2',
#            is_webserver = FALSE,
#            is_rm_intermediate_data = FALSE,
#            ms2_type = 'mgf',
#            thread = 4,
#            lib = 'zhuMetLib',
#            polarity = 'positive',
#            instrument = "SciexTripleTOF",
#            column = "hilic",
#            ce = "30",
#            method_lc = 'Amide23min',
#            is_rt_calibration = TRUE,
#            extension_step = '0',
#            dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
#            comp_group = c("pos_3_DDA", "QC"),
#            species = 'mmu',
#            quanti_pathway_method = 'mean',
#            p_cutoff = 0.05,
#            correct_p = TRUE,
#            fc_cutoff = 1,
#            is_check_data = TRUE,
#            is_anno_initial_seed = FALSE,
#            is_anno_mrn = FALSE,
#            is_credential = FALSE)

# interpretBiology(sample_file_pos = "data.csv",
#                  sample_info_file_pos = "sample.info.csv",
#                  table_identification_pair_file_pos = 'table3_identification_pair.csv',
#                  polarity = 'positive',
#                  path = '/home/zhouzw/Data_processing/20210522_debug/POS/',
#                  metdna_version = 'version2',
#                  comp_group = c('a', 'b'),
#                  uni_test = 't',
#                  correct_p = TRUE,
#                  p_cutoff = 0.05,
#                  fc_cutoff = 1,
#                  species = 'mmu',
#                  quanti_pathway_method = 'mean',
#                  organze_basis = 'kegg_id',
#                  extension_step = '0')

# path_output <- '/home/zhouzw/Data_processing/20210522_debug/POS/04_biology_intepretation/'
# sample <- readr::read_csv(file.path(path_output, 'data.csv'))
# sample_info <- readr::read_csv(file.path(path_output, 'sample.info.csv'), col_types = 'cc')
# table_identification_pair <- readr::read_csv(file.path(path_output, 'table3_identification_pair.csv'))
#
# runMetDNA2(ms1_file_pos = "data.csv",
#            sample_info_file_pos = "sample.info.csv",
#            path_pos ='/home/zhouzw/Data_processing/20210522_debug/05_20210723_200std_in_mice_liver_pos/',
#            metdna_version = 'version2',
#            is_webserver = FALSE,
#            is_rm_intermediate_data = FALSE,
#            ms2_type = 'mgf',
#            thread = 4,
#            lib = 'zhuMetLib',
#            polarity = 'positive',
#            instrument = "SciexTripleTOF",
#            column = "hilic",
#            ce = "30",
#            method_lc = 'Amide23min',
#            is_rt_calibration = TRUE,
#            extension_step = '0',
#            dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
#            comp_group = c("a", "b"),
#            species = 'mmu',
#            quanti_pathway_method = 'mean',
#            p_cutoff = 0.05,
#            correct_p = TRUE,
#            fc_cutoff = 1,
#            is_check_data = TRUE,
#            is_anno_initial_seed = FALSE,
#            is_anno_mrn = FALSE,
#            is_credential = FALSE
#            )

################################################################################
# library(MetDNA2)
# runMetDNA2(ms1_file_pos = "data.csv",
#            sample_info_file_pos = "sample.info.csv",
#            path_pos = '/home/zhouzw/Data_processing/20210522_debug/POS/',
#            metdna_version = 'version2',
#            is_webserver = FALSE,
#            is_rm_intermediate_data = FALSE,
#            ms2_type = 'mgf',
#            thread = 4,
#            lib = 'zhuMetLib',
#            polarity = 'positive',
#            instrument = "SciexTripleTOF",
#            column = "hilic",
#            ce = "30",
#            method_lc = 'Amide23min',
#            is_rt_calibration = TRUE,
#            extension_step = '0',
#            dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
#            comp_group = c("pos_3_DDA", "QC"),
#            species = 'hsa',
#            quanti_pathway_method = 'mean',
#            p_cutoff = 0.05,
#            correct_p = TRUE,
#            fc_cutoff = 1,
#            is_check_data = TRUE,
#            is_anno_initial_seed = FALSE,
#            is_anno_mrn = FALSE,
#            is_credential = FALSE)

################################################################################
# 20210727
#
# runMetDNA2(ms1_file_neg = "data.csv",
#            sample_info_file_neg = "sample.info.csv",
#            path_neg = "/home/zhouzw/Data_processing/20210522_debug/10_nist_urine_neg_emrn0_without_mz_segment/",
#            metdna_version = 'version2',
#            is_webserver = FALSE,
#            is_rm_intermediate_data = FALSE,
#            ms2_type = 'mgf',
#            thread = 3,
#            lib = 'zhuMetLib',
#            polarity = 'negative',
#            instrument = "SciexTripleTOF",
#            column = "hilic",
#            ce = "30",
#            method_lc = 'Amide23min',
#            is_rt_calibration = TRUE,
#            extension_step = '0',
#            dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
#            comp_group = c("neg_only_full_DDA", "QC"),
#            species = 'hsa',
#            quanti_pathway_method = 'mean',
#            p_cutoff = 0.05,
#            correct_p = TRUE,
#            fc_cutoff = 1,
#            is_check_data = TRUE,
#            is_anno_initial_seed = FALSE,
#            is_anno_mrn = FALSE,
#            is_credential = FALSE)
# interpretBiology(sample_file_neg = "data.csv",
#                  sample_info_file_neg = "sample.info.csv",
#                  table_identification_pair_file_neg = 'table3_identification_pair.csv',
#                  polarity = 'negative',
#                  path = '/home/zhouzw/Data_processing/20210522_debug/10_nist_urine_neg_emrn0_without_mz_segment/',
#                  metdna_version = 'version2',
#                  comp_group = c("neg_only_full_DDA", "QC"),
#                  uni_test = 't',
#                  correct_p = TRUE,
#                  p_cutoff = 0.05,
#                  fc_cutoff = 1,
#                  species = 'hsa',
#                  quanti_pathway_method = 'mean',
#                  organze_basis = 'kegg_id',
#                  extension_step = '0')

################################################################################
# load('/home/zhouzw/Data_processing/20210522_debug/10_nist_urine_neg_emrn0_without_mz_segment/00_annotation_table/00_intermediate_data/id_merge')
# load('/home/zhouzw/Data_processing/20210522_debug/10_nist_urine_neg_emrn0_without_mz_segment/00_annotation_table/00_intermediate_data/table_identification')
# load('/home/zhouzw/Data_processing/20210522_debug/10_nist_urine_neg_emrn0_without_mz_segment/01_result_initial_seed_annotation/00_intermediate_data/ms1_data')


################################################################################
# MetDNA2(ms1_file = 'data.csv',
#         sample_info_file = "sample.info.csv",
#         metdna_version = 'version2',
#         ms2_type = 'mgf',
#         path = '/home/zhouzw/Data_processing/20210522_debug/20210823_46std_test/',
#         thread = 3,
#         polarity = 'positive',
#         instrument = "SciexTripleTOF",
#         column = "hilic",
#         lib = 'zhuMetLib',
#         ce = "30",
#         method_lc = 'Amide23min',
#         extension_step = '2',
#         scoring_approach_recursive = 'dp',
#         is_check_data = TRUE,
#         is_anno_initial_seed = TRUE,
#         is_anno_mrn = TRUE,
#         is_credential = TRUE,
#         is_bio_interpret = FALSE,
#         is_cred_formula_filter = TRUE,
#         test_evaluation = '46STD')

################################################################################
#
# MetDNA2(ms1_file = 'data.csv',
#         sample_info_file = "sample.info.csv",
#         metdna_version = 'version2',
#         ms2_type = 'mgf',
#         path = '/home/zhouzw/Data_processing/20210830_s9fraction_incubation_on_46STDs/peak_area_filterd/pos/',
#         thread = 3,
#         polarity = 'positive',
#         instrument = "SciexTripleTOF",
#         column = "hilic",
#         lib = 'zhuMetLib',
#         ce = "30",
#         method_lc = 'Amide23min',
#         extension_step = '2',
#         scoring_approach_recursive = 'dp',
#         is_check_data = TRUE,
#         is_anno_initial_seed = TRUE,
#         is_anno_mrn = TRUE,
#         is_credential = TRUE,
#         is_bio_interpret = FALSE,
#         is_cred_formula_filter = TRUE,
#         test_evaluation = '46STD')

################################################################################
# debug for Yiyang Zhou 20211001 -----------------------------------------------

# path_dir <- '/home/zhouzw/Data_processing/01_debug/211001_zhouyiyang/1bd80c84979b3b2f1ef54eda9fc6afca/POS'
# load('/home/zhouzw/Data_processing/01_debug/211001_zhouyiyang/1bd80c84979b3b2f1ef54eda9fc6afca/POS/03_annotation_credential/00_intermediate_data/annot_all_before_credential.RData')
# load('/home/zhouzw/Data_processing/01_debug/211001_zhouyiyang/1bd80c84979b3b2f1ef54eda9fc6afca/POS/03_annotation_credential/00_intermediate_data/list_peak_group_annotation.RData')
# load('/home/zhouzw/Data_processing/01_debug/211001_zhouyiyang/1bd80c84979b3b2f1ef54eda9fc6afca/POS/03_annotation_credential/00_intermediate_data/record_conflict_peak_group.RData')
# record_conflict_peak_group

#
# runMetDNA2(
#   path_pos = "/home/zhouzw/Data_processing/01_debug/211001_zhouyiyang/1bd80c84979b3b2f1ef54eda9fc6afca/POS",
#   path_neg = "/home/zhouzw/Data_processing/01_debug/211001_zhouyiyang/1bd80c84979b3b2f1ef54eda9fc6afca/NEG",
#   metdna_version = "version2",
#   polarity = "positive",
#   instrument = "SciexTripleTOF",
#   column = "rp",
#   ce = "35,15",
#   method_lc = "Other",
#   correct_p = FALSE,
#   dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
#   extension_step = "0",
#   comp_group = c("w1", "w2"),
#   species = "hsa",
#   p_cutoff = 0.050000,
#   fc_cutoff = 1.000000,
#   is_rt_calibration = FALSE,
#   is_anno_initial_seed = FALSE,
#   is_anno_mrn = FALSE)

################################################################################
# debug for Junmiao 20211015 ---------------------------------------------------
#
# runMetDNA2(
#   path_pos = "/home/zhouzw/Data_processing/01_debug/211015_junmiao/d44d3cfb00c552a2275a17e5fd33243f/POS",
#   path_neg = "/home/zhouzw/Data_processing/01_debug/211015_junmiao/d44d3cfb00c552a2275a17e5fd33243f/NEG",
#   metdna_version = "version2",
#   polarity = "positive",
#   instrument = "SciexTripleTOF",
#   column = "rp",
#   ce = "35,15",
#   method_lc = "Other",
#   correct_p = FALSE,
#   dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
#   extension_step = "0",
#   comp_group = c("46", "NC"),
#   species = "rno",
#   p_cutoff = 0.050000,
#   fc_cutoff = 1.000000,
#   is_rt_calibration = FALSE,
#   is_anno_initial_seed = FALSE,
#   is_anno_mrn = FALSE)

# load('./')
# generateExportTable(annot_all = annot_all,
#                     direction = 'reverse',
#                     tolerance_rt_range = 30,
#                     lib = 'zhuMetLib',
#                     path = path,
#                     instrument = 'SciexTripleTOF',
#                     is_cred_formula_filter = FALSE,
#                     candidate_num = 5,
#                     mz_tol = 25,
#                     rt_tol = 30,
#                     ccs_tol = 3)

# peak_group_id_table$included_peak_group %>% stringr::str_extract(pattern = '(\\[\\d*M.+)|\\[M\\+|\\[M\\-')
# peak_group_id_table$included_peak_group %>% stringr::str_extract(pattern = '\\[\\d*M.+|\\[M\\+|\\[M\\-')

#
# runMetDNA2(
#   path_pos = "/home/zhouzw/Data_processing/01_debug/211015_junmiao/d44d3cfb00c552a2275a17e5fd33243f/POS",
#   path_neg = "/home/zhouzw/Data_processing/01_debug/211015_junmiao/d44d3cfb00c552a2275a17e5fd33243f/NEG",
#   metdna_version = "version2",
#   polarity = "positive",
#   instrument = "SciexTripleTOF",
#   column = "rp",
#   ce = "35,15",
#   method_lc = "Other",
#   correct_p = FALSE,
#   dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
#   extension_step = "0",
#   comp_group = c("46", "NC"),
#   species = "rno",
#   p_cutoff = 0.050000,
#   fc_cutoff = 1.000000,
#   is_check_data = TRUE,
#   is_rt_calibration = FALSE,
#   is_anno_initial_seed = FALSE,
#   is_anno_mrn = FALSE)

# checkQuality(ms1_file = "data.csv",
#              sample_info_file = "sample.info.csv",
#              ms2_type = "mgf",
#              path = "/home/zhouzw/Data_processing/01_debug/211015_junmiao/d44d3cfb00c552a2275a17e5fd33243f/POS")

################################################################################
# debug for Jiaqiang 20211015 ---------------------------------------------------

# runMetDNA2(
#   path_pos = "/home/zhouzw/Data_processing/01_debug/211015_jiaqiang/4fe1afad947fba72644b8c8e8067d0cb/POS",
#   path_neg = "/home/zhouzw/Data_processing/01_debug/211015_jiaqiang/4fe1afad947fba72644b8c8e8067d0cb/NEG/",
#   metdna_version = "version2",
#   polarity = "both",
#   instrument = "ThermoOrbitrap",
#   column = "hilic",
#   ce = "30",
#   method_lc = "Other",
#   correct_p = TRUE,
#   dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
#   extension_step = "0",
#   comp_group = c("HT", "PHT"),
#   species = "rno",
#   p_cutoff = 0.050000,
#   fc_cutoff = 1.000000,
#   is_rt_calibration = FALSE,
#   is_anno_initial_seed = FALSE,
#   is_anno_mrn = FALSE)


# path_dir <- "/home/zhouzw/Data_processing/01_debug/211015_jiaqiang/4fe1afad947fba72644b8c8e8067d0cb/POS"
# load(file.path(path_dir,
#                '03_annotation_credential',
#                "00_intermediate_data",
#                'annotation_initial.RData'))

################################################################################
# debug for CaiYuping ----------------------------------------------------------
# runMetDNA2(ms1_file_pos = "data.csv",
#            ms1_file_neg = "data.csv",
#            sample_info_file_pos = "sample.info.csv",
#            sample_info_file_neg = "sample.info.csv",
#            path_pos = "/home/zhouzw/Data_processing/01_debug/211018_YupingCai/RPLC_MetDNA2_posneg/POS",
#            path_neg = "/home/zhouzw/Data_processing/01_debug/211018_YupingCai/RPLC_MetDNA2_posneg/NEG",
#            metdna_version = 'version2',
#            is_webserver = FALSE,
#            is_rm_intermediate_data = FALSE,
#            ms2_type = 'mgf',
#            thread = 2,
#            lib = 'zhuMetLib_orbitrap',
#            polarity = 'both',
#            instrument = "ThermoExploris",
#            column = "rp",
#            ce = "SNCE20_30_40%",
#            method_lc = 'Other',
#            is_rt_calibration = FALSE,
#            extension_step = '0',
#            dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
#            comp_group = c("QC", "Subject"),
#            species = 'hsa',
#            quanti_pathway_method = 'mean',
#            p_cutoff = 0.05,
#            correct_p = TRUE,
#            fc_cutoff = 1,
#            is_plot_ms2 = FALSE,
#            seed_neighbor_match_plot = FALSE,
#            is_anno_initial_seed = FALSE,
#            is_anno_mrn = FALSE)


# peak_table_file <- 'peak_table.csv'
# annotation_initial_file <- 'annotation_initial.csv'
# ms2_data_file <- 'ms2_data.RData'
# path_dir <- '/home/zhouzw/Data_processing/01_debug/211018_YupingCai/RPLC_MetDNA2_posneg/POS/'
# polarity <- 'positive'
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
# ppm = 25 # ms1 tolerance
# acc = 25 # ms2 tolerance
# elements = "CHNOPS"
# num_formula_candidate = 3 # return top 3 candidate
# is_pred_formula_all = FALSE
# platform = 'linux'
# type_order = c('level1', 'level2', 'level3')
# tol_rt = 3


# load(file.path(path_dir, '03_annotation_credential', "00_intermediate_data", 'list_peak_group.RData'))
#
# ms2_data = raw_msms

# tol_mz = 10
# isotope_delta = 1.003355
# isotope_max_num = 4
# isotope_int_ratio_check = TRUE
# isotope_int_ratio_cutoff = 500
# monotonic_dec_check = FALSE
# monotonic_mz_cutoff = 800
# cutoff_ssc = NULL
# cutoff_ssc_int = 3000

# checkQuality(ms1_file = "data.csv",
#              sample_info_file = "sample.info.csv",
#              ms2_type = "mgf",
#              path = "/home/zhouzw/Data_processing/01_debug/211018_YupingCai/RPLC_MetDNA2_posneg/POS")

################################################################################
# debug for Junmiao ------------------------------------------------------------
# runMetDNA2(
#   path_pos = "/home/zhouzw/Data_processing/01_debug/211019_junmiao/2b5953e68a90b24857de39b1d0623966/POS/",
#   path_neg = "/home/zhouzw/Data_processing/01_debug/211019_junmiao/2b5953e68a90b24857de39b1d0623966/NEG/",
#   metdna_version = "version1",
#   polarity = "negative",
#   instrument = "SciexTripleTOF",
#   column = "rp",
#   ce = "35,15",
#   method_lc = "Other",
#   correct_p = FALSE,
#   dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
#   extension_step = "0",
#   comp_group = c("OA", "3A"),
#   species = "rno",
#   p_cutoff = 0.050000,
#   fc_cutoff = 1.000000,
#   is_rt_calibration = FALSE,
#   is_anno_initial_seed = FALSE,
#   is_anno_mrn = FALSE,
#   is_credential = FALSE)
#
# interpretBiology(sample_file_neg = 'data.csv',
#                  sample_info_file_neg = 'sample.info.csv',
#                  table_identification_pair_file_neg = 'table3_identification_pair.csv',
#                  polarity = 'negative',
#                  path = "/home/zhouzw/Data_processing/01_debug/211019_junmiao/2b5953e68a90b24857de39b1d0623966/NEG/",
#                  metdna_version = 'version1',
#                  comp_group = c("OA", "3A"),
#                  uni_test = 't',
#                  correct_p = FALSE,
#                  p_cutoff = 0.05,
#                  fc_cutoff = 1,
#                  species = 'rno',
#                  quanti_pathway_method = 'mean',
#                  organze_basis = 'kegg_id',
#                  extension_step = '0')

################################################################################
# load('/home/zhouzw/Data_processing/20211027_credential/list_peak_group.RData')
# load('/home/zhouzw/Data_processing/20211027_credential/ms2_data.RData')
#
# all_peak_groups <- pbapply::pblapply(list_peak_group, function(x){
#   MetDNA2::annotatePeakGroup(peak_group = x,
#                              ms2_data = raw_msms,
#                              polarity = 'positive',
#                              tol_mz = 25,
#                              is_ms2_check = FALSE,
#                              ms2_score_cutoff = 0.5,
#                              lib_adduct_nl = lib_adduct_nl,
#                              cutoff_ssc = 0.3,
#                              cutoff_ssc_int = 3000,
#                              is_rule_limitation = FALSE,
#                              cutoff_topN = 5)
#
# })
#
# table_annotation <- lapply(all_peak_groups, function(x){
#   x@peak_list_annotated
# }) %>% dplyr::bind_rows()
#
# table_annotation %>% dplyr::distinct(peak_name, .keep_all = TRUE)
#
# save(all_peak_groups, file = '/home/zhouzw/Data_processing/20211027_credential/all_peak_groups_211027.RData')

################################################################################
#
#
# interpretBiology(sample_file_pos = "data.csv",
#                  sample_info_file_pos = "sample.info.csv",
#                  table_identification_pair_file_pos = 'table3_identification_pair.csv',
#                  sample_file_neg = "data.csv",
#                  sample_info_file_neg = "sample.info.csv",
#                  table_identification_pair_file_neg = 'table3_identification_pair.csv',
#                  polarity = 'both',
#                  path = '/home/zhouzw/Data_processing/01_debug/211019_guangxueLiu/8ac819d1fabb516e84f6b724e6622b1f/',
#                  metdna_version = "version2",
#                  comp_group = c("control", "case"),
#                  uni_test = 't',
#                  correct_p = FALSE,
#                  p_cutoff = 0.05,
#                  fc_cutoff = 1,
#                  species = 'hsa',
#                  quanti_pathway_method = 'mean',
#                  organze_basis = 'kegg_id',
#                  extension_step = '0')
#
#
#
# runMetDNA2(
#   path_neg = '/home/zhouzw/Data_processing/20210522_debug/10_nist_urine_neg_emrn0_without_mz_segment',
#   metdna_version = "version2",
#   polarity = "negative",
#   instrument = "SciexTripleTOF",
#   column = "rp",
#   ce = "30",
#   method_lc = "Other",
#   correct_p = FALSE,
#   extension_step = "0",
#   comp_group = c("QC", "neg_only_full_DDA"),
#   species = "hsa",
#   p_cutoff = 0.050000,
#   fc_cutoff = 1.000000,
#   dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
#   is_rt_calibration = FALSE)


################################################################################
# 20211030 ZhangHaosong RP library
#
# runMetDNA2(
#   path_pos = '/home/zhouzw/Data_processing/20211028_zhulab_rp_library/test/',
#   metdna_version = "version2",
#   polarity = "positive",
#   instrument = "ThermoExploris",
#   lib = 'zhuRPLib',
#   column = "rp",
#   ce = "30",
#   method_lc = "zhulabRP",
#   correct_p = FALSE,
#   extension_step = "0",
#   comp_group = c("g1", "g2"),
#   species = "hsa",
#   p_cutoff = 0.050000,
#   fc_cutoff = 1.000000,
#   dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
#   is_bio_interpret = FALSE,
#   is_rt_calibration = TRUE,
#   is_webserver = FALSE)
#
#
# runMetDNA2(
#   path_neg = '/home/zhouzw/Data_processing/20211028_zhulab_rp_library/test2/',
#   metdna_version = "version2",
#   polarity = "negative",
#   instrument = "ThermoExploris",
#   lib = 'zhuRPLib',
#   column = "rp",
#   ce = "SNCE20_30_40%",
#   method_lc = "zhulabRP",
#   correct_p = FALSE,
#   extension_step = "0",
#   comp_group = c("LPS", "Ctrl"),
#   species = "hsa",
#   p_cutoff = 0.050000,
#   fc_cutoff = 1.000000,
#   dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
#   is_bio_interpret = FALSE,
#   is_rt_calibration = TRUE,
#   is_webserver = FALSE)
#
#
# authenticatePeakGroupCredential(list_peak_group = list_peak_group,
#                                 ms2_data = raw_msms,
#                                 path_dir = 'I:/00_projects/03_MetDNA2/00_data/20200805_targeted_annotation_evaluation_200STD/00_demo_data',
#                                 polarity = 'negative',
#                                 tol_mz = 10,
#                                 isotope_int_ratio_check = TRUE,
#                                 isotope_int_ratio_cutoff = 500,
#                                 is_ms2_check = TRUE,
#                                 ms2_score_cutoff = -1,
#                                 is_plot_pseudo_MS1 = TRUE)
#
# polarity = c('positive', 'negative'),
# tol_mz = 25,
# isotope_int_ratio_check = TRUE, # para annotateIsotope
# isotope_int_ratio_cutoff = 500,
# is_ms2_check = TRUE, # compare ms2 similarity of neutral loss
# ms2_score_cutoff = -1, # -1 represent not filter
# cutoff_ssc = 0.3,
# cutoff_ssc_int = 3000,
# is_rule_limitation = TRUE,
# cutoff_topN = 5,
# is_plot_pseudo_MS1 = TRUE,
# type_order = c('level1', 'level2', 'level3'),
#
# load('/home/zhouzw/Data_processing/20211028_zhulab_rp_library/test/03_annotation_credential/00_intermediate_data/list_peak_group.RData')
# load('/home/zhouzw/Data_processing/20211028_zhulab_rp_library/test/03_annotation_credential/ms2_data.RData')
# authenticatePeakGroupCredential(list_peak_group = list_peak_group,
#                                 ms2_data = raw_msms,
#                                 path_dir = '/home/zhouzw/Data_processing/20211028_zhulab_rp_library/test/',
#                                 thread = 2,
#                                 polarity = 'positive',
#                                 tol_mz = 25,
#                                 isotope_int_ratio_check = TRUE,
#                                 isotope_int_ratio_cutoff = 500,
#                                 is_ms2_check = TRUE,
#                                 ms2_score_cutoff = -1,
#                                 cutoff_ssc = 0.3,
#                                 cutoff_ssc_int = 3000,
#                                 is_rule_limitation = TRUE,
#                                 cutoff_topN = 5,
#                                 is_plot_pseudo_MS1 = TRUE,
#                                 type_order = c('level1', 'level2', 'level3'))


# list_peak_group_annotation <- pbapply::pblapply(seq_along(list_peak_group), function(i){
#   # mapProgressPrint(progress)
#   cat(i, ' ')
#   result <- annotatePeakGroup(peak_group = list_peak_group[[i]],
#                               ms2_data = raw_msms,
#                               polarity = 'positive',
#                               tol_mz = 25,
#                               isotope_int_ratio_check = 25,
#                               isotope_int_ratio_cutoff = 500,
#                               is_ms2_check = TRUE,
#                               ms2_score_cutoff = -1,
#                               cutoff_ssc = 0.3,
#                               cutoff_ssc_int = 3000,
#                               is_rule_limitation = TRUE,
#                               cutoff_topN = 5)
#
#   return(result)
#
# })

# result <- annotatePeakGroup(peak_group = list_peak_group[[i]],
#                             ms2_data = raw_msms,
#                             polarity = 'positive',
#                             tol_mz = 25,
#                             isotope_int_ratio_check = 25,
#                             isotope_int_ratio_cutoff = 500,
#                             is_ms2_check = TRUE,
#                             ms2_score_cutoff = -1,
#                             cutoff_ssc = 0.3,
#                             cutoff_ssc_int = 3000,
#                             is_rule_limitation = TRUE,
#                             cutoff_topN = 5)


# result_adduct <- annotateAdduct(peak_list = peak_list,
#                                   query_mz = query_mz,
#                                   query_rt = query_rt,
#                                   query_ccs = query_ccs,
#                                   query_peak_name = query_peak_name,
#                                   query_adduct = query_adduct,
#                                   polarity = polarity,
#                                   tol_mz = tol_mz,
#                                   cutoff_ssc = cutoff_ssc,
#                                   cutoff_ssc_int = cutoff_ssc_int,
#                                   is_rule_limitation = is_rule_limitation)



################################################################################
# 20211101 Haosong Zhang -------------------------------------------------------
#
# runMetDNA2(ms1_file_neg = "data.csv",
#            sample_info_file_neg = "sample.info.csv",
#            path_neg = '/home/zhouzw/Data_processing/20211028_zhulab_rp_library/test3/NEG/',
#            metdna_version = "version2",
#            polarity = "negative",
#            instrument = "ThermoExploris",
#            lib = 'zhuRPLib',
#            column = "rp",
#            ce = 'SNCE20_30_40%',
#            method_lc = "zhulabRP",
#            correct_p = FALSE,
#            extension_step = "0",
#            comp_group = c("LPS", "Ctrl"),
#            species = "hsa",
#            p_cutoff = 0.050000,
#            fc_cutoff = 1.000000,
#            dir_GenForm = "/share/home/zhuzhengjiang/share/bin",
#            is_bio_interpret = FALSE,
#            is_rt_calibration = TRUE,
#            is_webserver = FALSE)
#
# runMetDNA2(ms1_file_pos = "data.csv",
#            sample_info_file_pos = "sample.info.csv",
#            path_pos = "/home/zhouzw/Data_processing/01_debug/211101_zhanghaosong/POS/",
#            metdna_version = 'version2',
#            is_webserver = FALSE,
#            is_rm_intermediate_data = FALSE,
#            is_plot_ms2 = TRUE,
#            seed_neighbor_match_plot = TRUE,
#            ms2_type = 'msp',
#            thread = 2,
#            lib = 'zhuRPLib',
#            polarity = 'positive',
#            instrument = "ThermoExploris",
#            column = "rp",
#            ce = "SNCE20_30_40%",
#            method_lc = 'zhulabRP',
#            is_rt_calibration = TRUE,
#            extension_step = '0',
#            comp_group = c("Liver1","Liver2"),
#            species = "hsa",
#            p_cutoff = 0.050000,
#            fc_cutoff = 1.000000,
#            dir_GenForm = "/share/home/zhuzhengjiang/share/bin",
#            is_bio_interpret = FALSE)
#
# table_annotation <- readr::read_csv("/home/zhouzw/Data_processing/01_debug/211101_zhanghaosong/POS/01_result_initial_seed_annotation/ms2_match_annotation_result.csv")
# path_output <- '/home/zhouzw/Data_processing/01_debug/211101_zhanghaosong/POS/01_result_initial_seed_annotation'

################################################################################
# 20211112 Linfeng -------------------------------------------------------------
# load('/home/zhouzw/Data_processing/01_debug/211112_linfeng/49eca5614d883b5c3ba1a9605ce2105d/POS/01_result_initial_seed_annotation/00_intermediate_data/ms2')
#
# runMetDNA2(
#   path_pos = "/home/zhouzw/Data_processing/01_debug/211112_linfeng/49eca5614d883b5c3ba1a9605ce2105d/POS/",
#   path_neg = "/home/zhouzw/Data_processing/01_debug/211112_linfeng/49eca5614d883b5c3ba1a9605ce2105d/NEG/",
#   metdna_version = "version2",
#   polarity = "positive",
#   instrument = "ThermoOrbitrap",
#   column = "rp",
#   ce = "SNCE20_30_40%",
#   method_lc = "Other",
#   correct_p = FALSE,
#   extension_step = "0",
#   comp_group = c("CK", "D90"),
#   species = "dre",
#   p_cutoff = 0.050000,
#   fc_cutoff = 1.000000,
#   is_rt_calibration = FALSE)

# load('/home/zhouzw/Data_processing/01_debug/211112_linfeng/49eca5614d883b5c3ba1a9605ce2105d/POS/01_result_initial_seed_annotation/00_intermediate_data/ms1_result')
# load('/home/zhouzw/Data_processing/01_debug/211112_linfeng/49eca5614d883b5c3ba1a9605ce2105d/POS/01_result_initial_seed_annotation/00_intermediate_data/ms1_data')
#
# annotateInitialSeed(ms1_file = "data.csv",
#                     ms2_file = NULL,
#                     metdna_version = 'version2',
#                     ms2_type = 'mgf',
#                     path = "/home/zhouzw/Data_processing/01_debug/211112_linfeng/49eca5614d883b5c3ba1a9605ce2105d/POS/",
#
#                     polarity = 'positive',
#                     instrument = "ThermoOrbitrap",
#                     lib = 'zhuMetLib',
#                     column = 'rp',
#                     ce = '30',
#                     method_lc = 'Other',
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
#                     is_plot_ms2 = TRUE
#                     )
# path_output <- "/home/zhouzw/Data_processing/01_debug/211112_linfeng/49eca5614d883b5c3ba1a9605ce2105d/POS/01_result_initial_seed_annotation/test"
# dir.create(path_output, showWarnings = FALSE, recursive = TRUE)
#
# test <- combineMs1Ms2(ms1_data = ms1_data,
#               ms2_data = ms2_data,
#               ms2_type = ms2_type,
#               mz_tol_combine_ms1_ms2 = mz_tol_combine_ms1_ms2,
#               rt_tol_combine_ms1_ms2 = rt_tol_combine_ms1_ms2,
#               ccs_tol_combine_ms1_ms2 = ccs_tol_combine_ms1_ms2,
#               path = file.path(path_output, "00_intermediate_data"))
#
# ms2_data$info %>%
#   dplyr::filter(mz >= 472.31 & mz <= 472.33 & )

################################################################################
# 20211112 Linfeng -------------------------------------------------------------

# runMetDNA2(
#   path_pos = "/home/zhouzw/Data_processing/01_debug/211112_linfeng/9f657a0ccf6a0eb700abcf4e21093a20/POS/",
#   path_neg = "/home/zhouzw/Data_processing/01_debug/211112_linfeng/9f657a0ccf6a0eb700abcf4e21093a20/NEG/",
#   metdna_version = "version2",
#   polarity = "positive",
#   instrument = "ThermoOrbitrap",
#   column = "rp",
#   ce = "NCE30",
#   method_lc = "Other",
#   correct_p = FALSE,
#   extension_step = "0",
#   comp_group = c("CK", "D90"),
#   species = "dre",
#   p_cutoff = 0.050000,
#   fc_cutoff = 1.000000,
#   is_rt_calibration = FALSE,
#   dir_GenForm = "/share/home/zhuzhengjiang/share/bin",
#   is_anno_initial_seed = TRUE,
#   is_anno_mrn = FALSE,
#   is_credential = FALSE)
#
# interpretBiology(sample_file_pos = "data.csv",
#                  sample_info_file_pos = "sample.info.csv",
#                  table_identification_pair_file_pos = 'table3_identification_pair.csv',
#                  polarity = 'positive',
#                  path = "/home/zhouzw/Data_processing/01_debug/211112_linfeng/9f657a0ccf6a0eb700abcf4e21093a20/POS/",
#                  metdna_version = 'version2',
#                  comp_group = c("CK", "D90"),
#                  uni_test = 't',
#                  correct_p = FALSE,
#                  p_cutoff = 0.05,
#                  fc_cutoff = FALSE,
#                  species = 'dre',
#                  quanti_pathway_method = 'mean',
#                  organze_basis = 'kegg_id',
#                  extension_step = '0')

################################################################################
# runMetDNA2(ms1_file_pos = "data.csv",
#            sample_info_file_pos = "sample.info.csv",
#            path_pos = "/home/zhouzw/Data_processing/01_debug/211101_zhanghaosong/POS/",
#            metdna_version = 'version2',
#            is_webserver = FALSE,
#            is_rm_intermediate_data = FALSE,
#            is_plot_ms2 = TRUE,
#            seed_neighbor_match_plot = TRUE,
#            ms2_type = 'msp',
#            thread = 2,
#            lib = 'zhuRPLib',
#            polarity = 'positive',
#            instrument = "ThermoExploris",
#            column = "rp",
#            ce = "SNCE20_30_40%",
#            method_lc = 'zhulabRP',
#            is_rt_calibration = TRUE,
#            extension_step = '0',
#            comp_group = c("Liver1","Liver2"),
#            species = "hsa",
#            p_cutoff = 0.050000,
#            fc_cutoff = 1.000000,
#            dir_GenForm = "/share/home/zhuzhengjiang/share/bin",
#            is_bio_interpret = FALSE)

# load("/home/zhouzw/Data_processing/01_debug/211101_zhanghaosong/POS/01_result_initial_seed_annotation/00_intermediate_data/lib_meta")
# load("/home/zhouzw/Data_processing/01_debug/211101_zhanghaosong/POS/01_result_initial_seed_annotation/00_intermediate_data/rt_calibration_result")

################################################################################
# 20211217 Linfeng -------------------------------------------------------------
#
# load('/home/zhouzw/Data_processing/01_debug/211217_linfeng/557a32e88fc449fc50e198488e40642f/POS/01_result_initial_seed_annotation/00_intermediate_data/ms1_data')
# # load('/home/zhouzw/Data_processing/01_debug/211217_linfeng/557a32e88fc449fc50e198488e40642f/POS/01_result_initial_seed_annotation/00_intermediate_data/')
# path <- '/home/zhouzw/Data_processing/01_debug/211217_linfeng/557a32e88fc449fc50e198488e40642f/POS/'
# ms2_type <- 'mgf'
# file_list <- list.files(path)
# ms2_file <- grep(paste0('\\.(', ms2_type, ')$'), file_list, value = TRUE)
#
# ms2_data <- readMs2(ms2_file = file.path(path, ms2_file),
#                     ms2_type = ms2_type)
#
# ms2_data <- inteMs2(ms2_data = ms2_data,
#                     ms2_type = ms2_type,
#                     metdna_version = 'version2',
#                     instrument = "ThermoOrbitrap",
#                     is_include_precursor = TRUE,
#                     is_deisotope = FALSE,
#                     int_ms2_min_abs = 50,
#                     int_ms2_min_relative = 0.01,
#                     mz_range_ms2 = NULL)
#
# ms2 <- combineMs1Ms2(ms1_data = ms1_data,
#                      ms2_data = ms2_data,
#                      ms2_type = ms2_type,
#                      mz_tol_combine_ms1_ms2 = mz_tol_combine_ms1_ms2,
#                      rt_tol_combine_ms1_ms2 = rt_tol_combine_ms1_ms2,
#                      ccs_tol_combine_ms1_ms2 = ccs_tol_combine_ms1_ms2,
#                      path = file.path(path_output, "00_intermediate_data"))

################################################################################
# 20220115 ---------------------------------------------------------------------
# # test whether_use_predRT function -------------------------------------------
# MetDNA2(ms1_file = 'data.csv',
#         sample_info_file = "sample.info.csv",
#         metdna_version = 'version2',
#         ms2_type = 'mgf',
#         path = '/home/zhouzw/Data_processing/20220115_modify_RT_para/no_predRT/',
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
#         is_anno_mrn = TRUE,
#         is_credential = TRUE,
#         is_bio_interpret = FALSE,
#         whether_use_predRT = FALSE)


# MetDNA2(ms1_file = 'data.csv',
#         sample_info_file = "sample.info.csv",
#         metdna_version = 'version2',
#         ms2_type = 'mgf',
#         path = '/home/zhouzw/Data_processing/20220115_modify_RT_para/test_workflow/',
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
#         is_anno_initial_seed = TRUE,
#         is_anno_mrn = FALSE,
#         is_credential = FALSE,
#         is_bio_interpret = FALSE,
#         adduct_limit_specLib = c('[M+H]+', '[M]+'),
#         test_force_filtering_rt = 20,
#         whether_use_predRT = FALSE)


################################################################################
# test function of trainRandomForestWithCaret ----------------------------------
#
# data <- 'ms2_match_annotation_result.csv'
# path <- '/home/zhouzw/Data_processing/20220118_test_pred_rt_caret_package/nist_urine_pos/'
# inHouse_compound <- zhuMetlib$meta$compound %>% as.data.frame()
#
# data <- readInitialAnnotation(data = 'ms2_match_annotation_result.csv',
#                               direction = 'reverse',
#                               rt_filter = FALSE,
#                               inHouse_compound = inHouse_compound,
#                               instrument = 'SciexTripleTOF',
#                               path = file.path(path, "01_result_initial_seed_annotation"))
# # prefer_adduct <- 'M+H'
# # threads <- 3
# # data("md_inHouse_cpd", envir = environment())
# # use_default_md = TRUE
# # column <- 'hilic'
# # method_lc <- 'Amide23min'
# # md_inHouse_cpd <- md_zhumetlib
# # md_kegg_cpd <- md_mrn_emrn$version2
# md_inHouse_cpd <- md_zhumetlib
# md_emrn_cpd <- md_mrn_emrn$version2
#
# test1 <- predictRT(data = data,
#                   prefer_adduct = 'M+H',
#                   md_inHouse_cpd = md_inHouse_cpd,
#                   md_kegg_cpd = md_emrn_cpd,
#                   threads = 3,
#                   use_default_md = TRUE,
#                   column = 'hilic',
#                   method_lc = 'Amide23min',
#                   metdna_version = 'version2',
#                   use_pretrained_model = FALSE)

# test <- test1$KEGG.rt %>% tibble::rownames_to_column(var = 'id') %>% dplyr::as_tibble()
# test %>% dplyr::filter(id == 'KeggExd001791')
# test %>% dplyr::filter(id == 'KeggExd005706')
# test %>% dplyr::filter(id == 'KeggExd000923')
# test %>% dplyr::filter(id == 'KeggExd085779')
# test %>% dplyr::filter(id == 'KeggExd148705')
# test %>% dplyr::filter(id == 'KeggExd039980')
# test %>% dplyr::filter(id == 'KeggExd000970')
# test %>% dplyr::filter(id == 'KeggExd011330')
# test %>% dplyr::filter(id == 'KeggExd007445')
# test %>% dplyr::filter(id == 'KeggExd024929')

# test2 <- predictRT(data = data,
#                    prefer_adduct = 'M+H',
#                    md_inHouse_cpd = md_inHouse_cpd,
#                    md_kegg_cpd = md_emrn_cpd,
#                    threads = 3,
#                    use_default_md = TRUE,
#                    column = 'hilic',
#                    method_lc = 'Amide23min',
#                    metdna_version = 'version2',
#                    use_pretrained_model = TRUE)
#
# test <- test2$KEGG.rt %>% tibble::rownames_to_column(var = 'id') %>% dplyr::as_tibble()
# test %>% dplyr::filter(id == 'KeggExd001791')
# test %>% dplyr::filter(id == 'KeggExd005706')
# test %>% dplyr::filter(id == 'KeggExd000923')
# test %>% dplyr::filter(id == 'KeggExd085779')
# test %>% dplyr::filter(id == 'KeggExd148705')
# test %>% dplyr::filter(id == 'KeggExd039980')
# test %>% dplyr::filter(id == 'KeggExd000970')
# test %>% dplyr::filter(id == 'KeggExd011330')
# test %>% dplyr::filter(id == 'KeggExd007445')
# test %>% dplyr::filter(id == 'KeggExd024929')
#
#
# test3 <- predictRT(data = data,
#                    prefer_adduct = 'M+H',
#                    md_inHouse_cpd = md_inHouse_cpd,
#                    md_kegg_cpd = md_emrn_cpd,
#                    threads = 3,
#                    use_default_md = TRUE,
#                    column = 'hilic',
#                    method_lc = 'Amide23min',
#                    metdna_version = 'version1',
#                    use_pretrained_model = FALSE)
#
# test <- test3$KEGG.rt %>% tibble::rownames_to_column(var = 'id') %>% dplyr::as_tibble()
# test %>% dplyr::filter(id == 'KeggExd001791')
# test %>% dplyr::filter(id == 'KeggExd005706')
# test %>% dplyr::filter(id == 'KeggExd000923')
# test %>% dplyr::filter(id == 'KeggExd085779')
# test %>% dplyr::filter(id == 'KeggExd148705')
# test %>% dplyr::filter(id == 'KeggExd039980')
# test %>% dplyr::filter(id == 'KeggExd000970')
# test %>% dplyr::filter(id == 'KeggExd011330')
# test %>% dplyr::filter(id == 'KeggExd007445')
# test %>% dplyr::filter(id == 'KeggExd024929')

# load('/home/zhouzw/Data_processing/20220118_test_pred_rt_caret_package/test_remove/00_annotation_table/00_intermediate_data/table_identification')
# removeReplicateSeedAnnotationFromFinalTable(dir_path = '/home/zhouzw/Data_processing/20220118_test_pred_rt_caret_package/test_remove',
#                                             instrument = "SciexTripleTOF")
#
# table_identification_old <- table_identification
#
# load('/home/zhouzw/Data_processing/20220118_test_pred_rt_caret_package/test_remove/00_annotation_table/00_intermediate_data/table_identification')


################################################################################
# 20220118 test whole workflow -------------------------------------------------


# runMetDNA2(ms1_file_pos = "data.csv",
#            sample_info_file_pos = "sample.info.csv",
#            path_pos = "/home/zhouzw/Data_processing/20220118_test_pred_rt_caret_package/test_workflow/",
#            metdna_version = 'version2',
#            is_webserver = FALSE,
#            is_rm_intermediate_data = FALSE,
#            is_plot_ms2 = TRUE,
#            seed_neighbor_match_plot = TRUE,
#            ms2_type = 'mgf',
#            thread = 2,
#            lib = 'zhuMetLib',
#            polarity = 'positive',
#            instrument = "SciexTripleTOF",
#            column = "hilic",
#            ce = "30",
#            method_lc = 'Amide23min',
#            is_rt_calibration = TRUE,
#            extension_step = '2',
#            # dir_GenForm = "/share/home/zhuzhengjiang/share/bin",
#            comp_group = c("pos_3_DDA", "QC"),
#            species = 'hsa',
#            quanti_pathway_method = 'mean',
#            p_cutoff = 0.05,
#            correct_p = TRUE,
#            fc_cutoff = 1,
#            adduct_limit_specLib = c('[M+H]+', '[M]+'),
#            test_force_filtering_rt = 20,
#            use_pretrained_model = FALSE,
#            remove_conflict_seed_from_final_table = TRUE)
# load('/home/zhouzw/Data_processing/20220118_test_pred_rt_caret_package/test_workflow/01_result_initial_seed_annotation/00_intermediate_data/result_annotation')
#
# test3 <- lapply(result_annotation, function(x){
#   if (nrow(x@annotation_result) > 0) {
#     return(x@annotation_result)
#   } else {
#     return(NULL)
#   }
# })
#
# test3 %>% dplyr::bind_rows() %>% View()

################################################################################
# test v0.6.62 -----------------------------------------------------------------
# library(MetDNA2)
# runMetDNA2(ms1_file_pos = "data.csv",
#            sample_info_file_pos = "sample.info.csv",
#            path_pos = "/home/zhouzw/Data_processing/20220118_test_pred_rt_caret_package/test2/",
#            metdna_version = 'version2',
#            is_webserver = FALSE,
#            is_rm_intermediate_data = FALSE,
#            is_plot_ms2 = TRUE,
#            seed_neighbor_match_plot = FALSE,
#            is_plot_pseudo_MS1 = FALSE,
#            ms2_type = 'mgf',
#            thread = 2,
#            lib = 'zhuMetLib',
#            polarity = 'positive',
#            instrument = "SciexTripleTOF",
#            column = "hilic",
#            ce = "30",
#            method_lc = 'Amide23min',
#            is_rt_calibration = TRUE,
#            extension_step = '2',
#            candidate_num = 10,
#            dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
#            comp_group = c("pos_3_DDA", "QC"),
#            species = 'hsa',
#            quanti_pathway_method = 'mean',
#            p_cutoff = 0.05,
#            correct_p = TRUE,
#            fc_cutoff = 1,
#            adduct_limit_specLib = c('[M+H]+', '[M]+'),
#            test_force_filtering_rt = 20,
#            mz_ppm_thr = 400,
#            use_pretrained_model = FALSE,
#            remove_conflict_seed_from_final_table = TRUE,
#            is_check_data = FALSE,
#            is_anno_initial_seed = FALSE,
#            is_anno_mrn = TRUE,
#            is_credential = TRUE,
#            is_bio_interpret = FALSE)

################################################################################
# test v0.6.64 -----------------------------------------------------------------
# library(MetDNA2)
# runMetDNA2(ms1_file_pos = "data.csv",
#            sample_info_file_pos = "sample.info.csv",
#            path_pos = "/home/zhouzw/Data_processing/20220118_test_pred_rt_caret_package/test2/",
#            metdna_version = 'version2',
#            is_webserver = FALSE,
#            is_rm_intermediate_data = FALSE,
#            is_plot_ms2 = TRUE,
#            seed_neighbor_match_plot = FALSE,
#            is_plot_pseudo_MS1 = FALSE,
#            ms2_type = 'mgf',
#            thread = 2,
#            lib = 'zhuMetLib',
#            polarity = 'positive',
#            instrument = "SciexTripleTOF",
#            column = "hilic",
#            ce = "30",
#            method_lc = 'Amide23min',
#            is_rt_calibration = TRUE,
#            extension_step = '2',
#            candidate_num = 10,
#            dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
#            comp_group = c("pos_3_DDA", "QC"),
#            species = 'hsa',
#            quanti_pathway_method = 'mean',
#            p_cutoff = 0.05,
#            correct_p = TRUE,
#            fc_cutoff = 1,
#            adduct_limit_specLib = c('[M+H]+', '[M]+'),
#            test_force_filtering_rt = 20,
#            mz_tol = 15,
#            mz_ppm_thr = 300,
#            matched_frag_cutoff = 4,
#            whether_link_frag = TRUE,
#            use_pretrained_model = FALSE,
#            remove_conflict_seed_from_final_table = TRUE,
#            is_check_data = FALSE,
#            is_anno_initial_seed = FALSE,
#            is_anno_mrn = TRUE,
#            is_credential = TRUE,
#            is_bio_interpret = FALSE)

################################################################################
# #
# runMetDNA2(
#   path_pos = "/home/zhouzw/Data_processing/20220118_test_pred_rt_caret_package/testV1.0/",
#   metdna_version = "version2",
#   polarity = "positive",
#   instrument = "SciexTripleTOF",
#   column = "hilic",
#   ce = "30",
#   method_lc = "Amide23min",
#   correct_p = FALSE,
#   extension_step = "2",
#   comp_group = c("pos_3_DDA", "QC"),
#   species = "hsa",
#   dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
#   p_cutoff = 0.050000,
#   fc_cutoff = 1.000000,
#   is_rt_calibration = TRUE,
#   is_anno_initial_seed = FALSE)

################################################################################
# runMetDNA2(
#   path_pos = "/home/zhouzw/Data_processing/20220222_rp_method_test/POS/",
#   metdna_version = "version2",
#   polarity = "positive",
#   instrument = "SciexTripleTOF",
#   column = "hilic",
#   ce = "30",
#   method_lc = "Amide23min",
#   correct_p = FALSE,
#   extension_step = "2",
#   comp_group = c("pos_3_DDA", "QC"),
#   species = "hsa",
#   is_webserver = FALSE,
#   dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
#   p_cutoff = 0.050000,
#   fc_cutoff = 1.000000,
#   is_rt_calibration = TRUE,
#   is_anno_initial_seed = FALSE)

#
# runMetDNA2(
#   path_pos = "/home/zhouzw/Data_processing/20220222_rp_method_test/POS/",
#   metdna_version = "version2",
#   polarity = "positive",
#   instrument = "ThermoExploris",
#   column = "rp",
#   ce = "SNCE20_30_40%",
#   method_lc = "RP12min",
#   correct_p = FALSE,
#   extension_step = "2",
#   comp_group = c("QC", "plasma"),
#   species = "hsa",
#   dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
#   p_cutoff = 0.050000,
#   fc_cutoff = 1.000000,
#   is_rt_calibration = TRUE)
#
#
# runMetDNA2(
#   path_pos = "/home/zhouzw/Data_processing/20220222_rp_method_test/POS/",
#   metdna_version = "version2",
#   polarity = "positive",
#   instrument = "ThermoExploris",
#   column = "rp",
#   ce = "SNCE20_30_40%",
#   method_lc = "RP12min",
#   correct_p = FALSE,
#   extension_step = "2",
#   comp_group = c("QC", "plasma"),
#   species = "hsa",
#   is_webserver = FALSE,
#   dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
#   p_cutoff = 0.050000,
#   fc_cutoff = 1.000000,
#   is_rt_calibration = TRUE)

################################################################################
# 20220222 update rt_ref for RP12min -------------------------------------------
# lib_rt$zhulabRP_exp
# pos_table <- readr::read_csv('/home/zhouzw/Data_processing/20220222_rp_method_test/RT_recalibration_table_pos.csv')
# idx_pos <- match(pos_table$id.zhulab, lib_rt$zhulabRP_exp$id)
#
# neg_table <- readr::read_csv('/home/zhouzw/Data_processing/20220222_rp_method_test/RT_recalibration_table_neg.csv')
# idx_neg <- match(pos_table$id.zhulab, lib_rt$zhulabRP_exp$id)
#
# ids <- c(pos_table$id.zhulab, neg_table$id.zhulab) %>% unique()
# idx <- match(ids, lib_rt$zhulabRP_exp$id)
# temp <- data.frame(name = ids, rt = lib_rt$zhulabRP_exp$rt[idx], stringsAsFactors = FALSE)
# rt_ref$zhulabRP <- temp
#
# names(rt_ref)[4] <- 'RP12min'
#
# usethis::use_data(rt_ref, overwrite = TRUE)


################################################################################

# load('/home/zhouzw/Data_processing/01_debug/220301_ZhangdanXie/MetDNA2/NEG/00_annotation_table/00_intermediate_data/annot_all')
# generateExportTable(annot_all = annot_all,
#                     path = '/home/zhouzw/Data_processing/01_debug/220301_ZhangdanXie/MetDNA2/NEG/',
#                     direction = 'reverse',
#                     tolerance_rt_range = 30,
#                     path = 'reverse',
#                     instrument = 'ThermoExploris',
#                     is_cred_formula_filter = FALSE,
#                     is_anno_mrn = TRUE,
#                     is_credential = TRUE)


# generateExportTable(annot_all = annot_all,
#                     direction = 'reverse',
#                     tolerance_rt_range = 30,
#                     lib = 'zhuMetLib',
#                     path = '/home/zhouzw/Data_processing/01_debug/220301_ZhangdanXie/MetDNA2/NEG/',
#                     instrument = 'SciexTripleTOF',
#                     is_cred_formula_filter = FALSE,
#                     candidate_num = 10,
#                     mz_tol = 25,
#                     rt_tol = 30)


# /home/zhouzw/Data_processing/20220115_modify_RT_para/test_workflow/


# runMetDNA2(
#   path_pos = "/home/zhouzw/Data_processing/20220115_modify_RT_para/test_workflow/",
#   metdna_version = "version2",
#   polarity = "positive",
#   instrument = "SciexTripleTOF",
#   column = "hilic",
#   ce = "30",
#   method_lc = "Amide23min",
#   correct_p = FALSE,
#   extension_step = "2",
#   comp_group = c("pos_3_DDA", "QC"),
#   species = "hsa",
#   is_webserver = FALSE,
#   dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
#   p_cutoff = 0.050000,
#   fc_cutoff = 1.000000,
#   is_rt_calibration = TRUE)

# runMetDNA2(
#   path_pos = "/home/zhouzw/Data_processing/20220115_modify_RT_para/test_workflow/",
#   metdna_version = "version2",
#   polarity = "positive",
#   instrument = "SciexTripleTOF",
#   column = "hilic",
#   ce = "30",
#   method_lc = "Amide23min",
#   correct_p = FALSE,
#   extension_step = "2",
#   comp_group = c("pos_3_DDA", "QC"),
#   species = "hsa",
#   is_webserver = FALSE,
#   dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
#   p_cutoff = 0.050000,
#   fc_cutoff = 1.000000,
#   is_rt_calibration = TRUE, 
#   is_anno_mrn = FALSE,
#   is_rm_intermediate_data = FALSE)
# 
# 
# runMetDNA2(
#   path_pos = "/home/zhouzw/Data_processing/20220115_modify_RT_para/test_workflow3/",
#   metdna_version = "version2",
#   polarity = "positive",
#   instrument = "SciexTripleTOF",
#   column = "hilic",
#   ce = "30",
#   method_lc = "Amide23min",
#   correct_p = FALSE,
#   extension_step = "2",
#   comp_group = c("pos_3_DDA", "QC"),
#   species = "hsa",
#   is_webserver = FALSE,
#   dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
#   p_cutoff = 0.050000,
#   fc_cutoff = 1.000000,
#   tolerance_rt_range = 20,
#   is_rt_calibration = TRUE,
#   is_anno_mrn = FALSE,
#   is_credential = FALSE,
#   is_rm_intermediate_data = FALSE)


################################################################################
# 
# # library(MetDNA2)
# runMetDNA2(ms1_file_pos = "data.csv",
#            ms1_file_neg = "data.csv",
#            sample_info_file_pos = "sample.info.csv",
#            sample_info_file_neg = "sample.info.csv",
#            path_pos = "/home/zhouzw/Data_processing/01_debug/220322_zhanghaosong/POS",
#            # path_neg = "/home/zhouzw/Data_processing/01_debug/220322_zhanghaosong/NEG",
#            metdna_version = 'version2',
#            is_webserver = FALSE,
#            is_rm_intermediate_data = FALSE,
#            is_plot_ms2 = TRUE,
#            ms2_type = 'msp',
#            thread = 8,
#            lib = 'zhumetlib_orbitrap',
#            polarity = 'positive',
#            instrument = "ThermoExploris",
#            column = "rp",
#            ce = "SNCE20_30_40%",
#            method_lc = 'RP12min',
#            is_rt_calibration = TRUE,
#            extension_step = '0',
#            comp_group = c("urine"),
#            species = "hsa",
#            p_cutoff = 0.05,
#            fc_cutoff = 1.00,
#            dir_GenForm = '/home/zhouzw/software/genform-code-r8-trunk/src/GenForm/GenForm/bin/Release',
#            is_bio_interpret = F,
#            is_anno_mrn = T,
#            is_credential = T,
#            adduct_limit_specLib = c('[M+H]+','[M]+','[M-H2O+H]+'),
#            is_exported_report = T)

################################################################################
# 
# 
################################################################################
# 20220412 test ----------------------------------------------------------------
# 
# runMetDNA2(
#   path_pos = "D:/project/00_zhulab/01_metdna2/00_data/20220412_test_metdna2/",
#   metdna_version = "version2",
#   polarity = "positive",
#   instrument = "ThermoExploris",
#   column = "hilic",
#   ce = "30",
#   method_lc = "Amide23min",
#   correct_p = FALSE,
#   extension_step = "0",
#   comp_group = c("pos_3_DDA", "QC"),
#   species = "hsa",
#   dir_GenForm = "D:/project/packages/GenForm/",
#   platform = 'windows',
#   p_cutoff = 0.050000,
#   fc_cutoff = 1.000000,
#   is_rt_calibration = TRUE,
#   thread = 2,
#   is_webserver = FALSE,
#   is_rm_intermediate_data = FALSE,
#   is_plot_ms2 = TRUE,
#   is_anno_initial_seed = FALSE)

# 
# annotation_result = "ms2_match_annotation_result.csv"
# ms2_data = "ms2"
# prefer_adduct = "all"
# use_default_md = TRUE
# use_pretrained_model = FALSE
# metdna_version = 'version2'
# scoring_approach = 'dp'
# direction = 'reverse'
# instrument = 'ThermoExploris'
# lib = 'zhumetlib_orbitrap'
# column = "hilic"
# method_lc = 'Amide23min'
# excluded_adduct = NULL
# polarity = "positive"
# extension_step = '0'
# threads = 3
# path = '.'
# max_isotope = 4
# mz_tol = 15
# mz_ppm_thr = 400
# rt_tol1 = 3
# rt_tol2 = 30
# ccs_tol = 4
# cor_tol = 0
# int_tol = 500
# dp_tol = 0.5
# matched_frag_tol = 1
# whether_link_frag = FALSE
# max_step = 3
# score_cutoff = 0
# remain = FALSE
# remain_per = 0.5
# seed_neighbor_match_plot = TRUE
# candidate_num = 3
# # test_old_mrn = c('v0.3', 'v0.4'),
# test_adduct_version = 'version2'
# test_evaluation = 'No'
# whether_use_predRT = TRUE
# 
# 

################################################################################
# load('D:/project/00_zhulab/01_metdna2/00_data/20220412_test_metdna2/test/adduct_table')
# load('D:/project/00_zhulab/01_metdna2/00_data/20220412_test_metdna2/test/metabolic_network')
# load('D:/project/00_zhulab/01_metdna2/00_data/20220412_test_metdna2/test/metabolite')
# load('D:/project/00_zhulab/01_metdna2/00_data/20220412_test_metdna2/test/metabolite_ccs')
# load('D:/project/00_zhulab/01_metdna2/00_data/20220412_test_metdna2/test/ms2')
# load('D:/project/00_zhulab/01_metdna2/00_data/20220412_test_metdna2/test/peak_int')
# load('D:/project/00_zhulab/01_metdna2/00_data/20220412_test_metdna2/test/round')
# load('D:/project/00_zhulab/01_metdna2/00_data/20220412_test_metdna2/test/seed_idx')
# load('D:/project/00_zhulab/01_metdna2/00_data/20220412_test_metdna2/test/tags2')
# load('D:/project/00_zhulab/01_metdna2/00_data/20220412_test_metdna2/test/adduct_result')
# 
# max_isotope = 4
# polarity = "positive"
# instrument = 'ThermoExploris'
# scoring_approach = 'dp'
# rt_tol1 = 3 # rt_tol1 is absolute
# rt_tol2 = 30 # rt_tol2 is relative
# mz_tol = 15
# mz_ppm_thr = 0
# cor_tol = 0
# int_tol = 500 # for isotope annotation
# ccs_tol = 4 # for metaboliteAnnotation, %
# dp_tol = 0.5
# matched_frag_tol = 1
# whether_link_frag = TRUE # whehter consider propagation with matched fragment number
# max_step = 3
# mz_tol_ms2 = 25 # fragment match tolerance
# # cor.matrix,
# remove_index = NULL
# iso_annotation = TRUE
# add_annotation = TRUE
# met_annotation = TRUE
# threads = 3
# whether_use_predRT = TRUE

################################################################################