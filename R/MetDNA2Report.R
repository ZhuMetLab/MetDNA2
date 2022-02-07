################################################################################
# generateMetDNA2Report --------------------------------------------------------
# path <- '/home/zhouzw/Data_processing/20210330_MetDNA2_report/'
#
# library(dplyr)
# generateMetDNA2Report(sample_info = 'sample.info.csv',
#                       path = '/home/zhouzw/Data_processing/20210330_MetDNA2_report/',
#                       export_type = 'html')
#
# generateMetDNA2Report(sample_info = 'sample.info.csv',
#                       path = '/home/zhouzw/Data_processing/20210330_MetDNA2_report/',
#                       export_type = 'pdf')
#
# generateMetDNA2Report(sample_info = 'sample.info.csv',
#                       path = '/home/zhouzw/Data_processing/20210330_debug/200STD_mouse_liver_pos/',
#                       export_type = 'html',
#                       extension_step = '2',
#                       is_rt_calibration = TRUE,
#                       is_bio_interpret = FALSE)

#' @title generateMetDNA2Report
#' @author Zhiwei Zhou
#' @param sample_info
#' @param path
#' @param export_type
#' @export

setGeneric(name = 'generateMetDNA2Report',
           def = function(
             sample_info = 'sample.info.csv',
             path = '.',
             polarity = c('positive', 'negative', 'both'),
             export_type = c('html', 'pdf', 'all'),
             extension_step = c('0', '1', '2', '3', '4', '5', '6', '7', '8'),
             is_rt_calibration = TRUE,
             is_bio_interpret = TRUE
           ){

             export_type <- match.arg(export_type)

             if ('05_analysis_report' %in% list.files(path)) {
               unlink(file.path(path, '05_analysis_report'), recursive = TRUE, force = TRUE)
             }


             # generate data for positive or negative mode ---------------------
             if (polarity %in% c('positive', 'negative')) {
               path_output <- file.path(path, '05_analysis_report')
               dir.create(path_output, showWarnings = FALSE, recursive = TRUE)

               # copy template to directory
               if (extension_step == '0') {
                 rmarkdown::draft(file = file.path(path_output, '01_template'),
                                  template = "MetDNA2", package = "MetDNA2",
                                  create_dir = TRUE, edit = FALSE)
               } else {
                 rmarkdown::draft(file = file.path(path_output, '01_template'),
                                  template = "MetDNA2Extension", package = "MetDNA2",
                                  create_dir = TRUE, edit = FALSE)
               }

               # generate parameter list
               para_table <- generateParaList4Report(para_file = 'para_list.txt',
                                                     path = path,
                                                     showed_para = c('package_version', 'ms1_file', 'sample_info_file',
                                                                     'lib', 'polarity', 'instrument', 'column', 'ce', 'method_lc',
                                                                     'is_check_data', 'is_anno_initial_seed', 'is_anno_mrn', 'is_credential',
                                                                     'mz_tol', 'extension_step', 'formula_elements', 'comp_group', 'uni_test', 'species'))
               save(para_table,
                    file = file.path(path_output, '01_template', 'para_table.RData'),
                    version = 2)

               # generate sammary info
               file.copy(from = file.path(path, '01_result_initial_seed_annotation/00_intermediate_data', 'ms1_data'),
                         to = file.path(path_output, '01_template', 'ms1_data'), overwrite = TRUE, recursive = TRUE)

               load(file.path(path, '01_result_initial_seed_annotation/00_intermediate_data', 'ms1_data'))
               options(readr.num_columns = 0)
               sample_info <- readr::read_csv(file = file.path(path, 'sample.info.csv'))

               summary_info <- generateOverviewInfo4Report(ms1_data = ms1_data, sample_info = sample_info)
               save(summary_info, file = file.path(path_output, '01_template', 'summary_info.RData'), version = 2)

               # cp rt_calibration
               file.copy(from = file.path(path, '01_result_initial_seed_annotation/00_intermediate_data', 'rt_calibration_result'),
                         to = file.path(path_output, '01_template', 'rt_calibration_result'), overwrite = TRUE, recursive = TRUE)

               # generate annotation summary table
               load(file.path(path, '00_annotation_table/00_intermediate_data', 'table_identification'))
               annotation_summary <- generateAnnotationSummary4Report(result_annotation = table_identification)
               save(annotation_summary, file = file.path(path_output, '01_template', 'annotation_summary.RData'), version = 2)

               # generate recursive summary table
               load(file.path(path, '02_result_MRN_annotation/00_intermediate_data', 'tags2_after_redundancy_remove'))
               data("cpd_emrn", envir = environment())
               recursive_summary <- generateRecursiveSummary4Report(result_recursive = tags2_after_redundancy_remove,
                                                                    cpd_emrn = cpd_emrn)
               save(recursive_summary, file = file.path(path_output, '01_template', 'recursive_summary.RData'), version = 2)

               # generate credential summary table
               load(file.path(path, '03_annotation_credential/00_intermediate_data', 'annotation_initial.RData'))
               load(file.path(path, '03_annotation_credential/00_intermediate_data', 'list_peak_group_annotation_concised.RData'))
               options(readr.num_columns = 0)
               result_credential_initial <- annotation_initial
               result_credential_pg <- list_peak_group_annotation_concised
               result_credential <- readr::read_csv(file.path(path, '03_annotation_credential', 'annontation_credential_long.csv'))
               result_credential_filter <- readr::read_csv(file.path(path, '03_annotation_credential', 'annontation_credential_filter.csv'))

               credential_summary <- generateCredentialSummary4Report(result_credential_initial = result_credential_initial,
                                                                      result_credential = result_credential,
                                                                      result_credential_pg = result_credential_pg,
                                                                      result_credential_filter = result_credential_filter)

               save(credential_summary, file = file.path(path_output, '01_template', 'credential_summary.RData'), version = 2)


               # generate biolog interpretation summary table
               if (is_bio_interpret) {
                 file.copy(from = file.path(path, '04_biology_intepretation/00_intermediate_data', 'result_stat.RData'),
                           to = file.path(path_output, '01_template', 'result_stat.RData'), overwrite = TRUE, recursive = TRUE)
                 file.copy(from = file.path(path, '04_biology_intepretation/00_intermediate_data', 'result_pathway_enrichment.RData'),
                           to = file.path(path_output, '01_template', 'result_pathway_enrichment.RData'), overwrite = TRUE, recursive = TRUE)

                 load(file = file.path(path, '04_biology_intepretation/00_intermediate_data', 'result_stat.RData'))
                 load(file = file.path(path, '04_biology_intepretation/00_intermediate_data', 'result_pathway_enrichment.RData'))

                 bio_interp_summary <- generateBioInterpSummary4Report(result_stat_feature = result_stat,
                                                                       result_pathway = result_pathway_enrichment,
                                                                       p_cutoff = 0.05)

                 save(bio_interp_summary, file = file.path(path_output, '01_template', 'bio_interp_summary.RData'), version = 2)
               }

             }

             # generate data for both mode -------------------------------------
             if (polarity == 'both') {
               path_output <- file.path(path, 'BOTH', '05_analysis_report')
               dir.create(path_output, showWarnings = FALSE, recursive = TRUE)

               if (extension_step == '0') {
                 rmarkdown::draft(file = file.path(path_output, '01_template'),
                                  template = "MetDNA2BothPolarity", package = "MetDNA2",
                                  create_dir = TRUE, edit = FALSE)
               } else {
                 rmarkdown::draft(file = file.path(path_output, '01_template'),
                                  template = "MetDNA2ExtensionBothPolarity", package = "MetDNA2",
                                  create_dir = TRUE, edit = FALSE)
               }

               # generate parameter list
               para_table <- generateParaList4Report(para_file = 'para_list.txt',
                                                     path = file.path(path, 'BOTH'),
                                                     showed_para = c('package_version', 'ms1_file', 'sample_info_file',
                                                                     'lib', 'polarity', 'instrument', 'column', 'ce', 'method_lc',
                                                                     'is_check_data', 'is_anno_initial_seed', 'is_anno_mrn', 'is_credential',
                                                                     'mz_tol', 'extension_step', 'formula_elements', 'comp_group', 'uni_test', 'species'))
               save(para_table,
                    file = file.path(path_output, '01_template', 'para_table.RData'),
                    version = 2)

               # generate sammary info
               file.copy(from = file.path(path, 'POS', '01_result_initial_seed_annotation/00_intermediate_data', 'ms1_data'),
                         to = file.path(path_output, '01_template', 'ms1_data_pos'), overwrite = TRUE, recursive = TRUE)

               file.copy(from = file.path(path, 'NEG', '01_result_initial_seed_annotation/00_intermediate_data', 'ms1_data'),
                         to = file.path(path_output, '01_template', 'ms1_data_neg'), overwrite = TRUE, recursive = TRUE)

               file.copy(from = file.path(path, 'POS', 'sample.info.csv'),
                         to = file.path(path_output, '01_template', 'sample.info_pos.csv'), overwrite = TRUE, recursive = TRUE)

               file.copy(from = file.path(path, 'NEG', 'sample.info.csv'),
                         to = file.path(path_output, '01_template', 'sample.info_neg.csv'), overwrite = TRUE, recursive = TRUE)

               options(readr.num_columns = 0)
               load(file.path(path_output, '01_template', 'ms1_data_pos'))
               ms1_data_pos <- ms1_data
               sample_info_pos <- readr::read_csv(file = file.path(path_output, '01_template', 'sample.info_pos.csv'))
               summary_info_pos <- generateOverviewInfo4Report(ms1_data = ms1_data_pos, sample_info = sample_info_pos)
               save(summary_info_pos, file = file.path(path_output, '01_template', 'summary_info_pos.RData'), version = 2)
               rm('ms1_data');gc()

               load(file.path(path_output, '01_template', 'ms1_data_neg'))
               ms1_data_neg <- ms1_data
               sample_info_neg <- readr::read_csv(file = file.path(path_output, '01_template', 'sample.info_neg.csv'))
               summary_info_neg <- generateOverviewInfo4Report(ms1_data = ms1_data_neg, sample_info = sample_info_neg)
               save(summary_info_neg, file = file.path(path_output, '01_template', 'summary_info_neg.RData'), version = 2)
               rm('ms1_data');gc()

               # cp rt_calibration
               file.copy(from = file.path(path, 'POS', '01_result_initial_seed_annotation/00_intermediate_data', 'rt_calibration_result'),
                         to = file.path(path_output, '01_template', 'rt_calibration_result_pos'), overwrite = TRUE, recursive = TRUE)

               file.copy(from = file.path(path, 'NEG', '01_result_initial_seed_annotation/00_intermediate_data', 'rt_calibration_result'),
                         to = file.path(path_output, '01_template', 'rt_calibration_result_neg'), overwrite = TRUE, recursive = TRUE)


               # generate annotation summary table
               load(file.path(path, 'POS', '00_annotation_table/00_intermediate_data', 'table_identification'))
               annotation_summary_pos <- generateAnnotationSummary4Report(result_annotation = table_identification)
               save(annotation_summary_pos, file = file.path(path_output, '01_template', 'annotation_summary_pos.RData'), version = 2)
               rm('table_identification');gc()

               load(file.path(path, 'NEG', '00_annotation_table/00_intermediate_data', 'table_identification'))
               annotation_summary_neg <- generateAnnotationSummary4Report(result_annotation = table_identification)
               save(annotation_summary_neg, file = file.path(path_output, '01_template', 'annotation_summary_neg.RData'), version = 2)
               rm('table_identification');gc()

               # generate recursive summary table
               load(file.path(path, 'POS', '02_result_MRN_annotation/00_intermediate_data', 'tags2_after_redundancy_remove'))
               data("cpd_emrn", envir = environment())
               recursive_summary_pos <- generateRecursiveSummary4Report(result_recursive = tags2_after_redundancy_remove,
                                                                        cpd_emrn = cpd_emrn)
               save(recursive_summary_pos, file = file.path(path_output, '01_template', 'recursive_summary_pos.RData'), version = 2)
               rm('tags2_after_redundancy_remove');gc()

               load(file.path(path, 'NEG', '02_result_MRN_annotation/00_intermediate_data', 'tags2_after_redundancy_remove'))
               recursive_summary_neg <- generateRecursiveSummary4Report(result_recursive = tags2_after_redundancy_remove,
                                                                        cpd_emrn = cpd_emrn)
               save(recursive_summary_neg, file = file.path(path_output, '01_template', 'recursive_summary_neg.RData'), version = 2)
               rm('tags2_after_redundancy_remove', 'cpd_emrn');gc()

               # generate credential summary table
               load(file.path(path, 'POS', '03_annotation_credential/00_intermediate_data', 'annotation_initial.RData'))
               load(file.path(path, 'POS', '03_annotation_credential/00_intermediate_data', 'list_peak_group_annotation_concised.RData'))
               options(readr.num_columns = 0)
               result_credential_initial_pos <- annotation_initial
               result_credential_pg_pos <- list_peak_group_annotation_concised
               result_credential_pos <- readr::read_csv(file.path(path, 'POS', '03_annotation_credential', 'annontation_credential_long.csv'))
               result_credential_filter_pos <- readr::read_csv(file.path(path, 'POS', '03_annotation_credential', 'annontation_credential_filter.csv'))
               credential_summary_pos <- generateCredentialSummary4Report(result_credential_initial = result_credential_initial_pos,
                                                                          result_credential = result_credential_pos,
                                                                          result_credential_pg = result_credential_pg_pos,
                                                                          result_credential_filter = result_credential_filter_pos)
               save(credential_summary_pos, file = file.path(path_output, '01_template', 'credential_summary_pos.RData'), version = 2)
               rm('annotation_initial', 'list_peak_group_annotation_concised');gc()


               load(file.path(path, 'NEG', '03_annotation_credential/00_intermediate_data', 'annotation_initial.RData'))
               load(file.path(path, 'NEG', '03_annotation_credential/00_intermediate_data', 'list_peak_group_annotation_concised.RData'))
               options(readr.num_columns = 0)
               result_credential_initial_neg <- annotation_initial
               result_credential_pg_neg <- list_peak_group_annotation_concised
               result_credential_neg <- readr::read_csv(file.path(path, 'NEG', '03_annotation_credential', 'annontation_credential_long.csv'))
               result_credential_filter_neg <- readr::read_csv(file.path(path, 'NEG', '03_annotation_credential', 'annontation_credential_filter.csv'))
               credential_summary_neg <- generateCredentialSummary4Report(result_credential_initial = result_credential_initial_neg,
                                                                          result_credential = result_credential_neg,
                                                                          result_credential_pg = result_credential_pg_neg,
                                                                          result_credential_filter = result_credential_filter_neg)
               save(credential_summary_neg, file = file.path(path_output, '01_template', 'credential_summary_neg.RData'), version = 2)
               rm('annotation_initial', 'list_peak_group_annotation_concised');gc()

             }

             # render report
             if (export_type == "html"){
               rmarkdown::render(file.path(path_output, '01_template', "MetDNA2_template.Rmd"),
                                 output_format = rmarkdown::html_document(),
                                 params = list(is_rt_calibration = is_rt_calibration,
                                               is_bio_interpret = is_bio_interpret)
               )
             }

             if (export_type == "pdf"){
               rmarkdown::render(file.path(path_output, '01_template', "MetDNA2_template.Rmd"),
                                 output_format = rmarkdown::pdf_document(),
                                 params = list(is_rt_calibration = is_rt_calibration,
                                               is_bio_interpret = is_bio_interpret)
               )
             }

             if (export_type == "all"){
               rmarkdown::render(file.path(path_output, '01_template', "MetDNA2_template.Rmd"),
                                 output_format = rmarkdown::html_document(),
                                 params = list(is_rt_calibration = is_rt_calibration,
                                               is_bio_interpret = is_bio_interpret))
               rmarkdown::render(file.path(path_output, '01_template', "MetDNA2_template.Rmd"),
                                 output_format = rmarkdown::pdf_document(),
                                 params = list(is_rt_calibration = is_rt_calibration,
                                               is_bio_interpret = is_bio_interpret))
             }

             file.copy(from = file.path(path_output, '01_template', "MetDNA2_template.html"),
                       to = file.path(path_output, 'MetDNA2_analysis_report.html'),
                       overwrite = TRUE,
                       recursive = TRUE)

             # unlink(file.path(path_output, '01_template'), recursive = TRUE, force = TRUE)

           })


################################################################################
#   generateParaList4Report ----------------------------------------------------
#' @title generateParaList4Report
#' @author Zhiwei Zhou
#' @param para_file 'para_list.txt'
#' @param path '.'
#' @param showed_para c('package_version', 'ms1_file', 'ms2_file', 'sample_info_file',
#' 'lib', 'polarity', 'instrument', 'column', 'ce', 'method_lc',
#' is_check_data', 'is_anno_initial_seed', 'is_anno_mrn', 'is_credential'
#' mz_tol', 'extension_step', 'formula_elements', 'comp_group', 'uni_test', 'species')


setGeneric(name = 'generateParaList4Report',
           def = function(
             para_file = 'para_list.txt',
             path = '.',
             showed_para = c('package_version', 'ms1_file', 'ms2_file', 'sample_info_file',
                             'lib', 'polarity', 'instrument', 'column', 'ce', 'method_lc',
                             'is_check_data', 'is_anno_initial_seed', 'is_anno_mrn', 'is_credential',
                             'mz_tol', 'extension_step', 'formula_elements', 'comp_group', 'uni_test', 'species')
           ){
             options(readr.num_columns = 0)
             para_list <- readr::read_tsv(file = file.path(path, para_file))
             para_list <- para_list %>% dplyr::filter(para %in% showed_para)
             colnames(para_list) <- c('Parameter', 'Value')

             return(para_list)
           })



#   generateOverviewInfo4Report -------------------------------------------------------

#' @title generateOverviewInfo4Report
#' @author Zhiwei Zhou
#' @description generate overview of data for report

setGeneric(name = 'generateOverviewInfo4Report',
           def = function(
             ms1_data,
             sample_info
           ){
             feature_table <- ms1_data$info
             sample_table <- ms1_data$subject


             n_feature <- nrow(feature_table)
             range_mz <- range(feature_table$mz)
             range_rt <- range(feature_table$rt)
             n_sample <- nrow(sample_info)
             n_group <- length(unique(sample_info$group))

             summary_info <- list(n_feature = n_feature,
                                  range_mz = range_mz,
                                  range_rt = range_rt,
                                  n_sample = n_sample,
                                  n_group = n_group)

             return(summary_info)

           })


#   generateAnnotationSummary4Report --------------------------------------------------
#' @title generateAnnotationSummary4Report
#' @author Zhiwei Zhou
#' @description generate annotation summary for report
#' @param result_annotaion

setGeneric(name = 'generateAnnotationSummary4Report',
           def = function(
             result_annotation
           ){
             n_feature_id <- result_annotation %>%
               # dplyr::filter(!is.na(id)) %>%
               dplyr::distinct(peak_name, .keep_all = TRUE) %>%
               nrow()

             n_met_id <- result_annotation %>%
               # dplyr::filter(!is.na(id)) %>%
               # dplyr::distinct(inchikey1, .keep_all = TRUE) %>%
               nrow()

             summary_table <- result_annotation %>%
               dplyr::count(confidence_level) %>%
               dplyr::ungroup() %>%
               dplyr::filter(!is.na(confidence_level)) %>%
               dplyr::rename('Confidence level' = confidence_level,
                             'Annotation number' = n)

             n_feature_id_initial <- result_annotation %>%
               # dplyr::filter(!is.na(id)) %>%
               dplyr::filter(confidence_level %in% c('level1', 'level2')) %>%
               dplyr::distinct(peak_name, .keep_all = TRUE) %>%
               nrow()

             n_met_id_initial <- result_annotation %>%
               # dplyr::filter(!is.na(id)) %>%
               dplyr::filter(confidence_level %in% c('level1', 'level2')) %>%
               nrow()

             n_feature_id_initial_level1 <- result_annotation %>%
               # dplyr::filter(!is.na(id)) %>%
               dplyr::filter(confidence_level == 'level1') %>%
               dplyr::distinct(peak_name, .keep_all = TRUE) %>%
               nrow()

             n_met_id_initial_level1 <- result_annotation %>%
               # dplyr::filter(!is.na(id)) %>%
               dplyr::filter(confidence_level == 'level1') %>%
               nrow()

             n_feature_id_initial_level2 <- result_annotation %>%
               # dplyr::filter(!is.na(id)) %>%
               dplyr::filter(confidence_level == 'level2') %>%
               dplyr::distinct(peak_name, .keep_all = TRUE) %>%
               nrow()

             n_met_id_initial_level2 <- result_annotation %>%
               # dplyr::filter(!is.na(id)) %>%
               dplyr::filter(confidence_level == 'level2') %>%
               nrow()

             # n_feature_id_recursive <- result_annotation %>%
             #   dplyr::filter(!is.na(id)) %>%
             #   dplyr::filter(confidence_level %in% c('level3.1', 'level3.2')) %>%
             #   dplyr::distinct(peak_name, .keep_all = TRUE) %>%
             #   nrow()
             #
             # n_met_id_recursive <- result_annotation %>%
             #   dplyr::filter(!is.na(id)) %>%
             #   dplyr::filter(confidence_level %in% c('level3.1', 'level3.2')) %>%
             #   # dplyr::distinct(inchikey, .keep_all = TRUE) %>%
             #   nrow()

             # result_annotation %>% count(recursive_type)

             result_summary <- list(n_feature_id = n_feature_id,
                                    n_met_id = n_met_id,
                                    summary_table = summary_table,
                                    n_feature_id_initial = n_feature_id_initial,
                                    n_met_id_initial = n_met_id_initial,
                                    n_feature_id_initial_level1 = n_feature_id_initial_level1,
                                    n_met_id_initial_level1 = n_met_id_initial_level1,
                                    n_feature_id_initial_level2 = n_feature_id_initial_level2,
                                    n_met_id_initial_level2 = n_met_id_initial_level2
                                    # ,n_feature_id_recursive =  n_feature_id_recursive,
                                    # n_met_id_recursive = n_met_id_recursive
                                    )

             return(result_summary)

           })

#   generateRecursiveSummary4Report ---------------------------------------------------
#' @title generateRecursiveSummary4Report
#' @author Zhiwei Zhou
#' @param result_recursive
#' @param cpd_emrn

# load('/home/zhouzw/Data_processing/20210330_MetDNA2_report/02_result_MRN_annotation/00_intermediate_data/tags2_after_redundancy_remove')
# result_recursive <- tags2_after_redundancy_remove

setGeneric(name = 'generateRecursiveSummary4Report',
           def = function(
             result_recursive,
             cpd_emrn
           ){
             table_recursive <- trans2Matrix(result_recursive) %>%
               tibble::as_tibble() %>%
               readr::type_convert()

             cpd_type <- match(table_recursive$to, cpd_emrn$id) %>% cpd_emrn$type[.]
             table_recursive <- table_recursive %>%
               dplyr::mutate(cpd_type = cpd_type) %>%
               tidyr::replace_na(list(cpd_type = 'known_known'))

             n_annotation_recursive <- table_recursive %>% nrow()
             n_feature_recursive <- table_recursive %>%
               dplyr::distinct(name, .keep_all = TRUE) %>%
               nrow()
             n_round <- table_recursive %>%
               dplyr::arrange(dplyr::desc(level)) %>%
               dplyr::pull(level) %>%
               .[1]


             n_seed <- table_recursive %>% dplyr::filter(type == 'seed') %>% nrow()
             n_kegg_met <- table_recursive %>% dplyr::filter(cpd_type %in% c('known_known')) %>% nrow()
             n_kegg_met <- n_kegg_met - n_seed

             n_curated_unknown_met <- table_recursive %>% dplyr::filter(cpd_type %in% c('known_unknown', 'unknown_unknown')) %>% nrow()

             summary_recursive <- list(n_annotation_recursive = n_annotation_recursive,
                                       n_feature_recursive = n_feature_recursive,
                                       n_round = n_round,
                                       n_seed = n_seed,
                                       n_keeg_met = n_kegg_met,
                                       n_curated_unknown_met = n_curated_unknown_met,
                                       table_recursive = table_recursive)

             return(summary_recursive)
           })


#   generateCredentialSummary4Report -------------------------------------------
#' @title generateCredentialSummary4Report
#' @author Zhiwei Zhou
#' @param result_credential_initial
#' @param result_credential
#' @param result_credential_pg
#' @param result_credential_filter


#
# load('/home/zhouzw/Data_processing/20210330_MetDNA2_report/03_annotation_credential/00_intermediate_data/list_peak_group_annotation_concised.RData')
# load('/home/zhouzw/Data_processing/20210330_MetDNA2_report/03_annotation_credential/00_intermediate_data/record_formula_filter.RData')
# load('/home/zhouzw/Data_processing/20210330_MetDNA2_report/03_annotation_credential/00_intermediate_data/annotation_initial.RData')

# load('/home/zhouzw/Data_processing/20210330_MetDNA2_report/03_annotation_credential/00_intermediate_data/list_peak_group_formula.RData')
# result_credential_filter <- readr::read_csv('/home/zhouzw/Data_processing/20210330_MetDNA2_report/03_annotation_credential/annontation_credential_filter.csv')
# result_credential <- readr::read_csv('/home/zhouzw/Data_processing/20210330_MetDNA2_report/03_annotation_credential/annontation_credential_long.csv')


# result_credential_pg <- list_peak_group_annotation_concised
# result_credential_initial <- annotation_initial


setGeneric(name = 'generateCredentialSummary4Report',
           def = function(
             result_credential_initial,
             result_credential,
             result_credential_pg,
             result_credential_filter
           ){

             # - A total of xxx feature with xxx annotations out of xxx features with xxx annotations (initial seed annotation and recursive annotation) were credentialed.

             # overview
             n_feature_initial <- result_credential_initial$name %>% unique() %>% length()
             n_metabolite_initial <- result_credential_initial$id %>% unique() %>% length()
             n_annotation_initial <- result_credential_initial %>% nrow()

             n_feature_final <- result_credential$name %>% unique() %>% length()
             n_metabolite_final <- result_credential$id %>% unique() %>% length()
             n_annotation_final <- result_credential %>% nrow()

             summary_overview <- tibble::tibble(Features = c(n_feature_initial, n_feature_final),
                                                Metabolites = c(n_metabolite_initial, n_metabolite_final),
                                                Annotations = c(n_annotation_initial, n_annotation_final))

             rownames(summary_overview) <- c('Recursive annotation', 'Annotation credential')

             # - A total of xxx feature with xxx annotation out of xxx features with xxx annotations (initial seed annotation and recursive annotation) were credentialed.
             # - These credentialed putative annotations cover xxx features and xxx related annotations as figure 6.

             # credential 1 summary
             n_credential_feature_pg <- sapply(result_credential_pg, function(x){
               x@base_peak_name
             }) %>% unique() %>% length()

             temp_credential_filter <- result_credential_filter %>% dplyr::filter(type != 'type3') %>% nrow()
             n_credential_annotations <- nrow(result_credential_initial) - temp_credential_filter

             summary_credential_pg <- statTypeListPeakGroup(result_credential_pg)


             # - A total of xxx feature with xxx annotation out of xxx features with xxx annotations were credentialed.
             # - A total of xxx annotations were assigned as level4 due to the conflict with predicted formula

             # credential 2 summary
             n_filter_annotation_formula <- result_credential_filter %>% dplyr::filter(type == 'type3') %>% nrow()

             result_summary_credential <- list(n_feature_initial = n_feature_initial,
                                               n_annotation_initial = n_annotation_initial,
                                               n_feature_final = n_feature_final,
                                               n_annotation_final = n_annotation_final,
                                               summary_overview = summary_overview,
                                               n_credential_feature_pg = n_credential_feature_pg,
                                               n_credential_feature_annotations = n_credential_annotations,
                                               summary_credential_pg = summary_credential_pg,
                                               n_credential_filter_formula = n_filter_annotation_formula)

             return(result_summary_credential)
           })


#   generateBioInterpSummary4Report --------------------------------------------

#' @title generateBioInterpSummary4Report
#' @author Zhiwei Zhou
#' @param result_stat_feature
#' @param result_pathway
#' @param p_cutoff

# load('/home/zhouzw/Data_processing/20210402_aging_fly_report_demonstration/04_biology_intepretation/00_intermediate_data/result_pathway_enrichment.RData')
# load('/home/zhouzw/Data_processing/20210402_aging_fly_report_demonstration/04_biology_intepretation/00_intermediate_data/result_stat.RData')
#
# result_stat_feature <- result_stat
# result_pathway <- result_pathway_enrichment
#
# test <- funtcion(result_stat_feature,
#                  result_pathway,
#                  ...) {
#
# }

# generateBioInterpSummary4Report(result_stat_feature = result_stat_feature, result_pathway = result_pathway, p_cutoff = 0.05)

setGeneric(name = 'generateBioInterpSummary4Report',
           def = function(
             result_stat_feature,
             result_pathway,
             p_cutoff = 0.05,
             ...
           ){
             n_sig_feature <- result_stat_feature %>%
               dplyr::filter(p_values <= 0.05 & (fold_changes >= 1.25 | fold_changes <= 0.8)) %>%
               nrow()

             table_pathway_enrich <- result_pathway %>%
               dplyr::filter(p.value <= 0.05)

             n_sig_pathway <- nrow(table_pathway_enrich)

             result_summary_bioInterp <- list(n_sig_feature = n_sig_feature,
                                              n_sig_pathway = n_sig_pathway,
                                              table_pathway_enrich = table_pathway_enrich)

             return(result_summary_bioInterp)

           })

