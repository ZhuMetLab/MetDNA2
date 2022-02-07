################################################################################
# generateMetDNA2AnnotationResult ----------------------------------------------
#' @title generateMetDNA2AnnotationResult
#' @author Zhiwei Zhou
#' @param path '.'
#' @param thread Default: 4
#' @param is_cred_pg_filter
#' @param is_cred_formula_filter

# generateMetDNA2AnnotationResult(path = '/home/zhouzw/Data_processing/20201102_metdna2_pos_development',
#                                 is_cred_pg_filter = TRUE,
#                                 is_cred_formula_filter = TRUE)


# generateMetDNA2AnnotationResult(path = '/home/zhouzw/Data_processing/20201201_debug_mouse_liver_200std/02_POS/modified_package_result/',
#                                 thread = 4,
#                                 is_cred_pg_filter = FALSE,
#                                 is_cred_formula_filter = FALSE,
#                                 is_pred_formula_all = FALSE,
#                                 tolerance_rt_range = 30)
#
# generateMetDNA2AnnotationResult(path = '/home/zhouzw/Data_processing/20210315_tims_data_processing_metdna/',
#                                 thread = 4,
#                                 is_cred_pg_filter = TRUE,
#                                 is_cred_formula_filter = TRUE,
#                                 is_pred_formula_all = FALSE,
#                                 tolerance_rt_range = 30)
#
#
# generateMetDNA2AnnotationResult(path = '/home/zhouzw/Data_processing/20210522_debug/01_aging_fly_pos_version2/',
#                                 thread = 4,
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
# generateExportTable(annot_all = annot_all,
#                     direction = direction,
#                     tolerance_rt_range = tolerance_rt_range,
#                     path = path,
#                     instrument = instrument,
#                     is_cred_formula_filter = is_cred_formula_filter,
#                     ...)

setGeneric(name = 'generateMetDNA2AnnotationResult',
           def = function(
             path = '.',
             thread = 4,
             is_cred_pg_filter = TRUE,
             is_cred_formula_filter = TRUE,
             is_pred_formula_all = FALSE,
             direction = c('reverse', 'forward'),
             tolerance_rt_range = 30,
             instrument = c("SciexTripleTOF", "AgilentQTOF", "BrukerQTOF", "ThermoOrbitrap",
                            'ThermoExploris', "AgilentDTIMMS", "BrukerTIMS", 'WatersQTOF','WatersTWIMMS'),
             remove_conflict_seed_from_final_table = FALSE,
             test_rm_extra_anno_from_ini_seed = FALSE,
             ...
           ){
             path_output <- file.path(path, '00_annotation_table')

             if ('annot_all' %in% list.files(file.path(path_output, '00_intermediate_data'))) {
               cat('Note: annot_all exsited, and load the old RData\n')
               load(file.path(path_output, '00_intermediate_data', 'annot_all'))
             } else {

               # add credential result into annot_all
               if ('annot_all_before_credential.RData' %in% list.files(file.path(path, '03_annotation_credential/00_intermediate_data'))) {
                 load(file.path(path, '03_annotation_credential/00_intermediate_data', 'annot_all_before_credential.RData'))
                 load(file.path(path, '03_annotation_credential/00_intermediate_data', 'list_peak_group_annotation_concised.RData'))
                 load(file.path(path, '03_annotation_credential/00_intermediate_data', 'list_peak_group_formula.RData'))

                 result_credential_pg <- list_peak_group_annotation_concised
                 result_credential_formula <- list_peak_group_formula

                 annot_all <- addCredentialResult2MetDNA2AnnoClass(annot_all = annot_all,
                                                                   result_credential_pg = result_credential_pg,
                                                                   result_credential_formula = result_credential_formula,
                                                                   is_pred_formula_all = is_pred_formula_all)

                 dir.create(file.path(path_output, '00_intermediate_data'), showWarnings = FALSE, recursive = TRUE)
                 save(annot_all, file = file.path(path_output, '00_intermediate_data', 'annot_all'), version = 2)

               } else {
                 load(file.path(path, '01_result_initial_seed_annotation/00_intermediate_data', 'result_annotation'))
                 load(file.path(path, '02_result_MRN_annotation/00_intermediate_data', 'tags2_after_redundancy_remove'))
                 load(file.path(path, '03_annotation_credential/00_intermediate_data', 'list_peak_group_annotation_concised.RData'))
                 load(file.path(path, '03_annotation_credential/00_intermediate_data', 'list_peak_group_formula.RData'))

                 result_initial_annotation <- result_annotation
                 result_recursive_annotation <- tags2_after_redundancy_remove
                 names(result_recursive_annotation) <- names(result_initial_annotation)
                 result_credential_pg <- list_peak_group_annotation_concised
                 result_credential_formula <- list_peak_group_formula

                 annot_all <- generateMetDNA2AnnoClass(result_initial_annotation = result_initial_annotation)
                 annot_all <- addInitialSeedResult2MetDNA2AnnoClass(annot_all = annot_all,
                                                                    result_initial_annotation = result_initial_annotation)
                 annot_all <- addRecursiveResult2MetDNA2AnnoClass(annot_all = annot_all,
                                                                  result_recursive_annotation = result_recursive_annotation,
                                                                  ...)
                 annot_all <- addCredentialResult2MetDNA2AnnoClass(annot_all = annot_all,
                                                                   result_credential_pg = result_credential_pg,
                                                                   result_credential_formula = result_credential_formula,
                                                                   is_pred_formula_all = is_pred_formula_all)
               }

               rm(list = c('result_annotaion', 'tags2_after_redundancy_remove',
                           'list_peak_group_annotation_concised.RData', 'list_peak_group_formula.RData',
                           'result_initial_annotation', 'result_recursive_annotation',
                           'result_credential_pg', 'result_credential_formula'))
               gc()

               dir.create(file.path(path_output, '00_intermediate_data'), showWarnings = FALSE, recursive = TRUE)
               save(annot_all, file = file.path(path_output, '00_intermediate_data', 'annot_all'), version = 2)

             }


             # browser()
             generateExportTable(annot_all = annot_all,
                                 direction = direction,
                                 tolerance_rt_range = tolerance_rt_range,
                                 path = path,
                                 instrument = instrument,
                                 is_cred_formula_filter = is_cred_formula_filter,
                                 ...)

             if (remove_conflict_seed_from_final_table) {
               removeReplicateSeedAnnotationFromFinalTable(dir_path = path,
                                                           instrument = instrument)
             }

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


             if ('id_merge' %in% list.files(file.path(path_output, '00_intermediate_data'))) {
               cat('\n'); cat('Note: id_merge exsited, and load the old RData\n')
               load(file.path(path_output, '00_intermediate_data', 'id_merge'))
             } else {
               id_merge <- generateResultTable4AllAnnotation(annot_all = annot_all,
                                                             path = path_output,
                                                             thread = thread,
                                                             is_cred_pg_filter = FALSE,
                                                             is_cred_formula_filter = FALSE,
                                                             direction = direction,
                                                             tolerance_rt_range = tolerance_rt_range)

               dir.create(file.path(path_output, '00_intermediate_data'), showWarnings = FALSE, recursive = TRUE)
               save(id_merge, file = file.path(path_output, '00_intermediate_data', 'id_merge'), version = 2)
             }


             # if (instrument %in% c("AgilentDTIMMS", "BrukerTIMS")) {
             #   # export long annotation table
             #   dir.create(file.path(path), showWarnings = FALSE, recursive = TRUE)
             #   readr::write_csv(id_merge, file.path(path_output, 'id_merge_long.csv'))
             # } else {
             #   temp <- id_merge %>% dplyr::select(-c(ccs, ccs_error))
             #   dir.create(file.path(path), showWarnings = FALSE, recursive = TRUE)
             #   readr::write_csv(temp, file.path(path_output, 'id_merge_long.csv'))
             # }


             # id_merge_wide <- id_merge %>%
             #   dplyr::group_by(peak_name) %>%
             #   dplyr::arrange(confidence_level) %>%
             #   dplyr::summarise(source = paste(source, collapse = ';'),
             #                    id = paste(id, collapse = ';'),
             #                    id_zhulab= paste(id, collapse = ';'),
             #                    name = paste(name, collapse = ';'),
             #                    formula = paste(formula, collapse = ';'),
             #                    confidence_level = paste(confidence_level, collapse = ';'),
             #                    smiles = paste(smiles, collapse = ';'),
             #                    inchikey = paste(inchikey, collapse = ';'),
             #                    inchikey1 = paste(inchikey1, collapse = ';'),
             #                    adduct = paste(adduct, collapse = ';'),
             #                    isotope = paste(isotope, collapse = ';'),
             #                    mz_error = paste(mz_error, collapse = ';'),
             #                    rt_error_abs = paste(rt_error_abs, collapse = ';'),
             #                    rt_error_rela = paste(rt_error_rela, collapse = ';'),
             #                    ccs_error = paste(ccs_error, collapse = ';'),
             #                    ms2_score = paste(ms2_score, collapse = ';'),
             #                    cpd_type = paste(cpd_type, collapse = ';'),
             #                    from_peak = paste(from_peak, collapse = ';'),
             #                    from_annotation = paste(from_annotation, collapse = ';'),
             #                    recursive_type = paste(recursive_type, collapse = ';'),
             #                    round = paste(round, collapse = ';'),
             #                    as_seed = paste(as_seed, collapse = ';'),
             #                    as_seed_round = paste(as_seed_round, collapse = ';'),
             #                    total_score = paste(total_score, collapse = ';'),
             #                    is_credential_pg = paste(is_credential_pg, collapse = ';'),
             #                    is_credential_formula = paste(is_credential_formula, collapse = ';')) %>%
             #   dplyr::ungroup()
             #
             # id_merge_wide <- id_merge %>%
             #   dplyr::select(peak_name:ccs) %>%
             #   dplyr::distinct(peak_name, .keep_all = TRUE) %>%
             #   dplyr::left_join(id_merge_wide, by = 'peak_name')
             #
             # if (instrument %in% c("AgilentDTIMMS", "BrukerTIMS")) {
             #   # export long annotation table
             #   dir.create(file.path(path), showWarnings = FALSE, recursive = TRUE)
             #   readr::write_csv(id_merge_wide, file.path(path_output, 'id_merge_wide.csv'))
             # } else {
             #   temp <- id_merge_wide %>% dplyr::select(-c(ccs, ccs_error))
             #   dir.create(file.path(path), showWarnings = FALSE, recursive = TRUE)
             #   readr::write_csv(temp, file.path(path_output, 'id_merge_wide.csv'))
             # }

             # modify id_merge result
             # cat('\n'); cat('Export modified id_merge_result\n')
             # id_merge_modify <- modifyOutputResult(id_merge = id_merge,
             #                                       test_rm_extra_anno_from_ini_seed = test_rm_extra_anno_from_ini_seed)
             #
             # dir.create(file.path(path_output, '00_intermediate_data'), showWarnings = FALSE, recursive = TRUE)
             # save(id_merge_modify, file = file.path(path_output, '00_intermediate_data', 'id_merge_modify'), version = 2)
             #
             # if (instrument == 'IMMS') {
             #   dir.create(file.path(path), showWarnings = FALSE, recursive = TRUE)
             #   readr::write_csv(id_merge_modify, file.path(path_output, 'id_merge_modify_long.csv'))
             # } else {
             #   temp <- id_merge_modify %>% dplyr::select(-c(ccs, ccs_error))
             #   dir.create(file.path(path), showWarnings = FALSE, recursive = TRUE)
             #   readr::write_csv(temp, file.path(path_output, 'id_merge_modify_long.csv'))
             # }

             #
             # id_merge_modify_wide <- id_merge_modify %>%
             #   dplyr::group_by(peak_name) %>%
             #   dplyr::arrange(confidence_level) %>%
             #   dplyr::summarise(source = paste(source, collapse = ';'),
             #                    id = paste(id, collapse = ';'),
             #                    id_zhulab= paste(id, collapse = ';'),
             #                    name = paste(name, collapse = ';'),
             #                    formula = paste(formula, collapse = ';'),
             #                    confidence_level = paste(confidence_level, collapse = ';'),
             #                    smiles = paste(smiles, collapse = ';'),
             #                    inchikey = paste(inchikey, collapse = ';'),
             #                    inchikey1 = paste(inchikey1, collapse = ';'),
             #                    adduct = paste(adduct, collapse = ';'),
             #                    isotope = paste(isotope, collapse = ';'),
             #                    mz_error = paste(mz_error, collapse = ';'),
             #                    rt_error_abs = paste(rt_error_abs, collapse = ';'),
             #                    rt_error_rela = paste(rt_error_rela, collapse = ';'),
             #                    ccs_error = paste(ccs_error, collapse = ';'),
             #                    ms2_score = paste(ms2_score, collapse = ';'),
             #                    cpd_type = paste(cpd_type, collapse = ';'),
             #                    from_peak = paste(from_peak, collapse = ';'),
             #                    from_annotation = paste(from_annotation, collapse = ';'),
             #                    recursive_type = paste(recursive_type, collapse = ';'),
             #                    round = paste(round, collapse = ';'),
             #                    total_score = paste(total_score, collapse = ';'),
             #                    is_credential_pg = paste(is_credential_pg, collapse = ';'),
             #                    is_credential_formula = paste(is_credential_formula, collapse = ';'),
             #                    note = paste(note, collapse = ';')) %>%
             #   dplyr::ungroup() %>%
             #   dplyr::arrange(match(peak_name, unique(id_merge_modify$peak_name)))
             #
             # # keep one formula and adduct for level4
             # peak_level4 <- which(id_merge_modify$confidence_level == 'level4') %>% id_merge_modify$peak_name[.] %>% unique()
             # idx_peak_level4 <- match(peak_level4, id_merge_modify_wide$peak_name)
             #
             # id_merge_modify_wide[idx_peak_level4, ] <- id_merge_modify_wide[idx_peak_level4, ] %>%
             #   dplyr::mutate_at(c('source', 'id', 'id_zhulab', 'name', 'smiles', 'inchikey',
             #                      'inchikey1', 'isotope', 'mz_error', 'rt_error_abs', 'rt_error_rela', 'ccs_error', 'ms2_score',
             #                      'cpd_type', 'from_peak', 'from_annotation', 'recursive_type', 'round', 'total_score',
             #                      'is_credential_pg', 'is_credential_formula'),
             #                    function(x){NA})
             #
             # id_merge_modify_wide$formula[idx_peak_level4] <- match(peak_level4, id_merge_modify$peak_name) %>% id_merge_modify$formula[.]
             # id_merge_modify_wide$adduct[idx_peak_level4] <- match(peak_level4, id_merge_modify$peak_name) %>% id_merge_modify$adduct[.]
             # id_merge_modify_wide$confidence_level[idx_peak_level4] <- match(peak_level4, id_merge_modify$peak_name) %>% id_merge_modify$confidence_level[.]
             #
             # id_merge_modify_wide <- id_merge %>%
             #   dplyr::select(peak_name:ccs) %>%
             #   dplyr::distinct(peak_name, .keep_all = TRUE) %>%
             #   dplyr::left_join(id_merge_modify_wide, by = 'peak_name')
             #
             # if (instrument %in% c("AgilentDTIMMS", "BrukerTIMS")) {
             #   dir.create(file.path(path), showWarnings = FALSE, recursive = TRUE)
             #   readr::write_csv(id_merge_modify_wide, file.path(path_output, 'id_merge_modify_wide.csv'))
             # } else {
             #   temp <- id_merge_modify_wide %>% dplyr::select(-c(ccs, ccs_error))
             #   dir.create(file.path(path), showWarnings = FALSE, recursive = TRUE)
             #   readr::write_csv(temp, file.path(path_output, 'id_merge_modify_wide.csv'))
             # }
             #
             #
             # # export id_merge for service
             # id_merge_modify_service <- id_merge_modify_wide %>%
             #   dplyr::select(peak_name:id, name, formula, adduct, isotope, confidence_level, recursive_type, ms2_score, total_score)
             #
             # if (instrument %in% c("AgilentDTIMMS", "BrukerTIMS")) {
             #   dir.create(file.path(path), showWarnings = FALSE, recursive = TRUE)
             #   readr::write_csv(id_merge_modify_service, file.path(path_output, 'id_merge_modify_service.csv'))
             # } else {
             #   temp <- id_merge_modify_service %>% dplyr::select(-c(ccs))
             #   dir.create(file.path(path), showWarnings = FALSE, recursive = TRUE)
             #   readr::write_csv(temp, file.path(path_output, 'id_merge_modify_service.csv'))
             # }


             cat('\n');cat('The MetDNA2 annotation table has been generated.\n')

           })


#   modifyOutputRes---------------------------------------------------------
#' @title modifyOutputResult
#' @author Zhiwei Zhou
#' @description output id_merge result as level1, level2, level3.1, level3.2, level4
#' @param id_merge
#' @export

setGeneric(name = 'modifyOutputResult',
           def = function(
             id_merge,
             test_rm_extra_anno_from_ini_seed = FALSE
           ){
             # temp_peak_info <- id_merge %>%
             #   dplyr::distinct(peak_name, .keep_all = TRUE)
             peak_list <- unique(id_merge$peak_name)

             id_merge_modify <- pbapply::pblapply(seq_along(peak_list), function(i){
               temp_peak <- peak_list[i]
               temp_peak_result <- id_merge %>%
                 dplyr::filter(peak_name == temp_peak)

               if (is.na(temp_peak_result$id[1])) {
                 result <- temp_peak_result %>%
                   dplyr::select(peak_name:is_credential_formula) %>%
                   dplyr::mutate(note = NA)

                 return(result)
               }

               # keep the highest level annotation
               temp_highest_level <- temp_peak_result$confidence_level %>% sort() %>% .[1]
               temp_peak_result <- temp_peak_result %>% dplyr::filter(confidence_level == temp_highest_level)

               # keep annotation passing the peak grouping
               if (!any(temp_peak_result$is_credential_pg)) {
                 temp_peak_result <- temp_peak_result %>%
                   dplyr::filter(dplyr::row_number() == 1)

                 result <- temp_peak_result %>% dplyr::select(peak_name:rt)
                 temp_peak_result <- temp_peak_result %>%
                   dplyr::select(-c('peak_name', 'mz', 'rt', 'ccs')) %>%
                   dplyr::mutate_all(function(x){NA})

                 result <- result %>%
                   dplyr::bind_cols(temp_peak_result) %>%
                   dplyr::select(peak_name:is_credential_formula) %>%
                   dplyr::mutate(note = NA)

                 return(result)
               } else {
                 temp_peak_result <- temp_peak_result %>%
                   dplyr::filter(is_credential_pg)
               }

               # assign exact confidence level
               temp_peak_result <- temp_peak_result %>%
                 dplyr::mutate(confidence_level_new = dplyr::case_when(
                   !is_credential_formula ~ 'level4',
                   confidence_level == 'level1' ~ 'level1',
                   confidence_level == 'level2' ~ 'level2',
                   confidence_level == 'level3' & cpd_type == 'known_known' ~ 'level3.1',
                   confidence_level == 'level3' & cpd_type != 'known_known' ~ 'level3.2'
                 ))


               # if one feature only has annotation > level4, remove level 4 annotation
               # if one feature only has level4 annotation, added it into notes
               if (!all(temp_peak_result$confidence_level_new == 'level4')) {
                 result <- temp_peak_result %>%
                   dplyr::filter(confidence_level_new != 'level4') %>%
                   dplyr::mutate(confidence_level = confidence_level_new) %>%
                   dplyr::select(peak_name:is_credential_formula) %>%
                   dplyr::mutate(note = NA)

                 return(result)
               } else {
                 # add level 4 annotation result into note
                 result <- temp_peak_result %>%
                   dplyr::arrange(confidence_level) %>%
                   dplyr::mutate(note = paste0('id{', id, '}',
                                               'source{', source, '}',
                                               'name{', name, '}',
                                               'confidence_level{',confidence_level_new, '}')) %>%
                   dplyr::mutate(confidence_level = confidence_level_new,
                                 formula = best_pred_formula,
                                 adduct = best_pred_formula_adduct) %>%
                   dplyr::mutate_at(c('source', 'id', 'id_zhulab', 'name', 'smiles', 'inchikey',
                                      'inchikey1', 'isotope', 'mz_error', 'rt_error_abs', 'rt_error_rela', 'ccs_error', 'ms2_score',
                                      'cpd_type', 'from_peak', 'from_annotation', 'recursive_type', 'round', 'total_score',
                                      "as_seed", "as_seed_round", 'is_credential_pg', 'is_credential_formula'),
                                    function(x){NA}) %>%
                   dplyr::select(peak_name:is_credential_formula, note)

                 return(result)

               }


             })


             result_id_merge_modify <- id_merge_modify %>% dplyr::bind_rows()

             # remove extra annotation records if this feature has adduct or isotope annotations are directly from initial seeds
             if (test_rm_extra_anno_from_ini_seed) {
               peak_initial_annotation <- result_id_merge_modify %>%
                 dplyr::filter(confidence_level %in% c('level1', 'level2')) %>%
                 dplyr::pull(peak_name) %>%
                 unique()

               idx_initial_seed_adducts <- which((result_id_merge_modify$from_peak %in% peak_initial_annotation) &
                                                   result_id_merge_modify$round == 2 &
                                                   result_id_merge_modify$recursive_type %in% c('adductAnnotation', 'isotopeAnnotation'))
               temp_peaks <- result_id_merge_modify$peak_name[idx_initial_seed_adducts] %>% unique()
               idx_remove <- which(result_id_merge_modify$peak_name %in% temp_peaks)
               idx_remove <- setdiff(idx_remove, idx_initial_seed_adducts)
               if (length(idx_remove) > 0) {
                 result_id_merge_modify <- result_id_merge_modify[-idx_remove,]
               }
             }

             # remove annotation without formula prediction
             idx_remove <- which((result_id_merge_modify$confidence_level) == 'level4' & is.na(result_id_merge_modify$formula))

             if (length(idx_remove) > 0) {
               result_id_merge_modify <- result_id_merge_modify[-idx_remove,]
             }

             return(result_id_merge_modify)

           })

# generateMetDNA2AnnotationResultWithCredential --------------------------------

# #' @title generateMetDNA2AnnotationResultWithCredential
# #' @author Zhiwei ZHou
# #' @param path
# #' @param type_order
#
# # generateMetDNA2AnnotationResultWithCredential(path = path,
# #                                               type_order = c('level1', 'level2', 'level3'))
#
# setGeneric(name = 'generateMetDNA2AnnotationResultWithCredential',
#            def = function(
#              path = '.',
#              type_order = c('level1', 'level2', 'level3.1', 'level3.2', 'level4')
#            ){
#              cat('Generate MetDNA2 annotation result with credential\n\n')
#
#              path_output <- file.path(path, '00_annotation_table')
#
#              if (!all(c('id_merge', 'annot_all') %in% list.files(file.path(path_output, '00_intermediate_data')))) {
#                stop('Please run generateMetDNA2AnnotationResult first to generate annot_all and id_merge file\n')
#              }
#
#              # load data
#              load(file.path(path_output, '00_intermediate_data', 'annot_all'))
#              load(file.path(path_output, '00_intermediate_data', 'id_merge'))
#
#              # convert id_merge to list
#              temp_idx <- lapply(unique(id_merge$peak_name), function(x){
#                which(id_merge$peak_name == x)
#              })
#
#              names(temp_idx) <- unique(id_merge$peak_name)
#
#              id_merge_list <- vector(mode = 'list', length = length(temp_idx))
#              for (i in seq_along(temp_idx)) {
#                id_merge_list[[i]] <- id_merge[temp_idx[[i]],]
#              }
#
#              names(id_merge_list) <- names(temp_idx)
#
#
#
#
#              pg_annotation_all <- pbapply::pblapply(annot_all, function(x){
#                x@credential_annotation$peak_group
#              })
#
#              pg_annotation_all <- pg_annotation_all %>% dplyr::bind_rows()
#
#              # extract all annotations: peak group annotation, and formula prediction
#              #   for all reserved recursive annotations (pg_annotation & formula_prediction), retrieve all effective peak group
#              #      for all effective peak group, retrieve all annotations, and assign annotation result of base peaks.
#              #        the confidence level: level3, source: annotation_credential
#              id_effective <- id_merge %>%
#                dplyr::filter(!is.na(source)) %>%
#                dplyr::filter(is_credential_pg & is_credential_formula) %>%
#                dplyr::mutate(pg = paste(peak_name, adduct, sep = '_'))
#
#              pg_effective <- id_effective %>%
#                dplyr::distinct(peak_name, adduct) %>%
#                dplyr::mutate(pg = paste(peak_name, adduct, sep = '_')) %>%
#                dplyr::filter(pg %in% unique(pg_annotation_all$included_peak_group))
#
#              id_effective_with_credential <- pbapply::pblapply(seq_along(pg_effective$pg), function(i){
#                # cat(i, ' ')
#                temp_pg_basepeak <- pg_effective$peak_name[i]
#                temp_pg_basepeak_adduct <- pg_effective$adduct[i]
#                temp_pg_name <- pg_effective$pg[i]
#                temp_pg_annotation <- pg_annotation_all %>% dplyr::filter(included_peak_group == temp_pg_name)
#
#                temp_recursive_annotation <- id_effective %>%
#                  dplyr::filter(peak_name == temp_pg_basepeak &
#                                  adduct == temp_pg_basepeak_adduct) %>%
#                  dplyr::rename(peak_group = pg) %>%
#                  dplyr::mutate(insource_frag = NA) %>%
#                  dplyr::select(peak_name, source:isotope, insource_frag, peak_group, ms2_score, cpd_type, recursive_type) %>%
#                  dplyr::mutate_all(as.character)
#
#                temp_col <- colnames(temp_recursive_annotation)
#
#                # assign base peak ID inforamtion to its related peaks
#                result_pg_annotation <- lapply(seq_along(temp_recursive_annotation$id), function(j){
#                  # cat(j, ' ')
#                  temp_pg_annotation_adduct_nl <- temp_pg_annotation %>%
#                    dplyr::filter(type_annotation %in% c('adductAnnotation', 'neutralLossAnnotation')) %>%
#                    dplyr::select(peak_name, label, included_peak_group) %>%
#                    dplyr::rename(adduct = label,
#                                  peak_group = included_peak_group)
#
#                  temp_pg_annotation_isotope <- temp_pg_annotation %>%
#                    dplyr::filter(type_annotation %in% 'isotopeAnnotation') %>%
#                    dplyr::select(peak_name, label, included_peak_group) %>%
#                    dplyr::rename(isotope = label,
#                                  peak_group = included_peak_group)
#
#                  temp_pg_annotation_isf <- temp_pg_annotation %>%
#                    dplyr::filter(type_annotation %in% 'isfAnnotation') %>%
#                    dplyr::select(peak_name, label, included_peak_group) %>%
#                    dplyr::rename(insource_frag = label,
#                                  peak_group = included_peak_group)
#
#                  # for adductAnnotation and neutral loss annotation, assign isotope as [M]
#                  temp_pg_annotation_adduct_nl <- temp_pg_annotation_adduct_nl %>%
#                    dplyr::mutate(isotope = '[M]',
#                                  insource_frag = NA,
#                                  ms2_score = NA,
#                                  recursive_type = NA) %>%
#                    dplyr::mutate(source = 'annotation_credential',
#                                  id = temp_recursive_annotation$id[j],
#                                  id_zhulab = temp_recursive_annotation$id_zhulab[j],
#                                  name = temp_recursive_annotation$name[j],
#                                  formula = temp_recursive_annotation$formula[j],
#                                  confidence_level = 'level3',
#                                  smiles = temp_recursive_annotation$smiles[j],
#                                  inchikey = temp_recursive_annotation$inchikey[j],
#                                  inchikey1 = temp_recursive_annotation$inchikey1[j],
#                                  cpd_type = temp_recursive_annotation$cpd_type[j]) %>%
#                    dplyr::select(temp_col) %>%
#                    dplyr::mutate_all(as.character)
#
#
#                  # for isotopeAnnotation, assign adduct as NA
#                  temp_pg_annotation_isotope <- temp_pg_annotation_isotope %>%
#                    dplyr::mutate(adduct = NA,
#                                  insource_frag = NA,
#                                  ms2_score = NA,
#                                  recursive_type = NA) %>%
#                    dplyr::mutate(source = 'annotation_credential',
#                                  id = temp_recursive_annotation$id[j],
#                                  id_zhulab = temp_recursive_annotation$id_zhulab[j],
#                                  name = temp_recursive_annotation$name[j],
#                                  formula = temp_recursive_annotation$formula[j],
#                                  confidence_level = 'level3',
#                                  smiles = temp_recursive_annotation$smiles[j],
#                                  inchikey = temp_recursive_annotation$inchikey[j],
#                                  inchikey1 = temp_recursive_annotation$inchikey1[j],
#                                  cpd_type = temp_recursive_annotation$cpd_type[j]) %>%
#                    dplyr::select(temp_col) %>%
#                    dplyr::mutate_all(as.character)
#
#
#                  # for IsfAnnotation, assign adduct as NA, isotope as [M]
#                  temp_pg_annotation_isf <- temp_pg_annotation_isf %>%
#                    dplyr::mutate(isotope = '[M]',
#                                  adduct = NA,
#                                  ms2_score = NA,
#                                  recursive_type = NA) %>%
#                    dplyr::mutate(source = 'annotation_credential',
#                                  id = temp_recursive_annotation$id[j],
#                                  id_zhulab = temp_recursive_annotation$id_zhulab[j],
#                                  name = temp_recursive_annotation$name[j],
#                                  formula = temp_recursive_annotation$formula[j],
#                                  confidence_level = 'level3',
#                                  smiles = temp_recursive_annotation$smiles[j],
#                                  inchikey = temp_recursive_annotation$inchikey[j],
#                                  inchikey1 = temp_recursive_annotation$inchikey1[j],
#                                  cpd_type = temp_recursive_annotation$cpd_type[j]) %>%
#                    dplyr::select(temp_col) %>%
#                    dplyr::mutate_all(as.character)
#
#
#                  result_pg_annotation <- temp_pg_annotation_adduct_nl %>%
#                    dplyr::bind_rows(temp_pg_annotation_isotope) %>%
#                    dplyr::bind_rows(temp_pg_annotation_isf)
#
#                  return(result_pg_annotation)
#                })
#
#                result_pg_annotation <- result_pg_annotation %>% dplyr::bind_rows()
#
#                # merge all recursive annotation and pg annotation, and keep unique annotation
#                result <- temp_recursive_annotation %>%
#                  dplyr::bind_rows(result_pg_annotation) %>%
#                  dplyr::arrange(match(source, c('initial_seed', 'recursive_annotation', 'annotation_credential')),
#                                 confidence_level) %>%
#                  dplyr::distinct(peak_name, id, .keep_all = TRUE)
#
#                return(result)
#              })
#
#              id_effective_with_credential <- id_effective_with_credential %>% dplyr::bind_rows()
#
#              # browser()
#
#              peak_info <- id_merge %>% dplyr::distinct(peak_name, .keep_all = TRUE) %>% dplyr::select(peak_name:rt)
#              id_effective_with_credential <- peak_info %>% dplyr::left_join(id_effective_with_credential, by = 'peak_name')
#
#              dir.create(path_output, showWarnings = FALSE, recursive = TRUE)
#              readr::write_csv(id_effective_with_credential,
#                               path = file.path(path_output, 'id_effective_with_credential_result.csv'))
#
#              id_effective_with_credential_wide <- id_effective_with_credential %>%
#                dplyr::group_by(peak_name) %>%
#                dplyr::arrange(confidence_level) %>%
#                dplyr::summarise(source = paste(source, collapse = ';'),
#                                 id = paste(id, collapse = ';'),
#                                 id_zhulab= paste(id, collapse = ';'),
#                                 name = paste(name, collapse = ';'),
#                                 formula = paste(formula, collapse = ';'),
#                                 confidence_level = paste(confidence_level, collapse = ';'),
#                                 smiles = paste(smiles, collapse = ';'),
#                                 inchikey = paste(inchikey, collapse = ';'),
#                                 inchikey1 = paste(inchikey1, collapse = ';'),
#                                 adduct = paste(adduct, collapse = ';'),
#                                 isotope = paste(isotope, collapse = ';'),
#                                 insource_frag = paste(insource_frag, collapse = ';'),
#                                 peak_group = paste(peak_group, collapse = ';'),
#                                 ms2_score = paste(ms2_score, collapse = ';'),
#                                 cpd_type = paste(cpd_type, collapse = ';')) %>%
#                dplyr::ungroup()
#
#              id_effective_with_credential_wide <- peak_info %>% dplyr::left_join(id_effective_with_credential_wide, by = 'peak_name')
#
#              readr::write_csv(id_effective_with_credential_wide,
#                               path = file.path(path_output, 'id_effective_with_credential_result_wide.csv'))
#
#            })
#

# generateExportTable ----------------------------------------------------------

#' @title generateExportTable
#' @author Zhiwei Zhou
#' @param annot_all
#' @param path '.'
#' @param lib  c('zhuMetLib', 'zhuMetLib_orbitrap', 'fiehnHilicLib')
#' @param direction 'reverse', 'forward'
#' @param instrument c("SciexTripleTOF", "AgilentQTOF", "BrukerQTOF", "ThermoOrbitrap", 'ThermoExploris', "AgilentDTIMMS", "BrukerTIMS", 'WatersQTOF','WatersTWIMMS')
#' @param tolerance_rt_range 30
#' @param is_cred_formula_filter Default: FALSE
#' @param candidate_num 5
#' @export


# generateResultTable4AllAnnotation(annot_all = annot_all, path = path, lib = 'zhuMetLib', )
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

setGeneric(name = 'generateExportTable',
           function(annot_all,
                    path = '.',
                    direction = c('reverse', 'forward'),
                    lib =  c('zhuMetLib', 'zhuMetLib_orbitrap', 'fiehnHilicLib', 'zhuRPLib'),
                    instrument = c("SciexTripleTOF", "AgilentQTOF", "BrukerQTOF", "ThermoOrbitrap",
                                   'ThermoExploris', "AgilentDTIMMS", "BrukerTIMS", 'WatersQTOF','WatersTWIMMS'),
                    tolerance_rt_range = 30,
                    is_cred_formula_filter = FALSE,
                    candidate_num = 10,
                    test_evaluation = c('No', '200STD', '46STD'),
                    ...){

             # browser()
             path_output <- file.path(path, '00_annotation_table')
             load(file.path(path, '01_result_initial_seed_annotation/00_intermediate_data', 'ms1_data'))

             # extract credentialed peak group ---------------------------------
             list_credential_pg <- lapply(annot_all, function(x){
               x@credential_annotation$peak_group
             })

             peak_group_id_table <- list_credential_pg %>%
               dplyr::bind_rows() %>%
               dplyr::select(peak_group_id, included_peak_group) %>%
               dplyr::distinct(peak_group_id, .keep_all = TRUE)

             base_num <- peak_group_id_table$peak_group_id %>% stringr::str_extract(pattern = '\\d+') %>% as.numeric()
             base_peak <- peak_group_id_table$included_peak_group %>% stringr::str_split(pattern = '_\\[|_M', simplify = TRUE) %>% .[,1]
             # base_adduct <- peak_group_id_table$included_peak_group %>% stringr::str_extract(pattern = '\\[\\d*M.+|M\\+|M\\-')
             base_adduct <- peak_group_id_table$included_peak_group %>% stringr::str_extract(pattern = '(\\[\\d*M.+)')

             peak_all <- lapply(annot_all, function(x){
               x@peak_info %>% tibble::as_tibble()
             }) %>% dplyr::bind_rows()
             base_mz <- match(base_peak, peak_all$name) %>% peak_all$mz[.]
             base_neutral_mass <- mapply(function(x, y){
               convertMz2Adduct(base_mz = x, base_adduct = y, adduct_list = y)$exact_mass
             },
             x = base_mz,
             y = base_adduct)
             # convertMz2Adduct(base_mz = 220.1176, base_adduct = '[M+H]+', adduct_list = '[M-H]-')$exact_mass
             # convertMz2Adduct(base_mz = base_mz[4], base_adduct = base_adduct[4], adduct_list = base_adduct[4])$exact_mass
             # convertMz2Adduct(base_mz = 664.1144, base_adduct = '[M]+', adduct_list = '[M]+')$exact_mass

             peak_group_id_table <- peak_group_id_table %>%
               dplyr::mutate(base_peak = base_peak,
                             base_adduct = base_adduct,
                             base_num = base_num,
                             base_mz = base_mz,
                             base_neutral_mass = base_neutral_mass) %>%
               dplyr::arrange(base_num)

             dir.create(file.path(path_output, '00_intermediate_data'), showWarnings = FALSE, recursive = TRUE)
             save(peak_group_id_table, file = file.path(path_output, '00_intermediate_data', 'peak_group_id_table'), version = 2)



             cat('Export peak group table ...\n')
             list_peak_group <- generateExportTablePeakGroup(annot_all = annot_all, peak_group_id_table = peak_group_id_table)

             dir.create(file.path(path_output, '00_intermediate_data'), showWarnings = FALSE, recursive = TRUE)
             save(list_peak_group, file = file.path(path_output, '00_intermediate_data', 'list_peak_group'), version = 2)

             table_peak_group <- list_peak_group %>%
               dplyr::bind_rows() %>%
               dplyr::select(peak_name:ccs, isotope, adduct, peak_group_id, base_peak, rela_peaks) %>%
               dplyr::mutate(pg_num = as.numeric(stringr::str_extract(peak_group_id, pattern = '\\d+'))) %>%
               dplyr::arrange(pg_num) %>%
               dplyr::distinct(peak_name, peak_group_id, .keep_all = TRUE) %>%
               dplyr::select(-pg_num) %>%
               dplyr::select(peak_group_id, dplyr::everything()) %>%
               dplyr::mutate(mz = round(mz, 4),
                             rt = round(rt, 1)) %>%
               tidyr::replace_na(list('adduct' = ''))

             # replace base_peak and neutral mass
             table_peak_group$base_peak <- match(table_peak_group$peak_group_id, peak_group_id_table$peak_group_id) %>%
               peak_group_id_table$base_peak[.]
             table_peak_group$neutral_mass <- match(table_peak_group$peak_group_id, peak_group_id_table$peak_group_id) %>%
               peak_group_id_table$base_neutral_mass[.]
             table_peak_group <- table_peak_group %>%
               dplyr::select(peak_group_id:base_peak, neutral_mass, dplyr::everything()) %>%
               dplyr::mutate(neutral_mass = round(neutral_mass, 4)) %>%
               dplyr::rename(num_peaks = rela_peaks)

             if (!(instrument %in% c("AgilentDTIMMS", "BrukerTIMS", "WatersTWIMMS"))) {
               table_peak_group <- table_peak_group %>%
                 dplyr::select(-ccs)
             }

             readr::write_csv(table_peak_group, file = file.path(path_output, 'table2_peak_group.csv'))



             cat('Export identification table ...\n')
             list_identification <- generateExporTableIden(annot_all = annot_all,
                                                           peak_group_id_table = peak_group_id_table,
                                                           direction = direction,
                                                           tolerance_rt_range = tolerance_rt_range,
                                                           lib = lib,
                                                           test_evaluation = test_evaluation)

             dir.create(file.path(path_output, '00_intermediate_data'), showWarnings = FALSE, recursive = TRUE)
             save(list_identification, file = file.path(path_output, '00_intermediate_data', 'list_identification'), version = 2)

             table_identification <- list_identification %>%
               dplyr::bind_rows() %>%
               dplyr::select(peak_name:ccs, id:inchikey, isotope, adduct, mz_error, rt_error_abs, rt_error_rela, ccs_error,
                             ms2_score, recursive_type, peak_group_id, base_peak, rela_peaks, cons_formula_pred) %>%
               dplyr::rename(id_kegg = id,
                             iden_type = recursive_type)


             # retrieve rt errors for identification from recursive
             table_identification <- replaceRtError(table_identification = table_identification, path = path)
             table_identification <- calculateExportTableScore(table_identification = table_identification,
                                                               instrument = instrument,
                                                               ...)

             # table_identification <- calculateExportTableScore(table_identification = table_identification,
             #                                                   instrument = instrument, mz_tol = 25, rt_tol = 30, ccs_tol = 3)

             if (!(instrument %in% c("AgilentDTIMMS", "BrukerTIMS", "WatersTWIMMS"))) {
               table_identification <- table_identification %>%
                 dplyr::select(-c('ccs', 'ccs_error'))
             }

             if (is_cred_formula_filter) {
               table_identification <- table_identification %>%
                 dplyr::mutate(confidence_level = dplyr::case_when(cons_formula_pred~confidence_level,
                                                                   !cons_formula_pred~'level4'))
               #
               # table_identification <- table_identification %>%
               #   dplyr::filter(cons_formula_pred)
             }

             # keep top N candidates
             table_identification <- table_identification %>%
               dplyr::group_by(peak_name) %>%
               dplyr::arrange(dplyr::desc(total_score)) %>%
               dplyr::top_n(candidate_num, total_score) %>%
               dplyr::ungroup() %>%
               dplyr::rename(num_peaks = rela_peaks)


             # add stereo_isomers
             if (test_evaluation == 'No') {
               data("cpd_emrn", envir = environment())
               temp_idx <- table_identification$inchikey %>%
                 splitInchiKey() %>%
                 .[['inchikey1']] %>%
                 match(cpd_emrn$inchikey1)

               table_identification <- table_identification %>%
                 dplyr::mutate(stereo_isomer_id = cpd_emrn$id_kegg_synonyms[temp_idx],
                               stereo_isomer_name = cpd_emrn$name_synonyms[temp_idx])

               # replace id and name for in-silico compounds whose included in KEGG but not in MRN
               known_insilico_cpd_set <- cpd_emrn %>% dplyr::filter(!is.na(id_kegg) & !stringr::str_detect(id, 'C\\d+'))
               idx1 <- which(table_identification$id_kegg %in% known_insilico_cpd_set$id)
               idx2 <- match(table_identification$id_kegg[idx1], known_insilico_cpd_set$id)
               table_identification$id_kegg[idx1] <- known_insilico_cpd_set$id_kegg[idx2]
               table_identification$name[idx1] <- known_insilico_cpd_set$name[idx2]
             } else {
               table_identification <- table_identification %>%
                 dplyr::mutate(stereo_isomer_id = NA,
                               stereo_isomer_name = NA)
             }

             dir.create(file.path(path_output, '00_intermediate_data'), showWarnings = FALSE, recursive = TRUE)
             save(table_identification, file = file.path(path_output, '00_intermediate_data', 'table_identification'), version = 2)


             # generate wide identification table ------------------------------
             table_identification_wide <- table_identification %>%
               dplyr::group_by(peak_name) %>%
               dplyr::summarise(peak_name = peak_name[1],
                                mz = mz[1],
                                rt = rt[1],
                                id_kegg = paste(id_kegg, collapse = ';'),
                                id_zhulab = paste(id_zhulab, collapse = ';'),
                                name = paste(name, collapse = ';'),
                                formula = paste(formula, collapse = ';'),
                                confidence_level = paste(confidence_level, collapse = ';'),
                                smiles = paste(smiles, collapse = ';'),
                                inchikey = paste(inchikey, collapse = ';'),
                                isotope = paste(isotope, collapse = ';'),
                                adduct = paste(adduct, collapse = ';'),
                                total_score = paste(total_score, collapse = ';'),
                                mz_error = paste(mz_error, collapse = ';'),
                                rt_error_abs = paste(rt_error_abs, collapse = ';'),
                                rt_error_rela = paste(rt_error_rela, collapse = ';'),
                                ms2_score = paste(ms2_score, collapse = ';'),
                                iden_score = paste(iden_score, collapse = ';'),
                                iden_type = paste(iden_type, collapse = ';'),
                                peak_group_id = paste(peak_group_id, collapse = ';'),
                                base_peak = paste(base_peak, collapse = ';'),
                                num_peaks = paste(num_peaks, collapse = ';'),
                                cons_formula_pred = paste(cons_formula_pred, collapse = ';'),
                                stereo_isomer_id = paste(unique(stereo_isomer_id), collapse = ';'),
                                stereo_isomer_name = paste(unique(stereo_isomer_name), collapse = ';')) %>%
               dplyr::ungroup() %>%
               dplyr::arrange(mz)

             # join the quantification table
             quanti_table <- ms1_data$info %>%
               dplyr::bind_cols(ms1_data$subject) %>%
               dplyr::mutate(mz = round(mz, 4),
                             rt = round(rt, 1),
                             ccs = round(ccs, 1))

             if (!(instrument %in% c("AgilentDTIMMS", "BrukerTIMS", "WatersTWIMMS"))) {
               quanti_table <- quanti_table %>%
                 dplyr::select(-c('ccs'))
             }


             table_identification_wide_quanti <- table_identification_wide %>%
               dplyr::full_join(quanti_table, by = c('peak_name' = 'name')) %>%
               dplyr::rename(mz = mz.y,
                             rt = rt.y) %>%
               dplyr::select(-c('mz.x', 'rt.x')) %>%
               dplyr::select(peak_name, mz, rt, dplyr::everything())
             # %>%
             #   dplyr::mutate(
             #     dplyr::across(dplyr::everything(), ~tidyr::replace_na(.x, ''))
             #   )

             readr::write_csv(table_identification_wide_quanti, file = file.path(path_output, 'table1_identification.csv'))

             cat('Export unique index feature table ...\n')
             table_identification_unique <- table_identification
             temp_peak_name <- table_identification_unique$peak_name %>%
               rename4UniqueIndex()

             # add quanti table for table_identification_unique
             table_identification_unique <- table_identification %>%
               dplyr::left_join(quanti_table, by = c('peak_name' = 'name')) %>%
               dplyr::rename(mz = mz.y,
                             rt = rt.y) %>%
               dplyr::select(-c('mz.x', 'rt.x')) %>%
               dplyr::select(peak_name, mz, rt, dplyr::everything())

             table_identification_unique$peak_name <- temp_peak_name

             readr::write_csv(table_identification_unique, file = file.path(path_output, 'table3_identification_pair.csv'))


             cat('Result table export has done!')

           })


#   generateExportTablePeakGroup ----------------------------------------------
#' @title generateExportTablePeakGroup
#' @author Zhiwei Zhou
#' @param annot_all
#' @param peak_group_id_table
#' @export

setGeneric(name = 'generateExportTablePeakGroup',
           def = function(
             annot_all,
             peak_group_id_table
           ){
             list_credential_pg <- lapply(annot_all, function(x){
               x@credential_annotation$peak_group
             })

             temp_credential_pg_all <- list_credential_pg %>% dplyr::bind_rows()

             idx_pg <- mapply(function(x, y){
               which(temp_credential_pg_all$peak_group_id == x)
             },
             x = peak_group_id_table$peak_group_id,
             y = peak_group_id_table$base_peak)

             annot_cred_peaks <- lapply(seq_along(idx_pg), function(i){
               # cat(i, ' ')
               temp_idx <- idx_pg[[i]]

               template <- tibble::tibble(peak_name = as.character(),
                                          mz = numeric(),
                                          rt = numeric(),
                                          ccs = numeric(),
                                          source = character(),
                                          id = character(),
                                          id_zhulab = character(),
                                          name = character(),
                                          formula = character(),
                                          confidence_level = character(),
                                          smiles = character(),
                                          inchikey = character(),
                                          inchikey1 = character(),
                                          adduct = character(),
                                          isotope = character(),
                                          mz_error = numeric(),
                                          rt_error_abs = numeric(),
                                          rt_error_rela = numeric(),
                                          ccs_error = numeric(),
                                          ms2_score = numeric(),
                                          ms2_matched_frag = numeric(),
                                          cpd_type = character(),
                                          from_peak = character(),
                                          from_annotation = character(),
                                          recursive_type = character(),
                                          as_seed = logical(),
                                          as_seed_round = numeric(),
                                          round = numeric(),
                                          total_score = numeric(),
                                          peak_group_id = character(),
                                          base_peak = character(),
                                          rela_peaks = numeric(),
                                          cons_formula_pred = logical())

               if (length(temp_idx) == 0) {
                 result <- template
                 return(template)
               }

               temp_credential_peaks <- temp_credential_pg_all[temp_idx,] %>%
                 dplyr::select(peak_name:rt, type_annotation:annotated_peaks)

               temp_isotope <- temp_credential_peaks %>%
                 dplyr::filter(type_annotation == 'isotopeAnnotation')

               if (nrow(temp_isotope) > 0) {
                 temp_isotope <- temp_isotope %>%
                   # tidyr::pivot_wider(names_from = 'type_annotation', values_from = 'label') %>%
                   dplyr::rename(adduct = annotation,
                                 base_peak = included_peak_group,
                                 rela_peaks = annotated_peaks) %>%
                   dplyr::mutate(ccs = NA,
                                 source = 'annotation_credential',
                                 id = NA,
                                 id_zhulab = NA,
                                 name = NA,
                                 formula = NA,
                                 confidence_level = NA,
                                 smiles = NA,
                                 inchikey = NA,
                                 inchikey1 = NA,
                                 mz_error = NA,
                                 rt_error_abs = NA,
                                 rt_error_rela = NA,
                                 ccs_error = NA,
                                 ms2_score = NA,
                                 ms2_matched_frag = NA,
                                 cpd_type = NA,
                                 from_peak = NA,
                                 from_annotation = NA,
                                 recursive_type = NA,
                                 as_seed = NA,
                                 as_seed_round = NA,
                                 round = NA,
                                 total_score = NA,
                                 cons_formula_pred = NA) %>%
                   dplyr::select( "peak_name", "mz", "rt", "ccs", "source", "id", "id_zhulab", "name", "formula",
                                  "confidence_level", "smiles", "inchikey", "inchikey1", "adduct", "isotope", "mz_error",
                                  "rt_error_abs", "rt_error_rela", "ccs_error", "ms2_score", "ms2_matched_frag",
                                  "cpd_type", "from_peak", "from_annotation", "recursive_type", "as_seed", "as_seed_round", "round",
                                  "total_score", "peak_group_id", "base_peak", "rela_peaks","cons_formula_pred")
               } else {
                 temp_isotope <- template
               }

               temp_adduct <- temp_credential_peaks %>%
                 dplyr::filter(type_annotation != 'isotopeAnnotation')

               if (nrow(temp_adduct) > 0) {
                 temp_adduct <- temp_adduct %>%
                   dplyr::rename(adduct = annotation,
                                 base_peak = included_peak_group,
                                 rela_peaks = annotated_peaks) %>%
                   dplyr::mutate(ccs = NA,
                                 source = 'annotation_credential',
                                 id = NA,
                                 id_zhulab = NA,
                                 name = NA,
                                 formula = NA,
                                 confidence_level = NA,
                                 smiles = NA,
                                 inchikey = NA,
                                 inchikey1 = NA,
                                 mz_error = NA,
                                 rt_error_abs = NA,
                                 rt_error_rela = NA,
                                 ccs_error = NA,
                                 ms2_score = NA,
                                 ms2_matched_frag = NA,
                                 cpd_type = NA,
                                 from_peak = NA,
                                 from_annotation = NA,
                                 recursive_type = NA,
                                 as_seed = NA,
                                 as_seed_round = NA,
                                 round = NA,
                                 total_score = NA,
                                 cons_formula_pred = NA) %>%
                   dplyr::select( "peak_name", "mz", "rt", "ccs", "source", "id", "id_zhulab", "name", "formula",
                                  "confidence_level", "smiles", "inchikey", "inchikey1", "adduct", "isotope", "mz_error",
                                  "rt_error_abs", "rt_error_rela", "ccs_error", "ms2_score", "ms2_matched_frag",
                                  "cpd_type", "from_peak", "from_annotation", "recursive_type", "as_seed", "as_seed_round", "round",
                                  "total_score", "peak_group_id", "base_peak", "rela_peaks","cons_formula_pred")
               } else {
                 temp_adduct <- template
               }

               result <- temp_adduct %>%
                 dplyr::bind_rows(temp_isotope) %>%
                 dplyr::arrange(mz)

               return(result)

             })

             names(annot_cred_peaks) <- peak_group_id_table$peak_group_id
             return(annot_cred_peaks)

           })


#   gemerateExportTableIden ---------------------------------------------------

#' @title generateExporTableIden
#' @author Zhiwei Zhou
#' @param annot_all

setGeneric(name = 'generateExporTableIden',
           def = function(
             annot_all,
             peak_group_id_table,
             direction = c('reverse', 'forward'),
             tolerance_rt_range = 30,
             lib = 'zhuMetLib',
             test_evaluation = c('No', '200STD', '46STD'),
             ...
           ){
             # load some intermediate data
             if (lib %in% c("zhuMetLib", 'zhuMetLib_orbitrap')) {
               data("zhuMetlib", envir = environment())
               cpd_info <- zhuMetlib$meta$compound
             }

             if (lib == "fiehnHilicLib") {
               data("fiehnHilicLib", envir = environment())
               cpd_info <- fiehnHilicLib$meta$compound
             }

             if (lib == 'zhuRPLib') {
               data("zhuRPlib", envir = environment())
               cpd_info <- zhuRPlib$meta$compound
             }

             temp_colname_ms2 <- ifelse(direction == 'forward', 'ms2_score_forward', 'ms2_score_reverse')

             idx <- match(peak_group_id_table$base_peak, names(annot_all))

             annot_base_peaks <- lapply(seq_along(idx), function(i){
               # cat(i, ' ')
               temp_base_peak <- peak_group_id_table$base_peak[i]
               temp_base_peak_adduct <- peak_group_id_table$base_adduct[i]
               temp_peak_group_id <- peak_group_id_table$peak_group_id[i]

               temp_idx <- idx[i]
               temp_annot_all <- annot_all[[temp_idx]]

               peak_info <- tibble::tibble(peak_name = temp_annot_all@peak_info$name,
                                           mz = temp_annot_all@peak_info$mz,
                                           rt = temp_annot_all@peak_info$rt,
                                           ccs = temp_annot_all@peak_info$ccs)

               temp_initial_seed <- temp_annot_all@initial_seed_annotation

               # modify formats for initial_seed_annotation and recursive_annotation
               if (nrow(temp_initial_seed) > 0) {
                 temp_initial_seed <- temp_initial_seed %>%
                   # dplyr::select(id:adduct, mz_error:rt_error_rela, ms2_score_forward, cpd_type) %>%
                   dplyr::select(id:adduct, isotope, mz_error:ccs_error, temp_colname_ms2, ms2_matched_frag, cpd_type) %>%
                   dplyr::mutate(from_peak = as.character(NA),
                                 from_annotation = as.character(NA),
                                 round = as.numeric(NA),
                                 total_score = as.numeric(NA),
                                 recursive_type = 'MS1 + ExpLib Match',
                                 as_seed = NA,
                                 as_seed_round = as.numeric(NA),
                                 source = 'initial_seed',
                                 confidence_level = 'level2') %>%
                   dplyr::rename(ms2_score = temp_colname_ms2) %>%
                   dplyr::select(source, id, id_zhulab, name, formula, confidence_level, smiles, inchikey, inchikey1, adduct, isotope, mz_error,
                                 rt_error_abs, rt_error_rela, ccs_error, ms2_score, ms2_matched_frag, cpd_type, from_peak, from_annotation, recursive_type,
                                 as_seed, as_seed_round, round, total_score)


                 if (test_evaluation == 'No') {
                   # replace ID with KEGG ID
                   #   if the compound not included in KEGG, use ZhuMetLib ID
                   temp_kegg_id <- match(temp_initial_seed$id, cpd_info$id) %>% cpd_info$id_kegg[.]
                   # if (any(is.na(temp_kegg_id))) {
                   #   temp_idx <- which(is.na(temp_kegg_id))
                   #   temp_kegg_id[temp_idx] <- temp_initial_seed$id[temp_idx]
                   # }

                   temp_initial_seed$id <- temp_kegg_id
                 } else {
                   temp_kegg_id <- match(temp_initial_seed$id, cpd_info$id) %>% cpd_info$id[.]
                   temp_initial_seed$id <- temp_kegg_id
                 }



                 if (peak_info$rt <= 500) {
                   temp_rt_tolerance <- tolerance_rt_range
                 } else {
                   temp_rt_tolerance <- peak_info$rt*0.06
                 }

                 # idx_level1 <- which(!(is.na(temp_initial_seed$rt_error_abs)) & temp_initial_seed$rt_error_abs <= 30)
                 idx_level1 <- which(!(is.na(temp_initial_seed$rt_error_abs)) & temp_initial_seed$rt_error_abs <= temp_rt_tolerance)
                 temp_initial_seed$confidence_level[idx_level1] <- 'level1'
                 temp_initial_seed$recursive_type[idx_level1] <- 'MS1 + RT + ExpLib Match'

                 temp_initial_seed <- temp_initial_seed %>%
                   dplyr::arrange(confidence_level, dplyr::desc(ms2_score))
               } else {
                 temp_initial_seed <- tibble::tibble(source = character(),
                                                     id = character(),
                                                     id_zhulab= character(),
                                                     name = character(),
                                                     formula = character(),
                                                     confidence_level = character(),
                                                     smiles = character(),
                                                     inchikey = character(),
                                                     inchikey1 = character(),
                                                     adduct = character(),
                                                     isotope = character(),
                                                     mz_error = numeric(),
                                                     rt_error_abs = numeric(),
                                                     rt_error_rela = numeric(),
                                                     ccs_error = numeric(),
                                                     ms2_score = numeric(),
                                                     ms2_matched_frag = numeric(),
                                                     cpd_type = character(),
                                                     from_peak = character(),
                                                     from_annotation = character(),
                                                     recursive_type = character(),
                                                     as_seed = logical(),
                                                     as_seed_round = numeric(),
                                                     round = numeric(),
                                                     total_score = numeric())
               }


               temp_recursive <- temp_annot_all@recursive_annotation
               if (nrow(temp_recursive) > 0) {
                 temp_recursive <- temp_recursive %>%
                   dplyr::mutate(source = 'recursive_annotation',
                                 confidence_level = 'level3.1') %>%
                   dplyr::rename(recursive_type = annotation_type) %>%
                   dplyr::select(source, id, id_zhulab, name, formula, confidence_level, smiles, inchikey, inchikey1,
                                 adduct, isotope, mz_error, rt_error_abs, rt_error_rela, ccs_error, ms2_score, ms2_matched_frag,
                                 cpd_type, from_peak, from_annotation, recursive_type, as_seed, as_seed_round, round, total_score) %>%
                   dplyr::mutate(from_peak = as.character(from_peak),
                                 from_annotation = as.character(from_annotation),
                                 mz_error = as.numeric(mz_error),
                                 rt_error_abs = as.numeric(rt_error_abs),
                                 rt_error_rela = as.numeric(rt_error_rela),
                                 ccs_error = as.numeric(ccs_error),
                                 as_seed_round = as.numeric(as_seed_round),
                                 ms2_score = as.numeric(ms2_score),
                                 round = as.numeric(round),
                                 total_score = as.numeric(total_score)) %>%
                   dplyr::mutate(adduct = dplyr::case_when(
                     adduct == 'M' ~ '[M]+',
                     adduct == 'M+' ~ '[M]+',
                     adduct == 'M+H' ~ '[M+H]+',
                     adduct == 'M+NH4' ~ '[M+NH4]+',
                     adduct == 'M+Na' ~ '[M+Na]+',
                     adduct == "M-H+2Na" ~ "[M-H+2Na]+",
                     adduct == 'M-2H+3Na' ~ '[M-2H+3Na]+',
                     adduct == "M+K" ~ "[M+K]+",
                     adduct == 'M-H+2K' ~ '[M-H+2K]+',
                     adduct == "M-2H+3K" ~ "[M-2H+3K]+",
                     adduct == 'M+CH3CN+H' ~ '[M+CH3CN+H]+',
                     adduct == "M+CH3CN+Na" ~ "[M+CH3CN+Na]+",
                     adduct == '2M+H' ~ '[2M+H]+',
                     adduct == '2M+NH4' ~ '[2M+NH4]+',
                     adduct == "2M+K" ~ '[2M+K]+',
                     adduct == '2M+Na' ~ '[2M+Na]+',
                     adduct == "M+CH3COO+2H" ~ "[M+CH3COO+2H]+",
                     adduct == 'M+HCOO+2H' ~ '[M+HCOO+2H]+',
                     adduct == 'M+HCOO+H+K' ~ '[M+HCOO+H+K]+',
                     adduct == 'M+HCOO+H+Na' ~ '[M+HCOO+H+Na]+',
                     adduct == 'M-H2O+H' ~ '[M-H2O+H]+',
                     adduct == "M-2H2O+H" ~ "[M-2H2O+H]+",
                     adduct == "M-" ~ "[M]-",
                     adduct == 'M-H' ~ '[M-H]-',
                     adduct == 'M+Na-2H' ~ '[M+Na-2H]-',
                     adduct == 'M+K-2H' ~ '[M+K-2H]-',
                     adduct == 'M+NH4-2H' ~ '[M+NH4-2H]-',
                     adduct == '2M-H' ~ '[2M-H]-',
                     adduct == 'M+CH3COO' ~ '[M+CH3COO]-',
                     adduct == 'M+F' ~ '[M+F]-',
                     adduct == 'M+HCOO' ~ '[M+HCOO]-',
                     adduct == "2M+Na-2H" ~ "[2M+Na-2H]-",
                     adduct == "3M+Na-2H" ~ "[3M+Na-2H]-",
                     adduct == "M+2Na-3H" ~ "[M+2Na-3H]-",
                     adduct == "3M-H" ~ "[3M-H]-",
                     adduct == 'M+Cl' ~ '[M+Cl]-',
                     adduct == 'M+CH3CN-H' ~ '[M+CH3CN-H]-',
                     adduct == 'M+NH3+Cl' ~ '[M+NH3+Cl]-',
                     adduct == "M-2H" ~ "[M-2H]2-",
                     adduct == 'M-H2O-H' ~ '[M-H2O-H]-'))


                 temp_recursive <- temp_recursive %>%
                   dplyr::mutate(confidence_level = dplyr::case_when(cpd_type %in% c('known_known') ~ 'level3.1',
                                                                     cpd_type %in% c('known_unknown', 'unknown_unknown') ~ 'level3.2'))

               } else {
                 temp_recursive <- tibble::tibble(source = character(),
                                                  id = character(),
                                                  id_zhulab= character(),
                                                  name = character(),
                                                  formula = character(),
                                                  confidence_level = character(),
                                                  smiles = character(),
                                                  inchikey = character(),
                                                  inchikey1 = character(),
                                                  adduct = character(),
                                                  isotope = character(),
                                                  mz_error = numeric(),
                                                  rt_error_abs = numeric(),
                                                  rt_error_rela = numeric(),
                                                  ccs_error = numeric(),
                                                  ms2_score = numeric(),
                                                  ms2_matched_frag = numeric(),
                                                  cpd_type = character(),
                                                  from_peak = character(),
                                                  from_annotation = character(),
                                                  recursive_type = character(),
                                                  as_seed = logical(),
                                                  as_seed_round = numeric(),
                                                  round = numeric(),
                                                  total_score = numeric())
               }

               # add seed label for initial annotation
               initial_as_seed <- temp_initial_seed$id %>%
                 match(temp_recursive$id) %>%
                 temp_recursive$as_seed[.]

               initial_as_seed_round <- temp_initial_seed$id %>%
                 match(temp_recursive$id) %>%
                 temp_recursive$as_seed_round[.]

               temp_initial_seed$as_seed <- initial_as_seed
               temp_initial_seed$as_seed_round <- initial_as_seed_round

               # merge initial_seed_annotation and recursive annotation
               #   1. keep ids with same adduct
               #   2. keep the highest confidence level
               temp_merge1 <- temp_initial_seed %>%
                 dplyr::bind_rows(temp_recursive) %>%
                 dplyr::filter(adduct == temp_base_peak_adduct) %>%
                 # dplyr::distinct(inchikey1, .keep_all = TRUE) %>%
                 # dplyr::distinct(id, id_zhulab, .keep_all = TRUE) %>%
                 dplyr::arrange(confidence_level)

               # dereplicates for No_NA and NAs
               temp_merge1_1 <- temp_merge1 %>%
                 dplyr::filter(!is.na(id)) %>%
                 dplyr::distinct(id, .keep_all = TRUE)

               temp_merge1_2 <- temp_merge1 %>%
                 dplyr::filter(is.na(id)) %>%
                 dplyr::distinct(id_zhulab, .keep_all = TRUE)

               temp_merge1 <- temp_merge1_1 %>%
                 dplyr::bind_rows(temp_merge1_2) %>%
                 dplyr::arrange(confidence_level)

               # keep the highest confidence level
               # level 1 > level 2 > level 3.1/level3.2
               # Note: the level of level 3.1 and level 3.2 are considered same for candidate reservation after v0.6.62
               temp_confidence_level <- temp_merge1$confidence_level[1]

               if (temp_confidence_level %in% c('level1', 'level2')) {
                 temp_merge1 <- temp_merge1 %>%
                   dplyr::filter(confidence_level == temp_confidence_level)
               }

               temp_merge1 <- peak_info %>% dplyr::bind_cols(temp_merge1)

               # add credential peak group information
               temp_credential_pg <- temp_annot_all@credential_annotation$peak_group
               temp_idx <- paste(temp_merge1$peak_name, temp_merge1$adduct, sep = '_') %>% match(temp_credential_pg$included_peak_group)
               temp_merge2 <- temp_merge1 %>%
                 dplyr::mutate(peak_group_id = temp_credential_pg$peak_group_id[temp_idx],
                               base_peak = temp_credential_pg$included_peak_group[temp_idx],
                               rela_peaks = temp_credential_pg$annotated_peaks[temp_idx]) %>%
                 dplyr::filter(!is.na(peak_group_id))

               # add credential peak formula information
               temp_credential_formula <- temp_annot_all@credential_annotation$formula_prediction
               temp_merge3 <- temp_merge2 %>%
                 dplyr::mutate(cons_formula_pred = formula %in% temp_credential_formula$formula)

               return(temp_merge3)
             })

             names(annot_base_peaks) <- peak_group_id_table$peak_group_id

             return(annot_base_peaks)
           })


#   replaceRtError -------------------------------------------------------------

setGeneric(name = 'replaceRtError',
           function(
             table_identification,
             # item = c('level3.1', 'level3.2'),
             path = '.'
           ){
             load(file.path(path, '02_result_MRN_annotation/00_intermediate_data/rt_result'))

             temp_pred_rt_lib <- rt_result$KEGG.rt

             idx1 <- which(table_identification$confidence_level %in% c('level3.1', 'level3.2'))
             pred_rt <- match(table_identification$id_kegg[idx1], rownames(temp_pred_rt_lib)) %>% temp_pred_rt_lib$RT[.]
             exp_rt <- table_identification$rt[idx1]

             rt_error <- mapply(function(x, y){
               abs(x-y)/y*100
             },
             x = exp_rt,
             y = pred_rt)

             table_identification$rt_error_rela[idx1] <- rt_error
             table_identification$rt_error_abs[idx1] <- NA
             return(table_identification)
})



#   calculateExportTableScore --------------------------------------------------
#' @title calculateExportTableScore
#' @author Zhiwei Zhou
#' @param table_identification
#' @param mz_tol
#' @param rt_tol
#' @param ccs_tol
#' @param instrument

setGeneric(name = 'calculateExportTableScore',
           def = function(table_identification,
                          mz_tol = 25,
                          rt_tol = 30,
                          ccs_tol = 3,
                          instrument = c("SciexTripleTOF", "AgilentQTOF", "BrukerQTOF", "ThermoOrbitrap",
                                         'ThermoExploris', "AgilentDTIMMS", "BrukerTIMS", 'WatersQTOF','WatersTWIMMS'),
                          ...){


             if (instrument %in% c("AgilentDTIMMS", "BrukerTIMS", "WatersTWIMMS")) {
               weight_mz <- 0.1
               weight_rt <- weight_ccs <- 0.2
               weight_ms2 <- 0.5
             } else {
               weight_mz <- weight_rt <- 0.25
               weight_ccs <- 0
               weight_ms2 <- 0.5
             }

             table_identification <- table_identification %>%
               tidyr::replace_na(list('confidence_level' = 'level3.1'))

             score_iden <- mapply(function(x_mz, x_rt1, x_rt2, x_ccs, x_ms2, y){
               if (y %in% c('level1', 'level2')) {
                 result <- getIdenScore(mz_error = x_mz,
                                        rt_error = x_rt1,
                                        ccs_error = x_ccs,
                                        ms2_score = x_ms2,
                                        mz_tol = mz_tol,
                                        rt_tol = rt_tol,
                                        ccs_tol = ccs_tol,
                                        weight_mz = weight_mz,
                                        weight_rt = weight_rt,
                                        weight_ccs = weight_ccs,
                                        weight_ms2 = weight_ms2)
               } else {
                 result <- getIdenScore(mz_error = x_mz,
                                        rt_error = x_rt2,
                                        ccs_error = x_ccs,
                                        ms2_score = x_ms2,
                                        mz_tol = mz_tol,
                                        rt_tol = rt_tol,
                                        ccs_tol = ccs_tol,
                                        weight_mz = weight_mz,
                                        weight_rt = weight_rt,
                                        weight_ccs = weight_ccs,
                                        weight_ms2 = weight_ms2)
               }

               return(result)
             },
             x_mz = table_identification$mz_error,
             x_rt1 = table_identification$rt_error_abs,
             x_rt2 = table_identification$rt_error_rela,
             x_ccs = table_identification$ccs_error,
             x_ms2 = table_identification$ms2_score,
             y = table_identification$confidence_level)

             score_total <- mapply(function(x, y){
               switch(y,
                      'level1' = {
                        3 + x
                      },
                      'level2' = {
                        2 + x
                      },
                      'level3.1' = {
                        1 + x
                      },
                      'level3.2' = {
                        1 + x
                      },
                      'level4' = {
                        0 + x
                      })
             },
             x = score_iden,
             y = table_identification$confidence_level)

             table_identification <- table_identification %>%
               dplyr::mutate(iden_score = score_iden,
                             total_score = score_total) %>%
               dplyr::select(peak_name:adduct, total_score, mz_error:ms2_score, iden_score, dplyr::everything()) %>%
               dplyr::mutate(mz = round(mz, 4),
                             rt = round(rt, 1),
                             ccs = round(ccs, 1),
                             total_score = round(total_score, 2),
                             mz_error = round(mz_error, 1),
                             rt_error_abs = round(rt_error_abs, 1),
                             rt_error_rela = round(rt_error_rela, 1),
                             ccs_error = round(ccs_error, 1),
                             ms2_score = round(ms2_score, 4),
                             iden_score = round(iden_score, 2))

             return(table_identification)
           })

#     getIdenScore ---------------------------------------------------------------
setGeneric(name = 'getIdenScore',
           def = function(mz_error,
                          rt_error,
                          ccs_error = NULL,
                          ms2_score,
                          mz_tol = 25,
                          rt_tol = 30,
                          ccs_tol = 3,
                          weight_mz = 0.25,
                          weight_rt = 0.25,
                          weight_ccs = 0,
                          weight_ms2 = 0.5){

             if (is.na(mz_error)) {
               mz_error <- mz_tol
             }

             if (is.na(rt_error)) {
               rt_error <- rt_tol
             }

             if (is.na(ccs_error) | ccs_error < 0) {
               ccs_error <- ccs_tol
             }

             if (is.na(ms2_score)) {
               ms2_score <- 0
             }

             mz_score <- getLinerScore(delta = mz_error, tolerance = mz_tol)
             if (mz_score < 0) {
               mz_score <- 0
             }

            rt_score <- getLinerScore(delta = rt_error, tolerance = rt_tol)
            if (rt_score < 0) {
              rt_score <- 0
            }

            ccs_score <- getLinerScore(delta = ccs_error, tolerance = ccs_tol)
            if (ccs_score < 0) {
              ccs_score <- 0
            }

            score_iden <- weight_mz*mz_score + weight_rt*rt_score + weight_ccs*ccs_score + weight_ms2*ms2_score

            return(score_iden)
           })

#   rename4UniqueIndex ---------------------------------------------------------

setGeneric(name = 'rename4UniqueIndex',
           def = function(
             obj
           ){
             # uni_obj <- unique(obj)
             obj_rep <- table(obj) %>% .[.>=2] %>% names()
             for (x in obj_rep) {
               idx <- which(obj == x)
               obj[idx] <- paste(x, letters[1:length(idx)], sep = '_')
             }

             return(obj)
           })


# removeReplicateSeedAnnotationFromFinalTable ----------------------------------

#' @title removeReplicateSeedAnnotationFromFinalTable
#' @description remove replicated metabolite annotations between seed annotation and recursive annotation
#' @author Zhiwei Zhou
#' @param dir_path project path
#' @param instrument
#' @export

# removeReplicateSeedAnnotationFromFinalTable(dir_path = '/home/zhouzw/Data_processing/20220118_test_pred_rt_caret_package/test_remove',
#                                             instrument = "SciexTripleTOF")

removeReplicateSeedAnnotationFromFinalTable <- function(
  dir_path,
  instrument = c("SciexTripleTOF", "AgilentQTOF", "BrukerQTOF", "ThermoOrbitrap",
                 'ThermoExploris', "AgilentDTIMMS", "BrukerTIMS", 'WatersQTOF','WatersTWIMMS')
){
  # remove replicate IDs in initial seed
  load(file.path(dir_path, '00_annotation_table/00_intermediate_data/peak_group_id_table'))
  load(file.path(dir_path, '01_result_initial_seed_annotation/00_intermediate_data/ms1_data'))
  load(file.path(dir_path, '00_annotation_table/00_intermediate_data/list_peak_group'))
  load(file.path(dir_path, '00_annotation_table/00_intermediate_data/list_identification'))
  load(file.path(dir_path, '00_annotation_table/00_intermediate_data/list_peak_group'))
  load(file.path(dir_path, '00_annotation_table/00_intermediate_data/table_identification'))

  id_kegg_level12 <- table_identification %>%
    dplyr::filter(confidence_level %in% c('level1', 'level2')) %>%
    dplyr::filter(!is.na(id_kegg)) %>%
    dplyr::pull(id_kegg) %>%
    unique()

  # remove annotations replicated with level 1 and level 2 annotation
  # keep other annotations
  idx_remove <- which((table_identification$id_kegg %in% id_kegg_level12) &
                        (!(table_identification$confidence_level %in% c('level1', 'level2'))))

  if (length(idx_remove) > 0) {
    table_identification <- table_identification[-idx_remove, ]

    save(table_identification,
         file = file.path(dir_path, '00_annotation_table/00_intermediate_data/table_identification'),
         version = 2)
  }

  # modify table2
  id_peak_group_remove <- setdiff(names(list_peak_group), unique(table_identification$peak_group_id))
  idx_remove2 <- match(id_peak_group_remove, names(list_peak_group))

  if (length(idx_remove2) > 0) {
    list_peak_group <- list_peak_group[-idx_remove2]
    save(list_peak_group,
         file = file.path(dir_path, '00_annotation_table/00_intermediate_data/list_peak_group'),
         version = 2)
  }

  table_peak_group <- list_peak_group %>%
    dplyr::bind_rows() %>%
    dplyr::select(peak_name:ccs, isotope, adduct, peak_group_id, base_peak, rela_peaks) %>%
    dplyr::mutate(pg_num = as.numeric(stringr::str_extract(peak_group_id, pattern = '\\d+'))) %>%
    dplyr::arrange(pg_num) %>%
    dplyr::distinct(peak_name, peak_group_id, .keep_all = TRUE) %>%
    dplyr::select(-pg_num) %>%
    dplyr::select(peak_group_id, dplyr::everything()) %>%
    dplyr::mutate(mz = round(mz, 4),
                  rt = round(rt, 1)) %>%
    tidyr::replace_na(list('adduct' = ''))

  # replace base_peak and neutral mass
  table_peak_group$base_peak <- match(table_peak_group$peak_group_id, peak_group_id_table$peak_group_id) %>%
    peak_group_id_table$base_peak[.]
  table_peak_group$neutral_mass <- match(table_peak_group$peak_group_id, peak_group_id_table$peak_group_id) %>%
    peak_group_id_table$base_neutral_mass[.]
  table_peak_group <- table_peak_group %>%
    dplyr::select(peak_group_id:base_peak, neutral_mass, dplyr::everything()) %>%
    dplyr::mutate(neutral_mass = round(neutral_mass, 4)) %>%
    dplyr::rename(num_peaks = rela_peaks)

  if (!(instrument %in% c("AgilentDTIMMS", "BrukerTIMS", "WatersTWIMMS"))) {
    table_peak_group <- table_peak_group %>%
      dplyr::select(-ccs)
  }

  readr::write_csv(table_peak_group,
                   file = file.path(dir_path, '00_annotation_table', 'table2_peak_group.csv'))


  # modify table 1
  idx_remove3 <- match(id_peak_group_remove, names(list_identification))

  if (length(idx_remove3) > 0) {
    list_identification <- list_identification[-idx_remove3]
    save(list_identification,
         file = file.path(dir_path, '00_annotation_table/00_intermediate_data/list_identification'),
         version = 2)
  }


  table_identification_wide <- table_identification %>%
    dplyr::group_by(peak_name) %>%
    dplyr::summarise(peak_name = peak_name[1],
                     mz = mz[1],
                     rt = rt[1],
                     id_kegg = paste(id_kegg, collapse = ';'),
                     id_zhulab = paste(id_zhulab, collapse = ';'),
                     name = paste(name, collapse = ';'),
                     formula = paste(formula, collapse = ';'),
                     confidence_level = paste(confidence_level, collapse = ';'),
                     smiles = paste(smiles, collapse = ';'),
                     inchikey = paste(inchikey, collapse = ';'),
                     isotope = paste(isotope, collapse = ';'),
                     adduct = paste(adduct, collapse = ';'),
                     total_score = paste(total_score, collapse = ';'),
                     mz_error = paste(mz_error, collapse = ';'),
                     rt_error_abs = paste(rt_error_abs, collapse = ';'),
                     rt_error_rela = paste(rt_error_rela, collapse = ';'),
                     ms2_score = paste(ms2_score, collapse = ';'),
                     iden_score = paste(iden_score, collapse = ';'),
                     iden_type = paste(iden_type, collapse = ';'),
                     peak_group_id = paste(peak_group_id, collapse = ';'),
                     base_peak = paste(base_peak, collapse = ';'),
                     num_peaks = paste(num_peaks, collapse = ';'),
                     cons_formula_pred = paste(cons_formula_pred, collapse = ';'),
                     stereo_isomer_id = paste(unique(stereo_isomer_id), collapse = ';'),
                     stereo_isomer_name = paste(unique(stereo_isomer_name), collapse = ';')) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(mz)

  # join the quantification table
  quanti_table <- ms1_data$info %>%
    dplyr::bind_cols(ms1_data$subject) %>%
    dplyr::mutate(mz = round(mz, 4),
                  rt = round(rt, 1),
                  ccs = round(ccs, 1))

  if (!(instrument %in% c("AgilentDTIMMS", "BrukerTIMS", "WatersTWIMMS"))) {
    quanti_table <- quanti_table %>%
      dplyr::select(-c('ccs'))
  }

  table_identification_wide_quanti <- table_identification_wide %>%
    dplyr::full_join(quanti_table, by = c('peak_name' = 'name')) %>%
    dplyr::rename(mz = mz.y,
                  rt = rt.y) %>%
    dplyr::select(-c('mz.x', 'rt.x')) %>%
    dplyr::select(peak_name, mz, rt, dplyr::everything())

  readr::write_csv(table_identification_wide_quanti,
                   file = file.path(dir_path, '00_annotation_table', 'table1_identification.csv'))


  # modify table 3
  table_identification_unique <- table_identification
  temp_peak_name <- table_identification_unique$peak_name %>%
    rename4UniqueIndex()

  # add quanti table for table_identification_unique
  table_identification_unique <- table_identification %>%
    dplyr::left_join(quanti_table, by = c('peak_name' = 'name')) %>%
    dplyr::rename(mz = mz.y,
                  rt = rt.y) %>%
    dplyr::select(-c('mz.x', 'rt.x')) %>%
    dplyr::select(peak_name, mz, rt, dplyr::everything())

  table_identification_unique$peak_name <- temp_peak_name

  readr::write_csv(table_identification_unique,
                   file = file.path(dir_path, '00_annotation_table', 'table3_identification_pair.csv'))

}
