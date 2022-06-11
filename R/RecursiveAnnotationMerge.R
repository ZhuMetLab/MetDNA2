################################################################################
# mergeRecursiveAnnotation -----------------------------------------------------
#' @title mergeRecursiveAnnotation
#' @author Zhiwei Zhou
#' @param path '.'

# path <- '/home/zhouzw/Data_processing/20201214_metdna2_development'
# load(file.path(path, 'MS2_match_result/intermediate_data', 'result_annotation'))
# load(file.path(path, 'MRN_annotation_result/intermediate_data', 'tags2_after_redundancy_remove'))
#
# result_initial_annotation <- result_annotation
# result_recursive_annotation <- tags2_after_redundancy_remove
# names(result_recursive_annotation) <- names(result_initial_annotation)

# test <- mergeRecursiveAnnotation(path = path,
#                                  thread = 3,
#                                  direction = 'reverse',
#                                  tolerance_rt_range = 30)
#
# path <- '/home/zhouzw/Data_processing/20210720_modify_MetDNA2_for_200STD_validation/POS_200std_mouse_liver/'
# load('H:/00_projects/03_MetDNA2/06_files_for_package/library/20210223_update_emrn/cpd_emrn_210223.RData')
# load('H:/00_projects/03_MetDNA2/06_files_for_package/library/20210403_update_cpd_info/zhuMetlib_20210403.RData')

# temp <- mergeRecursiveAnnotation(path = path,
#                                  thread = 3,
#                                  lib = 'zhuMetLib',
#                                  direction = 'reverse',
#                                  tolerance_rt_range = 30,
#                                  use_redun_rm_result = TRUE)

setGeneric(name = 'mergeRecursiveAnnotation',
           def = function(
             path = '.',
             use_redun_rm_result = FALSE,
             # for skip functions
             is_anno_mrn = TRUE,
             is_credential = TRUE,
             ...
           ){
             if ('id_merge_before_credential.RData' %in%
                 list.files(file.path(path, '03_annotation_credential', "00_intermediate_data"))) {
               cat("Note: id_merge_before_credential existed. Skip the step of merging initial_seed annotation and recursive annotation results.\n")
               return(NULL)
             }

             load(file.path(path, '01_result_initial_seed_annotation/00_intermediate_data', 'result_annotation'))
             result_initial_annotation <- result_annotation

             if (use_redun_rm_result) {
               load(file.path(path, '02_result_MRN_annotation/00_intermediate_data',
                              'tags2_after_redundancy_remove'))
               result_recursive_annotation <- tags2_after_redundancy_remove
             } else {
               load(file.path(path, '02_result_MRN_annotation/00_intermediate_data',
                              'tags2_after_annotation'))
               result_recursive_annotation <- tags2_after_annotation
             }

             names(result_recursive_annotation) <- names(result_initial_annotation)
             rm(list = c('result_annotation', 'tags2_after_redundancy_remove', 'tags_after_annotation'));gc()

             # browser()
             cat('Merge initial_seed annotation and recursive annotation results\n')
             annot_all <- generateMetDNA2AnnoClass(result_initial_annotation = result_initial_annotation)

             annot_all <- addInitialSeedResult2MetDNA2AnnoClass(annot_all = annot_all,
                                                                result_initial_annotation = result_initial_annotation)
             
             if (is_anno_mrn) {
               annot_all <- addRecursiveResult2MetDNA2AnnoClass(annot_all = annot_all,
                                                                result_recursive_annotation = result_recursive_annotation,
                                                                ...)
             } else {
               cat('Skip recursive annotation results addition.\n')
             }

             dir.create(file.path(path, '03_annotation_credential', "00_intermediate_data"),
                        showWarnings = FALSE,
                        recursive = TRUE)

             save(annot_all,
                  file = file.path(path, '03_annotation_credential', "00_intermediate_data", 'annot_all_before_credential.RData'),
                  version = 2)

             # id_merge <- generateResultTable4AllAnnotation(annot_all = annot_all,
             #                                               is_cred_pg_filter = FALSE,
             #                                               is_cred_formula_filter = FALSE,
             #                                               ...)
             # 
             # # id_merge <- generateResultTable4AllAnnotation(annot_all = annot_all,
             # #                                               is_cred_pg_filter = FALSE,
             # #                                               is_cred_formula_filter = FALSE)
             # 
             # save(id_merge,
             #      file = file.path(path, '03_annotation_credential', "00_intermediate_data", 'id_merge_before_credential.RData'),
             #      version = 2)
             # 
             # return(id_merge)

             # for skip functions
             if (is_credential) {
               id_merge <- generateResultTable4AllAnnotation(annot_all = annot_all,
                                                             is_cred_pg_filter = FALSE,
                                                             is_cred_formula_filter = FALSE,
                                                             ...)
               
               # id_merge <- generateResultTable4AllAnnotation(annot_all = annot_all,
               #                                               is_cred_pg_filter = FALSE,
               #                                               is_cred_formula_filter = FALSE)
               
               save(id_merge,
                    file = file.path(path, '03_annotation_credential', "00_intermediate_data", 'id_merge_before_credential.RData'),
                    version = 2)
               
               return(id_merge)
             } else {
               cat('Skip id merge before credential.\n');return(NULL)
             }
             
           })


################################################################################
#   generateMetDNA2AnnoClass -----------------------------------------------------

#' @title generateMetDNA2AnnoClass
#' @description Convert Ms1DataInfor to AnnotationResult class
#' @author Zhiwei Zhou
#' @param ms1_info
#' @export

# path <- '/home/zhouzw/Data_processing/20201214_metdna2_development'
# load(file.path(path, '01_result_initial_seed_annotation/00_intermediate_data', 'result_annotation'))
# load(file.path(path, '02_result_MRN_annotation/intermediate_data', 'tags2_after_redundancy_remove'))
# load(file.path(path, '03_annotation_credential/00_intermediate_data', 'list_peak_group_annotation_concised.RData'))
# load(file.path(path, '03_annotation_credential/00_intermediate_data', 'list_peak_group_formula.RData'))
#
# result_initial_annotation <- result_annotation
# result_recursive_annotation <- tags2_after_redundancy_remove
# names(result_recursive_annotation) <- names(result_initial_annotation)
# result_credential_pg <- list_peak_group_annotation_concised
# result_credential_formula <- list_peak_group_formula
# annot_all <- generateMetDNA2AnnoClass(result_initial_annotation = result_initial_annotation)

setGeneric(name = 'generateMetDNA2AnnoClass',
           def = function(
             result_initial_annotation
           ){
             cat('Generate MetDNA2Annotation...\n')

             temp_initial_annotation <- tibble::tibble(id = character(),
                                                       id_zhulab = character(),
                                                       name = character(),
                                                       formula = character(),
                                                       smiles = character(),
                                                       inchikey = character(),
                                                       inchikey1 = character(),
                                                       adduct = character(),
                                                       isotope = character(),
                                                       mz_lib = numeric(),
                                                       rt_lib = numeric(),
                                                       ccs_lib = numeric(),
                                                       mz_error = numeric(),
                                                       # mz_score=mz_score,
                                                       rt_error_abs = numeric(),
                                                       rt_error_rela = numeric(),
                                                       # rt_score = numeric(),
                                                       ccs_error = numeric(),
                                                       # ccs_score = numeric(),
                                                       ms2_score_forward = numeric(),
                                                       ms2_score_reverse = numeric(),
                                                       ms2_matched_frag = numeric(),
                                                       cpd_type = character())

             temp_recursive_annotation <- tibble::tibble(id = character(),
                                                         id_zhulab = character(),
                                                         annotation_type = character(),
                                                         name = character(),
                                                         formula = character(),
                                                         smiles = character(),
                                                         inchikey = character(),
                                                         inchikey1 = character(),
                                                         adduct = character(),
                                                         isotope = character(),
                                                         from_peak = character(),
                                                         from_annotation = character(),
                                                         round = numeric(),
                                                         as_seed = logical(),
                                                         as_seed_round = numeric(),
                                                         mz_error = numeric(),
                                                         rt_error_abs = numeric(),
                                                         rt_error_rela = numeric(),
                                                         ccs_error = numeric(),
                                                         ms2_score = numeric(),
                                                         ms2_matched_frag = numeric(),
                                                         total_score = numeric(),
                                                         cpd_type = character())

             temp_credential_pg <- tibble::tibble(peak_name = character(),
                                                  mz = numeric(),
                                                  rt = numeric(),
                                                  ssc = numeric(),
                                                  int = numeric(),
                                                  mz_error = numeric(),
                                                  rt_error = numeric(),
                                                  ccs_error = numeric(),
                                                  int_ratio = numeric(),
                                                  ms2_score = numeric(),
                                                  type_annotation = character(),
                                                  isotope = character(),
                                                  annotation = character(),
                                                  included_peak_group = character(),
                                                  peak_group_id = character())

             # type_annotation, isotope, annotation,

             temp_credential_formula <- tibble::tibble(peak_name = character(),
                                                       adduct = character(),
                                                       peak_group = character(),
                                                       formula = character(),
                                                       score = numeric())

             temp_credential_annotation <- list('peak_group' = temp_credential_pg,
                                                'formula_prediction' = temp_credential_formula)

             temp_result <- vector(mode = 'list', length = length(result_initial_annotation))

             result <- pbapply::pblapply(seq_along(temp_result), function(i){
               peak_info <- result_initial_annotation[[i]]@peak_info


               new(Class = "MetDNA2Annotation",
                   peak_info = peak_info,
                   initial_seed_annotation = temp_initial_annotation,
                   recursive_annotation = temp_recursive_annotation,
                   credential_annotation = temp_credential_annotation)
             })

             names(result) <- names(result_initial_annotation)

             return(result)
           }
)





################################################################################
#   addInitialSeedResult2MetDNA2AnnoClass --------------------------------------
# load(file.path(path, 'MS2_match_result/intermediate_data', 'result_annotation'))
# load(file.path(path, 'MRN_annotation_result/intermediate_data', 'tags2_after_redundancy_remove'))
# load(file.path(path, 'annot_credential/00_intermediate_data', 'list_peak_group_annotation_concised.RData'))
# load(file.path(path, 'annot_credential/00_intermediate_data', 'list_peak_group_formula.RData'))
#
# result_initial_annotation <- result_annotation
# result_recursive_annotation <- tags2_after_redundancy_remove
# names(result_recursive_annotation) <- names(result_initial_annotation)
# result_credential_pg <- list_peak_group_annotation_concised
# result_credential_formula <- list_peak_group_formula
#
# annot_all <- addInitialSeedResult2MetDNA2AnnoClass(annot_all = annot_all,
#                                                    result_initial_annotation = result_initial_annotation)

setGeneric(name = 'addInitialSeedResult2MetDNA2AnnoClass',
           def = function(
             annot_all,
             result_initial_annotation
           ){
             cat('Add initial annotation result to MetDNA2Annotation\n')
             idx <- sapply(result_initial_annotation, function(x){nrow(x@annotation_result) > 0}) %>% which()

             temp_col <- annot_all[[1]]@initial_seed_annotation %>% colnames()
             pb <- utils::txtProgressBar(min = 0, max = length(idx), style = 3)
             for (i in seq_along(idx)) {
               utils::setTxtProgressBar(pb, i)
               temp_idx <- idx[[i]]
               temp_anno_initial <- result_initial_annotation[[temp_idx]]@annotation_result %>%
                 dplyr::select(id, name, formula, smiles, inchikey, inchikey1, adduct, mz_lib, rt_lib, ccs_lib,
                               mz_error, rt_error, ccs_error, msms_score_forward, msms_score_reverse, msms_matched_frag) %>%
                 dplyr::rename(rt_error_abs = rt_error,
                               ms2_score_forward = msms_score_forward,
                               ms2_score_reverse = msms_score_reverse,
                               ms2_matched_frag = msms_matched_frag) %>%
                 dplyr::mutate(id_zhulab = id,
                               rt_error_rela = NA,
                               cpd_type = 'known_known',
                               isotope = '[M]') %>%
                 dplyr::select(temp_col)

               annot_all[[temp_idx]]@initial_seed_annotation <- temp_anno_initial

             }

             return(annot_all)
           })


#   addRecursiveResult2MetDNA2AnnoClass ----------------------------------------------
#
# test <- addRecursiveResult2MetDNA2AnnoClass(annot_all = annot_all,
#                                             result_recursive_annotation = result_recursive_annotation, test_evaluation = '200STD')

setGeneric(name = 'addRecursiveResult2MetDNA2AnnoClass',
           def = function(
             annot_all,
             result_recursive_annotation,
             test_evaluation = c('No', '200STD', '46STD'),
             ...
           ){
             cat('\n');cat('Add recursive annotation result to MetDNA2Annotation\n')
             idx <- sapply(result_recursive_annotation, function(x){length(x@annotation) > 0}) %>% which()
             temp_col <- annot_all[[1]]@recursive_annotation %>% colnames()

             data("cpd_emrn", envir = environment())
             if (test_evaluation == '200STD') {
               rm(list = c('cpd_emrn'));gc()
               data("cpd_200stdExd", envir = environment())
               cpd_emrn <- cpd_200stdExd
               rm(list = c('cpd_200stdExd'));gc()
             }
             if (test_evaluation == '46STD') {
               rm(list = c('cpd_emrn'));gc()
               data("cpd_46stdExd", envir = environment())
               cpd_emrn <- cpd_46stdExd
               rm(list = c('cpd_46stdExd'));gc()
             }

             pb <- utils::txtProgressBar(min = 0, max = length(idx), style = 3)
             for (i in seq_along(idx)) {
               utils::setTxtProgressBar(pb, i)
               temp_idx <- idx[[i]]
               options(readr.num_columns = 0)
               temp_anno_recursive <- result_recursive_annotation[[temp_idx]]@annotation %>%
                 unlist() %>%
                 matrix(nrow = 19) %>%
                 t() %>%
                 tibble::as_tibble() %>%
                 readr::type_convert()
               colnames(temp_anno_recursive) <- c("type", 'From', 'From.peak', 'to', 'step', 'level', 'as.seed', 'as.seed.round',
                                                  'isotope', 'adduct', 'charge', 'formula', 'mz.error', 'rt.error', 'ccs.error',
                                                  'int.error', 'ms2.sim', 'nfrag', 'score')

               temp_anno_recursive <- temp_anno_recursive %>%
                 dplyr::select(type, to, adduct, isotope, From.peak, From, level, as.seed, as.seed.round,
                               mz.error, rt.error, ccs.error, ms2.sim, nfrag, score) %>%
                 dplyr::rename(id = to,
                               annotation_type = type,
                               from_peak = From.peak,
                               from_annotation = From,
                               round = level,
                               as_seed = as.seed,
                               as_seed_round = as.seed.round,
                               mz_error = mz.error,
                               rt_error = rt.error,
                               ccs_error = ccs.error,
                               ms2_score = ms2.sim,
                               ms2_matched_frag = nfrag,
                               total_score = score) %>%
                 dplyr::mutate(mz_error = as.numeric(mz_error),
                               rt_error = as.numeric(rt_error),
                               ccs_error = as.numeric(ccs_error),
                               ms2_score = as.numeric(ms2_score),
                               ms2_matched_frag = as.numeric(ms2_matched_frag),
                               total_score = as.numeric(total_score)) %>%
                 dplyr::mutate(id_zhulab = NA,
                               name = NA,
                               formula = NA,
                               smiles = NA,
                               inchikey = NA,
                               inchikey1 = NA,
                               cpd_type = NA,
                               rt_error_abs = NA,
                               rt_error_rela = NA) %>%
                 dplyr::select(temp_col, rt_error)

               idx_metannot <- which(temp_anno_recursive$annotation_type == 'metAnnotation')
               idx_adduct_iso <- which(temp_anno_recursive$annotation_type != 'metAnnotation')

               if (length(idx_metannot) > 0) {
                 temp_anno_recursive$rt_error_rela[idx_metannot] <- temp_anno_recursive$rt_error[idx_metannot]
               }

               if (length(idx_adduct_iso) > 0) {
                 temp_anno_recursive$rt_error_abs[idx_adduct_iso] <- temp_anno_recursive$rt_error[idx_adduct_iso]
               }

               temp_anno_recursive <- temp_anno_recursive %>% dplyr::select(-rt_error)

               # add structure information into recursive annotation
               temp_idx2 <- match(temp_anno_recursive$id, cpd_emrn$id)
               
               # for (i_na_check in seq_along(temp_idx2)) {
               #   if (is.na(temp_idx2[i_na_check])) {
               #     temp_idx2[i_na_check] = grep(pattern = temp_anno_recursive$id[i_na_check], x = cpd_emrn$id_kegg_synonyms)[1]
               #   }
               #   if (is.na(temp_idx2[i_na_check])) {
               #     temp_idx2[i_na_check] = grep(pattern = temp_anno_recursive$id[i_na_check], x = cpd_emrn$id_kegg)[1]
               #   }
               # } # 20220423 id_kegg redirection by Haosong Zhang
               
               temp_anno_recursive$name <- cpd_emrn$name[temp_idx2]
               temp_anno_recursive$formula <- cpd_emrn$formula[temp_idx2]
               temp_anno_recursive$smiles <- cpd_emrn$smiles[temp_idx2]
               temp_anno_recursive$inchikey <- cpd_emrn$inchikey[temp_idx2]
               temp_anno_recursive$inchikey1 <- cpd_emrn$inchikey1[temp_idx2]
               temp_anno_recursive$cpd_type <- cpd_emrn$type[temp_idx2]
               
               for (i_na_check in seq_along(temp_idx2)) {
                 if (is.na(temp_idx2[i_na_check])) {
                   data("lib_kegg", envir = environment())
                   temp_idx_kegg <- match(temp_anno_recursive$id[i_na_check], lib_kegg$id)
                   temp_anno_recursive$name[i_na_check] <- lib_kegg$name[temp_idx_kegg]
                   temp_anno_recursive$formula[i_na_check] <- lib_kegg$formula[temp_idx_kegg]
                   temp_anno_recursive$smiles[i_na_check] <- lib_kegg$smiles[temp_idx_kegg]
                   temp_anno_recursive$inchikey[i_na_check] <- lib_kegg$inchikey[temp_idx_kegg]
                   temp_anno_recursive$inchikey1[i_na_check] <- lib_kegg$inchikey1[temp_idx_kegg]
                   temp_anno_recursive$cpd_type[i_na_check] <- lib_kegg$type[temp_idx_kegg]
                   if (is.na(temp_idx_kegg)) {
                     temp_idx_emrn = grep(pattern = temp_anno_recursive$id[i_na_check], x = cpd_emrn$id_kegg_synonyms)[1]
                     if (is.na(temp_idx_emrn)) {
                       temp_idx_emrn = grep(pattern = temp_anno_recursive$id[i_na_check], x = cpd_emrn$id_kegg)[1]
                     }
                     temp_anno_recursive$name[i_na_check] <- cpd_emrn$name[temp_idx_emrn]
                     temp_anno_recursive$formula[i_na_check] <- cpd_emrn$formula[temp_idx_emrn]
                     temp_anno_recursive$smiles[i_na_check] <- cpd_emrn$smiles[temp_idx_emrn]
                     temp_anno_recursive$inchikey[i_na_check] <- cpd_emrn$inchikey[temp_idx_emrn]
                     temp_anno_recursive$inchikey1[i_na_check] <- cpd_emrn$inchikey1[temp_idx_emrn]
                     temp_anno_recursive$cpd_type[i_na_check] <- cpd_emrn$type[temp_idx_emrn]
                   }
                 }
               } # 20220427 id_kegg redirection from lib_kegg by Haosong Zhang

               annot_all[[temp_idx]]@recursive_annotation <- temp_anno_recursive

             }

             return(annot_all)
           })







#   addCredentialResult2MetDNA2AnnoClass ---------------------------------------


setGeneric(name = 'addCredentialResult2MetDNA2AnnoClass',
           def = function(
             annot_all,
             result_credential_pg,
             result_credential_formula,
             is_pred_formula_all = FALSE
           ){
             cat('\n');cat('Add credential result to MetDNA2Annotation\n')
             # result_credential_pg <- concisePeakGroup2Annotation(list_peak_group_annotation = list_peak_group_annotation_concised)
             list_credential_pg <- concisePeakGroup2Annotation(list_peak_group_annotation = result_credential_pg)

             result_credential_pg_ids <- list_credential_pg$peak_group_id
             result_credential_pg <- list_credential_pg$result_annotation_table

             result_credential_pg1 <- result_credential_pg %>%
               # dplyr::filter(stringr::str_detect(type_annotation, 'isotopeAnnotation')) %>%
               tidyr::separate_rows(ssc, int, mz_error, rt_error, ccs_error, int_ratio, ms2_score,
                                    isotope, type_annotation, annotation, peak_group, peak_group_id, peak_group_scale,
                                    sep = ';', convert = TRUE) %>%
               dplyr::filter(type_annotation == 'isotopeAnnotation') %>%
               dplyr::select(peak_name:ms2_score, type_annotation, isotope, annotation, peak_group, peak_group_id, peak_group_scale) %>%
               dplyr::rename(included_peak_group = peak_group,
                             annotated_peaks = peak_group_scale)

             result_credential_pg2 <- result_credential_pg %>%
               tidyr::separate_rows(ssc, int, mz_error, rt_error, ccs_error, int_ratio, ms2_score,
                                    isotope, type_annotation, annotation, peak_group, peak_group_id, peak_group_scale,
                                    sep = ';', convert = TRUE) %>%
               dplyr::filter(type_annotation != 'isotopeAnnotation') %>%
               dplyr::select(peak_name:ms2_score, type_annotation, isotope, annotation, peak_group, peak_group_id, peak_group_scale) %>%
               dplyr::rename(included_peak_group = peak_group,
                             annotated_peaks = peak_group_scale)

             result_credential_pg <- result_credential_pg2 %>%
               dplyr::bind_rows(result_credential_pg1) %>%
               dplyr::arrange(mz, peak_name)

             idx <- match(unique(result_credential_pg$peak_name), names(annot_all))
             pb <- utils::txtProgressBar(min = 0, max = length(idx), style = 3)
             for (i in seq_along(idx)) {
               utils::setTxtProgressBar(pb, i)
               temp_idx <- idx[i]
               temp_peak <- annot_all[[temp_idx]]@peak_info$name
               annot_all[[temp_idx]]@credential_annotation$peak_group <- result_credential_pg %>% dplyr::filter(peak_name == temp_peak)
             }

             # result_credential_formula <- generateFormulaResult(list_peak_group_formula = list_peak_group_formula, num_formula_candidate = 3)
             result_credential_formula <- generateFormulaResult(list_peak_group_formula = result_credential_formula, num_formula_candidate = 3)
             result_credential_formula <- result_credential_formula %>%
               tidyr::separate_rows(formula, score, sep = ';') %>%
               dplyr::mutate(peak_group = paste(name, adduct, sep = '_')) %>%
               dplyr::select(name, adduct, peak_group, formula, score) %>%
               dplyr::rename(peak_name = name)


             # if pred_formula_all,
             #   direct link formula result for peaks with annotation
             # else
             #   link formula result for all peaks within same peak group

             if (is_pred_formula_all) {
               idx <- match(unique(result_credential_formula$peak_name), names(annot_all))
               pb <- utils::txtProgressBar(min = 0, max = length(idx), style = 3)
               for (i in seq_along(idx)) {
                 utils::setTxtProgressBar(pb, i)
                 temp_idx <- idx[i]
                 temp_peak <- annot_all[[temp_idx]]@peak_info$name
                 annot_all[[temp_idx]]@credential_annotation$formula_prediction <- result_credential_formula %>% dplyr::filter(peak_name == temp_peak)
               }
             } else {

               # add base_peak formula result into all peaks
               temp_pg <- result_credential_pg %>% dplyr::filter(included_peak_group %in% unique(result_credential_formula$peak_group))
               idx <- match(unique(temp_pg$peak_name), names(annot_all))
               for (i in seq_along(idx)) {
                 utils::setTxtProgressBar(pb, i)
                 temp_idx <- idx[i]
                 temp_peak <- annot_all[[temp_idx]]@peak_info$name

                 # extract all related formula result of this peak
                 temp_pg_name <- temp_pg %>%
                   dplyr::filter(peak_name == temp_peak) %>%
                   dplyr::pull(included_peak_group) %>%
                   unique()

                 temp_formula_result <-  result_credential_formula %>% dplyr::filter(peak_group %in% temp_pg_name)

                 if (nrow(temp_formula_result) > 0) {
                   annot_all[[temp_idx]]@credential_annotation$formula_prediction <- temp_formula_result
                 } else {
                   next()
                 }

               }

             }

             return(annot_all)

           })




################################################################################
#   generateResultTable4AllAnnotation --------------------------------------------
#' @title generateResultTable4AllAnnotation
#' @description generate result table combining initial seed, recursive annotation, and annotation credential
#' @author Zhiwei Zhou
#' @param annot_all
#' @param path '.'
#' @param thread Default: 4
#' @param is_cred_pg_filter
#' @param is_cred_formula_filter
#' @export

# test <- generateResultTable4AllAnnotation(annot_all = annot_all,
#                                           path = '/home/zhouzw/Data_processing/20201113_rerun_200std_pos/modified_package_result',
#                                           lib = 'zhuMetLib',
#                                           direction = 'reverse',
#                                           thread = 3,
#                                           is_cred_pg_filter = FALSE,
#                                           is_cred_formula_filter = FALSE)
#
# test <- generateResultTable4AllAnnotation(annot_all = annot_all,
#                                           path = '/home/zhouzw/Data_processing/20210522_debug/01_aging_fly_pos/',
#                                           lib = 'zhuMetLib',
#                                           direction = 'reverse',
#                                           is_cred_pg_filter = FALSE,
#                                           is_cred_formula_filter = FALSE)

setGeneric(name = 'generateResultTable4AllAnnotation',
           def = function(
             annot_all,
             path = '.',
             lib = c('zhumetlib_qtof', 'zhumetlib_orbitrap', 'fiehnHilicLib'),
             thread = 4,
             direction = c('forward', 'reverse'),
             is_cred_pg_filter = TRUE,
             is_cred_formula_filter = TRUE,
             tolerance_rt_range = 30,
             test_evaluation = c('No', '200STD', '46STD')
           ){
             # direction <- match.arg(direction)
             cat('\n');cat('Generating merged annotation table of MetDNA2\n')
             temp_colname_ms2 <- ifelse(direction == 'forward', 'ms2_score_forward', 'ms2_score_reverse')

             tempFun <- function(i){
               # library(dplyr)
               # library(tibble)
               # library(tidyr)

               peak_info <- tibble::tibble(peak_name = as.character(annot_all[[i]]@peak_info$name),
                                           mz = annot_all[[i]]@peak_info$mz,
                                           rt = annot_all[[i]]@peak_info$rt,
                                           ccs = annot_all[[i]]@peak_info$ccs)

               # modify formats for initial_seed_annotation and recursive_annotation
               if (nrow(annot_all[[i]]@initial_seed_annotation) > 0) {
                 annot_initial_seed <- annot_all[[i]]@initial_seed_annotation %>%
                   # dplyr::select(id:adduct, mz_error:rt_error_rela, ms2_score_forward, cpd_type) %>%
                   dplyr::select(id:adduct, isotope, mz_error:ccs_error, temp_colname_ms2, ms2_matched_frag, cpd_type) %>%
                   dplyr::mutate(from_peak = as.character(NA),
                                 from_annotation = as.character(NA),
                                 round = NA,
                                 total_score = NA,
                                 recursive_type = NA,
                                 as_seed = NA,
                                 as_seed_round = NA,
                                 source = 'initial_seed',
                                 confidence_level = 'level2') %>%
                   dplyr::rename(ms2_score = temp_colname_ms2) %>%
                   dplyr::select(source, id, id_zhulab, name, formula, confidence_level, smiles, inchikey, inchikey1, adduct, isotope, mz_error,
                                 rt_error_abs, rt_error_rela, ccs_error, ms2_score, ms2_matched_frag, cpd_type, from_peak, from_annotation, recursive_type,
                                 as_seed, as_seed_round, round, total_score)

                 if (test_evaluation != '200STD') {
                   # replace ID with KEGG ID
                   #   if the compound not included in KEGG, use ZhuMetLib ID
                   temp_kegg_id <- match(annot_initial_seed$id, cpd_info$id) %>% cpd_info$id_kegg[.]
                   if (any(is.na(temp_kegg_id))) {
                     temp_idx <- which(is.na(temp_kegg_id))
                     temp_kegg_id[temp_idx] <- annot_initial_seed$id[temp_idx]
                   }

                   annot_initial_seed$id <- temp_kegg_id
                 }

                 if (test_evaluation != '46STD') {
                   # replace ID with KEGG ID
                   #   if the compound not included in KEGG, use ZhuMetLib ID
                   temp_kegg_id <- match(annot_initial_seed$id, cpd_info$id) %>% cpd_info$id_kegg[.]
                   if (any(is.na(temp_kegg_id))) {
                     temp_idx <- which(is.na(temp_kegg_id))
                     temp_kegg_id[temp_idx] <- annot_initial_seed$id[temp_idx]
                   }

                   annot_initial_seed$id <- temp_kegg_id
                 }

                 if (peak_info$rt <= 500) {
                   temp_rt_tolerance <- tolerance_rt_range
                 } else {
                   temp_rt_tolerance <- peak_info$rt*0.06
                 }

                 # idx_level1 <- which(!(is.na(annot_initial_seed$rt_error_abs)) & annot_initial_seed$rt_error_abs <= 30)
                 idx_level1 <- which(!(is.na(annot_initial_seed$rt_error_abs)) & annot_initial_seed$rt_error_abs <= temp_rt_tolerance)
                 annot_initial_seed$confidence_level[idx_level1] <- 'level1'
               } else {
                 annot_initial_seed <- tibble::tibble(source = character(),
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

               if (nrow(annot_all[[i]]@recursive_annotation) > 0) {
                 annot_recursive <- annot_all[[i]]@recursive_annotation %>%
                   dplyr::mutate(source = 'recursive_annotation',
                                 confidence_level = 'level3') %>%
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
               } else {
                 annot_recursive <- tibble::tibble(source = character(),
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
               initial_as_seed <- annot_initial_seed$id %>%
                 match(annot_recursive$id) %>%
                 annot_recursive$as_seed[.]

               initial_as_seed_round <- annot_initial_seed$id %>%
                 match(annot_recursive$id) %>%
                 annot_recursive$as_seed_round[.]

               annot_initial_seed$as_seed <- initial_as_seed
               annot_initial_seed$as_seed_round <- initial_as_seed_round

               # merge initial_seed_annotation and recursive annotation
               temp_annot_all <- annot_initial_seed %>%
                 dplyr::bind_rows(annot_recursive) %>%
                 dplyr::distinct(inchikey1, .keep_all = TRUE) %>%
                 dplyr::arrange(confidence_level)


               # Whether annotation is included in reserved peak_group
               if (nrow(annot_all[[i]]@credential_annotation$peak_group) > 0) {
                 # if the peak have annotation,
                 #    for each annotation,
                 #    check isotope, adduct whether appeared in the credential_annotation peak_group table
                 temp_result_pg <- annot_all[[i]]@credential_annotation$peak_group

                 if (nrow(temp_annot_all) > 0) {
                   temp_annot_all_isotope <- temp_annot_all$isotope
                   temp_annot_all_adduct <- temp_annot_all$adduct

                   annot_credential_pg <- sapply(seq_along(temp_annot_all_adduct), function(j){
                     if (temp_annot_all_isotope[j] != '[M]') {
                       result <- temp_annot_all_isotope[j] %in% temp_result_pg$label
                       return(result)
                     }

                     # for neutral loss ('[M-H2O+H]+', '[M-H2O-H]-' for annotation),
                     #   they are same as ISF actually
                     if (temp_annot_all_adduct[j] %in% c('[M-H2O+H]+', '[M-H2O-H]-')) {
                       result <- (temp_annot_all_adduct[j] %in% temp_result_pg$label) | (any(temp_result_pg$label == 'ISF'))
                       return(result)
                     }

                     result <- temp_annot_all_adduct[j] %in% temp_result_pg$label
                     return(result)
                   })

                   temp_annot_all <- temp_annot_all %>%
                     dplyr::mutate(is_credential_pg = annot_credential_pg)

                 } else {
                   temp_annot_all <- temp_annot_all %>%
                     dplyr::mutate(is_credential_pg = FALSE)
                 }

               } else {
                 temp_annot_all <- temp_annot_all %>% dplyr::mutate(is_credential_pg = FALSE)
               }


               # Whether annotation is included in top 3 formula
               if (nrow(annot_all[[i]]@credential_annotation$formula_prediction) > 0) {

                 temp_formula_result <- annot_all[[i]]@credential_annotation$formula_prediction %>%
                   dplyr::arrange(dplyr::desc(score))

                 # whether included in credential formula
                 annot_credential_formula <- temp_formula_result %>% dplyr::pull(formula)
                 annot_credential_formula <- temp_annot_all$formula %in% annot_credential_formula

                 # provide the top1 formula
                 temp_best_formula_result <- temp_formula_result %>%
                   dplyr::filter(dplyr::row_number() == 1)
                 temp_best_formula <- temp_best_formula_result$formula
                 temp_best_formula_adduct <- temp_best_formula_result$adduct


                 temp_annot_all <- temp_annot_all %>%
                   dplyr::mutate(is_credential_formula = annot_credential_formula,
                                 best_pred_formula = temp_best_formula,
                                 best_pred_formula_adduct = temp_best_formula_adduct)

               } else {
                 temp_annot_all <- temp_annot_all %>%
                   dplyr::mutate(is_credential_formula = FALSE,
                                 best_pred_formula = NA,
                                 best_pred_formula_adduct = NA)
               }

               # Impute NA for features with no annotation
               if (nrow(temp_annot_all) > 0) {
                 temp_annot_all <- peak_info %>% dplyr::bind_cols(temp_annot_all)
               } else {
                 temp <- matrix(rep(NA, ncol(temp_annot_all)), nrow = 1) %>% tibble::as_tibble()
                 colnames(temp) <- colnames(temp_annot_all)
                 temp_annot_all <- peak_info %>% dplyr::bind_cols(temp)
               }

               return(temp_annot_all)
             }

             library(parallel)
             # cl <- makeCluster(3L)
             cl <- makeCluster(thread)

             if (lib %in% c("zhumetlib_qtof", 'zhumetlib_orbitrap')) {
               data("zhuMetlib", envir = environment())
               cpd_info <- zhuMetlib$meta$compound
             }

             if (lib == "fiehnHilicLib") {
               data("fiehnHilicLib", envir = environment())
               cpd_info <- fiehnHilicLib$meta$compound
             }

             clusterExport(cl, c("tempFun", 'annot_all', 'cpd_info', 'tolerance_rt_range', 'temp_colname_ms2', '%>%', 'test_evaluation'),
                           envir = environment())
             system.time(
               id_merge <- parLapply(cl = cl,
                                     seq_along(annot_all),
                                     function(i) tempFun(i))
             )
             stopCluster(cl)


             # path_output <- file.path(path, 'annot_final')
             # dir.create(file.path(path, 'intermediate_data'), showWarnings = FALSE, recursive = TRUE)
             # save(id_merge, file = file.path(path, 'intermediate_data', 'id_merge_without_filter'), version = 2)

             # application of credential filters
             id_merge <- id_merge %>% dplyr::bind_rows()

             if (is_cred_pg_filter) {
               id_merge <- id_merge %>%
                 dplyr::filter(is.na(is_credential_pg) | is_credential_pg)
             }

             if (is_cred_formula_filter) {
               id_merge <- id_merge %>%
                 dplyr::filter(is.na(is_credential_formula) | is_credential_formula)
             }

             return(id_merge)

           })
