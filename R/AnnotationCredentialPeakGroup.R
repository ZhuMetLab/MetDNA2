################################################################################
# authenticatePeakGroupCredential ----------------------------------------------

#' @title authenticatePeakGroupCredential
#' @author Zhiwei Zhou
#' @description semi-targeted feature clustering (isotope, adduct, neutral loss, in-source fragments)
#' @param list_peak_group
#' @param ms2_data
#' @param path_dir '.'
#' @param polarity ionzation polarity, 'positive', 'negative'; Default: 'positive'
#' @param tol_mz mz tolerance. Default: 10 ppm
#' @param isotope_int_ratio_check whether check isotope intensity; Default: TRUE
#' @param isotope_int_ratio_cutoff isotope intensity ratio cutoff; Default: 500%
#' @param is_ms2_check whether compare ms2 of NL with base peak; Default: TRUE
#' @param ms2_score_cutoff Default: -1; # -1 represent not filter
#' @param is_plot_pseudo_MS1 whether output the pseudo MS1 spectrum.
#' @param cutoff_ssc
#' @export
#' @examples
#' load(system.file("tempdata", "list_peak_group_200805.RData", package="MetDNA2"))
#' load(system.file("tempdata", "raw_msms_200805.RData", package="MetDNA2"))
#' authenticatePeakGroupCredential(list_peak_group = list_peak_group,
#'                                 ms2_data = raw_msms,
#'                                 path_dir = 'I:/00_projects/03_MetDNA2/00_data/20200805_targeted_annotation_evaluation_200STD/00_demo_data',
#'                                 polarity = 'negative',
#'                                 tol_mz = 10,
#'                                 isotope_int_ratio_check = TRUE,
#'                                 isotope_int_ratio_cutoff = 500,
#'                                 is_ms2_check = TRUE,
#'                                 ms2_score_cutoff = -1,
#'                                 is_plot_pseudo_MS1 = TRUE)



# load('I:/00_projects/03_MetDNA2/00_data/20200710_semi_targeted_annotation_test/03_test/list_peak_group_200716.RData')
# load('I:/00_projects/03_MetDNA2/00_data/20200710_semi_targeted_annotation_test/03_test/raw_msms_200715.RData')

# polarity = 'negative'
# tol_mz = 10
# # para annotateIsotope
# isotope_int_ratio_check = TRUE
# isotope_int_ratio_cutoff = 500
# is_ms2_check = TRUE # compare ms2 similarity of neutral loss
# ms2_score_cutoff = -1 # -1 represent not filter

# path_dir <- ''
# load('./inst/tempdata/list_peak_group_200805.RData')
# load('./inst/tempdata/raw_msms_200805.RData')
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
# load('/home/zhouzw/Data_processing/20210319_different_instrument_test/QE/03_annotation_credential/00_intermediate_data/list_peak_group.RData')
# load('/home/zhouzw/Data_processing/20210319_different_instrument_test/QE/ms2_data.RData')
# ms2_data <- raw_msms
# path_dir <- '/home/zhouzw/Data_processing/20210319_different_instrument_test/QE/'
# polarity <- 'positive'
# tol_mz <- 25
# isotope_int_ratio_check = TRUE
# isotope_int_ratio_cutoff = 500
# is_ms2_check = TRUE # compare ms2 similarity of neutral loss
# ms2_score_cutoff = -1 # -1 represent not filter
# cutoff_ssc = 0.3
# cutoff_ssc_int = 3000
# is_rule_limitation = TRUE
# cutoff_topN = 5
# is_plot_pseudo_MS1 = TRUE
# thread = 4
# # type_order = c('seed', 'none_seed')
# type_order = c('level1', 'level2', 'level3')

setGeneric(name = 'authenticatePeakGroupCredential',
           def = function(
             list_peak_group,
             ms2_data,
             path_dir = '.',
             polarity = c('positive', 'negative'),
             tol_mz = 25,
             # para annotateIsotope
             isotope_int_ratio_check = TRUE,
             isotope_int_ratio_cutoff = 500,
             is_ms2_check = TRUE, # compare ms2 similarity of neutral loss
             ms2_score_cutoff = -1, # -1 represent not filter
             cutoff_ssc = 0.3,
             cutoff_ssc_int = 3000,
             is_rule_limitation = TRUE,
             cutoff_topN = 5,
             is_plot_pseudo_MS1 = TRUE,
             thread = 4,
             type_order = c('level1', 'level2', 'level3'), # the order of annotation credential
             ...
           ){
             match.arg(polarity)

             temp <- statNumListPeakGroup(list_peak_group)
             num_stat <- temp


             cat('Annotate isotopes, adducts, neutral loss and in-source fragments for each group\n')

             # progress <- mapProgress(n = length(list_peak_group))
             # list_peak_group_annotation <- purrr::map(seq_along(list_peak_group), function(i){

             # list_peak_group_annotation <- pbapply::pblapply(seq_along(list_peak_group), function(i){
             #   # mapProgressPrint(progress)
             #   cat(i, ' ')
             #   result <- annotatePeakGroup(peak_group = list_peak_group[[i]],
             #                               ms2_data = ms2_data,
             #                               polarity = polarity,
             #                               tol_mz = tol_mz,
             #                               isotope_int_ratio_check = isotope_int_ratio_check,
             #                               isotope_int_ratio_cutoff = isotope_int_ratio_cutoff,
             #                               is_ms2_check = is_ms2_check,
             #                               ms2_score_cutoff = ms2_score_cutoff,
             #                               cutoff_ssc = cutoff_ssc,
             #                               cutoff_ssc_int = cutoff_ssc_int,
             #                               is_rule_limitation = is_rule_limitation,
             #                               cutoff_topN = cutoff_topN)
             #
             #   return(result)
             #
             # })


             tempFun <- function(i){
               library(dplyr)
               library(SpectraTools)
               result <- annotatePeakGroup(peak_group = list_peak_group[[i]],
                                           ms2_data = ms2_data,
                                           polarity = polarity,
                                           tol_mz = tol_mz,
                                           isotope_int_ratio_check = isotope_int_ratio_check,
                                           isotope_int_ratio_cutoff = isotope_int_ratio_cutoff,
                                           is_ms2_check = is_ms2_check,
                                           ms2_score_cutoff = ms2_score_cutoff,
                                           cutoff_ssc = cutoff_ssc,
                                           cutoff_ssc_int = cutoff_ssc_int,
                                           is_rule_limitation = is_rule_limitation,
                                           cutoff_topN = cutoff_topN)

               return(result)
             }


             data('lib_adduct_nl', envir = environment())

             library(parallel)
             cl <- parallel::makeCluster(thread)
             parallel::clusterExport(cl, c('tempFun',
                                           'annotatePeakGroup',
                                           'initializePeakList',
                                           'annotateIsotope',
                                           'annotateAdduct',
                                           'annotateNeutralLoss',
                                           'annotateISF',
                                           'writePeakListAnnotation',
                                           'extractIsotopeNew',
                                           'mergePeakListAnnotation',
                                           'getMzRange',
                                           'simulateTheoIsotope',
                                           'calculateExactMass',
                                           'calculateMz',
                                           'transformMz',
                                           'convertMz2Adduct',
                                           'generateRuleList',
                                           'convertSpectraData',
                                           '%>%',
                                           'list_peak_group',
                                           'ms2_data',
                                           'polarity',
                                           'tol_mz',
                                           'isotope_int_ratio_check',
                                           'isotope_int_ratio_cutoff',
                                           'is_ms2_check',
                                           'ms2_score_cutoff',
                                           'cutoff_ssc',
                                           'cutoff_ssc_int',
                                           'is_rule_limitation',
                                           'cutoff_topN',
                                           'lib_adduct_nl'
                                           ),
                                     envir = environment())

             system.time(
               list_peak_group_annotation <- parallel::parLapply(cl = cl,
                                                                 seq_along(list_peak_group),
                                                                 function(i) {tempFun(i)}
               )
             )
             stopCluster(cl)

             names(list_peak_group_annotation) <- names(list_peak_group)

             dir.create(file.path(path_dir, '03_annotation_credential', "00_intermediate_data"),
                        showWarnings = FALSE,
                        recursive = TRUE)

             save(list_peak_group_annotation,
                  file = file.path(path_dir,
                                   '03_annotation_credential',
                                   "00_intermediate_data",
                                   'list_peak_group_annotation.RData'),
                  version = 2, compress = 'gzip')

             temp <- statNumListPeakGroup(list_peak_group_annotation)
             num_stat <- num_stat %>% dplyr::bind_rows(temp)

             temp <- statTypeListPeakGroup(list_peak_group_annotation)
             type_stat <- temp

             cat('\n\n');cat('Concise annotated peak groups ...\n')
             # reserve more confident peak group for one peak

             # 01 remove peak groups not meet the rule requirement
             cat('Remove peak_group filtered by rules\n')
             idx_rm <- sapply(list_peak_group_annotation, function(x){
               nrow(x@peak_list_annotated)
             })
             idx_rm <- which(idx_rm == 0)
             if (length(idx_rm) > 0) {
               list_peak_group_annotation_concised <- list_peak_group_annotation[-idx_rm]
               record_rule_filter_peak_group <- names(list_peak_group_annotation)[idx_rm]

               save(record_rule_filter_peak_group,
                    file = file.path(path_dir,
                                     '03_annotation_credential',
                                     "00_intermediate_data",
                                     'record_rule_filter_peak_group.RData'))
             } else {
               list_peak_group_annotation_concised <- list_peak_group_annotation
               record_rule_filter_peak_group <- NULL

               save(record_rule_filter_peak_group,
                    file = file.path(path_dir,
                                     '03_annotation_credential',
                                     "00_intermediate_data",
                                     'record_rule_filter_peak_group.RData'))
             }

             temp <- statNumListPeakGroup(list_peak_group_annotation_concised)
             num_stat <- num_stat %>% dplyr::bind_rows(temp)

             temp <- statTypeListPeakGroup(list_peak_group_annotation_concised)
             type_stat <- type_stat %>% dplyr::bind_rows(temp)

             # 02 remove conflict peak group
             cat('Remove conflict peak group\n')
             list_peak_group_annotation_concised <- concisePeak2PeakGroup(
               list_peak_group_annotation = list_peak_group_annotation_concised,
               type_order = type_order,
               path_dir = path_dir)

             temp <- statNumListPeakGroup(list_peak_group_annotation_concised)
             num_stat <- num_stat %>% dplyr::bind_rows(temp)

             temp <- statTypeListPeakGroup(list_peak_group_annotation_concised)
             type_stat <- type_stat %>% dplyr::bind_rows(temp)

             # 03 remove overlap peak group
             cat('Remove overlap peak group\n')
             list_peak_group_annotation_concised <- concisePeakGroup2PeakGroup(
               list_peak_group_annotation = list_peak_group_annotation_concised,
               type_order = type_order,
               path_dir = path_dir)

             temp <- statNumListPeakGroup(list_peak_group_annotation_concised)
             num_stat <- num_stat %>% dplyr::bind_rows(temp)

             temp <- statTypeListPeakGroup(list_peak_group_annotation_concised)
             type_stat <- type_stat %>% dplyr::bind_rows(temp)

             rownames(num_stat) <- c('input',
                                     'initial annotation',
                                     'rule removal',
                                     'conflict removal',
                                     'overlap removal')

             rownames(type_stat) <- c('initial annotation',
                                      'rule removal',
                                      'conflict removal',
                                      'overlap removal')

             cat('General statistics:')
             print(knitr::kable(num_stat)); cat('\n\n')

             cat('Annotation type summary:')
             print(knitr::kable(type_stat)); cat('\n\n')

             dir.create(file.path(path_dir,
                                  '03_annotation_credential',
                                  "00_intermediate_data"),
                        showWarnings = FALSE,
                        recursive = TRUE)

             save(list_peak_group_annotation_concised,
                  file = file.path(path_dir,
                                   '03_annotation_credential',
                                   "00_intermediate_data",
                                   'list_peak_group_annotation_concised.RData'))

             result_peak_group_concise <- concisePeakGroup2Annotation(
               list_peak_group_annotation = list_peak_group_annotation_concised)

             dir.create(file.path(path_dir,
                                  '03_annotation_credential',
                                  "00_intermediate_data"),
                        showWarnings = FALSE,
                        recursive = TRUE)

             save(result_peak_group_concise,
                  file = file.path(path_dir,
                                   '03_annotation_credential',
                                   "00_intermediate_data",
                                   'result_peak_group_concise.RData'))

             dir.create(file.path(path_dir,
                                  '03_annotation_credential'),
                        showWarnings = FALSE, recursive = TRUE)

             readr::write_csv(result_peak_group_concise$result_annotation_table,
                              path = file.path(path_dir,
                                               '03_annotation_credential',
                                               # '01_annotation_result',
                                               'semi_targeted_annotation_result.csv'))

             if (is_plot_pseudo_MS1) {
               cat('Plot pseudo MS1 spec...\n')
               dir.create(file.path(path_dir,
                                    '03_annotation_credential',
                                    '01_pseudo_MS1_spec'),
                          showWarnings = FALSE, recursive = TRUE)

               temp_fun <- function(i, list_peak_group_annotation, plotPseudoMs1Spec,
                                    `%>%`, path_dir) {

                 temp_plot <- list_peak_group_annotation[[i]] %>% plotPseudoMs1Spec()

                 temp_plot_name <- names(list_peak_group_annotation)[i] %>% paste0('.pdf')
                 ggplot2::ggsave(temp_plot,
                                 filename = file.path(path_dir,
                                                      '03_annotation_credential',
                                                      '01_pseudo_MS1_spec',
                                                      temp_plot_name),
                                 width = 10, height = 6)


               }

               idx <- sapply(list_peak_group_annotation, function(x){
                 nrow(x@peak_list_annotated)
               })
               idx <- which(idx > 0)
               temp_list_peak_group_annotation <- list_peak_group_annotation[idx]

               temp <- BiocParallel::bplapply(X = seq_along(temp_list_peak_group_annotation),
                                              FUN = temp_fun,
                                              BPPARAM = BiocParallel::SnowParam(workers = thread,
                                                                                progressbar = TRUE),
                                              list_peak_group_annotation = temp_list_peak_group_annotation,
                                              plotPseudoMs1Spec = plotPseudoMs1Spec,
                                              `%>%` = `%>%`,
                                              path_dir = path_dir)

               rm(list = 'temp');gc()

               # progress <- mapProgress(n = length(list_peak_group_annotation))
               # purrr::walk(seq_along(list_peak_group_annotation), function(i){
               #   # cat(i, ' ')
               #   mapProgressPrint(progress = progress)
               #
               #   temp_plot <- list_peak_group_annotation[[i]] %>% plotPseudoMs1Spec()
               #
               #   temp_plot_name <- names(list_peak_group_annotation)[i] %>% paste0('.pdf')
               #   ggplot2::ggsave(temp_plot,
               #                   filename = file.path(path_dir,
               #                                        'annot_credential',
               #                                        'pseudo_MS1_spec',
               #                                        temp_plot_name),
               #                   width = 10, height = 6)
               # })
             }



           })




################################################################################
# generate peak group ----------------------------------------------------------
setClass(Class = "PeakGroup",
         slots = list(base_peak_name = "character",
                      base_peak_mz = "numeric",
                      base_peak_rt = "numeric",
                      base_peak_ccs = 'numeric',
                      base_peak_adduct = 'character',
                      base_peak_seed = 'character',
                      group_size = 'numeric',
                      peak_list = 'data.frame',
                      peak_list_annotated = 'data.frame')
)


#' @title generatePeakGroup
#' @author Zhiwei Zhou
#' @param peak_table peak table (row--peaks; column--samples). First 3 column should be named as "name", "mz", "rt"
#' @param base_peak_name name of base peak
#' @param base_peak_adduct adduct format of base peak
#' @param base_peak_seed seed of base peak
#' @param tol_rt rt tolerance for peak grouping; Default: 3 (second)
#' @param cutoff_ssc sample-sample correlation; Default: 0
#' @export
#' @examples
#' peak_table <- read.csv(system.file("extdata", "peak_table_200STD_neg_200805.csv", package="MetDNA2"), stringsAsFactors = F)
#' test <- generatePeakGroup(peak_table = peak_table,
#'                           base_peak_name = 'M482T929',
#'                           base_peak_adduct = '[M-H]-',
#'                           tol_rt = 3,
#'                           cutoff_ssc = 0)

# peak_table <- readr::read_csv('./inst/extdata/peak_table_200STD_neg_200805.csv')
# peak_table_annotated <- readr::read_csv('./inst/extdata/peak_table_annotated_200STD_neg_200805.csv')
# base_peak_name <- 'M482T929'
# base_peak_adduct <- '[M-H]-'
# test <- generatePeakGroup(peak_table = peak_table,
#                           base_peak_name = 'M482T929',
#                           base_peak_adduct = '[M-H]-',
#                           tol_rt = 3,
#                           cutoff_ssc = 0)

setGeneric(name = 'generatePeakGroup',
           def = function(
             peak_table,
             base_peak_name,
             base_peak_adduct,
             base_peak_seed,
             tol_rt = 3,
             cutoff_ssc = 0 # cutoff of sample_sample_correlation
           ){
             base_peak <- peak_table %>%
               dplyr::filter(name == base_peak_name)

             base_rt <- base_peak$rt
             base_ccs <- base_peak$ccs

             peak_list <- peak_table %>%
               dplyr::filter(rt >= (base_rt-tol_rt) & rt <= (base_rt+tol_rt))


             ssc <- apply(peak_list, 1, function(x){
               ssc = cor(as.numeric(x[-c(1:4)]),
                         as.numeric(base_peak[-c(1:4)]),
                         method = 'pearson')
             })

             int <- apply(peak_list, 1, function(x){
               ssc = mean(as.numeric(x[-c(1:4)]))
             })

             # ssc <- apply(peak_list, 1, function(x){
             #   ssc = cor(as.numeric(x[-c(1:3)]),
             #             as.numeric(base_peak[-c(1:3)]),
             #             method = 'pearson')
             # })
             #
             # int <- apply(peak_list, 1, function(x){
             #   ssc = mean(as.numeric(x[-c(1:3)]))
             # })

             peak_list <- peak_list %>%
               dplyr::mutate(mz = as.numeric(mz),
                             rt = as.numeric(rt),
                             ccs = as.numeric(ccs)) %>%
               dplyr::mutate(ssc = ssc,
                             int = int) %>%
               dplyr::select(name, mz, rt, ccs, ssc, int) %>%
               dplyr::rename(peak_name = name)

             if (cutoff_ssc > 0) {
               peak_list <- peak_list %>%
                 dplyr::filter(abs(ssc) >= cutoff_ssc)
             }

             result_peak_group <- new('PeakGroup',
                                      base_peak_name = base_peak_name,
                                      base_peak_mz = as.numeric(base_peak$mz),
                                      base_peak_rt = as.numeric(base_peak$rt),
                                      base_peak_ccs = as.numeric(base_peak$ccs),
                                      base_peak_seed = base_peak_seed,
                                      base_peak_adduct = base_peak_adduct,
                                      group_size = nrow(peak_list),
                                      peak_list = peak_list)

             return(result_peak_group)
           })


setMethod(f = "show",
          signature = "PeakGroup",
          definition = function(object){
            # cat("-----------Meta information------------\n")
            cat("Base peak name:", object@base_peak_name, "\n")
            cat("m/z:", object@base_peak_mz, "\n")
            cat("RT:", object@base_peak_rt, "\n")
            cat("CCS:", object@base_peak_ccs, "\n")
            cat("Adduct:", object@base_peak_adduct, "\n")
            cat("Seed:", object@base_peak_seed, "\n")
            cat("Number of peak in group:", nrow(object@peak_list), "\n")
            cat("Number of annotated peak in group:", nrow(object@peak_list_annotated), "\n")
          }
)

# setMethod(f = "showPeakGroup",
#           signature = "PeakGroup",
#           definition = function(object){
#             cat("Base peak info: Name", object@base_peak_name,
#                 'Adduct', object@base_peak_adduct, "\n")
#
#             cat(object@peak_group)
#
#           }
# )






# Semi-annotatin ---------------------------------------------------------------
# setClass(Class = 'PeakGroupAnnotation',
#          slots = list(peak_list_annotated = 'data.frame'),
#          contains = 'PeakGroup')

setGeneric(name = 'initializePeakList',
           def = function(
             peak_list
           ){
             peak_list <- peak_list %>%
               dplyr::mutate(mz_error = NA,
                             rt_error = NA,
                             ccs_error = NA,
                             int_ratio = NA,
                             ms2_score = NA,
                             type = NA,
                             annotation = NA,
                             query_peak = NA)

             return(peak_list)
           })

setGeneric(name = 'writePeakListAnnotation',
           def = function(pl_annotation_old,
                          pl_annotation_new){
             if (missing(pl_annotation_new)) {
               stop('Please input pl_annotation_new')
             }

             if (missing(pl_annotation_old)) {
               pl_annotation_old <- tibble::tibble(peak_name = character(),
                                                   mz = numeric(),
                                                   rt = numeric(),
                                                   ccs = numeric(),
                                                   ssc = numeric(),
                                                   int = numeric(),
                                                   mz_error = numeric(),
                                                   rt_error = numeric(),
                                                   ccs_error = numeric(),
                                                   int_ratio = numeric(),
                                                   ms2_score = numeric(),
                                                   type = character(),
                                                   annotation = character(),
                                                   query_peak = character())
             }

             result <- pl_annotation_old %>% dplyr::bind_rows(pl_annotation_new)
             return(result)
           })

setGeneric(name = 'extractIsotopeNew',
           def = function(
             pl_annotation
           ){
             peak_name_isotope <- pl_annotation %>%
               dplyr::filter(type == 'isotopeAnnotation') %>%
               dplyr::pull(peak_name)

             peak_name_other <- pl_annotation %>%
               dplyr::filter(type != 'isotopeAnnotation')

             peak_name_other <- peak_name_other  %>%
               dplyr::filter(!(peak_name %in% peak_name_isotope)) %>%
               dplyr::pull(peak_name)

             # # when these pair existed: M+_[M+H]+ & M-_[M-H]-
             # #   M+ and M- will not used for annotated its isotopes
             # if (any(c('M+', 'M-') %in% peak_name_other$annotation)) {
             #
             #   if (any(c('[M+H]+', '[M-H]-') %in% peak_name_other$annotation)) {
             #     peak_name_other <- peak_name_other  %>%
             #       dplyr::filter(!(peak_name %in% peak_name_isotope)) %>%
             #       dplyr::filter(!(annotation %in% c('M+', 'M-'))) %>% # avoid [M+H]+ annotation
             #       dplyr::pull(peak_name)
             #   } else {
             #     peak_name_other <- peak_name_other  %>%
             #       dplyr::filter(!(peak_name %in% peak_name_isotope)) %>%
             #       dplyr::pull(peak_name)
             #   }
             #
             # } else {
             #   peak_name_other <- peak_name_other  %>%
             #     dplyr::filter(!(peak_name %in% peak_name_isotope)) %>%
             #     dplyr::pull(peak_name)
             # }


             if (length(peak_name_other) == 0) {
               return(NULL)
             }

             peak_name_other
           })





# annotateIsotope --------------------------------------------------------------
# example L0076 neg in 200STD
# load('I:/00_projects/03_MetDNA2/00_data/20200710_semi_targeted_annotation_test/03_test/list_peak_group_200716.RData')
# peak_group <- list_peak_group$`M482T929_[M-H]-`
# peak_list <- peak_group@peak_list %>% initializePeakList()
# query_peak_name <- peak_group@base_peak_name
# query_mz <- peak_group@base_peak_mz
# query_rt <- peak_group@base_peak_rt
# query_adduct <- peak_group@base_peak_adduct
# annotateIsotope(peak_list = peak_list,
#                 query_mz = query_mz,
#                 query_rt = query_rt,
#                 query_peak_name = query_peak_name,
#                 tol_mz = 10)

# peak_group <- list_peak_group$`M204T364_[M-H]-`
# peak_list <- peak_group@peak_list %>% initializePeakList()
# query_peak_name <- 'M159T364'
# query_mz <- 159.0322
# query_rt <- 364.199
# annotateIsotope(peak_list = peak_list,
#                 query_mz = query_mz,
#                 query_rt = query_rt,
#                 query_peak_name = query_peak_name,
#                 tol_mz = 10,
#                 isotope_int_ratio_check = TRUE,
#                 isotope_int_ratio_cutoff = 500)
#
# annotateIsotope(peak_list = peak_list,
#                 query_mz = 148.0381,
#                 query_rt = 65,
#                 query_ccs = 149.9,
#                 query_peak_name = "M148T65C150",
#                 tol_mz = 10)


setGeneric(name = 'annotateIsotope',
           def = function(
             peak_list,
             query_mz,
             query_rt,
             query_ccs,
             query_peak_name,
             # query_adduct,
             tol_mz = 10,
             isotope_delta = 1.003355,
             isotope_max_num = 4,
             isotope_int_ratio_check = TRUE,
             isotope_int_ratio_cutoff = 500,
             monotonic_dec_check = FALSE,
             monotonic_mz_cutoff = 800,
             cutoff_ssc = NULL,
             cutoff_ssc_int = 3000,
             ...
           ){

             # browser()
             mz_isotope_list <- query_mz + seq(0, isotope_max_num-1)*isotope_delta
             isotope_matrix <- getMzRange(mz = mz_isotope_list,
                                          ppm = 25,
                                          mz_ppm_thr = 400)

             # keep peaks with ssc larger than cutoff only
             if (length(cutoff_ssc) > 0) {
               peak_list <- peak_list %>%
                 dplyr::filter((ssc >= cutoff_ssc) | (int <= cutoff_ssc_int))
             }

             result_isotope <- lapply(seq_along(mz_isotope_list), function(i){
               mz_min <- isotope_matrix[i,1]
               mz_max <- isotope_matrix[i,2]

               if (i==1) {
                 label <- "[M]"
               } else {
                 label <- paste0("[M+", i-1, "]")
               }

               result <- peak_list %>%
                 dplyr::filter(mz >= mz_min & mz<= mz_max) %>%
                 dplyr::mutate(mz_error = as.numeric(mz - mz_isotope_list[i]),
                               rt_error = as.numeric(rt - query_rt),
                               ccs_error = as.numeric(ccs - query_ccs),
                               type = 'isotopeAnnotation',
                               annotation = label,
                               query_peak = query_peak_name)
             })

             result_isotope <- result_isotope %>% dplyr::bind_rows()

             # If multiple isotope has been annotated,
             #    reserve peak with minimum mz error
             result_isotope <- result_isotope %>%
               dplyr::arrange(annotation, abs(mz_error), abs(rt_error), abs(ccs_error)) %>%
               dplyr::distinct(annotation, .keep_all = TRUE) %>%
               dplyr::mutate(int_ratio = int/int[1])

             if (isotope_int_ratio_check) {
               simulate_isotope <- simulateTheoIsotope(mz = query_mz,
                                                       isotope_max_num = isotope_max_num)

               simulate_isotope <- simulate_isotope %>%
                 dplyr::filter(label %in% result_isotope$annotation) %>%
                 dplyr::select(int_theo, label)

               temp_result_isotope <- result_isotope %>%
                 dplyr::left_join(simulate_isotope, by = c('annotation' = 'label'))

               result_isotope <- temp_result_isotope %>%
                 dplyr::rowwise() %>%
                 dplyr::mutate(int_error = abs(int_ratio-int_theo)/int_theo*100) %>%
                 dplyr::filter(int_error <= isotope_int_ratio_cutoff) %>%
                 dplyr::select(-c('int_theo', 'int_error')) %>%
                 tibble::as_tibble()

               # purrr::map2(result_isotope$int_ratio,
               #             simulate_isotope$int_theo,
               #             function(x, y){
               #
               #             })

             }

             # If monotonic_dec_check, the M+1 need to small the M.
             # Note: it is not suitable for molecules with halogen elements!
             if (monotonic_dec_check) {
               if (query_mz <= monotonic_mz_cutoff) {
                 idx <- sapply(seq_along(result_isotope$int_ratio),
                               function(i){
                                 if (i==1) return(TRUE)

                                 if (result_isotope$int_ratio[i] <
                                     result_isotope$int_ratio[i-1]) {
                                   return(TRUE)
                                 } else {
                                   return(FALSE)
                                 }
                               })

                 result_isotope <- result_isotope[idx,]
               }
             }

             return(result_isotope)
           })



setGeneric(name = 'getMzRange',
           def = function(mz,
                          ppm = 10,
                          mz_ppm_thr = 500){
             result <- sapply(mz, function(x) {
               if (x >= mz_ppm_thr) {
                 x * (1 + c(-1, 1) * ppm * 1e-6)
               } else {
                 temp1 <- x + mz_ppm_thr * c(-1, 1) * ppm * 1e-6
               }
             })

             t(result)
           }
)


# simulateTheoIsotope(mz = 481.9783, isotope_max_num = 4)
setGeneric(name = 'simulateTheoIsotope',
           function(
             mz,
             isotope_max_num = 4
           ){
             # calculate simulated carbon number with alkane (CnH2n+2)
             num_carbon <- (mz - 1.0078*2) %/% 14.0156
             simulate_alkane_formula <- paste0('C', num_carbon, 'H', 2*num_carbon+2)

             options(readr.num_columns = 0)
             simulate_isotope <- Rdisop::getMolecule(simulate_alkane_formula, z=1) %>%
               Rdisop::getIsotope(index = seq(isotope_max_num)) %>%
               t() %>%
               tibble::as_tibble() %>%
               dplyr::rename(mz = V1, int_theo = V2)

             if (isotope_max_num > 1) {
               label <- c('[M]', paste0("[M+", seq(isotope_max_num-1), "]"))
             } else {
               label <- '[M]'
             }

             simulate_isotope <- simulate_isotope %>%
               dplyr::mutate(label = label) %>%
               dplyr::mutate(int_theo = int_theo/int_theo[1])

             return(simulate_isotope)

           })


# annotateAdduct ---------------------------------------------------------------

# # example L0076 neg in 200STD
# load('./lib_adduct_nl_200714.RData')
# load('I:/00_projects/03_MetDNA2/00_data/20200710_semi_targeted_annotation_test/03_test/list_peak_group_200716.RData')
#
# peak_group <- list_peak_group$`M482T929_[M-H]-`
# peak_list <- peak_group@peak_list %>% initializePeakList()
# query_peak_name <- peak_group@base_peak_name
# query_mz <- peak_group@base_peak_mz
# query_rt <- peak_group@base_peak_rt
# query_adduct <- peak_group@base_peak_adduct
# polarity <- 'negative'
# annotateAdduct(peak_list = peak_list,
#                query_mz = query_mz,
#                query_rt = query_rt,
#                query_peak_name = query_peak_name,
#                query_adduct = query_adduct,
#                polarity = polarity,
#                tol_mz = 10)

# peak_group <- list_peak_group_annotation$`M116T562_[M+H]+`
# peak_group <- list_peak_group$`M184T910_[M]+`
# peak_list <- peak_group@peak_list %>% initializePeakList()
# query_peak_name <- peak_group@base_peak_name
# query_mz <- peak_group@base_peak_mz
# query_rt <- peak_group@base_peak_rt
# query_adduct <- peak_group@base_peak_adduct
# polarity <- 'positive'
# annotateAdduct(peak_list = peak_list,
#                query_mz = query_mz,
#                query_rt = query_rt,
#                query_peak_name = query_peak_name,
#                query_adduct = query_adduct,
#                polarity = polarity,
#                tol_mz = 25,
#                cutoff_ssc = 0.3,
#                cutoff_ssc_int = 3000)

setGeneric(name = 'annotateAdduct',
           def = function(
             peak_list,
             query_mz,
             query_rt,
             query_ccs,
             query_peak_name,
             query_adduct,
             polarity = c('positive', 'negative'),
             tol_mz = 25,
             cutoff_ssc = NULL,
             cutoff_ssc_int = 3000,
             is_rule_limitation = TRUE,
             # lib_adduct_nl = lib_adduct_nl,
             ...
           ){
             # browser()

             adduct_list <- convertMz2Adduct(base_mz = query_mz,
                                             base_adduct = query_adduct,
                                             type = 'adduct',
                                             polarity = polarity)

             mz_adduct_list <- adduct_list$mz
             name_adduct_list <- adduct_list$adduct

             adduct_matrix <- getMzRange(mz = mz_adduct_list,
                                         ppm = tol_mz,
                                         mz_ppm_thr = 400)

             # keep peaks with ssc larger than cutoff only (except int less than cutoff)
             if (length(cutoff_ssc) > 0) {
               peak_list <- peak_list %>%
                 dplyr::filter((ssc >= cutoff_ssc) | (int <= cutoff_ssc_int))
             }

             result_adduct <- lapply(seq_along(mz_adduct_list), function(i){
               mz_min <- adduct_matrix[i,1]
               mz_max <- adduct_matrix[i,2]

               result <- peak_list %>%
                 dplyr::filter(mz >= mz_min & mz<= mz_max) %>%
                 dplyr::mutate(mz_error = as.numeric(mz - mz_adduct_list[i]),
                               rt_error = as.numeric(rt - query_rt),
                               ccs_error = as.numeric(ccs - query_ccs),
                               type = 'adductAnnotation',
                               annotation = name_adduct_list[i],
                               query_peak = query_peak_name)
             })

             result_adduct <- result_adduct %>% dplyr::bind_rows()

             # If multiple adduct has been annotated,
             #    reserve peak with minimum mz error
             result_adduct <- result_adduct %>%
               dplyr::arrange(annotation, abs(mz_error), abs(rt_error), abs(ccs_error)) %>%
               dplyr::distinct(annotation, .keep_all = TRUE)

             # Check rule requirements for adducts
             if (is_rule_limitation) {
               rule_list <- generateRuleList(polarity = polarity)

               idx <- which(result_adduct$annotation %in% rule_list$rule_adduct)
               if (length(idx) > 0) {
                 idx_eff <- sapply(seq_along(idx), function(i){
                   temp_rule <- rule_list %>%
                     dplyr::filter(rule_adduct == result_adduct$annotation[idx[i]]) %>%
                     dplyr::pull(rule)

                   result <- all(temp_rule %in% result_adduct$annotation)
                 })

                 if (!all(idx_eff)) {
                   adduct_rm <- idx[!idx_eff] %>%
                     result_adduct$annotation[.]

                   result_adduct <- result_adduct %>%
                     dplyr::filter(annotation != adduct_rm)
                 }
               }

             }

             return(result_adduct)

           })


# lib_adduct_pos <- readr::read_csv('./lib_adduct_pos_200713_manual.csv')
# lib_adduct_neg <- readr::read_csv('./lib_adduct_neg_200713_manual.csv')
# lib_nl_pos <- readr::read_csv('./lib_nl_pos_200713_manual.csv')
# lib_nl_neg <- readr::read_csv('./lib_nl_neg_200713_manual.csv')
#
# lib_pos <- lib_adduct_pos %>% bind_rows(lib_nl_pos)
# lib_neg <- lib_adduct_neg %>% bind_rows(lib_nl_neg)
#
# lib_adduct_nl <- list(positive = lib_pos,
#                       negative = lib_neg)
#
# save(lib_adduct_nl, file = './lib_adduct_nl_200714.RData', version = 2)


# calculateExactMass(formula = "C2H5OH")

setGeneric(name = 'calculateExactMass',
           def = function(
             formula
           ){
             molecule <- Rdisop::getMolecule(formula)
             # getFormula(molecule)
             Rdisop::getMass(molecule)
           })


# exact_mass <- 180.0634
# adduct <- '[M-H]-'
# delta_mz <- -1.0073
# calculateMz(exact_mass = 180.0634,
#             adduct = '[M-H]-',
#             delta_mz = -1.0073)

# calculateMz(exact_mass = 180.0634,
#             adduct = '[2M-H]-',
#             delta_mz = -1.0073)
#
# calculateMz(exact_mass = 180.0634,
#             adduct = '[M-2H]2-',
#             delta_mz = -1.0073)


setGeneric(name = 'calculateMz',
           def = function(
             exact_mass,
             adduct,
             delta_mz,
             nmol = NULL,
             ncharge = NULL
           ){

             if (length(nmol) == 0) {
               if (stringr::str_detect(adduct, pattern = '2M')) {
                 mz <- exact_mass*2 + delta_mz
               } else if (stringr::str_detect(adduct, pattern = '3M')) {
                 mz <- exact_mass*3 + delta_mz
               } else {
                 mz <- exact_mass + delta_mz
               }
             } else {
               mz <- exact_mass*nmol + delta_mz
             }


             if (length(ncharge) == 0) {
               if (stringr::str_detect(adduct, pattern = '\\]2\\-|\\]2\\+')) {
                 mz <- mz/2
               } else if (stringr::str_detect(adduct, pattern = '\\]3\\-|\\]3\\+')) {
                 mz <- mz/3
               } else {
                 mz
               }
             } else {
               mz <- mz/ncharge
             }

             mz
           })


# exact_mass <- 180.0634
# adduct <- '[M-H]-'
# delta_mz <- -1.0073
# transformMz(exact_mass = 180.0634, type = 'adduct', polarity = 'positive')
# transformMz(exact_mass = 180.0634, type = 'nl', polarity = 'positive')
# transformMz(exact_mass = 180.0634, adduct_list = c('[M-H]-', '[M+H]+'))

setGeneric(name = 'transformMz',
           function(
             exact_mass,
             formula = NULL,
             adduct_list = NULL,
             type = c('adduct', 'nl'),
             polarity = c('positive', 'negative'),
             # lib_adduct_nl = lib_adduct_nl,
             ...
           ){

             if (all(is.null(exact_mass), is.null(formula))) {
               stop('Please input exact_mass or formula.')
             }

             if (!is.null(formula)) {
               exact_mass <- calculateExactMass(formula)
             }

             if (is.null(adduct_list)) {
               lib <- switch(polarity,
                             'positive' = {
                               lib_adduct_nl$positive
                             },
                             'negative' = {
                               lib_adduct_nl$negative
                             }
               )

               lib <- switch(type,
                             'adduct' = {
                               lib %>%
                                 dplyr::filter(type == 'Adduct') %>%
                                 dplyr::filter(credential == 'Yes')
                             },
                             'nl' = {
                               lib %>%
                                 dplyr::filter(type == 'NeutralLoss') %>%
                                 dplyr::filter(credential == 'Yes')
                             })
             } else {
               lib <- lib_adduct_nl$positive %>%
                 dplyr::bind_rows(lib_adduct_nl$negative)

               if (!all(adduct_list %in% lib$adduct)) {
                 stop('Sorry, not all adduct included in the adduct list\n')
               }

               lib <- lib %>%
                 dplyr::filter(adduct %in% adduct_list) %>%
                 dplyr::arrange(match(adduct, adduct_list))
             }


             result_mz <- sapply(seq_along(lib$adduct), function(i){
               calculateMz(exact_mass = exact_mass,
                           adduct = lib$adduct[i],
                           delta_mz = lib$delta_mz[i])
             })

             result <- tibble::tibble(exact_mass = exact_mass,
                                      adduct = lib$adduct,
                                      mz = result_mz)


             return(result)

           })




# convertMz2Adduct(base_mz = 181.0707,
#                  base_adduct = '[M+H]+',
#                  adduct_list = NULL,
#                  type = 'adduct',
#                  polarity = 'positive')

# convertMz2Adduct(base_mz = 181.0707,
#                  base_adduct = '[M+H]+',
#                  adduct_list = NULL,
#                  type = 'nl',
#                  polarity = 'negative')

# convertMz2Adduct(base_mz = 181.0707,
#                  base_adduct = '[M+H]+',
#                  adduct_list = c('[M-H2O+H]+', '[M+Na]+'),
#                  type = 'nl',
#                  polarity = 'negative')

# convertMz2Adduct(base_mz = 143.0347,
#                  base_adduct = '[2M-H]-',
#                  adduct_list = c('[M-H2O+H]+', '[M+Na]+'))

# convertMz2Adduct(base_mz = 143.0347,
#                  base_adduct = '[3M-H]-',
#                  adduct_list = c('[M-H2O+H]+', '[M+Na]+'))
# convertMz2Adduct(base_mz = 664.1144, base_adduct = '[M]+', adduct_list = '[M]+')$exact_mass

setGeneric(name = 'convertMz2Adduct',
           def = function(
             base_mz,
             base_adduct,
             # lib_adduct_nl = lib_adduct_nl,
             ...
             # adduct_list = NULL,
             # type = c('adduct', 'nl'),
             # polarity = c('positive', 'negative')
           ){

             lib <- lib_adduct_nl$positive %>%
               dplyr::bind_rows(lib_adduct_nl$negative)

             if (!(base_adduct %in% lib$adduct)) {
               stop('Sorry, base_adduct is not included\n')
             }

             temp_delta_mz <- lib %>%
               dplyr::filter(adduct == base_adduct) %>%
               dplyr::pull(delta_mz)


             if (stringr::str_detect(base_adduct, pattern = '2M')) {
               temp_exact_mass <- (base_mz - temp_delta_mz)/2
             } else if (stringr::str_detect(base_adduct, pattern = '3M')) {
               temp_exact_mass <- (base_mz - temp_delta_mz)/3
             } else {
               temp_exact_mass <- base_mz - temp_delta_mz
             }

             # temp_exact_mass <- base_mz - temp_delta_mz

             result <- transformMz(exact_mass = temp_exact_mass,
                                   ...)

             # if initial seed base_adduct is '[M]+', '[M-2H]-', the adduct would be added for credential
             if ((base_adduct %in% c('[M]+', '[M-2H]-')) & !(base_adduct %in% result$adduct)) {
               temp_result <- transformMz(exact_mass = temp_exact_mass, adduct_list = base_adduct)
               result <- temp_result %>% dplyr::bind_rows(result)
               return(result)
             }

             return(result)

           })


setGeneric(name = 'generateRuleList',
           def = function(
             polarity = c('positive', 'negative')
           ){
             lib <- switch(polarity,
                           'positive' = {
                             lib_adduct_nl$positive
                           },
                           'negative' = {
                             lib_adduct_nl$negative
                           }
             )

             result <- lib %>%
               dplyr::filter(!is.na(rule_limitation)) %>%
               dplyr::select(adduct, rule_limitation) %>%
               tidyr::separate_rows(rule_limitation, sep = ';') %>%
               dplyr::rename(rule_adduct = adduct,
                             rule = rule_limitation)

             return(result)
           })


# annotateNeutralLoss ----------------------------------------------------------
# # example L0076 neg in 200STD
# # load('./lib_adduct_nl_200714.RData')
# load('I:/00_projects/03_MetDNA2/00_data/20200710_semi_targeted_annotation_test/03_test/list_peak_group_200716.RData')
# load('I:/00_projects/03_MetDNA2/00_data/20200710_semi_targeted_annotation_test/03_test/raw_msms_200715.RData')
#
# peak_group <- list_peak_group$`M482T929_[M-H]-`
# peak_list <- peak_group@peak_list %>% initializePeakList()
# query_peak_name <- peak_group@base_peak_name
# query_mz <- peak_group@base_peak_mz
# query_rt <- peak_group@base_peak_rt
# query_adduct <- peak_group@base_peak_adduct
# peak_list_ms2 <- match(peak_list$peak_name, names(raw_msms)) %>%
#   raw_msms[.] %>%
#   .[!is.na(names(.))]
# polarity <- 'negative'
# annotateNeutralLoss(peak_list = peak_list,
#                     peak_list_ms2 = peak_list_ms2,
#                     query_mz = query_mz,
#                     query_rt = query_rt,
#                     query_peak_name = query_peak_name,
#                     query_adduct = query_adduct,
#                     polarity = polarity,
#                     tol_mz = 10,
#                     is_ms2_check = TRUE)

# # example L0177 neg in 200STD
# load('I:/00_projects/03_MetDNA2/00_data/20200710_semi_targeted_annotation_test/03_test/list_peak_group_200716.RData')
# load('I:/00_projects/03_MetDNA2/00_data/20200710_semi_targeted_annotation_test/03_test/raw_msms_200715.RData')
# peak_group <- list_peak_group$`M303T808_[M-H]-`
#
# peak_list <- peak_group@peak_list %>% initializePeakList()
# query_peak_name <- peak_group@base_peak_name
# query_mz <- peak_group@base_peak_mz
# query_rt <- peak_group@base_peak_rt
# query_adduct <- peak_group@base_peak_adduct
# peak_list_ms2 <- match(peak_list$peak_name, names(raw_msms)) %>%
#   raw_msms[.] %>%
#   .[!is.na(names(.))]
#
# polarity <- 'negative'
# annotateNeutralLoss(peak_list = peak_list,
#                     peak_list_ms2 = peak_list_ms2,
#                     query_mz = query_mz,
#                     query_rt = query_rt,
#                     query_peak_name = query_peak_name,
#                     query_adduct = query_adduct,
#                     polarity = polarity,
#                     tol_mz = 10,
#                     is_ms2_check = TRUE)
#
# peak_group <- list_peak_group_annotation$`M298T163_[M+H]+`
# peak_list <- peak_group@peak_list %>% initializePeakList()
# query_peak_name <- peak_group@base_peak_name
# query_mz <- peak_group@base_peak_mz
# query_rt <- peak_group@base_peak_rt
# query_adduct <- peak_group@base_peak_adduct
# peak_list_ms2 <- match(peak_list$peak_name, names(raw_msms)) %>%
#   raw_msms[.] %>%
#   .[!is.na(names(.))]
#
# polarity <- 'positive'
# annotateNeutralLoss(peak_list = peak_list,
#                     peak_list_ms2 = peak_list_ms2,
#                     query_mz = query_mz,
#                     query_rt = query_rt,
#                     query_peak_name = query_peak_name,
#                     query_adduct = query_adduct,
#                     polarity = polarity,
#                     tol_mz = 25,
#                     is_ms2_check = TRUE,
#                     is_rule_limitation = TRUE)



setGeneric(name = 'annotateNeutralLoss',
           def = function(
             peak_list,
             peak_list_ms2,
             query_peak_name,
             query_mz,
             query_rt,
             query_ccs,
             query_adduct,
             polarity = c('positive', 'negative'),
             tol_mz = 25,
             is_ms2_check = FALSE,
             ms2_score_cutoff = -1,
             cutoff_ssc = NULL,
             cutoff_ssc_int = 3000,
             is_rule_limitation = TRUE,
             # lib_adduct_nl = lib_adduct_nl,
             ...
           ){
             nl_list <- convertMz2Adduct(base_mz = query_mz,
                                         base_adduct = query_adduct,
                                         type = 'nl',
                                         polarity = polarity,
                                         lib_adduct_nl = lib_adduct_nl)

             mz_nl_list <- nl_list$mz
             name_nl_list <- nl_list$adduct

             nl_matrix <- getMzRange(mz = mz_nl_list,
                                     ppm = tol_mz,
                                     mz_ppm_thr = 400)

             # keep peaks with ssc larger than cutoff only (except intensity less than cutoff)
             if (length(cutoff_ssc) > 0) {
               peak_list <- peak_list %>%
                 dplyr::filter((ssc >= cutoff_ssc) | (int <= cutoff_ssc_int))
             }

             result_nl <- lapply(seq_along(mz_nl_list), function(i){
               mz_min <- nl_matrix[i,1]
               mz_max <- nl_matrix[i,2]

               result <- peak_list %>%
                 dplyr::filter(mz >= mz_min & mz<= mz_max) %>%
                 dplyr::mutate(mz_error = as.numeric(mz - mz_nl_list[i]),
                               rt_error = as.numeric(rt - query_rt),
                               ccs_error = as.numeric(ccs - query_ccs),
                               type = 'neutralLossAnnotation',
                               annotation = name_nl_list[i],
                               query_peak = query_peak_name)
             })

             result_nl <- result_nl %>% dplyr::bind_rows()

             # If multiple nl has been annotated,
             #    reserve peak with minimum mz error
             result_nl <- result_nl %>%
               dplyr::arrange(annotation, abs(mz_error), abs(rt_error), abs(ccs_error)) %>%
               dplyr::distinct(annotation, .keep_all = TRUE)

             if (is_ms2_check) {
               if (missing(peak_list_ms2)) {
                 stop('Please input msms check')
               }

               if (nrow(result_nl) == 0) {
                 return(result_nl)
               }

               # if the annotation don't has ms/ms, return result
               temp_idx <- which(names(peak_list_ms2) == query_peak_name)

               if (length(temp_idx) == 0) {
                 return(result_nl)
               } else {
                 base_ms2 <- peak_list_ms2[[temp_idx]]
               }


               result_ms2_score <- lapply(seq_along(result_nl$peak_name),
                                          function(i){
                                            # cat(i, ' ')
                                            nl_mz <- result_nl$mz[i]

                                            # purify ms2 spectra and remove fragments >= NL precursor
                                            base_ms2$spec <- purifyMs2(spec = base_ms2$spec,
                                                                       is_include_precursor = TRUE,
                                                                       mz_range_ms2 = c(0, nl_mz+0.1),
                                                                       is_deisotope = FALSE,
                                                                       int_ms2_min_abs = 30,
                                                                       int_ms2_min_relative = 0.01,
                                                                       mz_precursor = nl_mz,
                                                                       ppm_precursor_filter = 10)

                                            # if no fragment mz less than nl mz, return NA
                                            if (length(base_ms2$spec) == 0) {
                                              return(NA)
                                            }

                                            temp_nl_name <- result_nl$peak_name[i]
                                            temp_nl_ms2 <- try(which(names(peak_list_ms2) == temp_nl_name) %>%
                                                                 peak_list_ms2[[.]], silent = TRUE)

                                            # if no nl ms2, return NA
                                            if (class(temp_nl_ms2) == 'try-error') {
                                              return(NA)
                                            }

                                            # purify ms2 spectra and remove fragments >= NL precursor
                                            temp_nl_ms2$spec <- purifyMs2(spec = temp_nl_ms2$spec,
                                                                          is_include_precursor = TRUE,
                                                                          mz_range_ms2 = c(0, nl_mz+0.1),
                                                                          is_deisotope = FALSE,
                                                                          int_ms2_min_abs = 30,
                                                                          int_ms2_min_relative = 0.01,
                                                                          mz_precursor = nl_mz,
                                                                          ppm_precursor_filter = 10)

                                            # if no fragment mz less than nl mz, return NA
                                            if (length(temp_nl_ms2$spec) == 0) {
                                              return(NA)
                                            }

                                            # modify the msms format for SpectraTools
                                            base_ms2 <- convertSpectraData(ms2_data = base_ms2)
                                            temp_nl_ms2 <- convertSpectraData(ms2_data = temp_nl_ms2)


                                            matchParam <- SpectraTools::MatchParam(ppm = 10,
                                                                                   methodScore = 'dp',
                                                                                   methodMatch = 'direct',
                                                                                   weightIntensity = 1,
                                                                                   weightMZ = 0,
                                                                                   cutoff = 0,
                                                                                   includePrecursor = TRUE,
                                                                                   intensityExpNormed = TRUE,
                                                                                   intensityLibNormed = TRUE,
                                                                                   tuneLibSpectra = FALSE)

                                            result <- try(SpectraTools::MatchSpectra(base_ms2,
                                                                                     temp_nl_ms2,
                                                                                     matchParam),
                                                          silent = TRUE)

                                            if ((class(result) == 'try-error') | (length(result) == 0)) {
                                              return(NA)
                                            }

                                            result <- result@info %>% tibble::as_tibble() %>% dplyr::pull(scoreReverse)

                                          })

               result_nl <- result_nl %>%
                 dplyr::mutate(ms2_score = unlist(result_ms2_score)) %>%
                 dplyr::filter((ms2_score >= ms2_score_cutoff) | is.na(ms2_score))

               # Check rule requirements for adducts
               if (is_rule_limitation) {
                 rule_list <- generateRuleList(polarity = polarity)

                 idx <- which(result_nl$annotation %in% rule_list$rule_adduct)
                 if (length(idx) > 0) {
                   idx_eff <- sapply(seq_along(idx), function(i){
                     temp_rule <- rule_list %>%
                       dplyr::filter(rule_adduct == result_nl$annotation[idx[i]]) %>%
                       dplyr::pull(rule)

                     result <- all(temp_rule %in% result_nl$annotation)
                   })

                   if (!all(idx_eff)) {
                     nl_rm <- idx[!idx_eff] %>%
                       result_nl$annotation[.]

                     result_nl <- result_nl %>%
                       dplyr::filter(annotation != nl_rm)
                   }
                 }

               }


             }

             return(result_nl)

           })


#' @title convertSpectraData
#' @param ms2_data
#' @importClassesFrom SpectraTools 'SpectraData'
#' @export

setGeneric(name = 'convertSpectraData',
           def = function(
             ms2_data
           ){
             options(readr.num_columns = 0)
             temp_info <- ms2_data$info %>%
               dplyr::rename(name = NAME,
                             mz = PRECURSORMZ) %>%
               dplyr::select(name:mz) %>%
               readr::type_convert()

             temp_ms2_data <- ms2_data$spec

             result <- new('SpectraData',
                           info = temp_info,
                           spectra = list(temp_ms2_data))

             return(result)
           })







# annotateISF ------------------------------------------------------------------

# # example L0076 neg in 200STD
# # load('./lib_adduct_nl_200714.RData')
# load('I:/00_projects/03_MetDNA2/00_data/20200710_semi_targeted_annotation_test/03_test/list_peak_group_200716.RData')
# load('I:/00_projects/03_MetDNA2/00_data/20200710_semi_targeted_annotation_test/03_test/raw_msms_200715.RData')
# peak_group <- list_peak_group$`M482T929_[M-H]-`
# peak_list <- peak_group@peak_list %>% initializePeakList()
# query_peak_name <- peak_group@base_peak_name
# query_mz <- peak_group@base_peak_mz
# query_rt <- peak_group@base_peak_rt
# query_adduct <- peak_group@base_peak_adduct
# peak_list_ms2 <- match(peak_list$peak_name, names(raw_msms)) %>%
#   raw_msms[.] %>%
#   .[!is.na(names(.))]
# polarity <- 'negative'
# annotateISF(peak_list = peak_list,
#             peak_list_ms2 = peak_list_ms2,
#             query_peak_name = query_peak_name,
#             query_mz = query_mz,
#             query_rt = query_rt,
#             tol_mz = 10)

# peak_group <- list_peak_group_annotation$`M112T341_2_[M+H]+`
# peak_list <- peak_group@peak_list %>% initializePeakList()
# query_peak_name <- peak_group@base_peak_name
# query_mz <- peak_group@base_peak_mz
# query_rt <- peak_group@base_peak_rt
# query_adduct <- peak_group@base_peak_adduct
# peak_list_ms2 <- match(peak_list$peak_name, names(raw_msms)) %>%
#   raw_msms[.] %>%
#   .[!is.na(names(.))]
# polarity <- 'positive'
# annotateISF(peak_list = peak_list,
#             peak_list_ms2 = peak_list_ms2,
#             query_peak_name = query_peak_name,
#             query_mz = query_mz,
#             query_rt = query_rt,
#             tol_mz = 10)

setGeneric(name = 'annotateISF',
           function(
             peak_list,
             peak_list_ms2,
             query_peak_name,
             query_mz,
             query_rt,
             query_ccs,
             tol_mz = 10,
             cutoff_ssc = NULL,
             cutoff_ssc_int = 3000,
             cutoff_topN = 5,
             ...
           ){

             # if the annotation don't has ms/ms, return result
             temp_idx <- which(names(peak_list_ms2) == query_peak_name)

             if (length(temp_idx) == 0) {
               return(NULL)
             } else {
               base_ms2 <- peak_list_ms2[[temp_idx]]
             }


             # base_ms2 <- which(names(peak_list_ms2) == query_peak_name) %>% peak_list_ms2[[.]]

             base_ms2 <- purifyMs2(spec = base_ms2$spec,
                                   is_include_precursor = FALSE,
                                   is_deisotope = FALSE,
                                   mz_range_ms2 = c(0, query_mz),
                                   int_ms2_min_abs = 30,
                                   int_ms2_min_relative = 0.01,
                                   mz_precursor = query_mz,
                                   ppm_precursor_filter = 10)

             if (length(base_ms2) == 0) {
               return(NULL)
             }


             # limit topN fragments in base peak for ISF annotation
             # mz_isf_list <- base_ms2[,'mz']
             mz_isf_list <- order(base_ms2[,'intensity'], decreasing = TRUE) %>%
               base_ms2[.,,drop = FALSE] %>%
               .[,'mz']

             if (length(mz_isf_list) >= cutoff_topN) {
               mz_isf_list <- mz_isf_list[1:cutoff_topN]
             }

             isf_matrix <- getMzRange(mz = mz_isf_list,
                                      ppm = 25,
                                      mz_ppm_thr = 400)

             # keep peaks with ssc larger than cutoff only (except intensity less than cutoff)
             if (length(cutoff_ssc) > 0) {
               peak_list <- peak_list %>%
                 dplyr::filter((ssc >= cutoff_ssc) | (int <= cutoff_ssc_int))
             }

             result_isf <- lapply(seq_along(mz_isf_list), function(i){
               mz_min <- isf_matrix[i,1]
               mz_max <- isf_matrix[i,2]

               result <- peak_list %>%
                 dplyr::filter(mz >= mz_min & mz<= mz_max) %>%
                 dplyr::mutate(mz_error = as.numeric(mz - mz_isf_list[i]),
                               rt_error = as.numeric(rt - query_rt),
                               ccs_error = as.numeric(ccs - query_ccs),
                               type = 'isfAnnotation',
                               annotation = 'ISF',
                               query_peak = query_peak_name)
             })

             result_isf <- result_isf %>% dplyr::bind_rows()

             # If multiple isf has been annotated,
             #    reserve peak with minimum mz error
             result_isf <- result_isf %>%
               dplyr::arrange(annotation, abs(mz_error), abs(rt_error), abs(ccs_error)) %>%
               dplyr::distinct(peak_name, .keep_all = TRUE)

             return(result_isf)
           })





# mergePeakListAnnotation ------------------------------------------------------
# mergePeakListAnnotation(pl_annotation = peak_list_annotation) %>% View()
setGeneric(name = 'mergePeakListAnnotation',
           def = function(
             pl_annotation
           ){

             # when these pair existed: M+_[M+H]+ & M-_[M-H]- & their isotope annotation overlaped
             #   rm isotope annotations from M+ abd M-

             if (all(c('M+', '[M+H]+') %in% pl_annotation$annotation) |
                 all(c('M-', '[M-H]-') %in% pl_annotation$annotation)) {

               temp_peak1 <- which(pl_annotation$annotation == '[M+H]+' | pl_annotation$annotation == '[M-H]-') %>% pl_annotation$peak_name[.]
               temp_peak2 <- which(pl_annotation$annotation == 'M+' | pl_annotation$annotation == 'M-') %>% pl_annotation$peak_name[.]

               temp_peak1_list <- pl_annotation %>% dplyr::filter(query_peak == temp_peak1 & type == 'isotopeAnnotation')
               temp_peak2_list <- pl_annotation %>% dplyr::filter(query_peak == temp_peak2 & type == 'isotopeAnnotation')

               is_overlap <- sum(temp_peak2_list$peak_name %in% temp_peak1_list$peak_name)

               if (is_overlap > 0) {
                 idx_rm <- which((pl_annotation$query_peak == temp_peak2) & (pl_annotation$annotation %in% c('[M+1]', '[M+2]', '[M+3]')))
                 pl_annotation <- pl_annotation[-idx_rm,]
               }

             }

             isotope_list <- pl_annotation %>% dplyr::filter(type == 'isotopeAnnotation')

             # if one peak were annotated as [M] & [M+1]/[M+2], remove reduancy
             #    [M+2] > [M+1] > [M]

             isotope_list <- isotope_list %>%
               dplyr::arrange(peak_name, desc(annotation)) %>%
               dplyr::distinct(peak_name, .keep_all = TRUE)


             other_list <- pl_annotation %>% dplyr::filter(type != 'isotopeAnnotation')

             # if the peak was annotated as [M+1], [M+2], [M+3],
             #    then remove other annotations of this peak
             list_confilict_annotation <- isotope_list %>%
               dplyr::filter(annotation != '[M]') %>%
               dplyr::pull(peak_name)

             # Other conflict was deplicated according to m/z error
             other_list <- other_list %>%
               dplyr::filter(!(peak_name %in% list_confilict_annotation)) %>%
               dplyr::arrange(abs(mz_error)) %>%
               dplyr::distinct(peak_name, .keep_all = TRUE)

             # convert to wider table
             result_list <- isotope_list %>%
               tidyr::pivot_wider(names_from = type, values_from = annotation)

             other_list <- other_list %>%
               tidyr::pivot_wider(names_from = type, values_from = annotation)

             temp_list <- other_list %>%
               dplyr::select(peak_name, dplyr::contains('Annotation'))

             result_list <- result_list %>%
               dplyr::left_join(temp_list, by = 'peak_name') %>%
               dplyr::arrange(mz)


             # replace mz_error, rt_error int_ratio, ms2_score & query_peak with correct annotation
             idx <- match(temp_list$peak_name, result_list$peak_name)
             result_list$mz_error[idx] <- other_list$mz_error
             result_list$rt_error[idx] <- other_list$rt_error
             result_list$ccs_error[idx] <- other_list$ccs_error
             result_list$int_ratio[idx] <- other_list$int_ratio
             result_list$ms2_score[idx] <- other_list$ms2_score
             result_list$query_peak[idx] <- other_list$query_peak


             # impute format
             idx <- c('isotopeAnnotation', 'adductAnnotation', 'neutralLossAnnotation', 'isfAnnotation') %>%
               match(., colnames(result_list)) %>%
               is.na() %>%
               which()

             if (length(idx)>0) {
               imputate_list <- matrix(nrow = nrow(result_list), ncol = length(idx))
               colnames(imputate_list) <- c('isotopeAnnotation', 'adductAnnotation', 'neutralLossAnnotation', 'isfAnnotation')[idx]

               options(readr.num_columns = 0)
               imputate_list <- imputate_list %>% tibble::as_tibble()

               result_list <- result_list %>% dplyr::bind_cols(imputate_list)
             }

             result_list <- result_list %>%
               dplyr::select(peak_name:query_peak, isotopeAnnotation, adductAnnotation,
                             neutralLossAnnotation, isfAnnotation)

             return(result_list)
           })



# annotatePeakGroup ------------------------------------------------------------

#' @title annotatePeakGroup
#' @author Zhiwei Zhou
#' @description annotate peak group with isotopes, adducts, neutral loss and in-source fragments
#' @param peak_group a object of peak_group
#' @param polarity
#' @param tol_mz Default: 10 ppm
#' @param is_ms2_check whether compare ms2 similarity of neutral loss to base peak ms2. Default: TRUE
#' @param ms2_score_cutoff Default: -1; # -1 represent not filter
#' @param cutoff_ssc Whether apply ssc filter in annotation. Default: 0.3
#' @param cutoff_ssc_int If the peak intensity less than this value, cutoff_ssc_int is invalid.
#' @param is_rule_limitation Whether apply rule limitation in adductAnnotation and nlAnnotation (e.g. '[2M+H]+' needs '[M+H]+'). Default: TRUE
#' @param cutoff_topN Top N fragments in base peak ms2 used for ISF annotation. Default: 5
#' @export
#' @example
#' load(system.file("tempdata", "list_peak_group_200805.RData", package="MetDNA2"))
#' load(system.file("tempdata", "raw_msms_200805.RData", package="MetDNA2"))
#' peak_group_L0076 <- list_peak_group$`M482T929_[M-H]-`
#' test_L0076 <- annotatePeakGroup(peak_group = peak_group_L0076,
#'                                 ms2_data = raw_msms,
#'                                 polarity = 'negative',
#'                                 tol_mz = 10,
#'                                 is_ms2_check = TRUE,
#'                                 ms2_score_cutoff = -1)


# load('./inst/tempdata/list_peak_group_200805.RData')
# load('./inst/tempdata/raw_msms_200805.RData')
# peak_group_L0076 <- list_peak_group$`M482T929_[M-H]-`
#
# test_L0076 <- annotatePeakGroup(peak_group = peak_group_L0076,
#                                 ms2_data = raw_msms,
#                                 polarity = 'negative',
#                                 tol_mz = 10,
#                                 is_ms2_check = TRUE,
#                                 ms2_score_cutoff = -1)


# peak_group_L0177 <- list_peak_group$`M303T808_[M-H]-`
#
# test_L0177 <- annotatePeakGroup(peak_group = peak_group_L0177,
#                                 ms2_data = raw_msms,
#                                 polarity = 'negative',
#                                 tol_mz = 10,
#                                 is_ms2_check = TRUE,
#                                 ms2_score_cutoff = -1)


# peak_group_L0306 <- list_peak_group$`M175T635_[M-H]-`
#
# test_L0306 <- annotatePeakGroup(peak_group = peak_group_L0306,
#                                 ms2_data = raw_msms,
#                                 polarity = 'negative',
#                                 tol_mz = 10,
#                                 is_ms2_check = TRUE,
#                                 ms2_score_cutoff = -1)


# peak_group_L0319 <- list_peak_group$`M218T716_[M-H]-`
#
# test_L0319 <- annotatePeakGroup(peak_group = peak_group_L0319,
#                                 ms2_data = raw_msms,
#                                 polarity = 'negative',
#                                 tol_mz = 10,
#                                 is_ms2_check = TRUE,
#                                 ms2_score_cutoff = -1)


# peak_group_S0079 <- list_peak_group$`M188T302_[M-H]-`
#
# test_S0079 <- annotatePeakGroup(peak_group = peak_group_S0079,
#                                 ms2_data = raw_msms,
#                                 polarity = 'negative',
#                                 tol_mz = 10,
#                                 is_ms2_check = TRUE,
#                                 ms2_score_cutoff = -1)


# peak_group_S0080 <- list_peak_group$`M168T185_[M-H]-`
#
# test_S0080 <- annotatePeakGroup(peak_group = peak_group_S0080,
#                                 ms2_data = raw_msms,
#                                 polarity = 'negative',
#                                 tol_mz = 10,
#                                 is_ms2_check = TRUE,
#                                 ms2_score_cutoff = -1)


# peak_group <- list_peak_group$`M428T300_[M+H]+`
#
# test_S0080 <- annotatePeakGroup(peak_group = peak_group_S0080,
#                                 ms2_data = raw_msms,
#                                 polarity = 'negative',
#                                 tol_mz = 10,
#                                 is_ms2_check = TRUE,
#                                 ms2_score_cutoff = -1)


# peak_group <- list_peak_group$`M145T551_[2M+H]+`
# peak_group <- list_peak_group$`M428T300_[M+H]+`
# peak_group <- list_peak_group$`M184T910_[M]+`
# ms2_data <- raw_msms
# polarity <- 'positive'
# tol_mz <- 25
# is_ms2_check <- TRUE
# ms2_score_cutoff <- -1
# cutoff_ssc <- 0.3
# cutoff_ssc_int <- 3000
# is_rule_limitation <- TRUE
# cutoff_topN <- 5

setGeneric(name = 'annotatePeakGroup',
           def = function(
             peak_group,
             ms2_data,
             polarity = c('positive', 'negative'),
             tol_mz = tol_mz,
             is_ms2_check = TRUE, # compare ms2 similarity of neutral loss
             ms2_score_cutoff = -1, # -1 represent not filter
             lib_adduct_nl = lib_adduct_nl,
             cutoff_ssc = 0.3,
             cutoff_ssc_int = 3000,
             is_rule_limitation = TRUE,
             cutoff_topN = 5,
             ...
           ){
             # match.arg(polarity)

             if (missing(peak_group)) {
               stop('Please input peak group!\n')
             }

             if (missing(ms2_data)) {
               stop('Please input ms2 data!\n')
             }

             # browser()

             peak_list <- peak_group@peak_list %>% initializePeakList()
             query_mz <- peak_group@base_peak_mz
             query_rt <- peak_group@base_peak_rt
             query_ccs <- peak_group@base_peak_ccs
             query_adduct <- peak_group@base_peak_adduct
             query_peak_name <- peak_group@base_peak_name
             peak_list_ms2 <- match(peak_list$peak_name, names(ms2_data)) %>%
               ms2_data[.] %>%
               .[!is.na(names(.))]


             # result_isotope <- annotateIsotope(peak_list = peak_list,
             #                                   query_mz = query_mz,
             #                                   query_rt = query_rt,
             #                                   query_ccs = query_ccs,
             #                                   query_peak_name = query_peak_name,
             #                                   tol_mz = tol_mz,
             #                                   cutoff_ssc = NULL)
             # cat('iso\t')
             result_isotope <- annotateIsotope(peak_list = peak_list,
                                               query_mz = query_mz,
                                               query_rt = query_rt,
                                               query_ccs = query_ccs,
                                               query_peak_name = query_peak_name,
                                               tol_mz = tol_mz,
                                               cutoff_ssc = NULL,
                                               ...)

             # cat('addu\t')
             result_adduct <- annotateAdduct(peak_list = peak_list,
                                             query_mz = query_mz,
                                             query_rt = query_rt,
                                             query_ccs = query_ccs,
                                             query_peak_name = query_peak_name,
                                             query_adduct = query_adduct,
                                             polarity = polarity,
                                             tol_mz = tol_mz,
                                             cutoff_ssc = cutoff_ssc,
                                             cutoff_ssc_int = cutoff_ssc_int,
                                             is_rule_limitation = is_rule_limitation)
             # cat('nl\t')
             result_nl <- annotateNeutralLoss(peak_list = peak_list,
                                              peak_list_ms2 = peak_list_ms2,
                                              query_mz = query_mz,
                                              query_rt = query_rt,
                                              query_ccs = query_ccs,
                                              query_peak_name = query_peak_name,
                                              query_adduct = query_adduct,
                                              polarity = polarity,
                                              tol_mz = tol_mz,
                                              is_ms2_check = is_ms2_check,
                                              ms2_score_cutoff = ms2_score_cutoff,
                                              cutoff_ssc = cutoff_ssc,
                                              cutoff_ssc_int = cutoff_ssc_int,
                                              is_rule_limitation = is_rule_limitation)


             # cat('isf\t')
             result_isf <- annotateISF(peak_list = peak_list,
                                       peak_list_ms2 = peak_list_ms2,
                                       query_peak_name = query_peak_name,
                                       query_mz = query_mz,
                                       query_rt = query_rt,
                                       query_ccs = query_ccs,
                                       tol_mz = tol_mz,
                                       cutoff_ssc = NULL,
                                       cutoff_topN = cutoff_topN)

             peak_list_annotation <- writePeakListAnnotation(pl_annotation_new = result_isotope) %>%
               writePeakListAnnotation(pl_annotation_new = result_adduct) %>%
               writePeakListAnnotation(pl_annotation_new = result_nl) %>%
               writePeakListAnnotation(pl_annotation_new = result_isf)

             # annotate isotopes for adducts, neutral loss, and in-source fragments
             # cat('PURR\t')
             temp_new <- extractIsotopeNew(peak_list_annotation)

             if (length(temp_new) > 0) {
               # result_isotope_new <- purrr::map(seq_along(temp_new),
               result_isotope_new <- lapply(seq_along(temp_new),
                                            function(i){
                                              temp <- temp_new[i]
                                              idx <- which(peak_list$peak_name == temp)
                                              # temp_data <- peak_list %>% dplyr::filter(peak_name == temp_new[i])

                                              temp_data <- peak_list[idx,]

                                              annotateIsotope(peak_list = peak_list,
                                                              query_mz = temp_data$mz,
                                                              query_rt = temp_data$rt,
                                                              query_ccs = temp_data$ccs,
                                                              query_peak_name = temp_new[i],
                                                              tol_mz = tol_mz,
                                                              cutoff_ssc = NULL)

                                            }) %>% dplyr::bind_rows()

               # cat('PURR2\t')
               peak_list_annotation <- writePeakListAnnotation(pl_annotation_old = peak_list_annotation,
                                                               pl_annotation_new = result_isotope_new)

             }

             # if the base_peak don't meet the rule,
             #    remove all annotations of this peak group
             is_rule_valid <- query_adduct %in% peak_list_annotation$annotation

             if (!is_rule_valid){
               peak_list_merged <- tibble::tibble(peak_name = character(),
                                                  mz = numeric(),
                                                  rt = numeric(),
                                                  ccs = numeric(),
                                                  ssc = numeric(),
                                                  int = numeric(),
                                                  mz_error = numeric(),
                                                  rt_error = numeric(),
                                                  ccs_error = numeric(),
                                                  int_ratio = numeric(),
                                                  ms2_score = numeric(),
                                                  query_peak = character(),
                                                  isotopeAnnotation = character(),
                                                  adductAnnotation = character(),
                                                  neutralLossAnnotation = character(),
                                                  isfAnnotation = character())

               peak_group@peak_list_annotated <- peak_list_merged

               return(peak_group)

             }

             # cat('merge\t')
             # remove redunancy of annotation
             peak_list_merged <- mergePeakListAnnotation(pl_annotation = peak_list_annotation)

             peak_group@peak_list_annotated <- peak_list_merged
             # cat('\n')
             return(peak_group)
           })


################################################################################
# concise peak groups ----------------------------------------------------------
# concisePeak2PeakGroup --------------------------------------------------------
#' @title concisePeak2PeakGroup
#' @author Zhiwei Zhou
#' @description concise conflict peak groups (one peak - multiple peak groups(i.e., adduct), M167T834_[M-H]- & M167T834_[M-H2O-H]-)
#' @param list_peak_group_annotation
#' @export
#' @examples
#' load(system.file("tempdata", "list_peak_group_annotation_200805.RData", package="MetDNA2"))
#' test <- concisePeak2PeakGroup(list_peak_group_annotation = list_peak_group_annotation)

# load('./inst/tempdata/list_peak_group_annotation_200805.RData')
# test <- concisePeak2PeakGroup(list_peak_group_annotation = list_peak_group_annotation)

setGeneric(name = 'concisePeak2PeakGroup',
           function(
             list_peak_group_annotation,
             path_dir,
             type_order = c('level1', 'level2', 'level3'),
             ...
           ){
             temp_stat <- purrr::map(seq_along(list_peak_group_annotation), function(i){
               tibble::tibble(name = list_peak_group_annotation[[i]]@base_peak_name,
                              mz = list_peak_group_annotation[[i]]@base_peak_mz,
                              rt = list_peak_group_annotation[[i]]@base_peak_rt,
                              ccs = list_peak_group_annotation[[i]]@base_peak_ccs,
                              adduct = list_peak_group_annotation[[i]]@base_peak_adduct,
                              type = list_peak_group_annotation[[i]]@base_peak_seed,
                              num_peak_group = nrow(list_peak_group_annotation[[i]]@peak_list),
                              num_peak_group_anno = nrow(list_peak_group_annotation[[i]]@peak_list_annotated))
             }) %>% dplyr::bind_rows()

             peak_conflict <- temp_stat %>%
               dplyr::count(name) %>%
               dplyr::filter(n>1) %>%
               dplyr::pull(name)

             peak_group_conflict_excluded <- purrr::map(seq_along(peak_conflict),
                                                        function(i){
                                                          temp_data <- temp_stat %>%
                                                            dplyr::filter(name == peak_conflict[i]) %>%
                                                            dplyr::arrange(match(type, type_order),
                                                                           desc(num_peak_group_anno)) %>%
                                                            dplyr::mutate(label = paste(name, adduct, sep = '_'))

                                                          temp_data_reserve <- temp_data %>%
                                                            dplyr::filter(type %in% type[1]) %>%
                                                            dplyr::filter(num_peak_group_anno == max(num_peak_group_anno)) %>%
                                                            dplyr::pull(label)

                                                          result <- temp_data %>%
                                                            dplyr::filter(!(label %in% temp_data_reserve)) %>%
                                                            dplyr::pull(label)

                                                          return(result)
                                                        }) %>% unlist()

             # if no conflict peak is removed, return result_list_peak_group
             if (length(peak_group_conflict_excluded) != 0) {
               idx <- match(peak_group_conflict_excluded, names(list_peak_group_annotation))
               result_list_peak_group <- list_peak_group_annotation[-idx]

               record_conflict_peak_group <- peak_group_conflict_excluded
               save(record_conflict_peak_group,
                    file = file.path(path_dir,
                                     '03_annotation_credential',
                                     "00_intermediate_data",
                                     'record_conflict_peak_group.RData'))

               return(result_list_peak_group)

             } else {
               result_list_peak_group <- list_peak_group_annotation
               record_conflict_peak_group <- peak_group_conflict_excluded
               save(record_conflict_peak_group,
                    file = file.path(path_dir,
                                     '03_annotation_credential',
                                     "00_intermediate_data",
                                     'record_conflict_peak_group.RData'))

               return(result_list_peak_group)

             }


           })


# concisePeakGroup2PeakGroup ---------------------------------------------------

#' @title concisePeakGroup2PeakGroup
#' @author Zhiwei Zhou
#' @description concise overlap peak groups (Peak group - peak group, e.g. M175T752_[M-H]-, ISF of M175T752 (M132T752_[M+NH4-2H]-, M132T752_[M-H]-), Isotope of M175T752 (M177T752_[M-H]-))
#' @param list_peak_group_annotation
#' @export
#' @examples
#' load(system.file("tempdata", "list_peak_group_annotation_concised_p2pg_200805.RData", package="MetDNA2"))
#' test <- concisePeak2PeakGroup(list_peak_group_annotation = list_peak_group_annotation_concised)

# load('./inst/tempdata/list_peak_group_annotation_concised_p2pg_200805.RData')
# test <- concisePeakGroup2PeakGroup(list_peak_group_annotation = list_peak_group_annotation_concised)

setGeneric(name = 'concisePeakGroup2PeakGroup',
           def = function(
             list_peak_group_annotation,
             # type_order = c('seed', 'none_seed'),
             type_order = c('level1', 'level2', 'level3'),
             path_dir,
             ...
           ){

             stat_peak_groups <- purrr::map(seq_along(list_peak_group_annotation), function(i){
               tibble::tibble(name = list_peak_group_annotation[[i]]@base_peak_name,
                              mz = list_peak_group_annotation[[i]]@base_peak_mz,
                              rt = list_peak_group_annotation[[i]]@base_peak_rt,
                              ccs = list_peak_group_annotation[[i]]@base_peak_ccs,
                              adduct = list_peak_group_annotation[[i]]@base_peak_adduct,
                              seed = list_peak_group_annotation[[i]]@base_peak_seed,
                              num_peak_group = nrow(list_peak_group_annotation[[i]]@peak_list),
                              num_peak_group_anno = nrow(list_peak_group_annotation[[i]]@peak_list_annotated))

             }) %>%
               dplyr::bind_rows() %>%
               dplyr::arrange(match(seed, type_order), # credential order
                              dplyr::desc(num_peak_group_anno))


             list_unconcised_peak_group <- paste(stat_peak_groups$name,
                                                 stat_peak_groups$adduct, sep = '_')
             list_concised_peak_group <- vector('list',
                                                length = length(list_unconcised_peak_group))
             record_overlap_peak_group <- vector('list',
                                                 length = length(list_unconcised_peak_group))
             i <- 1
             while (length(list_unconcised_peak_group) > 0) {
               # cat(i, ' ')

               list_unconcised_peak_name <- list_unconcised_peak_group %>%
                 sapply(., function(x){
                   stringr::str_split(x, pattern = '_\\[|_\\M')[[1]][1]
                 }) %>% unname()

               # select largest peak group
               base_peak_group_name <- list_unconcised_peak_group[1]
               base_peak_name <- (base_peak_group_name %>% stringr::str_split(pattern = '_\\[|_\\M'))[[1]][1]
               base_peak_group_peak_list <- list_peak_group_annotation[[base_peak_group_name]]@peak_list_annotated
               base_peak_adduct <- list_peak_group_annotation[[base_peak_group_name]]@base_peak_adduct

               # if dimer dectected and moner existed,
               #   skip the base_peak_adduct
               if (stringr::str_detect(base_peak_adduct, '2M|3M')) {
                 temp <- stat_peak_groups %>%
                   dplyr::filter(paste(name, adduct, sep = '_') %in% list_unconcised_peak_group) %>%
                   # dplyr::filter(name %in% list_unconcised_peak_name) %>%
                   dplyr::filter(name %in% base_peak_group_peak_list$peak_name) %>%
                   dplyr::filter(name != base_peak_name)

                 if (nrow(temp) > 0) {
                   if (sum(!stringr::str_detect(temp$adduct, '2M|3M')) > 0) {
                     list_unconcised_peak_group <- c(list_unconcised_peak_group[-1], base_peak_group_name)
                     next()
                   }
                 }

               }

               overlap_peak_group_name <- stat_peak_groups %>%
                 dplyr::filter(paste(name, adduct, sep = '_') %in% list_unconcised_peak_group) %>%
                 # dplyr::filter(name %in% list_unconcised_peak_name) %>%
                 dplyr::filter(name %in% base_peak_group_peak_list$peak_name) %>%
                 dplyr::filter(name != base_peak_name) %>%
                 dplyr::mutate(label = paste(name, adduct, sep = '_')) %>%
                 dplyr::pull(label)

               # if overlap peak group existed,
               #    remove overlap peak groups & add to list_concised_peak_group
               #  else
               #    add to list_concised_peak_group

               if (length(overlap_peak_group_name) > 0) {
                 # remove overlaped peak groups
                 idx <- match(overlap_peak_group_name, names(list_peak_group_annotation))
                 list_peak_group_annotation <- list_peak_group_annotation[-idx]

                 # modify concised and unconcised peak groups
                 list_concised_peak_group[[i]] <- c(overlap_peak_group_name, base_peak_group_name)
                 record_overlap_peak_group[[i]] <- paste(overlap_peak_group_name, base_peak_group_name, sep = '@')
               } else {
                 list_concised_peak_group[[i]] <- c(base_peak_group_name)
                 record_overlap_peak_group[[i]] <- NULL
               }

               idx <- list_concised_peak_group[[i]] %>%
                 match(., list_unconcised_peak_group)
               list_unconcised_peak_group <- list_unconcised_peak_group[-idx]

               i <- i + 1
             }

             save(record_overlap_peak_group,
                  file = file.path(path_dir,
                                   '03_annotation_credential',
                                   "00_intermediate_data",
                                   'record_overlap_peak_group.RData'))

             return(list_peak_group_annotation)
           })




# concisePeakGroup2Annotation --------------------------------------------------
#' @title concisePeakGroup2Annotation
#' @author Zhiwei Zhou
#' @description convert peak groups to annotation table
#' @param list_peak_group_annotation
#' @export
#' @examples
#' load(system.file("tempdata", "list_peak_group_annotation_concised_pg2pg_200805.RData", package="MetDNA2"))
#' test <- concisePeakGroup2Annotation(list_peak_group_annotation_concised)

# load('./inst/tempdata/list_peak_group_annotation_concised_pg2pg_200805.RData')
# test <- concisePeakGroup2Annotation(list_peak_group_annotation_concised)

setGeneric(name = 'concisePeakGroup2Annotation',
           def = function(
             list_peak_group_annotation
           ){
             # progress <- MetDNA2::mapProgress(n = length(list_peak_group_annotation))
             peak_group_id <- tibble::tibble(id = seq(length(list_peak_group_annotation)),
                                             name = names(list_peak_group_annotation)) %>%
               dplyr::mutate(id = paste0('[', id, ']'))

             peak_list_annotated <- purrr::map(seq_along(list_peak_group_annotation), function(i){
               # MetDNA2::mapProgressPrint(progress)
               temp_peak_group <- list_peak_group_annotation[[i]]

               temp_peak_list_annotated <- temp_peak_group@peak_list_annotated %>%
                 dplyr::mutate(peak_group = paste(temp_peak_group@base_peak_name,
                                                  temp_peak_group@base_peak_adduct,
                                                  sep = '_'),
                               scale_peak_group = nrow(temp_peak_group@peak_list_annotated))

               # temp_idx <- which(temp_peak_list_annotated$isotopeAnnotation %in% c('[M+1]', '[M+2]', '[M+3]', '[M+4]'))
               # temp_peak_list_annotated$query_peak

               return(temp_peak_list_annotated)
             })

             peak_list_annotated <- peak_list_annotated %>% dplyr::bind_rows()
             peak_list_annotated <- peak_list_annotated %>%
               dplyr::arrange(mz, peak_group, desc(scale_peak_group)) %>%
               dplyr::left_join(peak_group_id, by = c('peak_group' = 'name')) %>%
               dplyr::rename(peak_group_id = id)


             peak_list_annotated_long <- peak_list_annotated %>%
               tidyr::pivot_longer(cols = c('adductAnnotation',
                                            'neutralLossAnnotation',
                                            'isfAnnotation'),
                                   names_to = 'adduct_nl_isf',
                                   values_to = 'label') %>%
               dplyr::filter(!is.na(label) | isotopeAnnotation %in% c('[M+1]',
                                                                      '[M+2]',
                                                                      '[M+3]',
                                                                      '[M+4]'))

             peak_list_annotated_long_isotope <- peak_list_annotated_long %>%
               dplyr::filter(isotopeAnnotation %in% c('[M+1]',
                                                      '[M+2]',
                                                      '[M+3]',
                                                      '[M+4]')) %>%
               dplyr::distinct(peak_name,
                               isotopeAnnotation,
                               peak_group,
                               .keep_all = TRUE) %>%
               dplyr::mutate(adduct_nl_isf = 'isotopeAnnotation')

             peak_list_annotated_long_adduct_nl_isf <- peak_list_annotated_long %>%
               dplyr::filter(!is.na(label))


             isotope_adduct_label <- mapply(function(x,y){
               result <- peak_list_annotated_long_adduct_nl_isf %>%
                 dplyr::filter(peak_name == x & peak_group == y) %>%
                 dplyr::pull(label)

               # if the isotopes were annotated by the isotope of ISF, adduct, NL
               #    trace the initial annotation by recursive way
               temp_x <- x
               while(length(result) == 0) {
                 # cat(length(result))
                 temp_x <- peak_list_annotated_long_isotope %>%
                   dplyr::filter(peak_name == temp_x & peak_group == y) %>%
                   dplyr::pull(query_peak)

                 result <- peak_list_annotated_long_adduct_nl_isf %>%
                   dplyr::filter(peak_name == temp_x & peak_group == y) %>%
                   dplyr::pull(label)
               }

               return(result)
             },
             x = peak_list_annotated_long_isotope$query_peak,
             y = peak_list_annotated_long_isotope$peak_group)

             peak_list_annotated_long_isotope$label <- isotope_adduct_label

             result_annotation_table <- peak_list_annotated_long_adduct_nl_isf %>%
               dplyr::bind_rows(peak_list_annotated_long_isotope) %>%
               dplyr::arrange(mz, peak_name, desc(scale_peak_group)) %>%
               dplyr::group_by(peak_name) %>%
               dplyr::summarise(mz = mz[1],
                                rt = rt[1],
                                ssc = paste(ssc, collapse = ';'),
                                int = paste(int, collapse = ';'),
                                mz_error = paste(mz_error, collapse = ';'),
                                rt_error = paste(rt_error, collapse = ';'),
                                ccs_error = paste(ccs_error, collapse = ';'),
                                int_ratio = paste(int_ratio, collapse = ';'),
                                ms2_score = paste(ms2_score, collapse = ';'),
                                isotope = paste(isotopeAnnotation, collapse = ';'),
                                type_annotation = paste(adduct_nl_isf, collapse = ';'),
                                annotation = paste(label, collapse = ';'),
                                peak_group = paste(peak_group, collapse = ';'),
                                peak_group_id = paste(peak_group_id, collapse = ';'),
                                peak_group_scale = paste(scale_peak_group, collapse = ';')) %>%
               dplyr::ungroup() %>%
               dplyr::arrange(mz, peak_name)

             result <- list(peak_group_id = peak_group_id,
                            result_annotation_table = result_annotation_table)

             return(result)

           })





################################################################################
# statistic peak group ---------------------------------------------------------

#' @title statNumListPeakGroup
#' @author Zhiwei Zhou
#' @description statistics number of peak groups
#' @param list_pak_group
#' @export
#' @examples
#' load(system.file("tempdata", "list_peak_group_annotation_200805.RData", package="MetDNA2"))
#' statNumListPeakGroup(list_peak_group = list_peak_group_annotation)


# load('./inst/tempdata/list_peak_group_200805.RData')
# statNumListPeakGroup(list_peak_group = list_peak_group)

# load('./inst/tempdata/list_peak_group_annotation_200805.RData')
# statNumListPeakGroup(list_peak_group = list_peak_group_annotation)


setGeneric(name = 'statNumListPeakGroup',
           def = function(
             list_peak_group
           ){
             num_peak <- sapply(list_peak_group, function(x){
               x@base_peak_name
             }) %>%
               unname() %>%
               unique() %>%
               length()

             num_peak_group <- sapply(list_peak_group, function(x){
               paste0(x@base_peak_name, '_', x@base_peak_adduct)
             }) %>%
               unname() %>%
               unique() %>%
               length()

             num_annotated_peak <- sapply(list_peak_group, function(x){
               nrow(x@peak_list_annotated)
             }) %>%
               unname() %>%
               sum()

             result <- tibble::tibble('No. of peak' = num_peak,
                                      'No. of peak group' = num_peak_group,
                                      'No. of annotated peak' = num_annotated_peak)

             return(result)

           })


#' @title statTypeListPeakGroup
#' @author Zhiwei Zhou
#' @description
#' @param list_peak_group
#' @export
#' @examples

# load('./inst/tempdata/list_peak_group_annotation_200805.RData')
# statTypeListPeakGroup(list_peak_group = list_peak_group_annotation)

setGeneric(name = 'statTypeListPeakGroup',
           function(
             list_peak_group
           ){
             table_annotated <- purrr::map(list_peak_group, function(x){
               x@peak_list_annotated %>%
                 dplyr::mutate(peak_group = paste0(x@base_peak_name,
                                                   '_',
                                                   x@base_peak_adduct))
             }) %>% dplyr::bind_rows()

             if (nrow(table_annotated) == 0) {
               result <- tibble::tibble('No. of unique peak' = 0,
                                        'No. of peak' = 0,
                                        'No. of isotope' = 0,
                                        'No. of adduct' = 0,
                                        'No. of NL' = 0,
                                        'No. of ISF' = 0)

               return(result)
             }

             table_annotated <- table_annotated %>%
               tidyr::pivot_longer(cols = c('adductAnnotation',
                                            'neutralLossAnnotation',
                                            'isfAnnotation'),
                                   names_to = 'adduct_nl_isf',
                                   values_to = 'label') %>%
               dplyr::filter(!is.na(label) | isotopeAnnotation %in% c('[M+1]',
                                                                      '[M+2]',
                                                                      '[M+3]',
                                                                      '[M+4]'))

             temp_isotope <- table_annotated %>%
               dplyr::filter(isotopeAnnotation %in% c('[M+1]',
                                                      '[M+2]',
                                                      '[M+3]',
                                                      '[M+4]')) %>%
               dplyr::distinct(peak_name,
                               isotopeAnnotation,
                               peak_group,
                               .keep_all = TRUE)

             temp_adduct_nl_isf <- table_annotated %>%
               dplyr::filter(!is.na(label))


             num_unique_peak <- table_annotated %>%
               dplyr::distinct(peak_name) %>%
               dplyr::count() %>%
               dplyr::pull(n)

             num_peak <- nrow(temp_isotope) + nrow(temp_adduct_nl_isf)
             num_isotope <- nrow(temp_isotope)
             num_adduct <- temp_adduct_nl_isf %>%
               dplyr::filter(adduct_nl_isf == 'adductAnnotation') %>% nrow()
             num_nl <- temp_adduct_nl_isf %>%
               dplyr::filter(adduct_nl_isf == 'neutralLossAnnotation') %>% nrow()
             num_isf <- temp_adduct_nl_isf %>%
               dplyr::filter(adduct_nl_isf == 'isfAnnotation') %>% nrow()

             result <- tibble::tibble('No. of unique peak' = num_unique_peak,
                                      'No. of peak' = num_peak,
                                      'No. of isotope' = num_isotope,
                                      'No. of adduct' = num_adduct,
                                      'No. of NL' = num_nl,
                                      'No. of ISF' = num_isf)

             return(result)
           })


#' @title convertSemiTargetedResult2LongTable
#' @author Zhiwei Zhou
#' @description convert wide semi targeted clustering result to long table
#' @param result_semi_table
#' @export
#' @examples
#' result_semi_table <- readr::read_csv('I:/00_projects/03_MetDNA2/00_data/20200805_targeted_annotation_evaluation_200STD/01_targeted_annotation/200STD_neg_reverse/annot_credential/semi_targeted_annotation_result.csv')
#' test <- convertSemiTargetedResult2LongTable(result_semi_table)

# result_semi_table <- readr::read_csv('I:/00_projects/03_MetDNA2/00_data/20200805_targeted_annotation_evaluation_200STD/01_targeted_annotation/200STD_neg_reverse/annot_credential/semi_targeted_annotation_result.csv')
#
# test <- convertSemiTargetedResult2LongTable(result_semi_table)

setGeneric(name = 'convertSemiTargetedResult2LongTable',
           def = function(
             result_semi_table
           ){
             result <- result_semi_table %>%
               tidyr::separate_rows(ssc,
                                    int,
                                    mz_error,
                                    rt_error,
                                    ccs_error,
                                    int_ratio,
                                    ms2_score,
                                    isotope,
                                    type_annotation,
                                    annotation,
                                    peak_group,
                                    peak_group_scale,
                                    sep = ';') %>%
               readr::type_convert()

             return(result)
           })

################################################################################
# refine annotation ------------------------------------------------------------

# #' @title refineAnnotation
# #' @author Zhiwei Zhou
# #' @param annotation_initial
# #' @param list_peak_group_annotation_concised
# #' @export
#
# load('./inst/tempdata/list_peak_group_annotation_concised_pg2pg_200805.RData')
#
# annotation_initial <- readr::read_csv('I:/00_projects/03_MetDNA2/00_data/20200805_targeted_annotation_evaluation_200STD/01_targeted_annotation/200STD_neg_reverse/initial_id_result.csv')
#
# test <- refineAnnotation(annotation_initial = annotation_initial,
#                          list_peak_group_annotation_concised = list_peak_group_annotation_concised)
#
# setGeneric(name = 'refineAnnotation',
#            def = function(
#              annotation_initial,
#              list_peak_group_annotation_concised
#            ){
#              concised_table <- purrr::map(seq_along(list_peak_group_annotation_concised),
#                                           function(i){
#                                             result <- tibble::tibble(peak_name_concised = list_peak_group_annotation_concised[[i]]@base_peak_name,
#                                                                      adduct_concised = list_peak_group_annotation_concised[[i]]@base_peak_adduct)
#
#                                           }) %>%
#                dplyr::bind_rows() %>%
#                dplyr::mutate(label = paste(peak_name_concised, adduct_concised, sep = '_'))
#
#              annotation_concise <- annotation_initial %>%
#                dplyr::mutate(label = paste(name, adduct, sep = '_')) %>%
#                dplyr::left_join(concised_table, by = 'label') %>%
#                dplyr::filter(!is.na(peak_name_concised)) %>%
#                dplyr::select(-c('peak_name_concised', 'adduct_concised'))
#
#              num_peak_initial <- unique(annotation_initial$name) %>% length()
#              num_id_initial <- unique(annotation_initial$id) %>% length()
#              num_annotation_initial <- nrow(annotation_initial)
#
#              num_peak_concise <- unique(annotation_concise$name) %>% length()
#              num_id_concise <- unique(annotation_concise$id) %>% length()
#              num_annotation_concise <- nrow(annotation_concise)
#
#              num_stat <- tibble::tibble('No. peak' = c(num_peak_initial,
#                                                        num_peak_concise),
#                                         'No. metabolite' = c(num_id_initial,
#                                                              num_id_concise),
#                                         'No. annotation (peak-met)' = c(num_annotation_initial,
#                                                                         num_annotation_concise))
#
#              cat('Summary of ID refinement\n')
#              rownames(num_stat) <- c('Initial ID', 'Credential ID')
#              print(knitr::kable(num_stat))
#
#              return(annotation_concise)
#
#            })
#
#
