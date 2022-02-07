################################################################################
# Plots in RecursiveAnnotation -------------------------------------------------
#   plotRtCalibration ---------------------------------------------------------
#' @title plotRtCalibration
#' @author Zhiwei Zhou
#' @param result_rt_calibration
#' @param path '.'
#' @export

setGeneric(name = 'plotRtCalibration',
           def = function(
             result_rt_calibration,
             rt_recalibration_model,
             path = '.'
           ){

             cat("\n")

             pdf(file = file.path(path, 'plot_rt_calibration.pdf'),
                 width = 6,
                 height = 6)

             color <- c("black", "#F8766D", "#D89000", "#A3A500", "#39B600", "#00BF7D", "#00BFC4", "#00B0F6", "#9590FF", "#E76BF3", "#FF62BC", "gray", "#00FF00", "#0066CC", "#ff99bb", "#8833ff", "#22ddff", "#ffcc55", "#bbff33", "#00eeaa", "#5599ff", "#bb4400", "#008833", "black")

             plot(x = result_rt_calibration$ref.rt,
                  y = result_rt_calibration$exp.rt,
                  pch = 19,
                  cex = 1.25,
                  col=color[1:nrow(result_rt_calibration)],
                  xlab = "Reference RT of RTQC (s)",
                  ylab = "Experimental RT of RTQC (s)",
                  main = "Calibration Model")

             lines(lowess(rt_recalibration_model),lty=2,lwd=1)

             r2 <- lm(rt_recalibration_model$fitted~result_rt_calibration$exp.rt)
             r2 <- round(summary(r2)$r.squared, digits = 4)
             RMSE <- round(sqrt(mean((summary(rt_recalibration_model)$residuals)^2)),
                           digits = 4)

             text(x = min(result_rt_calibration$ref.rt)-10,
                  y = 0.95*max(result_rt_calibration$exp.rt),
                  labels = paste("R-Square:", r2),
                  pos = 4,
                  cex = 1)

             text(x = min(result_rt_calibration$ref.rt)-10,
                  y = 0.88*max(result_rt_calibration$exp.rt),
                  labels = paste("RMSE:", RMSE),
                  pos = 4,
                  cex = 1)

             legend("bottomright",
                    legend = rownames(result_rt_calibration),
                    fill = color[1:nrow(result_rt_calibration)],
                    border = NA,
                    bty = "n")

             dev.off()

           })

#   plotIdMs2 ------------------------------------------------------------------
#' @title plotIdMs2
#' @description generate ms2 mirror plot
#' @author Zhiwei Zhou
#' @param obj_ms2 object of SpectraData
#' @param obj_spec
#' @export
#' @examples

setGeneric(name = 'plotIdMs2',
           def = function(
             obj_ms2 = NULL,
             obj_spec = NULL
           ){
             # mode <- match(mode)

             if (length(obj_ms2) > 0) {
               obj_spec <- obj_ms2@matchedFragments[[1]]
             }

             if (length(obj_spec) == 0) {
               stop('Please input obj_spec')
             }

             temp_spec <- obj_spec %>%
               # dplyr::mutate(int = tidyr::replace_na(int, 0)) %>%
               dplyr::mutate(int_lib = intensity/max(intensity),
                             int_exp = intensityExp/max(intensityExp)) %>%
               dplyr::rename(mz_lib = mz,
                             mz_exp = mzExp) %>%
               dplyr::select(mz_lib, int_lib, mz_exp, int_exp, fragPrecursor) %>%
               dplyr::mutate(label = dplyr::case_when(
                 int_lib > 0 & int_exp > 0 ~ 'matched',
                 !(int_lib > 0 & int_exp > 0) ~ 'unmatched'
               ))

             # switch (mode,
             #         'SpecLibMatch' = {
             #
             #         },
             #         'NeighborMatch' = {
             #
             #         }
             # )

             temp_spec1 <- temp_spec %>%
               dplyr::select(mz_lib, int_lib, fragPrecursor, label) %>%
               dplyr::rename(mz = mz_lib, int = int_lib) %>%
               # dplyr::mutate(int = 0-int,
               #               tag = 'library')
               dplyr::mutate(int = 0-int,
                             tag = dplyr::case_when(
                               label == 'unmatched' ~ 'frag_unmatch',
                               label == 'matched' ~ 'library'
                             ))

             temp_spec2 <- temp_spec %>%
               dplyr::select(mz_exp, int_exp, fragPrecursor, label) %>%
               dplyr::rename(mz = mz_exp, int = int_exp) %>%
               # dplyr::mutate(tag = 'experiment')
               dplyr::mutate(tag = dplyr::case_when(
                               label == 'unmatched' ~ 'frag_unmatch',
                               label == 'matched' ~ 'experiment'
                             ))


             temp_data <- temp_spec1 %>%
               dplyr::bind_rows(temp_spec2)


             temp_plot <- ggplot2::ggplot(temp_data) +
               ggplot2::geom_segment(ggplot2::aes(x = mz, xend = mz,
                                                  y = 0, yend = int,
                                                  colour = tag)) +
               ggplot2::geom_point(ggplot2::aes(x = mz,
                                                y = int,
                                                shape = label,
                                                colour = tag))+
               ggplot2::scale_colour_manual(values = c(
                 'experiment' = 'dodgerblue',
                 'library' = 'tomato',
                 'frag_unmatch' = 'gray'
               )) +
               ggplot2::scale_shape_manual(values = c(
                 'matched' = 16,
                 'unmatched' = 4
               )) +
               ggplot2::geom_hline(yintercept = 0) +
               ggplot2::xlim(0.95*min(temp_data$mz),
                             1.05*max(temp_data$mz)) +
               ggplot2::xlab('m/z') +
               ggplot2::ylab('Relative intensity') +
               ZZWTheme() +
               # ggplot2::ggtitle(label = paste0(obj_ms2@info$name, '{',
               #                                 round(obj_ms2@info$scoreForward, 4), '}')) +
               ggplot2::theme(legend.position = c(0.8, 0.75),
                              title = ggplot2::element_text(vjust = 0.5))


             return(temp_plot)

           })

#   plotIdShiftMs2 ---------------------------------------------------------------
#' @title plotIdShiftMs2
#' @description generate ms2 mirror plot for shift match (bonanza, hybrid, gnps)
#' @author Zhiwei Zhou
#' @param obj_ms2 object of SpectraData
#' @export
#' @examples

# load('/home/zhouzw/Data_processing/20210224_metdna2_development_test/obj_ms2_cpd1_for_metdna2.RData')
# load('/home/zhouzw/Data_processing/20210224_metdna2_development_test/obj_ms2_cpd2_for_metdna2.RData')
# score_dp <- runSpecMatch(obj_ms2_cpd1 = obj_ms2_cpd1, obj_ms2_cpd2 = obj_ms2_cpd2, mz_tol_ms2 = 35, scoring_approach = 'dp')
# score_bonanza <- runSpecMatch(obj_ms2_cpd1 = obj_ms2_cpd1, obj_ms2_cpd2 = obj_ms2_cpd2, mz_tol_ms2 = 35, scoring_approach = 'bonanza')
# score_hybrid <- runSpecMatch(obj_ms2_cpd1 = obj_ms2_cpd1, obj_ms2_cpd2 = obj_ms2_cpd2, mz_tol_ms2 = 35, scoring_approach = 'hybrid')
# score_gnps <- runSpecMatch(obj_ms2_cpd1 = obj_ms2_cpd1, obj_ms2_cpd2 = obj_ms2_cpd2, mz_tol_ms2 = 35, scoring_approach = 'gnps')
#
# plotIdShiftMs2(obj_ms2 = score_bonanza)
# plotIdShiftMs2(obj_ms2 = score_hybrid)

setGeneric(name = 'plotIdShiftMs2',
           def = function(
             obj_ms2 = NULL,
             obj_spec_frag_match,
             obj_spec_frag_nl
           ){
             if (length(obj_ms2) > 0) {
               obj_spec_frag_match <- obj_ms2@matchedFragments[[1]]
               obj_spec_frag_nl <- obj_ms2@nlFragments[[1]]
             }

             if (length(obj_spec_frag_match) == 0 & length(obj_spec_frag_nl) == 0) {
               stop('Please input obj_spec_frag_match and obj_spec_frag_nl')
             }


             obj_spec_frag_match <- obj_spec_frag_match %>%
               dplyr::mutate(type = 'exact_match')
             obj_spec_frag_nl <- obj_spec_frag_nl %>%
               dplyr::mutate(type = 'nl_match')

             temp_spec_frag_match <- obj_spec_frag_match %>%
               dplyr::rename(mz_lib = mz,
                             mz_exp = mzExp,
                             int_lib = intensity,
                             int_exp = intensityExp) %>%
               dplyr::select(mz_lib, int_lib, mz_exp, int_exp, fragPrecursor, type) %>%
               dplyr::mutate(label = dplyr::case_when(
                 int_lib > 0 & int_exp > 0 ~ 'matched',
                 !(int_lib > 0 & int_exp > 0) ~ 'unmatched'
               ))

             temp_spec_frag_nl <- obj_spec_frag_nl %>%
               dplyr::rename(mz_lib = mz,
                             mz_exp = mzExp,
                             int_lib = intensity,
                             int_exp = intensityExp) %>%
               dplyr::select(mz_lib, int_lib, mz_exp, int_exp, fragPrecursor, type) %>%
               dplyr::filter(int_lib > 0 & int_exp > 0) %>%
               dplyr::mutate(mz_lib = dplyr::case_when(type == 'exact_match' ~ mz_lib,
                                                       type == 'nl_match' ~ mz_exp)) %>%
               dplyr::mutate(label = dplyr::case_when(
                 int_lib > 0 & int_exp > 0 ~ 'matched',
                 !(int_lib > 0 & int_exp > 0) ~ 'unmatched'
               ))

             temp_spec <- temp_spec_frag_match %>%
               dplyr::bind_rows(temp_spec_frag_nl)

             temp_spec1 <- temp_spec %>%
               dplyr::select(mz_lib, int_lib, fragPrecursor, type, label) %>%
               dplyr::mutate(int_lib = int_lib/max(int_lib)) %>%
               dplyr::rename(mz = mz_lib, int = int_lib) %>%
               dplyr::mutate(int = 0-int,
                             tag = dplyr::case_when(
                               label == 'unmatched' ~ 'frag_unmatch',
                               label == 'matched' & type == 'exact_match' ~ 'library',
                               label == 'matched' & type == 'nl_match' ~ 'library_shift'
                             ))

             temp_spec2 <- temp_spec %>%
               dplyr::select(mz_exp, int_exp, fragPrecursor, type, label) %>%
               dplyr::mutate(int_exp = int_exp/max(int_exp)) %>%
               dplyr::rename(mz = mz_exp, int = int_exp) %>%
               dplyr::mutate(tag = dplyr::case_when(
                 label == 'unmatched' ~ 'frag_unmatch',
                 label == 'matched' ~ 'experiment'
               ))


             temp_data <- temp_spec1 %>%
               dplyr::bind_rows(temp_spec2)

             temp_plot <- ggplot2::ggplot(temp_data) +
               ggplot2::geom_segment(ggplot2::aes(x = mz, xend = mz,
                                                  y = 0, yend = int,
                                                  colour = tag)) +
               ggplot2::geom_point(ggplot2::aes(x = mz,
                                                y = int,
                                                shape = label,
                                                colour = tag))+
               ggplot2::scale_colour_manual(values = c(
                 'experiment' = 'dodgerblue',
                 'library' = 'tomato',
                 'frag_unmatch' = 'gray',
                 'library_shift' = 'orange'
               )) +
               ggplot2::scale_shape_manual(values = c(
                 'matched' = 16,
                 'unmatched' = 4
               )) +
               ggplot2::geom_hline(yintercept = 0) +
               ggplot2::xlim(0.95*min(temp_data$mz),
                             1.05*max(temp_data$mz)) +
               ggplot2::ylim(-1, 1) +
               ggplot2::xlab('m/z') +
               ggplot2::ylab('Relative intensity') +
               ZZWTheme() +
               # ggplot2::ggtitle(label = paste0(obj_ms2@info$name, '{',
               #                                 round(obj_ms2@info$scoreForward, 4), '}')) +
               ggplot2::theme(legend.position = c(0.8, 0.75),
                              title = ggplot2::element_text(vjust = 0.5))


             return(temp_plot)
           })


#   plotNeighbor ---------------------------------------------------------------

#' @title plotNeighbor
#' @description export neighbor match plot during recursive annotation
#' @author Zhiwei Zhou
#' @param tags2
#' @param path '.'
#' @param ms2 ms2 RData
#' @param kegg.compound
#' @export

# path <- '/home/zhouzw/Data_processing/20210309_test_metdna_v0.40/aging_fly_POS/'
# load(file.path(path, "02_result_MRN_annotation", '00_intermediate_data', 'tags2_after_redundancy_remove'))
# load(file.path(path, '01_result_initial_seed_annotation', '00_intermediate_data', 'ms2'))
# load(file.path(path, "02_result_MRN_annotation", '00_intermediate_data', 'matchParam_recursive'))
# tags2 <- tags2_after_redundancy_remove
# # kegg.compound <- cpd_emrn
#
# plotNeighbor(tags2 = tags2,
#              is.include.precursor = TRUE,
#              path = file.path(path, "02_result_MRN_annotation"),
#              ms2 = ms2,
#              kegg.compound = cpd_emrn,
#              matchParam = matchParam)


setGeneric(name = "plotNeighbor",
           def = function(tags2,
                          path = ".",
                          ms2 = ms2,
                          kegg.compound = kegg.compound,
                          matchParam = matchParam,
                          ...){

             ms2.name <- unname(unlist(lapply(ms2, function(x){
               x[[1]][1,1]
             })))

             index <- which(showTags2(tags2 = tags2, slot = "annotation.len") > 0)
             if(length(index) == 0) return("No annotaion")
             # tags2 <- tags2[index]
             annotation.type <- showTags2(tags2 = tags2, slot = "annotation.type")

             index <- unname(which(unlist(lapply(annotation.type, function(x){
               any(unlist(x) == "metAnnotation")
             }))))

             annotation.type <- annotation.type[index]

             peak.index <- as.numeric(names(annotation.type))
             annotation.index <- lapply(annotation.type, function(x){
               which(x == "metAnnotation")
             })

             pb <- utils::txtProgressBar(min = 0, max = length(peak.index), style = 3)
             for(i in 1:length(peak.index)){
               # cat(match(i, peak.index)); cat(" ")
               # cat('i', i); cat(" ")
               utils::setTxtProgressBar(pb, i)
               temp.neighbor <- tags2[[peak.index[i]]]
               temp.annotation.index <- annotation.index[[i]]
               neighbor.peak.name <- temp.neighbor@name
               neighbor.mz <- temp.neighbor@mz

               path_output <- file.path(path, '01_surrogate_ms2_spec_plot')

               for(j in temp.annotation.index){
                 # cat('j', j); cat(" ")
                 temp.annotation <- temp.neighbor@annotation[[j]]
                 neighbor.id <- temp.annotation$to
                 neighbor.compound.name <- kegg.compound$name[match(neighbor.id, kegg.compound$id)]
                 neighbor.dp <- temp.annotation$ms2.sim
                 neighbor_info <- data.frame(name = neighbor.peak.name,
                                             mz = neighbor.mz,
                                             cpd_id = neighbor.id,
                                             cpd_name = neighbor.compound.name,
                                             stringsAsFactors = FALSE)
                 neighbor_ms2 <- ms2[[match(neighbor.peak.name, ms2.name)]][[2]]

                 obj_neighbor_ms2 <- new('SpectraData',
                                         info = neighbor_info,
                                         spectra = list(neighbor_ms2))

                 seed.peak.name <- temp.annotation$From.peak
                 seed.mz <- tags2[[match(seed.peak.name, showTags2(tags2, slot = "name"))]]@mz
                 seed.id <- temp.annotation$From
                 seed.compound.name <- kegg.compound$name[match(seed.id, kegg.compound$id)]
                 seed_info <- data.frame(name = seed.peak.name,
                                         mz = seed.mz,
                                         cpd_id = seed.id,
                                         cpd_name = seed.compound.name,
                                         stringsAsFactors = FALSE)

                 seed_ms2 <- ms2[[match(seed.peak.name, ms2.name)]][[2]]
                 obj_seed_ms2 <- new('SpectraData',
                                     info = seed_info,
                                     spectra = list(seed_ms2))

                 if (is.null(seed_ms2) | is.null(neighbor_ms2)){
                   # x <- c(x, seed.peak.name)
                   next()
                 } else {

                   # matchParam <- SpectraTools::MatchParam(ppm = 25,
                   #                                        methodScore = 'dp',
                   #                                        methodMatch = 'direct',
                   #                                        weightIntensity = 1,
                   #                                        weightMZ = 0,
                   #                                        cutoff = 0,
                   #                                        includePrecursor = TRUE,
                   #                                        intensityExpNormed = TRUE,
                   #                                        intensityLibNormed = TRUE,
                   #                                        tuneLibSpectra = FALSE)

                   temp_match_result <- try(SpectraTools::MatchSpectra(obj_neighbor_ms2,
                                                                       obj_seed_ms2,
                                                                       matchParam),
                                            silent = TRUE)

                   # for only precursor matched error
                   if (length(temp_match_result@matchedFragments$`#1`) == 0) {
                     next()
                   }

                   if (length(temp_match_result) == 0) {
                     next()
                   }

                   if (class(temp_match_result) == 'try-error') {
                     next()
                   }

                   if (matchParam@methodScore == 'dp') {
                     suppressMessages(
                       temp_plot <- plotIdMs2(temp_match_result) +
                         ggplot2::scale_colour_manual(
                           name = 'attribute',
                           labels= c(paste0('neigbhor ', '(', neighbor_info$name, ')'),
                                     'Unmatched fragments',
                                     paste0('seed ', '(', seed_info$name, ')')),
                           values = c(
                             'experiment' = 'dodgerblue',
                             'library' = 'tomato',
                             'frag_unmatch' = 'gray'
                           )
                         ) +
                         ggplot2::ggtitle(label = paste0(seed_info$cpd_id, ' ---> ', neighbor_info$cpd_id,
                                                         ' (DP: ', round(neighbor.dp, 4), ')'))
                     )
                   } else {
                     suppressMessages(
                       temp_plot <- plotIdShiftMs2(temp_match_result) +
                         ggplot2::scale_colour_manual(
                           name = 'attribute',
                           labels= c(paste0('neigbhor ', '(', neighbor_info$name, ')'),
                                     'unmatched fragments',
                                     paste0('seed ', '(', seed_info$name, ')'),
                                     paste0('seed shifted', '(', seed_info$name, ')')),
                           values = c(
                             'experiment' = 'dodgerblue',
                             'library' = 'tomato',
                             'frag_unmatch' = 'gray',
                             'library_shift' = 'orange'
                           )
                         ) +
                         ggplot2::ggtitle(label = paste0(seed_info$cpd_id, ' ---> ', neighbor_info$cpd_id,
                                                         ' (', matchParam@methodScore, ': ', round(neighbor.dp, 4), ')'))
                     )
                   }


                   dir.create(file.path(path_output, neighbor_info$name), showWarnings = FALSE, recursive = TRUE)
                   plot_name <- paste0(seed_info$name, '_', neighbor_info$name,
                                       '(', seed_info$cpd_id, '_', neighbor_info$cpd_id, ').pdf')
                   temp_name <- file.path(path_output, neighbor_info$name, plot_name)
                   ggplot2::ggsave(temp_plot,
                                   filename = temp_name,
                                   width = 10, height = 6)

                 }
               }
             }

           })


################################################################################
# Plots in AnnotationCredential ------------------------------------------------
#   plotPseudoMs1Spec ----------------------------------------------------------
#' @title plotPseudoMs1Spec
#' @description generate Pseudo MS1 spectrum for annotated peak group
#' @author Zhiwei Zhou
#' @param peak_group
#' @export
#' @examples
#' load(system.file("tempdata", "list_peak_group_annotation_200805.RData", package="MetDNA2"))
#' list_peak_group_annotation[[35]] %>% plotPseudoMs1Spec()

# load('./inst/tempdata/list_peak_group_annotation_200805.RData')
# list_peak_group_annotation[[35]] %>% plotPseudoMs1Spec()

# load('I:/00_projects/03_MetDNA2/00_data/20200710_semi_targeted_annotation_test/03_test/list_peak_group_200716.RData')
# load('I:/00_projects/03_MetDNA2/00_data/20200710_semi_targeted_annotation_test/03_test/raw_msms_200715.RData')
# peak_group_L0076 <- list_peak_group$`M482T929_[M-H]-`
#
# test_L0076 <- annotatePeakGroup(peak_group = peak_group_L0076,
#                                 ms2_data = raw_msms,
#                                 polarity = 'negative',
#                                 tol_mz = 10,
#                                 is_ms2_check = TRUE,
#                                 ms2_score_cutoff = -1)
# #
# save(test_L0076, file = './peak_group_test_L0076_200716.RData', version = 2)
#
# plotPseudoMs1Spec(peak_group = test_L0076)

# peak_group_L0177 <- list_peak_group$`M303T808_[M-H]-`
#
# test_L0177 <- annotatePeakGroup(peak_group = peak_group_L0177,
#                                 ms2_data = raw_msms,
#                                 polarity = 'negative',
#                                 tol_mz = 10,
#                                 is_ms2_check = TRUE,
#                                 ms2_score_cutoff = -1)
#
# plotPseudoMs1Spec(peak_group = test_L0177) + theme(legend.position = 'none')

#
# peak_group_L0306 <- list_peak_group$`M175T635_[M-H]-`
#
# test_L0306 <- annotatePeakGroup(peak_group = peak_group_L0306,
#                                 ms2_data = raw_msms,
#                                 polarity = 'negative',
#                                 tol_mz = 10,
#                                 is_ms2_check = TRUE,
#                                 ms2_score_cutoff = -1)
#
# plotPseudoMs1Spec(peak_group = test_L0306) + theme(legend.position = 'none')

setGeneric(name = 'plotPseudoMs1Spec',
           def = function(
             peak_group
           ){
             if (nrow(peak_group@peak_list_annotated) == 0) {
               stop('No annotated MS1 peaks\n')
             }
             temp_data <- peak_group@peak_list_annotated %>%
               tidyr::pivot_longer(cols = c('isotopeAnnotation', 'adductAnnotation',
                                            'neutralLossAnnotation', 'isfAnnotation'),
                                   names_to = 'type',
                                   values_to = 'label') %>%
               dplyr::filter(!is.na(label)) %>%
               dplyr::arrange(peak_name,
                              match(type, c('adductAnnotation',
                                            'neutralLossAnnotation',
                                            'isfAnnotation',
                                            'isotopeAnnotation'))) %>%
               dplyr::distinct(peak_name, .keep_all = TRUE) %>%
               dplyr::mutate(int = int/max(int))


             base_peak_mz <- peak_group@base_peak_mz
             base_peak_int <- temp_data %>%
               dplyr::filter(peak_name == peak_group@base_peak_name) %>%
               dplyr::pull(int)
             base_peak_adduct <- peak_group@base_peak_adduct


             temp_plot <- ggplot2::ggplot(temp_data) +
               ggplot2::geom_segment(ggplot2::aes(x = mz, xend = mz,
                                                  y = 0, yend = int,
                                                  colour = type)) +
               ggplot2::geom_point(ggplot2::aes(x = base_peak_mz, y = base_peak_int))+
               ggplot2::scale_colour_manual(values = c(
                 'adductAnnotation' = 'tomato',
                 'neutralLossAnnotation' = 'dodgerblue',
                 'isfAnnotation' = 'orange',
                 'isotopeAnnotation' = 'limegreen'
               )) +
               ggplot2::geom_hline(yintercept = 0) +
               ggplot2::xlim(0.95*min(temp_data$mz),
                             1.05*max(temp_data$mz)) +
               ggplot2::ylim(0, 1.05) +
               ggplot2::xlab('m/z') +
               ggplot2::ylab('relative intensity') +
               ZZWTheme() +
               ggplot2::theme(legend.position = c(0.8, 0.8))

             # add labels in the spectrum
             temp_adduct <- temp_data %>%
               dplyr::filter(type %in% c('adductAnnotation'))

             if (nrow(temp_adduct)>0) {
               temp_plot <- temp_plot +
                 ggplot2::annotate('text',
                                   x = temp_adduct$mz,
                                   y = temp_adduct$int + 0.04,
                                   label = temp_adduct$label)
             }

             temp_nl <- temp_data %>%
               dplyr::filter(type %in% c('neutralLossAnnotation'))

             if (nrow(temp_nl)>0) {

               # if protonation/deprotpnation base peak
               #    plot loss formula
               #   else only add label
               if (base_peak_adduct %in% c('[M+H]+', '[M-H]-')) {
                 temp_label <- temp_nl$label %>%
                   stringr::str_extract(pattern = 'M\\-(\\w+)') %>%
                   stringr::str_replace(pattern = 'M', replacement = '')

                 temp_nl <- temp_nl %>%
                   dplyr::select(peak_name:int, type, label) %>%
                   dplyr::mutate(xend = base_peak_mz)

                 temp_plot <- temp_plot +
                   ggplot2::geom_segment(ggplot2::aes(x = mz,
                                                      xend = xend,
                                                      y = 0.8*int,
                                                      yend = 0.8*int),
                                         linetype = 'dashed',
                                         data = temp_nl) +
                   ggplot2::annotate('text',
                                     x = (base_peak_mz - temp_nl$mz)/2 + temp_nl$mz,
                                     y = temp_nl$int*0.8 + 0.04,
                                     label = temp_label)

                 # temp_plot <- temp_plot +
                 #   ggplot2::geom_segment(aes(x = temp_nl$mz,
                 #                             xend = rep(base_peak_mz, nrow(temp_nl)),
                 #                             y = 0.8*temp_nl$int,
                 #                             yend = 0.8*temp_nl$int),
                 #                         linetype = 'dashed') +
                 #   ggplot2::annotate('text',
                 #                     x = (base_peak_mz - temp_nl$mz)/2 + temp_nl$mz,
                 #                     y = temp_nl$int*0.8 + 0.04,
                 #                     label = temp_label)
               } else {
                 temp_plot <- temp_plot +
                   ggplot2::annotate('text',
                                     x = temp_nl$mz,
                                     y = temp_nl$int + 0.04,
                                     label = temp_nl$label)
               }

             }

             temp_isf <- temp_data %>%
               dplyr::filter(type %in% c('isfAnnotation'))

             if (nrow(temp_adduct)>0) {
               temp_plot <- temp_plot +
                 ggplot2::annotate('text',
                                   x = temp_isf$mz,
                                   y = temp_isf$int + 0.04,
                                   label = temp_isf$label)

             }


             return(temp_plot)

           })




################################################################################
# Plots in BiologyIntepretation ------------------------------------------------
#   plotVolanco ----------------------------------------------------------------
#' @title plotVolcano
#' @author Zhiwei Zhou
#' @param raw_data 1st column: feature name; 2nd column: p_values; 3rd column: fold_changes
#' @param p_cutoff the p cutoff of significant; All p-values would be -log10 transformed in the plot. Default: 0.05
#' @param fc_cutoff the fold change cutoff of significant: Default: 1.25
#' @param log_tran log_transformation method for fold change, "log2" or "log10". Default: log2
#' @export
#' @examples
#' data('stat_result')
#' plotVolcano(raw_data = stat_result, p_cutoff = 0.05, fc_cutoff = 1.25, log_tran = 'log2')

# plotVolcano(test)


setGeneric(name = 'plotVolcano',
           def = function(
             raw_data,
             p_cutoff = 0.05,
             fc_cutoff = 1.25,
             log_tran = c('log2', 'log10')
           ){
             log_tran <- match.arg(log_tran)
             colnames(raw_data) <- c('feature', 'p_values', 'fold_changes')

             switch (log_tran,
                     'log2' = {
                       log_p <- -log10(raw_data$p_values)
                       log_fc <- log2(raw_data$fold_changes)
                     },
                     'log10' = {
                       log_p <- -log10(raw_data$p_values)
                       log_fc <- log10(raw_data$fold_changes)
                     }
             )

             volcano_result <- raw_data %>%
               dplyr::mutate(log_p = log_p,
                             log_fc = log_fc) %>%
               dplyr::mutate(is_sig=dplyr::case_when(
                 (p_values <= p_cutoff) & (fold_changes >= fc_cutoff) ~ 'Increase',
                 (p_values <= p_cutoff) & (fold_changes <= (1/fc_cutoff)) ~ 'Decrease',
                 (p_values > p_cutoff) | (fold_changes > (1/fc_cutoff)) | (fold_changes > fc_cutoff) ~ 'Not_significant'
               ))

             volcano_plot <- ggplot2::ggplot(data = volcano_result) +
               ggplot2::geom_point(ggplot2::aes(x = log_fc, y = log_p, colour = is_sig)) +
               ggplot2::geom_hline(yintercept = -log10(p_cutoff), linetype=2) +
               ggplot2::geom_vline(xintercept = ifelse(log_tran=='log2', log2(fc_cutoff), log10(fc_cutoff)),
                                   linetype = 2) +
               ggplot2::geom_vline(xintercept = ifelse(log_tran=='log2', -log2(fc_cutoff), -log10(fc_cutoff)),
                                   linetype = 2) +
               ggplot2::scale_color_manual(name = 'Group',
                                           values = c('Increase' = 'tomato',
                                                      'Decrease' = 'dodgerblue',
                                                      'Not_significant' = 'gray'),
                                           label = c('Increase' = 'Increase',
                                                     'Decrease' = 'Decrease',
                                                     'Not_significant' = 'Not significant')) +
               ggplot2::xlab(paste0(log_tran, '(Fold-change)')) +
               ggplot2::ylab(paste0('-log10(P-value)')) +
               ZZWTheme() +
               ggplot2::theme(legend.position = c(0.8, 0.8))

             return(volcano_plot)

           }
)



#   plotPathwayScatter ---------------------------------------------------------
#' @title plotPathwayScatter
#' @author Zhiwei Zhou
#' @param result_pathway_enrichment
#' @param p_cutoff
#' @param p_adjust
#' @export

setGeneric(name = 'plotPathwayScatter',
           function(
             result_pathway_enrichment,
             p_cutoff = 0.05,
             p_adjust = TRUE
           ){
             temp_data <- result_pathway_enrichment %>%
               tibble::rownames_to_column(var = 'pathway_name') %>%
               tidyr::separate(pathway_name, into = c('pathway_name', 'pathway_id'), sep = ';') %>%
               tidyr::separate(pathway_name, into = c('pathway_name', 'temp'), sep = ' - ') %>%
               dplyr::select(-temp) %>%
               # tidyr::unite(pathway_name, pathway_id, col = 'pathway_name', sep = ';') %>%
               dplyr::mutate(pcg_overlap = (Overlap/Pathway.length)*100)

             if (p_adjust) {
               temp_data <- temp_data %>%
                 dplyr::mutate(label = dplyr::case_when(q.value <= p_cutoff ~ 'Significant',
                                                        q.value > p_cutoff ~ 'Not_significant')) %>%
                 dplyr::mutate(log_p = -log10(q.value))
             } else {
               temp_data <- temp_data %>%
                 dplyr::mutate(label = dplyr::case_when(p.value <= p_cutoff ~ 'Significant',
                                                        p.value > p_cutoff ~ 'Not_significant')) %>%
                 dplyr::mutate(log_p = -log10(p.value))
             }

             idx <- which(temp_data$log_p >= -log10(0.05))

             if (length(idx) > 0) {
               temp_plot <- ggplot2::ggplot(temp_data) +
                 ggplot2::geom_point(ggplot2::aes(x = pcg_overlap,
                                                  y = log_p,
                                                  size = Pathway.length,
                                                  colour = label)) +
                 ggplot2::scale_size(range = c(1, 15)) +
                 ggplot2::scale_colour_manual(values = c('Significant' = 'tomato',
                                                         'Not_significant' = 'gray')) +
                 ggplot2::geom_hline(yintercept = -log10(p_cutoff),
                                     linetype = 'dashed') +
                 ggplot2::xlab('Overlap (%)') +
                 ggplot2::ylab('-log10(P-value)') +
                 ggplot2::annotate(geom = 'text',
                                   x = temp_data$pcg_overlap[idx],
                                   y = temp_data$log_p[idx],
                                   label = temp_data$pathway_name[idx]) +
                 ZZWTheme() +
                 ggplot2::theme(legend.position = c(0.15, 0.7))
             } else {
               temp_plot <- ggplot2::ggplot(temp_data) +
                 ggplot2::geom_point(ggplot2::aes(x = pcg_overlap,
                                                  y = log_p,
                                                  size = Pathway.length,
                                                  colour = label)) +
                 ggplot2::scale_size(range = c(1, 15)) +
                 ggplot2::scale_colour_manual(values = c('Significant' = 'tomato',
                                                         'Not_significant' = 'gray')) +
                 ggplot2::geom_hline(yintercept = -log10(p_cutoff),
                                     linetype = 'dashed') +
                 ggplot2::xlab('Overlap (%)') +
                 ggplot2::ylab('-log10(P-value)') +
                 # ggplot2::annotate(geom = 'text',
                 #                   x = temp_data$pcg_overlap[idx],
                 #                   y = temp_data$log_p[idx],
                 #                   label = temp_data$pathway_name[idx]) +
                 ZZWTheme() +
                 ggplot2::theme(legend.position = c(0.15, 0.7))
             }


             return(temp_plot)

           })



#   plotPathwayBar -------------------------------------------------------------
#' @title plotPathwayBar
#' @author Zhiwei Zhou
#' \email{zhouzw@@sioc.ac.cn}
#' @param result_pathway_enrichment
#' @param p_cutoff
#' @param p_adjust
#' @export

setGeneric(name = 'plotPathwayBar',
           def = function(
             result_pathway_enrichment,
             p_cutoff = 0.05,
             p_adjust = TRUE
           ){
             temp_data <- result_pathway_enrichment %>%
               tibble::rownames_to_column(var = 'pathway_name') %>%
               tidyr::separate(pathway_name, into = c('pathway_name', 'pathway_id'), sep = ';') %>%
               tidyr::separate(pathway_name, into = c('pathway_name', 'temp'), sep = ' - ') %>%
               dplyr::select(-temp) %>%
               tidyr::unite(pathway_name, pathway_id, col = 'pathway_name', sep = ';') %>%
               dplyr::mutate(pcg_overlap = (Overlap/Pathway.length)*100)


             if (p_adjust) {
               temp_data <- temp_data %>%
                 dplyr::mutate(label = dplyr::case_when(q.value <= p_cutoff ~ 'Significant',
                                                        q.value > p_cutoff ~ 'Not_significant')) %>%
                 dplyr::mutate(log_p = -log10(q.value))
             } else {
               temp_data <- temp_data %>%
                 dplyr::mutate(label = dplyr::case_when(p.value <= p_cutoff ~ 'Significant',
                                                        p.value > p_cutoff ~ 'Not_significant')) %>%
                 dplyr::mutate(log_p = -log10(p.value))
             }

             temp_data <- temp_data %>% dplyr::arrange(dplyr::desc(log_p))

             temp_data$pathway_name <- factor(x = seq(length(temp_data$pathway_name)),
                                              levels = rev(seq(length(temp_data$pathway_name))),
                                              labels = rev(temp_data$pathway_name))

             temp_plot <- ggplot2::ggplot(temp_data) +
               ggplot2::geom_bar(ggplot2::aes(x = pathway_name, y = log_p, fill = label), stat = 'identity') +
               ggplot2::scale_fill_manual(values = c('Significant' = 'tomato',
                                                     'Not_significant' = 'gray')) +
               ggplot2::ggtitle(label = 'Pathway enrichment analysis') +
               ggplot2::xlab('Pathway') +
               ggplot2::ylab('-Log10(adjusted P-value)') +
               ggplot2::coord_flip() +
               ZZWTheme(type = 'classic') +
               ggplot2::theme(legend.position = 'none',
                              axis.text.y = ggplot2::element_text(angle = 0))

             return(temp_plot)

           })


#   plotPathwayHeatmap ---------------------------------------------------------
#' @title plotPathwayHeatmap
#' @author Zhiwei Zhou
#' \email{zhouzw@@sioc.ac.cn}
#' @param quanti_pathway_data
#' @param sample_info
#' @param comp_group
#' @export

setGeneric(name = 'plotPathwayHeatmap',
           def = function(
             quanti_pathway_data,
             sample_info,
             comp_group = c("W30", "W03")
           ){

             temp_data <- quanti_pathway_data %>%
               tidyr::separate(pathway_name, into = c('pathway_name', 'pathway_id'), sep = ';') %>%
               tidyr::separate(pathway_name, into = c('pathway_name', 'temp'), sep = ' - ') %>%
               dplyr::select(-c('temp', 'pathway_id')) %>%
               tibble::column_to_rownames(var = 'pathway_name') %>%
               as.matrix()

             sample_info <- sample_info %>%
               dplyr::filter(group %in% comp_group) %>%
               dplyr::arrange(group, sample.name)

             temp_data <- match(sample_info$sample.name, colnames(temp_data)) %>% temp_data[,.]
             annotation_col <- data.frame(
               Group = factor(sample_info$group)
             )
             rownames(annotation_col) <- colnames(temp_data)

             # Specify colors
             temp <- c("tomato", "dodgerblue")
             names(temp) <- comp_group
             ann_colors <- list(
               Group = temp
             )

             temp_plot <- pheatmap::pheatmap(temp_data,
                                             scale = "row",
                                             cluster_cols = TRUE,
                                             cluster_rows = TRUE,
                                             show_rownames = TRUE,
                                             show_colnames = FALSE,
                                             # annotation_row = annotation_row,
                                             annotation_col = annotation_col,
                                             annotation_colors = ann_colors,
                                             annotation_legend = TRUE,
                                             clustering_distance_rows = "euclidean",
                                             clustering_distance_cols = "euclidean",
                                             clustering_method = 'complete',
                                             border_color = NA,
                                             display_numbers = FALSE,
                                             number_color = 'black',
                                             number_format = '%.2f',
                                             color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

             return(temp_plot)

           })



#   plotPathwayBoxplot ---------------------------------------------------------
#' @title plotPathwayBoxplot
#' @param quanti_pathway_data
#' @param sample_info
#' @param comp_group
#' @param path '.'
#' @export


setGeneric(name = 'plotPathwayBoxplot',
           def = function(
             quanti_pathway_data,
             sample_info,
             comp_group = c("W30", "W03"),
             path = '.'
           ){
             quanti_pathway_data <- quanti_pathway_data %>%
               tidyr::separate(pathway_name, into = c('pathway_name', 'pathway_id'), sep = ';') %>%
               tidyr::separate(pathway_name, into = c('pathway_name', 'temp'), sep = ' - ') %>%
               dplyr::select(-c('temp', 'pathway_id'))

             path_output <- file.path(path, 'Boxplot')
             dir.create(path_output, showWarnings = FALSE, recursive = TRUE)

             quanti_pathway_data2 <- quanti_pathway_data %>%
               tidyr::pivot_longer(cols = -pathway_name, names_to = 'sample_name')

             temp_group <- match(quanti_pathway_data2$sample_name, sample_info$sample.name) %>%
               sample_info$group[.]
             quanti_pathway_data2 <- quanti_pathway_data2 %>%
               dplyr::mutate(group = temp_group) %>%
               dplyr::select(pathway_name, sample_name, group, value)

             unique_pathway <- quanti_pathway_data2$pathway_name %>% unique()

             progress <- mapProgress(n = length(unique_pathway))
             purrr::walk(seq_along(unique_pathway), function(i){
               mapProgressPrint(progress = progress)
               temp_pathway_name <- unique_pathway[i]
               plot_name <- temp_pathway_name %>%
                 make.names() %>%
                 stringr::str_replace(pattern = '\\.+', replacement = ' ') %>%
                 paste0('.pdf')

               temp_data <- quanti_pathway_data2 %>%
                 dplyr::filter(pathway_name == temp_pathway_name)

               assign_comparison <- list(comp_group)

               temp_plot <- ggplot2::ggplot(temp_data, ggplot2::aes(x = group, y = value)) +
                 ggplot2::geom_boxplot(ggplot2::aes(fill = group)) +
                 ggplot2::geom_point(colour = 'black') +
                 ggpubr::stat_compare_means(ggplot2::aes(x = group, y = value),
                                            comparisons = assign_comparison) +
                 ggplot2::scale_fill_manual(name = 'Group',
                                            values = c('dodgerblue', 'tomato')) +
                 ggplot2::ylab('Dysregulated score') +
                 ggplot2::xlab('Group') +
                 ggplot2::ggtitle(temp_pathway_name) +
                 ZZWTheme() +
                 ggplot2::theme(legend.position = c(0.9, 0.9))

               ggplot2::ggsave(temp_plot,
                               filename = file.path(path_output, plot_name),
                               width = 6, height = 6)
             })


           })



#   plotPathwayHeatmapMet ------------------------------------------------------
#' @title plotPathwayHeatmapMet
#' @author Zhiwei Zhou
#' \email{zhouzw@@sioc.ac.cn}
#' @param quanti_met_data
#' @param quanti_pathway_data
#' @param sample_info
#' @param lib_pathway
#' @param comp_group
#' @param species
#' @param path '.'
#' @export


setGeneric(name = 'plotPathwayHeatmapMet',
           def = function(
             quanti_met_data,
             quanti_pathway_data,
             sample_info,
             lib_pathway,
             comp_group = c("W30", "W03"),
             species = c("hsa", "mmu", "rno", "bta", "gga",
                         "dre", "dme", "cel", "sce", "osa",
                         "ath", "smm", "pfa", "tbr", "eco",
                         "bsu", "ppu", "sau", "tma", "syf", "mln"),
             path = '.',
             ...
           ){

             temp_pathway <- which(names(lib_pathway) == species) %>% lib_pathway[[.]]
             temp_cpd_name_synonyms <- quanti_met_data %>%
               tidyr::separate_rows('cpd_name_synonyms', sep = ';') %>%
               dplyr::pull(cpd_name_synonyms)
             temp_metabolite_quant_data <- quanti_met_data %>%
               tidyr::separate_rows('id_synonyms', sep = ';') %>%
               dplyr::mutate(cpd_name_synonyms = temp_cpd_name_synonyms)
             temp_pathway_name <- names(temp_pathway)
             # temp_pathway_name <- temp_pathway %>%
             #   names() %>%
             #   stringr::str_split(pattern = ';', n = 2) %>%
             #   sapply(., function(x)x[1]) %>%
             #   stringr::str_split(pattern = ' - ', n = 2) %>%
             #   sapply(., function(x)x[1])

             unique_pathway <- quanti_pathway_data$pathway_name
             idx <- match(unique_pathway, temp_pathway_name)
             temp_pathway <- temp_pathway[idx]

             temp_quanti_pathway_data <- lapply(seq_along(temp_pathway), function(i){
               temp_met <- temp_pathway[[i]]
               temp_metabolite_quant_data <- temp_metabolite_quant_data %>%
                 dplyr::filter(id_synonyms %in% temp_met)
             })

             names(temp_quanti_pathway_data) <- names(temp_pathway)


             path_output <- file.path(path, 'Heatmap')
             dir.create(path_output, showWarnings = FALSE, recursive = TRUE)

             unique_pathway <- unique_pathway %>%
               stringr::str_split(pattern = ';', n = 2) %>%
               sapply(., function(x)x[1]) %>%
               stringr::str_split(pattern = ' - ', n = 2) %>%
               sapply(., function(x)x[1])

             progress <- mapProgress(n = length(temp_quanti_pathway_data))
             purrr::walk(seq_along(temp_quanti_pathway_data), function(i){
               mapProgressPrint(progress = progress)
               temp_pathway_name <- unique_pathway[i]
               plot_name <- temp_pathway_name %>%
                 make.names() %>%
                 stringr::str_replace(pattern = '\\.+', replacement = ' ') %>%
                 paste0('.pdf')

               temp_data <- temp_quanti_pathway_data[[i]] %>%
                 dplyr::select(-c('inchikey1', 'peak_name', 'id', 'cpd_name',
                                  'id_synonyms', 'confidence_level')) %>%
                 tibble::column_to_rownames(var = 'cpd_name_synonyms') %>%
                 as.matrix()

               sample_info <- sample_info %>%
                 dplyr::filter(group %in% comp_group) %>%
                 dplyr::arrange(group, sample.name)

               temp_data <- match(sample_info$sample.name, colnames(temp_data)) %>% temp_data[,.]
               annotation_col <- data.frame(
                 Group = factor(sample_info$group)
               )
               rownames(annotation_col) <- colnames(temp_data)

               # Specify colors
               temp <- c("tomato", "dodgerblue")
               names(temp) <- comp_group
               ann_colors <- list(
                 Group = temp
               )

               temp_plot <- pheatmap::pheatmap(temp_data,
                                               scale = "row",
                                               cluster_cols = TRUE,
                                               cluster_rows = TRUE,
                                               show_rownames = TRUE,
                                               show_colnames = FALSE,
                                               # annotation_row = annotation_row,
                                               annotation_col = annotation_col,
                                               annotation_colors = ann_colors,
                                               annotation_legend = TRUE,
                                               clustering_distance_rows = "euclidean",
                                               clustering_distance_cols = "euclidean",
                                               clustering_method = 'complete',
                                               border_color = NA,
                                               display_numbers = FALSE,
                                               number_color = 'black',
                                               number_format = '%.2f',
                                               color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

               ggplot2::ggsave(temp_plot,
                               filename = file.path(path_output, plot_name),
                               width = 6, height = 6)

             })

           })

################################################################################
# Other functions --------------------------------------------------------------
#   ZZWTheme ------------------------------------------------------------------
#' @title ZZWTheme
#' @author Zhiwei Zhou
#' \email{zhouzw@@sioc.ac.cn}
#' @description ggplot theme
#' @param type: 'common', 'classic'
#' @examples
#' mydata <- data.frame(
#'  Lebal  = c("linerange1","linerange2","linerange3","linerange4","linerange5"),
#'  xstart = c(3.5,7,12,16,20),
#'  ymin   = c(2.5,6.5,3,4.5,3.8),
#'  ymax   = c(7.5,9.5,9,13.5,4.2),
#'  class  = c("A","A","A","C","C")
#' )
#' test <- ggplot(mydata) +
#'  geom_linerange(aes(x = xstart, ymin = ymin , ymax = ymax , colour = class) , size = 1.5)
#' test
#' mytheme <- ZZWTheme()
#' test + mytheme
#' mytheme2 <- ZZWTheme('classic')
#' test +
#'   mytheme2 +
#'   theme(axis.text.y = ggplot2::element_text(size = 10, angle = 0,
#'                                             colour = 'black'))


setGeneric(name = 'ZZWTheme',
           def = function(
             type = c('common', 'classic')
           ){

             type <- match.arg(type)
             switch (type,
                     'common' = {
                       result <- ggplot2::theme_bw() +
                         ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.4,
                                                                           face = "bold",
                                                                           size = 14),
                                        panel.grid.major = ggplot2::element_blank(),
                                        panel.grid.minor = ggplot2::element_blank(),
                                        panel.background = ggplot2::element_blank(),
                                        panel.border = ggplot2::element_rect(fill = NA, colour = 'black'),
                                        legend.background = ggplot2::element_blank(),
                                        axis.text.x = ggplot2::element_text(size = 10, colour = 'black'),
                                        axis.text.y = ggplot2::element_text(size = 10, angle = 90,
                                                                            colour = 'black'),
                                        axis.title.x = ggplot2::element_text(size = 12, colour = 'black'),
                                        axis.title.y = ggplot2::element_text(size = 12, colour = 'black'),
                                        axis.ticks = ggplot2::element_line(colour = 'black'))
                     },
                     'classic' = {
                       result <- ggplot2::theme_classic() +
                         ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.4,
                                                                           face = "bold",
                                                                           size = 14),
                                        panel.grid.major = ggplot2::element_blank(),
                                        panel.grid.minor = ggplot2::element_blank(),
                                        panel.background = ggplot2::element_rect(fill = NA, colour = NA),
                                        # panel.border = ggplot2::element_rect(fill = NA, colour = 'black'),
                                        legend.background = ggplot2::element_blank(),
                                        axis.line = ggplot2::element_line(colour = 'black'),
                                        axis.text.x = ggplot2::element_text(size = 10, colour = 'black'),
                                        axis.text.y = ggplot2::element_text(size = 10, angle = 90,
                                                                            colour = 'black'),
                                        axis.title.x = ggplot2::element_text(size = 12, colour = 'black'),
                                        axis.title.y = ggplot2::element_text(size = 12, colour = 'black'),
                                        axis.ticks = ggplot2::element_line(colour = 'black'))
                     }
             )



             return(result)
           }
)


