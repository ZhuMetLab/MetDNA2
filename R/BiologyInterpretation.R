################################################################################
# intepretBiology ------------------------------------------------------------

#' @title intepretBiology
#' @author Zhiwei Zhou
#' \email{zhouzw@@sioc.ac.cn}
#' @description a wrapper for MetDNA2
#' @param ms1_data
#' @param result_annotation
#' @param sample_info
#' @param path '.'
#' @param comp_group the group for comparison, which needed included in sample information
#' @param uni_test univariate test, including "t","wilcox","anova". Default: 't'
#' @param correct_p whether correct P value
#' @param species "hsa", "mmu", "rno", "bta", "gga", "dre", "dme", "cel", "sce", "osa", "ath", "smm", "pfa", "tbr", "eco", "bsu", "ppu", "sau", "tma", "syf", "mln"
#' @param quanti_pathway_method the mehtod of pathway quantification, including 'mean', 'sum', 'median'. Default: 'mean'
#' @param extension_step c('0', '1', '2', '3', '4', '5', '6', '7', '8')
#' @export
#' @examples
#'

# sample_file = 'data.csv'
# sample_info_file = 'sample.info.csv'
# table_identification_pair_file = 'table_identification_pair.csv'
# path = '/home/zhouzw/Data_processing/20210522_debug/test'
# metdna_version = 'version2'
# comp_group = c("W30", "W03")
# uni_test = "t"
# correct_p = TRUE
# p_cutoff = 0.05
# fc_cutoff = 1.25
# species = 'dme'
# quanti_pathway_method = 'mean'
# organze_basis = 'kegg_id'
# extension_step = '0'
#
# interpretBiology(sample_file_pos = 'data.csv',
#                  sample_info_file_pos = 'sample.info.csv',
#                  table_identification_pair_file_pos = 'table3_identification_pair.csv',
#                  sample_file_neg = 'data.csv',
#                  sample_info_file_neg = 'sample.info.csv',
#                  table_identification_pair_file_neg = 'table3_identification_pair.csv',
#                  path = '/home/zhouzw/Data_processing/20210522_debug/web_bug/907bd62c2668fe3e2a035ca9ab12a1c3/',
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



setGeneric(name = 'interpretBiology',
           def = function(
             sample_file_pos = 'data.csv',
             sample_info_file_pos = 'sample.info.csv',
             table_identification_pair_file_pos = 'table3_identification_pair.csv',
             sample_file_neg = 'data.csv',
             sample_info_file_neg = 'sample.info.csv',
             table_identification_pair_file_neg = 'table3_identification_pair.csv',
             path = '.',
             polarity = c('positive', 'negative', 'both'),
             metdna_version = c('version2', 'version1'),
             comp_group = c("W30", "W03"),
             uni_test = c("t","wilcox","anova"),
             correct_p = TRUE,
             p_cutoff = 0.05,
             fc_cutoff = 1.25,
             species = c("hsa", "mmu", "rno", "bta", "gga",
                         "dre", "dme", "cel", "sce", "osa",
                         "ath", "smm", "pfa", "tbr", "eco",
                         "bsu", "ppu", "sau", "tma", "syf", "mln"),
             quanti_pathway_method = c('mean', 'sum', 'median'),
             organze_basis = c('kegg_id', 'all'),
             extension_step = c('0', '1', '2', '3', '4', '5', '6', '7', '8'),
             ...
           ){
             library(MetBioInterpretation)

             message('Start biology intepretation ...\n')
             uni_test <- match.arg(uni_test)
             species <- match.arg(species)
             quanti_pathway_method <- match.arg(quanti_pathway_method)
             metdna_version <- match.arg(metdna_version)
             organze_basis <- match.arg(organze_basis)

             # browser()

             # positive mode only ------------------------
             if (polarity == 'positive') {
               # browser()
               sample_file <- sample_file_pos
               sample_info_file <- sample_info_file_pos
               table_identification_pair_file <- table_identification_pair_file_pos

               # copy files into the working directory
               if (!file.exists(file.path(path, sample_file))) {
                 stop('The file ', sample_file, ' NOT found in the directory, please check it!\n')
               }
               if (!file.exists(file.path(path, sample_info_file))) {
                 stop('The file ', sample_info_file, ' NOT found in the directory, please check it!\n')
               }
               if (!file.exists(file.path(path, '00_annotation_table', table_identification_pair_file))) {
                 stop('The file ', table_identification_pair_file, ' NOT found in the directory, please check it!\n')
               }



               path_output <- file.path(path, '04_biology_intepretation')
               dir.create(path_output, showWarnings = FALSE, recursive = TRUE)
               file.copy(from = file.path(path, sample_file),
                         to = file.path(path_output),
                         copy.date = TRUE,
                         overwrite = TRUE)
               file.copy(from = file.path(path, sample_info_file),
                         to = file.path(path_output),
                         copy.date = TRUE,
                         overwrite = TRUE)
               file.copy(from = file.path(path, '00_annotation_table', table_identification_pair_file),
                         to = file.path(path_output),
                         copy.date = TRUE,
                         overwrite = TRUE)


               options(readr.num_columns = FALSE)
               sample <- readr::read_csv(file.path(path_output, sample_file))
               sample_info <- readr::read_csv(file.path(path_output, sample_info_file), col_types = 'cc')
               table_identification_pair <- readr::read_csv(file.path(path_output, table_identification_pair_file))

               if (all(c('name', 'mz', 'rt', 'ccs') %in% colnames(sample))) {
                 sample <- sample %>%
                   tibble::column_to_rownames(var = 'name') %>%
                   dplyr::select(-c('mz', 'rt', 'ccs'))
               } else {
                 sample <- sample %>%
                   tibble::column_to_rownames(var = 'name') %>%
                   dplyr::select(-c('mz', 'rt'))
               }

               temp_error <- sample_info %>%
                 dplyr::filter(group %in% comp_group) %>%
                 dplyr::count(group)
               if (any(temp_error$n < 3)) {
                 cat('Note: one group has <3 samples, the biological sample interpreation will be skip\n')
                 return(NULL)
               }

               cat('Start univariate analsis ...\n')
               result_stat <- MetBioInterpretation::doStatAnalysis(sample = sample,
                                                                   sample_info = sample_info,
                                                                   comp_group = comp_group,
                                                                   uni_test = uni_test,
                                                                   correct = correct_p,
                                                                   p_cutoff = p_cutoff,
                                                                   fc_cutoff = fc_cutoff,
                                                                   by_what = 'mean')

               dir.create(file.path(path_output, '00_intermediate_data'),
                          showWarnings = FALSE, recursive = TRUE)
               save(result_stat,
                    file = file.path(path_output, '00_intermediate_data', 'result_stat.RData'),
                    version = 2)
               readr::write_csv(result_stat,
                                path = file.path(path_output, 'statistics_result.csv'))

               cat('Plot vocano plot...\n')
               temp_plot <- MetBioInterpretation::plotVolcano(raw_data = result_stat,
                                                              p_cutoff = p_cutoff,
                                                              fc_cutoff = fc_cutoff,
                                                              log_tran = 'log2')

               ggplot2::ggsave(temp_plot,
                               filename = file.path(path_output, 'volcano_plot.pdf'),
                               width = 6, height = 6)


               cat('\n\n');cat('Perform pathway enrichment analysis ...\n')


               if (metdna_version == 'version1') {
                 # extract significant feature
                 feature_significant <- result_stat %>%
                   dplyr::filter(p_values <= p_cutoff & ((fold_changes >= fc_cutoff) | (fold_changes <= 1/fc_cutoff))) %>%
                   dplyr::pull(name)

                 load(file.path(path, '02_result_MRN_annotation/00_intermediate_data/tags2_after_redundancy_remove'))

                 names(tags2_after_redundancy_remove) <- sapply(tags2_after_redundancy_remove, function(x){
                   x@name
                 })

                 idx <- match(feature_significant, names(tags2_after_redundancy_remove))
                 temp_result <- convertTags2ResultMatrix(tags2 = tags2_after_redundancy_remove[idx], candidate.num = 5)
                 met_sig <- temp_result$to %>% stringr::str_split(pattern = ';') %>% unlist()
                 met_sig <- (!is.na(met_sig)) %>% which() %>% met_sig[.] %>% unique()

               } else {
                 # extract significant feature
                 feature_significant <- result_stat %>%
                   dplyr::filter(p_values <= p_cutoff & ((fold_changes >= fc_cutoff) | (fold_changes <= 1/fc_cutoff))) %>%
                   dplyr::pull(name)

                 result_annotation_sig <- table_identification_pair %>%
                   dplyr::mutate(temp_peak_name = stringr::str_replace(peak_name, '_[a-z]+', '')) %>%
                   dplyr::filter(temp_peak_name %in% feature_significant) %>%
                   dplyr::select(-temp_peak_name)

                 if (organze_basis == 'kegg_id') {
                   met_sig <- result_annotation_sig %>%
                     dplyr::pull(id_kegg) %>%
                     unique() %>%
                     .[!is.na(.)]
                 } else {
                   met_sig <- result_annotation_sig %>%
                     tidyr::separate_rows(stereo_isomer_id) %>%
                     dplyr::pull(stereo_isomer_id) %>%
                     unique() %>%
                     .[!is.na(.)]
                 }
               }

               if (length(met_sig) > 0) {
                 result_pathway_enrichment <- MetBioInterpretation::enrichPathway(met_sig = met_sig,
                                                                                  species = species,
                                                                                  lib_pathway = lib_pathway,
                                                                                  test_method = 'hypergeometric')

                 if (length(result_pathway_enrichment) > 0) {
                   result_pathway_enrichment <- result_pathway_enrichment %>% tibble::rownames_to_column(var = 'pathway')

                   dir.create(file.path(path_output, '01_pathway_enrichment'),
                              showWarnings = FALSE, recursive = TRUE)
                   save(result_pathway_enrichment,
                        file = file.path(path_output, '00_intermediate_data', 'result_pathway_enrichment.RData'),
                        version = 2)
                   readr::write_csv(result_pathway_enrichment,
                                    path = file.path(path_output,
                                                     '01_pathway_enrichment',
                                                     'pathway_enrichment_analysis.csv'))

                   temp_plot <- MetBioInterpretation::plotPathwayScatter(result_pathway_enrichment = result_pathway_enrichment, p_cutoff = 0.05, p_adjust = FALSE)
                   ggplot2::ggsave(temp_plot,
                                   filename = file.path(path_output,
                                                        '01_pathway_enrichment',
                                                        'pathway_enrichment_overview.pdf'),
                                   width = 6, height = 6)

                   temp_plot <- MetBioInterpretation::plotPathwayBar(result_pathway_enrichment = result_pathway_enrichment, p_cutoff = 0.05, p_adjust = FALSE)
                   ggplot2::ggsave(temp_plot,
                                   filename = file.path(path_output,
                                                        '01_pathway_enrichment',
                                                        'pathway_enrichment_MSEA.pdf'),
                                   width = 6, height = 6)
                 }
               } else {
                 cat('Note: No significant metabolites, skip remain analysis, including pathway enrichment, pathway quantative analysis, and class enrichment\n')
                 return(NULL)
               }

               cat('\n\n');cat('Perform pathway quantitative analysis ...\n')
               dir.create(file.path(path_output, '02_pathway_quantitative_analysis'),
                          showWarnings = FALSE, recursive = TRUE)

               if (metdna_version == 'version1') {
                 # extract significant feature from MRN intermediate data
                 load(file.path(path, '02_result_MRN_annotation/00_intermediate_data/id_result_redun_rm'))

                 id_result_redun_rm <- lapply(unique(id_result_redun_rm$name), function(x){
                   temp_idx <- which(id_result_redun_rm$name == x)
                   temp_data <- id_result_redun_rm[temp_idx,,drop = FALSE]
                   if(nrow(temp_data) > 5) temp_data <- temp_data[1:5,]
                   temp_data
                 })

                 id_result_redun_rm <- do.call(rbind, id_result_redun_rm) %>% tibble::as_tibble()
                 id_result_redun_rm <- id_result_redun_rm %>%
                   dplyr::filter(isotope == '[M]') %>%
                   dplyr::filter(name %in% feature_significant)

                 temp_idx <- match(id_result_redun_rm$name, rownames(sample))
                 id_result_redun_rm <- id_result_redun_rm %>%
                   dplyr::select(name, mz, rt, to, Confidence) %>%
                   dplyr::bind_cols(sample[temp_idx,]) %>%
                   dplyr::rename(peak_name = name,
                                 id_kegg = to,
                                 confidence_level = Confidence)

                 data("cpd_emrn", package = 'MetDNA2')
                 data("zhuMetlib", package = 'MetDNA2')
                 temp_idx2 <- match(id_result_redun_rm$id_kegg, cpd_emrn$id)
                 id_result_redun_rm <- id_result_redun_rm %>%
                   dplyr::mutate(name = cpd_emrn$name[temp_idx2],
                                 inchikey = cpd_emrn$inchikey[temp_idx2],
                                 stereo_isomer_id = cpd_emrn$id_kegg_synonyms[temp_idx2],
                                 stereo_isomer_name = cpd_emrn$name_synonyms[temp_idx2]) %>%
                   dplyr::select(peak_name:confidence_level, name, inchikey, stereo_isomer_id, stereo_isomer_name, dplyr::everything())

                temp_idx3 <- id_result_redun_rm$name %>% is.na() %>% which()
                temp_idx4 <- match(id_result_redun_rm$id_kegg[temp_idx3], zhuMetlib$meta$compound$id_kegg)
                id_result_redun_rm$name[temp_idx3] <- zhuMetlib$meta$compound$name[temp_idx4]
                id_result_redun_rm$inchikey[temp_idx3] <- zhuMetlib$meta$compound$inchikey[temp_idx4]
                id_result_redun_rm$peak_name <- id_result_redun_rm$peak_name %>% rename4UniqueIndex()
                readr::write_csv(id_result_redun_rm, file = file.path(path_output, "table3_identification_pair_v1.csv"))

                MetBioInterpretation::selectQuantiMetabolite(table_identification_pair = "table3_identification_pair_v1.csv",
                                                             sample_info = sample_info_file,
                                                             comp_group = comp_group,
                                                             path = path_output,
                                                             is_scale = TRUE,
                                                             scale_method = "pareto",
                                                             represent_met = 'median',
                                                             organze_basis = organze_basis)


               } else {
                 options(readr.num_columns = FALSE)
                 # table_identification_pair <- readr::read_csv(file.path(path_output, table_identification_pair_file))
                 temp_result <- table_identification_pair %>%
                   dplyr::mutate(temp_peak_name = stringr::str_replace(peak_name, '_[a-z]+', '')) %>%
                   dplyr::filter(temp_peak_name %in% feature_significant) %>%
                   dplyr::select(-temp_peak_name)
                 readr::write_csv(temp_result, file = file.path(path_output, "table3_identification_pair_v2.csv"))

                 MetBioInterpretation::selectQuantiMetabolite(table_identification_pair = "table3_identification_pair_v2.csv",
                                                              sample_info = sample_info_file,
                                                              comp_group = comp_group,
                                                              path = path_output,
                                                              is_scale = TRUE,
                                                              scale_method = "pareto",
                                                              represent_met = 'median',
                                                              organze_basis = organze_basis)
               }


               temp <- try(MetBioInterpretation::quantiPathway(quanti_table = 'pathway_metabolite_quantitative_result.csv',
                                                               sample_info = sample_info_file,
                                                               comp_group = comp_group,
                                                               is_scale = FALSE, # peak area has been scaled in selectQuantiMetabolite
                                                               correct_p = FALSE,
                                                               uni_test = uni_test,
                                                               path = path_output,
                                                               by_what = 'mean',
                                                               species = species,
                                                               is_consider_stereo_isomer = ifelse(organze_basis == 'kegg_id', FALSE, TRUE)),
                           silent = TRUE)

               if (class(temp) == 'try-error') {
                 cat('Note: quantiPathway is error. Pay attention, please!\n')
                 return(NULL)
               }

               cat('\n\n');cat('Export pathway quantitative analysis results...\n')
               file.copy(from = file.path(path_output, 'pathway_quantitative_result.csv'),
                         to = file.path(path_output, '02_pathway_quantitative_analysis'),
                         overwrite = TRUE, copy.date = TRUE)
               file.copy(from = file.path(path_output, 'pathway_metabolite_quantitative_result.csv'),
                         to = file.path(path_output, '02_pathway_quantitative_analysis'),
                         overwrite = TRUE, copy.date = TRUE)

               options(readr.num_columns = 0)
               quanti_pathway_data <- readr::read_csv(file.path(path_output,
                                                                '02_pathway_quantitative_analysis',
                                                                'pathway_quantitative_result.csv'))
               quanti_met_data <- readr::read_csv(file.path(path_output,
                                                            '02_pathway_quantitative_analysis',
                                                            'pathway_metabolite_quantitative_result.csv'))

               temp_plot <- MetBioInterpretation::plotPathwayHeatmap(quanti_pathway_data = quanti_pathway_data,
                                                                     sample_info = sample_info,
                                                                     comp_group = comp_group)
               ggplot2::ggsave(temp_plot,
                               filename = file.path(path_output,
                                                    '02_pathway_quantitative_analysis',
                                                    'pathway_heatmap.pdf'))

               cat('\n');cat('Plot boxplot for each pathway ...\n')

               MetBioInterpretation::plotPathwayBoxplot(quanti_pathway_data = quanti_pathway_data,
                                                        sample_info = sample_info,
                                                        comp_group = comp_group,
                                                        path = file.path(path_output, '02_pathway_quantitative_analysis'))

               cat('\n');cat('Plot metabolite heatmap for each pathway ...\n')

               MetBioInterpretation::plotPathwayHeatmapMet(quanti_met_data = quanti_met_data,
                                                           quanti_pathway_data = quanti_pathway_data,
                                                           sample_info = sample_info,
                                                           comp_group = comp_group,
                                                           species = species,
                                                           is_consider_stereo_isomer = ifelse(organze_basis == 'kegg_id', FALSE, TRUE),
                                                           path = file.path(path_output, '02_pathway_quantitative_analysis'))


               # cat('\n');cat('Pathway analysis has done ...\n')
               cat('\n\n');cat('Enrich chemical class analysis...\n')
               dir.create(file.path(path_output, '03_class_enrichment'),
                          showWarnings = FALSE, recursive = TRUE)
               table_class_enrich <- MetBioInterpretation::generateClassEnrichTable(table_identification_pair = table_identification_pair,
                                                                                    result_stat = result_stat)

               purrr::walk(c('superclass', 'class', 'subclass'), function(x){
                 # cat(x, ' ')
                 switch(x,
                        'superclass' = {
                          temp <- table_class_enrich %>%
                            dplyr::select(peak_name, p_values, fold_changes, superclass) %>%
                            dplyr::rename('compound_name' = peak_name,
                                          'pvalue' = p_values,
                                          'foldchange' = fold_changes,
                                          'class' = superclass) %>%
                            dplyr::filter(!is.na(class)) %>%
                            dplyr::arrange(pvalue) %>%
                            dplyr::mutate(order = seq(dplyr::n())) %>%
                            dplyr::select(compound_name, order, dplyr::everything())
                        },
                        'class' = {
                          temp <- table_class_enrich %>%
                            dplyr::select(peak_name, p_values, fold_changes, class) %>%
                            dplyr::rename('compound_name' = peak_name,
                                          'pvalue' = p_values,
                                          'foldchange' = fold_changes,
                                          'class' = class) %>%
                            dplyr::filter(!is.na(class)) %>%
                            dplyr::arrange(pvalue) %>%
                            dplyr::mutate(order = seq(dplyr::n())) %>%
                            dplyr::select(compound_name, order, dplyr::everything())
                        },
                        'subclass' = {
                          temp <- table_class_enrich %>%
                            dplyr::select(peak_name, p_values, fold_changes, subclass) %>%
                            dplyr::rename('compound_name' = peak_name,
                                          'pvalue' = p_values,
                                          'foldchange' = fold_changes,
                                          'class' = subclass) %>%
                            dplyr::filter(!is.na(class)) %>%
                            dplyr::arrange(pvalue) %>%
                            dplyr::mutate(order = seq(dplyr::n())) %>%
                            dplyr::select(compound_name, order, dplyr::everything())
                        })

                 readr::write_csv(temp,
                                  file = file.path(path_output, '03_class_enrichment', paste0(x, '_table_class_enrichment.csv')))

                 enrichment_result <- try(MetBioInterpretation::enrichClass(table_class_file = paste0(x, '_table_class_enrichment.csv'),
                                                                            path = file.path(path_output, '03_class_enrichment')),
                                          silent = TRUE)

                 if (class(enrichment_result) == 'try-error') {
                   cat('Note: class enrichment has some error. Pay attention, please!\n')
                   return(NULL)
                 }

                 readr::write_csv(enrichment_result,
                                  file = file.path(path_output, '03_class_enrichment', paste0(x, '_result_class_enrichment.csv')))

                 if (sum(enrichment_result$pvalues <= 0.05) > 0) {
                   try(expr = {
                     temp_plot <- MetBioInterpretation::plotClassScater(enrichment_result = enrichment_result)
                     dir.create(file.path(path_output, '03_class_enrichment', 'plot'), showWarnings = FALSE, recursive = TRUE)
                     ggplot2::ggsave(temp_plot,
                                     filename = file.path(path_output, '03_class_enrichment', 'plot', paste0(x, '_class_enrichment.pdf')),
                                     width = 6, height = 6)
                   }, silent = TRUE)
                 } else {
                   cat('Note: No class enriched, skip this step!\n')
                 }

               })
             }

             # negative mode only ------------------------
             if (polarity == 'negative') {
               # browser()
               sample_file <- sample_file_neg
               sample_info_file <- sample_info_file_neg
               table_identification_pair_file <- table_identification_pair_file_neg

               # copy files into the working directory
               if (!file.exists(file.path(path, sample_file))) {
                 stop('The file ', sample_file, ' NOT found in the directory, please check it!\n')
               }
               if (!file.exists(file.path(path, sample_info_file))) {
                 stop('The file ', sample_info_file, ' NOT found in the directory, please check it!\n')
               }
               if (!file.exists(file.path(path, '00_annotation_table', table_identification_pair_file))) {
                 stop('The file ', table_identification_pair_file, ' NOT found in the directory, please check it!\n')
               }

               path_output <- file.path(path, '04_biology_intepretation')
               dir.create(path_output, showWarnings = FALSE, recursive = TRUE)
               file.copy(from = file.path(path, sample_file),
                         to = file.path(path_output),
                         copy.date = TRUE,
                         overwrite = TRUE)
               file.copy(from = file.path(path, sample_info_file),
                         to = file.path(path_output),
                         copy.date = TRUE,
                         overwrite = TRUE)
               file.copy(from = file.path(path, '00_annotation_table', table_identification_pair_file),
                         to = file.path(path_output),
                         copy.date = TRUE,
                         overwrite = TRUE)

               options(readr.num_columns = FALSE)
               sample <- readr::read_csv(file.path(path_output, sample_file))
               sample_info <- readr::read_csv(file.path(path_output, sample_info_file), col_types = c('cc'))
               table_identification_pair <- readr::read_csv(file.path(path_output, table_identification_pair_file))

               if (all(c('name', 'mz', 'rt', 'ccs') %in% colnames(sample))) {
                 sample <- sample %>%
                   tibble::column_to_rownames(var = 'name') %>%
                   dplyr::select(-c('mz', 'rt', 'ccs'))
               } else {
                 sample <- sample %>%
                   tibble::column_to_rownames(var = 'name') %>%
                   dplyr::select(-c('mz', 'rt'))
               }

               temp_error <- sample_info %>%
                 dplyr::filter(group %in% comp_group) %>%
                 dplyr::count(group)
               if (any(temp_error$n < 3)) {
                 cat('Note: one group has <3 samples, the biological sample interpreation will be skip\n')
                 return(NULL)
               }

               cat('Start univariate analsis ...\n')
               result_stat <- MetBioInterpretation::doStatAnalysis(sample = sample,
                                                                   sample_info = sample_info,
                                                                   comp_group = comp_group,
                                                                   uni_test = uni_test,
                                                                   correct = correct_p,
                                                                   p_cutoff = p_cutoff,
                                                                   fc_cutoff = fc_cutoff,
                                                                   by_what = 'mean')

               dir.create(file.path(path_output, '00_intermediate_data'),
                          showWarnings = FALSE, recursive = TRUE)
               save(result_stat,
                    file = file.path(path_output, '00_intermediate_data', 'result_stat.RData'),
                    version = 2)
               readr::write_csv(result_stat,
                                path = file.path(path_output, 'statistics_result.csv'))

               cat('Plot vocano plot...\n')
               temp_plot <- MetBioInterpretation::plotVolcano(raw_data = result_stat,
                                                              p_cutoff = p_cutoff,
                                                              fc_cutoff = fc_cutoff,
                                                              log_tran = 'log2')

               ggplot2::ggsave(temp_plot,
                               filename = file.path(path_output, 'volcano_plot.pdf'),
                               width = 6, height = 6)


               cat('\n\n');cat('Perform pathway enrichment analysis ...\n')


               if (metdna_version == 'version1') {
                 # extract significant feature
                 feature_significant <- result_stat %>%
                   dplyr::filter(p_values <= p_cutoff & ((fold_changes >= fc_cutoff) | (fold_changes <= 1/fc_cutoff))) %>%
                   dplyr::pull(name)

                 load(file.path(path, '02_result_MRN_annotation/00_intermediate_data/tags2_after_redundancy_remove'))

                 names(tags2_after_redundancy_remove) <- sapply(tags2_after_redundancy_remove, function(x){
                   x@name
                 })

                 idx <- match(feature_significant, names(tags2_after_redundancy_remove))
                 temp_result <- convertTags2ResultMatrix(tags2 = tags2_after_redundancy_remove[idx], candidate.num = 5)
                 met_sig <- temp_result$to %>% stringr::str_split(pattern = ';') %>% unlist()
                 met_sig <- (!is.na(met_sig)) %>% which() %>% met_sig[.] %>% unique()

               } else {
                 # extract significant feature
                 feature_significant <- result_stat %>%
                   dplyr::filter(p_values <= p_cutoff & ((fold_changes >= fc_cutoff) | (fold_changes <= 1/fc_cutoff))) %>%
                   dplyr::pull(name)

                 result_annotation_sig <- table_identification_pair %>%
                   dplyr::mutate(temp_peak_name = stringr::str_replace(peak_name, '_[a-z]+', '')) %>%
                   dplyr::filter(temp_peak_name %in% feature_significant) %>%
                   dplyr::select(-temp_peak_name)

                 if (organze_basis == 'kegg_id') {
                   met_sig <- result_annotation_sig %>%
                     dplyr::pull(id_kegg) %>%
                     unique() %>%
                     .[!is.na(.)]
                 } else {
                   met_sig <- result_annotation_sig %>%
                     tidyr::separate_rows(stereo_isomer_id) %>%
                     dplyr::pull(stereo_isomer_id) %>%
                     unique() %>%
                     .[!is.na(.)]
                 }
               }

               if (length(met_sig) > 0) {
                 result_pathway_enrichment <- MetBioInterpretation::enrichPathway(met_sig = met_sig,
                                                                                  species = species,
                                                                                  lib_pathway = lib_pathway,
                                                                                  test_method = 'hypergeometric')

                 if (length(result_pathway_enrichment) > 0) {
                   result_pathway_enrichment <- result_pathway_enrichment %>% tibble::rownames_to_column(var = 'pathway')

                   dir.create(file.path(path_output, '01_pathway_enrichment'),
                              showWarnings = FALSE, recursive = TRUE)
                   save(result_pathway_enrichment,
                        file = file.path(path_output, '00_intermediate_data', 'result_pathway_enrichment.RData'),
                        version = 2)
                   readr::write_csv(result_pathway_enrichment,
                                    path = file.path(path_output,
                                                     '01_pathway_enrichment',
                                                     'pathway_enrichment_analysis.csv'))

                   temp_plot <- MetBioInterpretation::plotPathwayScatter(result_pathway_enrichment = result_pathway_enrichment, p_cutoff = 0.05, p_adjust = FALSE)
                   ggplot2::ggsave(temp_plot,
                                   filename = file.path(path_output,
                                                        '01_pathway_enrichment',
                                                        'pathway_enrichment_overview.pdf'),
                                   width = 6, height = 6)

                   temp_plot <- MetBioInterpretation::plotPathwayBar(result_pathway_enrichment = result_pathway_enrichment, p_cutoff = 0.05, p_adjust = FALSE)
                   ggplot2::ggsave(temp_plot,
                                   filename = file.path(path_output,
                                                        '01_pathway_enrichment',
                                                        'pathway_enrichment_MSEA.pdf'),
                                   width = 6, height = 6)
                 }
               } else {
                 cat('Note: No significant metabolites, skip remain analysis, including pathway enrichment, pathway quantative analysis, and class enrichment\n')
                 return(NULL)
               }

               cat('\n\n');cat('Perform pathway quantitative analysis ...\n')
               dir.create(file.path(path_output, '02_pathway_quantitative_analysis'),
                          showWarnings = FALSE, recursive = TRUE)

               if (metdna_version == 'version1') {
                 # extract significant feature from MRN intermediate data
                 load(file.path(path, '02_result_MRN_annotation/00_intermediate_data/id_result_redun_rm'))

                 id_result_redun_rm <- lapply(unique(id_result_redun_rm$name), function(x){
                   temp_idx <- which(id_result_redun_rm$name == x)
                   temp_data <- id_result_redun_rm[temp_idx,,drop = FALSE]
                   if(nrow(temp_data) > 5) temp_data <- temp_data[1:5,]
                   temp_data
                 })

                 id_result_redun_rm <- do.call(rbind, id_result_redun_rm) %>% tibble::as_tibble()
                 id_result_redun_rm <- id_result_redun_rm %>%
                   dplyr::filter(isotope == '[M]') %>%
                   dplyr::filter(name %in% feature_significant)

                 temp_idx <- match(id_result_redun_rm$name, rownames(sample))
                 id_result_redun_rm <- id_result_redun_rm %>%
                   dplyr::select(name, mz, rt, to, Confidence) %>%
                   dplyr::bind_cols(sample[temp_idx,]) %>%
                   dplyr::rename(peak_name = name,
                                 id_kegg = to,
                                 confidence_level = Confidence)

                 data("cpd_emrn", package = 'MetDNA2')
                 data("zhuMetlib", package = 'MetDNA2')
                 temp_idx2 <- match(id_result_redun_rm$id_kegg, cpd_emrn$id)
                 id_result_redun_rm <- id_result_redun_rm %>%
                   dplyr::mutate(name = cpd_emrn$name[temp_idx2],
                                 inchikey = cpd_emrn$inchikey[temp_idx2],
                                 stereo_isomer_id = cpd_emrn$id_kegg_synonyms[temp_idx2],
                                 stereo_isomer_name = cpd_emrn$name_synonyms[temp_idx2]) %>%
                   dplyr::select(peak_name:confidence_level, name, inchikey, stereo_isomer_id, stereo_isomer_name, dplyr::everything())

                 temp_idx3 <- id_result_redun_rm$name %>% is.na() %>% which()
                 temp_idx4 <- match(id_result_redun_rm$id_kegg[temp_idx3], zhuMetlib$meta$compound$id_kegg)
                 id_result_redun_rm$name[temp_idx3] <- zhuMetlib$meta$compound$name[temp_idx4]
                 id_result_redun_rm$inchikey[temp_idx3] <- zhuMetlib$meta$compound$inchikey[temp_idx4]
                 id_result_redun_rm$peak_name <- id_result_redun_rm$peak_name %>% rename4UniqueIndex()
                 readr::write_csv(id_result_redun_rm, file = file.path(path_output, "table3_identification_pair_v1.csv"))

                 MetBioInterpretation::selectQuantiMetabolite(table_identification_pair = "table3_identification_pair_v1.csv",
                                                              sample_info = sample_info_file,
                                                              comp_group = comp_group,
                                                              path = path_output,
                                                              is_scale = TRUE,
                                                              scale_method = "pareto",
                                                              represent_met = 'median',
                                                              organze_basis = organze_basis)


               } else {
                 options(readr.num_columns = FALSE)
                 # table_identification_pair <- readr::read_csv(file.path(path_output, table_identification_pair_file))
                 temp_result <- table_identification_pair %>%
                   dplyr::mutate(temp_peak_name = stringr::str_replace(peak_name, '_[a-z]+', '')) %>%
                   dplyr::filter(temp_peak_name %in% feature_significant) %>%
                   dplyr::select(-temp_peak_name)
                 readr::write_csv(temp_result, file = file.path(path_output, "table3_identification_pair_v2.csv"))

                 MetBioInterpretation::selectQuantiMetabolite(table_identification_pair = "table3_identification_pair_v2.csv",
                                                              sample_info = sample_info_file,
                                                              comp_group = comp_group,
                                                              path = path_output,
                                                              is_scale = TRUE,
                                                              scale_method = "pareto",
                                                              represent_met = 'median',
                                                              organze_basis = organze_basis)
               }


               temp <- try(MetBioInterpretation::quantiPathway(quanti_table = 'pathway_metabolite_quantitative_result.csv',
                                                               sample_info = sample_info_file,
                                                               comp_group = comp_group,
                                                               is_scale = FALSE, # peak area has been scaled in selectQuantiMetabolite
                                                               correct_p = FALSE,
                                                               uni_test = uni_test,
                                                               path = path_output,
                                                               by_what = 'mean',
                                                               species = species,
                                                               is_consider_stereo_isomer = ifelse(organze_basis == 'kegg_id', FALSE, TRUE)),
                           silent = TRUE)

               if (class(temp) == 'try-error') {
                 cat('Note: quantiPathway is error. Pay attention, please!\n')
                 return(NULL)
               }

               cat('\n\n');cat('Export pathway quantitative analysis results...\n')
               file.copy(from = file.path(path_output, 'pathway_quantitative_result.csv'),
                         to = file.path(path_output, '02_pathway_quantitative_analysis'),
                         overwrite = TRUE, copy.date = TRUE)
               file.copy(from = file.path(path_output, 'pathway_metabolite_quantitative_result.csv'),
                         to = file.path(path_output, '02_pathway_quantitative_analysis'),
                         overwrite = TRUE, copy.date = TRUE)

               options(readr.num_columns = 0)
               quanti_pathway_data <- readr::read_csv(file.path(path_output,
                                                                '02_pathway_quantitative_analysis',
                                                                'pathway_quantitative_result.csv'))
               quanti_met_data <- readr::read_csv(file.path(path_output,
                                                            '02_pathway_quantitative_analysis',
                                                            'pathway_metabolite_quantitative_result.csv'))

               temp_plot <- MetBioInterpretation::plotPathwayHeatmap(quanti_pathway_data = quanti_pathway_data,
                                                                     sample_info = sample_info,
                                                                     comp_group = comp_group)
               ggplot2::ggsave(temp_plot,
                               filename = file.path(path_output,
                                                    '02_pathway_quantitative_analysis',
                                                    'pathway_heatmap.pdf'))

               cat('\n');cat('Plot boxplot for each pathway ...\n')

               MetBioInterpretation::plotPathwayBoxplot(quanti_pathway_data = quanti_pathway_data,
                                                        sample_info = sample_info,
                                                        comp_group = comp_group,
                                                        path = file.path(path_output, '02_pathway_quantitative_analysis'))

               cat('\n');cat('Plot metabolite heatmap for each pathway ...\n')

               MetBioInterpretation::plotPathwayHeatmapMet(quanti_met_data = quanti_met_data,
                                                           quanti_pathway_data = quanti_pathway_data,
                                                           sample_info = sample_info,
                                                           comp_group = comp_group,
                                                           species = species,
                                                           is_consider_stereo_isomer = ifelse(organze_basis == 'kegg_id', FALSE, TRUE),
                                                           path = file.path(path_output, '02_pathway_quantitative_analysis'))


               # cat('\n');cat('Pathway analysis has done ...\n')
               cat('\n\n');cat('Enrich chemical class analysis...\n')
               dir.create(file.path(path_output, '03_class_enrichment'),
                          showWarnings = FALSE, recursive = TRUE)
               table_class_enrich <- MetBioInterpretation::generateClassEnrichTable(table_identification_pair = table_identification_pair,
                                                                                    result_stat = result_stat)

               purrr::walk(c('superclass', 'class', 'subclass'), function(x){
                 switch(x,
                        'superclass' = {
                          temp <- table_class_enrich %>%
                            dplyr::select(peak_name, p_values, fold_changes, superclass) %>%
                            dplyr::rename('compound_name' = peak_name,
                                          'pvalue' = p_values,
                                          'foldchange' = fold_changes,
                                          'class' = superclass) %>%
                            dplyr::filter(!is.na(class)) %>%
                            dplyr::arrange(pvalue) %>%
                            dplyr::mutate(order = seq(dplyr::n())) %>%
                            dplyr::select(compound_name, order, dplyr::everything())
                        },
                        'class' = {
                          temp <- table_class_enrich %>%
                            dplyr::select(peak_name, p_values, fold_changes, class) %>%
                            dplyr::rename('compound_name' = peak_name,
                                          'pvalue' = p_values,
                                          'foldchange' = fold_changes,
                                          'class' = class) %>%
                            dplyr::filter(!is.na(class)) %>%
                            dplyr::arrange(pvalue) %>%
                            dplyr::mutate(order = seq(dplyr::n())) %>%
                            dplyr::select(compound_name, order, dplyr::everything())
                        },
                        'subclass' = {
                          temp <- table_class_enrich %>%
                            dplyr::select(peak_name, p_values, fold_changes, subclass) %>%
                            dplyr::rename('compound_name' = peak_name,
                                          'pvalue' = p_values,
                                          'foldchange' = fold_changes,
                                          'class' = subclass) %>%
                            dplyr::filter(!is.na(class)) %>%
                            dplyr::arrange(pvalue) %>%
                            dplyr::mutate(order = seq(dplyr::n())) %>%
                            dplyr::select(compound_name, order, dplyr::everything())
                        })

                 readr::write_csv(temp,
                                  file = file.path(path_output, '03_class_enrichment', paste0(x, '_table_class_enrichment.csv')))

                 enrichment_result <- try(MetBioInterpretation::enrichClass(table_class_file = paste0(x, '_table_class_enrichment.csv'),
                                                                            path = file.path(path_output, '03_class_enrichment')),
                                          silent = TRUE)

                 if (class(enrichment_result) == 'try-error') {
                   cat('Note: class enrichment has some error. Pay attention, please!\n')
                   return(NULL)
                 }

                 readr::write_csv(enrichment_result,
                                  file = file.path(path_output, '03_class_enrichment', paste0(x, '_result_class_enrichment.csv')))

                 temp_plot <- MetBioInterpretation::plotClassScater(enrichment_result = enrichment_result)
                 dir.create(file.path(path_output, '03_class_enrichment', 'plot'), showWarnings = FALSE, recursive = TRUE)
                 ggplot2::ggsave(temp_plot,
                                 filename = file.path(path_output, '03_class_enrichment', 'plot', paste0(x, '_class_enrichment.pdf')),
                                 width = 6, height = 6)

               })
             }

             # both polarity -----------------------------
             if (polarity == 'both') {

               # browser()
               # copy files into the working directory
               if (!file.exists(file.path(path, 'POS', sample_file_pos))) {
                 stop('The file ', sample_file_pos, ' NOT found in the directory, please check it!\n')
               }
               if (!file.exists(file.path(path, 'POS', sample_info_file_pos))) {
                 stop('The file ', sample_info_file, ' NOT found in the directory, please check it!\n')
               }
               if (!file.exists(file.path(path, 'POS', '00_annotation_table', table_identification_pair_file_pos))) {
                 stop('The file ', table_identification_pair_file_pos, ' NOT found in the directory, please check it!\n')
               }

               if (!file.exists(file.path(path, 'NEG', sample_file_pos))) {
                 stop('The file ', sample_file_pos, ' NOT found in the directory, please check it!\n')
               }
               if (!file.exists(file.path(path, 'NEG', sample_info_file_pos))) {
                 stop('The file', sample_info_file, ' NOT found in the directory, please check it!\n')
               }
               if (!file.exists(file.path(path, 'NEG', '00_annotation_table', table_identification_pair_file_pos))) {
                 stop('The file', table_identification_pair_file_pos, ' NOT found in the directory, please check it!\n')
               }


               # sample_file <- sample_file_neg
               # sample_info_file <- sample_info_file_neg
               # path <- path_neg
               # table_identification_pair_file <- table_identification_pair_file_neg

               # copy files
               path_output <- file.path(path, 'BOTH', '04_biology_intepretation')
               dir.create(path_output, showWarnings = FALSE, recursive = TRUE)
               file.copy(from = file.path(path, 'POS', sample_file_pos),
                         to = file.path(path_output, gsub('\\.csv', '_pos.csv', sample_file_pos)),
                         copy.date = TRUE,
                         overwrite = TRUE)
               file.copy(from = file.path(path, 'POS', sample_info_file_pos),
                         to = file.path(path_output, gsub('\\.csv', '_pos.csv', sample_info_file_pos)),
                         copy.date = TRUE,
                         overwrite = TRUE)
               file.copy(from = file.path(path, 'POS', '00_annotation_table', table_identification_pair_file_pos),
                         to = file.path(path_output, gsub('\\.csv', '_pos.csv', table_identification_pair_file_pos)),
                         copy.date = TRUE,
                         overwrite = TRUE)

               file.copy(from = file.path(path, 'NEG', sample_file_neg),
                         to = file.path(path_output, gsub('\\.csv', '_neg.csv', sample_file_neg)),
                         copy.date = TRUE,
                         overwrite = TRUE)
               file.copy(from = file.path(path, 'NEG', sample_info_file_neg),
                         to = file.path(path_output, gsub('\\.csv', '_neg.csv', sample_info_file_neg)),
                         copy.date = TRUE,
                         overwrite = TRUE)
               file.copy(from = file.path(path, 'NEG', '00_annotation_table', table_identification_pair_file_neg),
                         to = file.path(path_output, gsub('\\.csv', '_neg.csv', table_identification_pair_file_neg)),
                         copy.date = TRUE,
                         overwrite = TRUE)

               options(readr.num_columns = FALSE)
               sample_pos <- readr::read_csv(file.path(path_output, gsub('\\.csv', '_pos.csv', sample_file_pos)))
               sample_info_pos <- readr::read_csv(file.path(path_output, gsub('\\.csv', '_pos.csv', sample_info_file_pos)), col_types = c('cc'))
               table_identification_pair_pos <- readr::read_csv(file.path(path_output, gsub('\\.csv', '_pos.csv', table_identification_pair_file_pos)))

               sample_neg <- readr::read_csv(file.path(path_output, gsub('\\.csv', '_neg.csv', sample_file_neg)))
               sample_info_neg <- readr::read_csv(file.path(path_output, gsub('\\.csv', '_neg.csv', sample_info_file_neg)), col_types = c('cc'))
               table_identification_pair_neg <- readr::read_csv(file.path(path_output, gsub('\\.csv', '_neg.csv', table_identification_pair_file_neg)))

               if (!all(sort(sample_info_pos$sample.name) == sort(sample_info_neg$sample.name))) {
                 stop('Sample names are different in sample info between positive and negative modes. Please check!\n')
               }

               sample_neg <- match(colnames(sample_pos), colnames(sample_neg)) %>% sample_neg[,.]

               # merge pos and neg
               sample_pos$name <- paste(sample_pos$name, 'POS', sep = '_')
               table_identification_pair_pos$peak_name <- paste(table_identification_pair_pos$peak_name, 'POS', sep = '_')
               sample_neg$name <- paste(sample_neg$name, 'NEG', sep = '_')
               table_identification_pair_neg$peak_name <- paste(table_identification_pair_neg$peak_name, 'NEG', sep = '_')

               sample <- sample_pos %>% dplyr::bind_rows(sample_neg)
               sample_info <- sample_info_pos

               table_identification_pair <- table_identification_pair_pos %>% dplyr::bind_rows(table_identification_pair_neg)
               # export merged table_identification_pair_file for quantitative analysis
               readr::write_csv(table_identification_pair,
                                file = file.path(path_output, gsub('\\.csv', '_pos_and_neg.csv', "table_identification_pair.csv")))

               if (all(c('name', 'mz', 'rt', 'ccs') %in% colnames(sample))) {
                 sample <- sample %>%
                   tibble::column_to_rownames(var = 'name') %>%
                   dplyr::select(-c('mz', 'rt', 'ccs'))
               } else {
                 sample <- sample %>%
                   tibble::column_to_rownames(var = 'name') %>%
                   dplyr::select(-c('mz', 'rt'))
               }

               temp_error <- sample_info %>%
                 dplyr::filter(group %in% comp_group) %>%
                 dplyr::count(group)
               if (any(temp_error$n < 3)) {
                 cat('Note: one group has <3 samples, the biological sample interpreation will be skip\n')
                 return(NULL)
               }

               cat('Start univariate analsis ...\n')
               result_stat <- MetBioInterpretation::doStatAnalysis(sample = sample,
                                                                   sample_info = sample_info,
                                                                   comp_group = comp_group,
                                                                   uni_test = uni_test,
                                                                   correct = correct_p,
                                                                   p_cutoff = p_cutoff,
                                                                   fc_cutoff = fc_cutoff,
                                                                   by_what = 'mean')

               dir.create(file.path(path_output, '00_intermediate_data'),
                          showWarnings = FALSE, recursive = TRUE)
               save(result_stat,
                    file = file.path(path_output, '00_intermediate_data', 'result_stat.RData'),
                    version = 2)
               readr::write_csv(result_stat,
                                path = file.path(path_output, 'statistics_result.csv'))

               cat('Plot vocano plot...\n')
               temp_plot <- MetBioInterpretation::plotVolcano(raw_data = result_stat,
                                                              p_cutoff = p_cutoff,
                                                              fc_cutoff = fc_cutoff,
                                                              log_tran = 'log2')

               ggplot2::ggsave(temp_plot,
                               filename = file.path(path_output, 'volcano_plot.pdf'),
                               width = 6, height = 6)


               cat('\n\n');cat('Perform pathway enrichment analysis ...\n')

               if (metdna_version == 'version1') {
                 # extract significant feature
                 feature_significant <- result_stat %>%
                   dplyr::filter(p_values <= p_cutoff & ((fold_changes >= fc_cutoff) | (fold_changes <= 1/fc_cutoff))) %>%
                   dplyr::pull(name)

                 load(file.path(path, 'POS/02_result_MRN_annotation/00_intermediate_data/tags2_after_redundancy_remove'))
                 tags2_after_redundancy_remove_pos <- tags2_after_redundancy_remove;rm(tags2_after_redundancy_remove);gc()

                 load(file.path(path, 'NEG/02_result_MRN_annotation/00_intermediate_data/tags2_after_redundancy_remove'))
                 tags2_after_redundancy_remove_neg <- tags2_after_redundancy_remove;rm(tags2_after_redundancy_remove);gc()

                 names(tags2_after_redundancy_remove_pos) <- sapply(tags2_after_redundancy_remove_pos, function(x){
                   paste(x@name,'POS', sep = '_')
                 })

                 names(tags2_after_redundancy_remove_neg) <- sapply(tags2_after_redundancy_remove_neg, function(x){
                   paste(x@name,'NEG', sep = '_')
                 })

                 idx_pos <- match(feature_significant, names(tags2_after_redundancy_remove_pos))
                 idx_pos <- which(!is.na(idx_pos)) %>% idx_pos[.]
                 temp_result_pos <- convertTags2ResultMatrix(tags2 = tags2_after_redundancy_remove_pos[idx_pos], candidate.num = 5)

                 idx_neg <- match(feature_significant, names(tags2_after_redundancy_remove_neg))
                 idx_neg <- which(!is.na(idx_neg)) %>% idx_neg[.]
                 temp_result_neg <- convertTags2ResultMatrix(tags2 = tags2_after_redundancy_remove_neg[idx_neg], candidate.num = 5)

                 met_sig <- c(temp_result_pos$to, temp_result_neg$to) %>% stringr::str_split(pattern = ';') %>% unlist()
                 met_sig <- (!is.na(met_sig)) %>% which() %>% met_sig[.] %>% unique()

               } else {
                 # extract significant feature
                 feature_significant <- result_stat %>%
                   dplyr::filter(p_values <= p_cutoff & ((fold_changes >= fc_cutoff) | (fold_changes <= 1/fc_cutoff))) %>%
                   dplyr::pull(name)

                 result_annotation_sig <- table_identification_pair %>%
                   dplyr::mutate(temp_peak_name = stringr::str_replace(peak_name, '_[a-z]+', '')) %>%
                   dplyr::filter(temp_peak_name %in% feature_significant) %>%
                   dplyr::select(-temp_peak_name)

                 if (organze_basis == 'kegg_id') {
                   met_sig <- result_annotation_sig %>%
                     dplyr::pull(id_kegg) %>%
                     unique() %>%
                     .[!is.na(.)]
                 } else {
                   met_sig <- result_annotation_sig %>%
                     tidyr::separate_rows(stereo_isomer_id) %>%
                     dplyr::pull(stereo_isomer_id) %>%
                     unique() %>%
                     .[!is.na(.)]
                 }
               }

               if (length(met_sig) > 0) {
                 result_pathway_enrichment <- MetBioInterpretation::enrichPathway(met_sig = met_sig,
                                                                                  species = species,
                                                                                  lib_pathway = lib_pathway,
                                                                                  test_method = 'hypergeometric')

                 if (length(result_pathway_enrichment) > 0) {
                   result_pathway_enrichment <- result_pathway_enrichment %>% tibble::rownames_to_column(var = 'pathway')

                   dir.create(file.path(path_output, '01_pathway_enrichment'),
                              showWarnings = FALSE, recursive = TRUE)
                   save(result_pathway_enrichment,
                        file = file.path(path_output, '00_intermediate_data', 'result_pathway_enrichment.RData'),
                        version = 2)
                   readr::write_csv(result_pathway_enrichment,
                                    path = file.path(path_output,
                                                     '01_pathway_enrichment',
                                                     'pathway_enrichment_analysis.csv'))

                   temp_plot <- MetBioInterpretation::plotPathwayScatter(result_pathway_enrichment = result_pathway_enrichment, p_cutoff = 0.05, p_adjust = FALSE)
                   ggplot2::ggsave(temp_plot,
                                   filename = file.path(path_output,
                                                        '01_pathway_enrichment',
                                                        'pathway_enrichment_overview.pdf'),
                                   width = 6, height = 6)

                   temp_plot <- MetBioInterpretation::plotPathwayBar(result_pathway_enrichment = result_pathway_enrichment, p_cutoff = 0.05, p_adjust = FALSE)
                   ggplot2::ggsave(temp_plot,
                                   filename = file.path(path_output,
                                                        '01_pathway_enrichment',
                                                        'pathway_enrichment_MSEA.pdf'),
                                   width = 6, height = 6)
                 }
               } else {
                 cat('Note: No significant metabolites, skip remain analysis, including pathway enrichment, pathway quantative analysis, and class enrichment\n')
                 return(NULL)
               }


               cat('\n\n');cat('Perform pathway quantitative analysis ...\n')
               dir.create(file.path(path_output, '02_pathway_quantitative_analysis'),
                          showWarnings = FALSE, recursive = TRUE)

               if (metdna_version == 'version1') {
                 # extract significant feature from MRN intermediate data
                 load(file.path(path, 'POS/02_result_MRN_annotation/00_intermediate_data/id_result_redun_rm'))
                 id_result_redun_rm_pos <- id_result_redun_rm; rm(id_result_redun_rm); gc()

                 id_result_redun_rm_pos <- lapply(unique(id_result_redun_rm_pos$name), function(x){
                   temp_idx <- which(id_result_redun_rm_pos$name == x)
                   temp_data <- id_result_redun_rm_pos[temp_idx,,drop = FALSE]
                   if(nrow(temp_data) > 5) temp_data <- temp_data[1:5,]
                   temp_data
                 })

                 id_result_redun_rm_pos <- do.call(rbind, id_result_redun_rm_pos) %>% tibble::as_tibble()
                 id_result_redun_rm_pos <- id_result_redun_rm_pos %>%
                   dplyr::mutate(name = paste(name, 'POS', sep = '_')) %>%
                   dplyr::filter(isotope == '[M]') %>%
                   dplyr::filter(name %in% feature_significant)

                 temp_idx <- match(id_result_redun_rm_pos$name, rownames(sample))
                 id_result_redun_rm_pos <- id_result_redun_rm_pos %>%
                   dplyr::select(name, mz, rt, to, Confidence) %>%
                   dplyr::bind_cols(sample[temp_idx,]) %>%
                   dplyr::rename(peak_name = name,
                                 id_kegg = to,
                                 confidence_level = Confidence)

                 data("cpd_emrn", package = 'MetDNA2')
                 data("zhuMetlib", package = 'MetDNA2')
                 temp_idx2 <- match(id_result_redun_rm_pos$id_kegg, cpd_emrn$id)
                 id_result_redun_rm_pos <- id_result_redun_rm_pos %>%
                   dplyr::mutate(name = cpd_emrn$name[temp_idx2],
                                 inchikey = cpd_emrn$inchikey[temp_idx2],
                                 stereo_isomer_id = cpd_emrn$id_kegg_synonyms[temp_idx2],
                                 stereo_isomer_name = cpd_emrn$name_synonyms[temp_idx2]) %>%
                   dplyr::select(peak_name:confidence_level, name, inchikey, stereo_isomer_id, stereo_isomer_name, dplyr::everything())

                 temp_idx3 <- id_result_redun_rm_pos$name %>% is.na() %>% which()
                 temp_idx4 <- match(id_result_redun_rm_pos$id_kegg[temp_idx3], zhuMetlib$meta$compound$id_kegg)
                 id_result_redun_rm_pos$name[temp_idx3] <- zhuMetlib$meta$compound$name[temp_idx4]
                 id_result_redun_rm_pos$inchikey[temp_idx3] <- zhuMetlib$meta$compound$inchikey[temp_idx4]
                 id_result_redun_rm_pos$peak_name <- id_result_redun_rm_pos$peak_name %>% rename4UniqueIndex()



                 # extract significant feature from MRN intermediate data
                 load(file.path(path, 'NEG/02_result_MRN_annotation/00_intermediate_data/id_result_redun_rm'))
                 id_result_redun_rm_neg <- id_result_redun_rm; rm(id_result_redun_rm); gc()

                 id_result_redun_rm_neg <- lapply(unique(id_result_redun_rm_neg$name), function(x){
                   temp_idx <- which(id_result_redun_rm_neg$name == x)
                   temp_data <- id_result_redun_rm_neg[temp_idx,,drop = FALSE]
                   if(nrow(temp_data) > 5) temp_data <- temp_data[1:5,]
                   temp_data
                 })

                 id_result_redun_rm_neg <- do.call(rbind, id_result_redun_rm_neg) %>% tibble::as_tibble()
                 id_result_redun_rm_neg <- id_result_redun_rm_neg %>%
                   dplyr::mutate(name = paste(name, 'NEG', sep = '_')) %>%
                   dplyr::filter(isotope == '[M]') %>%
                   dplyr::filter(name %in% feature_significant)

                 temp_idx <- match(id_result_redun_rm_neg$name, rownames(sample))
                 id_result_redun_rm_neg <- id_result_redun_rm_neg %>%
                   dplyr::select(name, mz, rt, to, Confidence) %>%
                   dplyr::bind_cols(sample[temp_idx,]) %>%
                   dplyr::rename(peak_name = name,
                                 id_kegg = to,
                                 confidence_level = Confidence)

                 # data("cpd_emrn", package = 'MetDNA2')
                 # data("zhuMetlib", package = 'MetDNA2')
                 temp_idx2 <- match(id_result_redun_rm_neg$id_kegg, cpd_emrn$id)
                 id_result_redun_rm_neg <- id_result_redun_rm_neg %>%
                   dplyr::mutate(name = cpd_emrn$name[temp_idx2],
                                 inchikey = cpd_emrn$inchikey[temp_idx2],
                                 stereo_isomer_id = cpd_emrn$id_kegg_synonyms[temp_idx2],
                                 stereo_isomer_name = cpd_emrn$name_synonyms[temp_idx2]) %>%
                   dplyr::select(peak_name:confidence_level, name, inchikey, stereo_isomer_id, stereo_isomer_name, dplyr::everything())

                 temp_idx3 <- id_result_redun_rm_neg$name %>% is.na() %>% which()
                 temp_idx4 <- match(id_result_redun_rm_neg$id_kegg[temp_idx3], zhuMetlib$meta$compound$id_kegg)
                 id_result_redun_rm_neg$name[temp_idx3] <- zhuMetlib$meta$compound$name[temp_idx4]
                 id_result_redun_rm_neg$inchikey[temp_idx3] <- zhuMetlib$meta$compound$inchikey[temp_idx4]
                 id_result_redun_rm_neg$peak_name <- id_result_redun_rm_neg$peak_name %>% rename4UniqueIndex()


                 id_result_redun_rm <- id_result_redun_rm_pos %>% dplyr::bind_rows(id_result_redun_rm_neg)

                 readr::write_csv(id_result_redun_rm, file = file.path(path_output, "table3_identification_pair_v1_pos_and_neg.csv"))

                 MetBioInterpretation::selectQuantiMetabolite(table_identification_pair = "table3_identification_pair_v1_pos_and_neg.csv",
                                                              sample_info = gsub('\\.csv', '_pos.csv', sample_info_file_pos),
                                                              comp_group = comp_group,
                                                              path = path_output,
                                                              is_scale = TRUE,
                                                              scale_method = "pareto",
                                                              represent_met = 'median',
                                                              organze_basis = organze_basis)


               } else {
                 options(readr.num_columns = FALSE)
                 # table_identification_pair <- readr::read_csv(file.path(path_output, table_identification_pair_file))
                 temp_result <- table_identification_pair %>%
                   dplyr::mutate(temp_peak_name = stringr::str_replace(peak_name, '_[a-z]+', '')) %>%
                   dplyr::filter(temp_peak_name %in% feature_significant) %>%
                   dplyr::select(-temp_peak_name)
                 readr::write_csv(temp_result, file = file.path(path_output, "table3_identification_pair_v2_pos_and_neg.csv"))

                 MetBioInterpretation::selectQuantiMetabolite(table_identification_pair = "table3_identification_pair_v2_pos_and_neg.csv",
                                                              sample_info = gsub('\\.csv', '_pos.csv', sample_info_file_pos),
                                                              comp_group = comp_group,
                                                              path = path_output,
                                                              is_scale = TRUE,
                                                              scale_method = "pareto",
                                                              represent_met = 'median',
                                                              organze_basis = organze_basis)
               }

               # cat('\n\n');cat('Perform pathway quantitative analysis ...\n')
               # dir.create(file.path(path_output, '02_pathway_quantitative_analysis'),
               #            showWarnings = FALSE, recursive = TRUE)
               #
               # MetBioInterpretation::selectQuantiMetabolite(table_identification_pair = gsub('\\.csv', '_pos_and_neg.csv', "table_identification_pair.csv"),
               #                                              sample_info = gsub('\\.csv', '_pos.csv', sample_info_file_pos),
               #                                              comp_group = comp_group,
               #                                              path = path_output,
               #                                              is_scale = TRUE,
               #                                              scale_method = "pareto",
               #                                              organze_basis = organze_basis)

               temp <- try(MetBioInterpretation::quantiPathway(quanti_table = 'pathway_metabolite_quantitative_result.csv',
                                                               sample_info = gsub('\\.csv', '_pos.csv', sample_info_file_pos),
                                                               comp_group = comp_group,
                                                               is_scale = FALSE, # peak area has been scaled in selectQuantiMetabolite
                                                               correct_p = FALSE,
                                                               uni_test = uni_test,
                                                               path = path_output,
                                                               by_what = 'mean',
                                                               species = species,
                                                               is_consider_stereo_isomer = ifelse(organze_basis == 'kegg_id', FALSE, TRUE)),
                           silent = TRUE)

               if (class(temp) == 'try-error') {
                 cat('Note: quantiPathway is error. Pay attention, please!\n')
                 return(NULL)
               }

               cat('\n\n');cat('Export pathway quantitative analysis results...\n')
               file.copy(from = file.path(path_output, 'pathway_quantitative_result.csv'),
                         to = file.path(path_output, '02_pathway_quantitative_analysis'),
                         overwrite = TRUE, copy.date = TRUE)
               file.copy(from = file.path(path_output, 'pathway_metabolite_quantitative_result.csv'),
                         to = file.path(path_output, '02_pathway_quantitative_analysis'),
                         overwrite = TRUE, copy.date = TRUE)

               options(readr.num_columns = 0)
               quanti_pathway_data <- readr::read_csv(file.path(path_output,
                                                                '02_pathway_quantitative_analysis',
                                                                'pathway_quantitative_result.csv'))
               quanti_met_data <- readr::read_csv(file.path(path_output,
                                                            '02_pathway_quantitative_analysis',
                                                            'pathway_metabolite_quantitative_result.csv'))

               temp_plot <- MetBioInterpretation::plotPathwayHeatmap(quanti_pathway_data = quanti_pathway_data,
                                                                     sample_info = sample_info,
                                                                     comp_group = comp_group)
               ggplot2::ggsave(temp_plot,
                               filename = file.path(path_output,
                                                    '02_pathway_quantitative_analysis',
                                                    'pathway_heatmap.pdf'))

               cat('\n');cat('Plot boxplot for each pathway ...\n')

               MetBioInterpretation::plotPathwayBoxplot(quanti_pathway_data = quanti_pathway_data,
                                                        sample_info = sample_info,
                                                        comp_group = comp_group,
                                                        path = file.path(path_output, '02_pathway_quantitative_analysis'))

               cat('\n');cat('Plot metabolite heatmap for each pathway ...\n')

               MetBioInterpretation::plotPathwayHeatmapMet(quanti_met_data = quanti_met_data,
                                                           quanti_pathway_data = quanti_pathway_data,
                                                           sample_info = sample_info,
                                                           comp_group = comp_group,
                                                           species = species,
                                                           is_consider_stereo_isomer = ifelse(organze_basis == 'kegg_id', FALSE, TRUE),
                                                           path = file.path(path_output, '02_pathway_quantitative_analysis'))


               # cat('\n');cat('Pathway analysis has done ...\n')
               cat('\n\n');cat('Enrich chemical class analysis...\n')
               dir.create(file.path(path_output, '03_class_enrichment'),
                          showWarnings = FALSE, recursive = TRUE)
               table_class_enrich <- MetBioInterpretation::generateClassEnrichTable(table_identification_pair = table_identification_pair,
                                                                                    result_stat = result_stat)

               purrr::walk(c('superclass', 'class', 'subclass'), function(x){
                 switch(x,
                        'superclass' = {
                          temp <- table_class_enrich %>%
                            dplyr::select(peak_name, p_values, fold_changes, superclass) %>%
                            dplyr::rename('compound_name' = peak_name,
                                          'pvalue' = p_values,
                                          'foldchange' = fold_changes,
                                          'class' = superclass) %>%
                            dplyr::filter(!is.na(class)) %>%
                            dplyr::arrange(pvalue) %>%
                            dplyr::mutate(order = seq(dplyr::n())) %>%
                            dplyr::select(compound_name, order, dplyr::everything())
                        },
                        'class' = {
                          temp <- table_class_enrich %>%
                            dplyr::select(peak_name, p_values, fold_changes, class) %>%
                            dplyr::rename('compound_name' = peak_name,
                                          'pvalue' = p_values,
                                          'foldchange' = fold_changes,
                                          'class' = class) %>%
                            dplyr::filter(!is.na(class)) %>%
                            dplyr::arrange(pvalue) %>%
                            dplyr::mutate(order = seq(dplyr::n())) %>%
                            dplyr::select(compound_name, order, dplyr::everything())
                        },
                        'subclass' = {
                          temp <- table_class_enrich %>%
                            dplyr::select(peak_name, p_values, fold_changes, subclass) %>%
                            dplyr::rename('compound_name' = peak_name,
                                          'pvalue' = p_values,
                                          'foldchange' = fold_changes,
                                          'class' = subclass) %>%
                            dplyr::filter(!is.na(class)) %>%
                            dplyr::arrange(pvalue) %>%
                            dplyr::mutate(order = seq(dplyr::n())) %>%
                            dplyr::select(compound_name, order, dplyr::everything())
                        })

                 readr::write_csv(temp,
                                  file = file.path(path_output, '03_class_enrichment', paste0(x, '_table_class_enrichment.csv')))


                 enrichment_result <- try(MetBioInterpretation::enrichClass(table_class_file = paste0(x, '_table_class_enrichment.csv'),
                                                                            path = file.path(path_output, '03_class_enrichment')),
                                          silent = TRUE)

                 if (class(enrichment_result) == 'try-error') {
                   cat('Note: class enrichment has some error. Pay attention, please!\n')
                   return(NULL)
                 }

                 readr::write_csv(enrichment_result,
                                  file = file.path(path_output, '03_class_enrichment', paste0(x, '_result_class_enrichment.csv')))

                 temp_plot <- MetBioInterpretation::plotClassScater(enrichment_result = enrichment_result)
                 dir.create(file.path(path_output, '03_class_enrichment', 'plot'), showWarnings = FALSE, recursive = TRUE)
                 ggplot2::ggsave(temp_plot,
                                 filename = file.path(path_output, '03_class_enrichment', 'plot', paste0(x, '_class_enrichment.pdf')),
                                 width = 6, height = 6)

               })
             }

           })

