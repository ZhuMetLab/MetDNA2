#' ################################################################################
#' # enrichPathway ----------------------------------------------------------------
#'
#' #' @title enrichPathway
#' #' @author Zhiwei Zhou
#' #' \email{zhouzw@@sioc.ac.cn}
#' #' @param met_sig significant metabolite ID
#' #' @param species "hsa", "mmu", "rno", "bta", "gga", "dre", "dme", "cel", "sce", "osa", "ath", "smm", "pfa", "tbr", "eco", "bsu", "ppu", "sau", "tma", "syf", "mln"
#' #' @param test_method "hypergenometric" or "fisher"
#' #' @export
#'
#' # load('I:/00_projects/01_AllCCS/00 data/200420_mouse_liver_pathway_to_related_metabolites/result_stat_all_200420.RData')
#' # load('I:/00_projects/01_AllCCS/00 data/200414_mose_liver_extend_db_processing_update/03_mz_ccs_msms_match/01_annotation_result/msfinder_process/msms_score_msfinder_200419_2.RData')
#' #
#' # result_annotation <- readr::read_csv('H:/00_projects/03_MetDNA2/00_data/20210303_biologInterpretation_development/20210125_fly_pos_emrn0/00_annotation_table/id_merge_modify_long.csv')
#'
#' #
#' # load('H:/00_projects/03_MetDNA2/00_data/20210303_biologInterpretation_development/20210125_fly_pos_emrn0/00_annotation_table/00_intermediate_data/id_merge_modify')
#' # list_id <- result_annotation %>%
#' #   dplyr::filter(!is.na(id)) %>%
#' #   dplyr::filter(recursive_type!='isotopeAnnotation' | is.na(recursive_type))
#' #
#' # sig_feature <- stat_result %>%
#' #   filter(p_values <= 0.05) %>%
#' #   pull(name)
#' #
#' # sig_id_result <- list_id %>%
#' #   filter(peak_name %in% sig_feature)
#' #
#' # inchikey1_sig <- sig_id_result %>% pull(inchikey1)
#' # met_sig <- mapInChIKey1ToKEGG(inchikey1 = inchikey1_sig, lib_cpd = cpd_emrn)
#' #
#' # met_sig_unique <- met_sig %>% stringr::str_split(pattern = ';') %>% unlist() %>% unique() %>% .[!is.na(.)]
#' #
#' # test <- enrichPathway(met_sig = met_sig,
#' #                       species = 'dme',
#' #                       test_method = 'hypergeometric')
#'
#'
#'
#'
#' setGeneric(name = 'enrichPathway',
#'            def = function(
#'              met_sig,
#'              lib_pathway,
#'              species = c("hsa", "mmu", "rno", "bta", "gga",
#'                          "dre", "dme", "cel", "sce", "osa",
#'                          "ath", "smm", "pfa", "tbr", "eco",
#'                          "bsu", "ppu", "sau", "tma", "syf", "mln"),
#'              test_method = c("hypergeometric", "fisher")
#'            ){
#'
#'              metabolite_id <- unique(met_sig)
#'              species <- match.arg(species)
#'              test_method <- match.arg(test_method)
#'
#'              if (missing(lib_pathway)) {
#'                data('lib_pathway', envir =  environment())
#'              }
#'
#'
#'              M <- which(names(lib_pathway) == species) %>% lib_pathway[[.]]
#'
#'              metabolite_id <- unname(metabolite_id)
#'
#'              metabolite_id <- metabolite_id[which(metabolite_id %in% unique(unlist(M)))]
#'              if(length(metabolite_id) == 0) return(NULL)
#'
#'              ALL <- unname(unique(unlist(M))) # all metabolite id
#'              # ALL <- as.character(as.matrix(ALL))
#'              SIG <- as.character(as.matrix(metabolite_id)) # all significant metabolite id
#'
#'              num_all <- length(ALL) # number of all metabolites
#'              num_sig <- length(SIG) # number of significant metabolites
#'
#'              Lall0 <- unlist(lapply(M, function(x) {
#'                length(intersect(x, ALL))
#'              }))
#'
#'              Lall <- Lall0 # list of all metabolites for each pathway
#'
#'              Lsig <- unlist(lapply(M, function(x) {
#'                length(intersect(x, SIG))
#'              })) # list of significant metabolites for each pathway
#'
#'              # remove.idx <- which(unname(Lall0, length) == 0)
#'              remove.idx <- which(unname(Lsig, length) == 0)
#'
#'              if(length(remove.idx) > 0){
#'                Lall <- Lall0[-remove.idx]
#'                Lsig <- Lsig[-remove.idx]
#'              }
#'
#'              # error handling
#'              if (length(Lall) < 2){
#'                # stop function
#'                return(NULL)
#'              }
#'
#'
#'              P <- NaN
#'              system.time(for (i in 1:length(Lall)){
#'                # ------------------------------------
#'                #Generating 2?~2 table
#'                # -------------------------------------
#'                a1 <- Lsig[i] # significant and including pathway
#'                a2 <- Lall[i] - Lsig[i] # not significant and including pathway
#'                a3 <- length(SIG) - a1 # significant and not including pathway
#'                a4 <- (length(ALL) - length(SIG))-a2 # not significant and not including pathway
#'
#'                if(test_method == "hypergeometric"){
#'                  P[i] <- phyper(q = a1 - 1,
#'                                 m = Lall[i],
#'                                 n = num_all - Lall[i],
#'                                 k = num_sig,
#'                                 lower.tail = FALSE)
#'                } else {
#'                  tab <- t(matrix(c(a1,a2,a3,a4),2))
#'
#'                  # Fisher's exact test
#'
#'                  check <- tryCatch({
#'                    resfish <- fisher.test(tab, alternative="greater")
#'                  }, error = function(e){
#'                    NA
#'                  })
#'
#'                  if(class(check) != "htest"){
#'                    P[i] <- 1
#'                  } else {
#'                    resfish <- fisher.test(tab, alternative="greater")
#'                    P[i] <- resfish$p.value
#'                  }
#'                }
#'              })
#'
#'              # -----------------------
#'              # q-value
#'
#'              Q <- p.adjust(P, method="BH")
#'              FDR <- p.adjust(P, method="fdr")
#'
#'              # --------------------------------------------------------
#'              #significant metabolites for metabolite set
#'              # --------------------------------------------------------
#'              # LES <- NaN
#'              # for (i in 1:ncol(Lsig)){
#'              #   les <- SIG[Lsig[,i]==1]
#'              #   LES[i] <- list(les)
#'              #
#'              # }
#'              # names(LES) <- colnames(Lsig)
#'
#'              # Result ----------------------------------------------------------
#'              PQ <- cbind(P,Q, FDR)
#'              rownames(PQ) <- colnames(Lsig)
#'              colnames(PQ) <- c("p.value","q.value", "FDR")
#'
#'              info <- lapply(M[-remove.idx], function(module) {
#'                overlap.number <- length(intersect(module, metabolite_id))
#'                pathway.number <- length(module)
#'                c(pathway.number, overlap.number)
#'              })
#'
#'              info <- do.call(rbind, info)
#'              colnames(info) <- c("Pathway.length", "Overlap")
#'              info <- data.frame(PQ, info, stringsAsFactors = FALSE)
#'              info <- info[order(info[,1]),]
#'
#'              # return
#'              info <- info
#'
#'              return(info)
#'
#'
#'            }
#' )
#'
#' #   mapInChIKey1ToKEGG ---------------------------------------------------------
#' #' @title mapInChIKey1AndIds
#' #' @author Zhiwei Zhou
#' #' @param items
#' #' @param lib_cpd
#' #' @param mode 'inchikey1_kegg', 'kegg_inchikey1', 'inchikey1_kegg_name'
#' # load('H:/00_projects/03_MetDNA2/00_data/20210303_biologInterpretation_development/cpd_emrn_210223.RData')
#' # lib_cpd <- cpd_emrn
#'
#' setGeneric(name = 'mapInChIKey1AndIds',
#'            def = function(items,
#'                           lib_cpd,
#'                           mode = c('inchikey1_kegg',
#'                                    'kegg_inchikey1',
#'                                    'inchikey1_kegg_name'),
#'                           ...){
#'              mode <- match.arg(mode)
#'
#'              if (missing(lib_cpd)) {
#'                data('cpd_emrn', envir = environment())
#'                lib_cpd <- cpd_emrn
#'                rm(cpd_emrn);gc()
#'              }
#'
#'              if (mode == 'inchikey1_kegg') {
#'                idx <- match(items, lib_cpd$inchikey1)
#'                result_ids <- lib_cpd[idx,] %>% dplyr::pull(id_kegg_synonyms)
#'              }
#'
#'              if (mode == 'inchikey1_kegg_name') {
#'                idx <- match(items, lib_cpd$inchikey1)
#'                result_ids <- lib_cpd[idx,] %>% dplyr::pull(name_synonyms)
#'              }
#'
#'              if (mode == 'kegg_inchikey1') {
#'                idx <- match(items, lib_cpd$id)
#'                result_ids <- lib_cpd[idx,] %>% dplyr::pull(inchikey1)
#'              }
#'
#'              return(result_ids)
#'            })
#'
#'
#'
#' ################################################################################
#' # quantiPathway ----------------------------------------------------------------
#'
#' #' @title quantiPathway
#' #' @author Zhiwei Zhou
#' #' @description generate pathway quantification table
#' #' @param sample_pos
#' #' @param sample_neg
#' #' @param sample_info
#' #' @param result_annotation_pos
#' #' @param result_annotation_neg
#' #' @param result_stat
#' #' @param comp_group
#' #' @param metabolic_network
#' #' @param path Default: '.'
#' #' @param polarity 'positive', 'negative', 'both'. Default: 'positive'
#' #' @param is_scale
#' #' @param scale_method
#' #' @param species
#' #' @param represent_method 'mean', 'sum', 'median'. Default: 'mean'
#' #' @return "pathway_metabolite_quantitative_result.csv", "pathway_quantitative_result.csv",
#' #' @export
#' #' @examples
#'
#' setGeneric(name = 'quantiPathway',
#'            def = function(
#'              sample,
#'              sample_info,
#'              result_annotation,
#'              metabolic_network,
#'              lib_pathway,
#'              comp_group = c("W30", "W03"),
#'              uni_test = c("t","wilcox","anova"),
#'              correct_p = TRUE,
#'              by_what = c("median","mean"), # foldchange
#'              path = '.',
#'              # polarity = c('positive', 'negative', 'both'),
#'              is_scale = TRUE,
#'              scale_method = c("pareto", "auto", "z_score"),
#'              species = c("hsa", "mmu", "rno", "bta", "gga",
#'                          "dre", "dme", "cel", "sce", "osa",
#'                          "ath", "smm", "pfa", "tbr", "eco",
#'                          "bsu", "ppu", "sau", "tma", "syf", "mln"),
#'              represent_method = c('mean', 'sum', 'median'),
#'              ...
#'            ){
#'              # browser()
#'              # polarity <- match.arg(polarity)
#'              scale_method <- match.arg(scale_method)
#'              represent_method <- match.arg(represent_method)
#'
#'              # data('lib_pathway', envir =  environment())
#'              # data('cpd_emrn', envir =  environment())
#'
#'              if (missing(comp_group)){
#'                warning("You don't set group, so use the default group in sample_info.\n")
#'                comp_group <- unique(sample_info$group)
#'              }
#'
#'
#'              if (!all(comp_group %in% unique(sample_info$group))){
#'                idx <- which(!is.element(comp_group, unique(sample_info[,2])))
#'                stop(paste(comp_group[idx], collapse = " and "),
#'                     ifelse(length(idx) > 1, " are ", " is "),
#'                     "not in your sample_info.")
#'              }
#'
#'              sample_info <- sample_info %>%
#'                dplyr::filter(group %in% comp_group) %>%
#'                dplyr::arrange(group, sample.name)
#'
#'              # # retrieve sample data from ms1_data, and merge sample_pos and sample_neg
#'              # if (polarity == "positive" | polarity == "both"){
#'              #   temp <- paste(rownames(sample_pos), 'POS', sep = '_')
#'              #   idx <- match(sample_info$sample.name, colnames(sample_pos))
#'              #   sample_pos <- sample_pos[,idx]
#'              #   rownames(sample_pos) <- temp
#'              # }
#'              #
#'              # if (polarity == "negative" | polarity == "both"){
#'              #   temp <- paste(rownames(sample_neg), 'NEG', sep = '_')
#'              #   idx <- match(sample_info$sample.name, colnames(sample_neg))
#'              #   sample_neg <- sample_neg[,idx]
#'              #   rownames(sample_neg) <- temp
#'              # }
#'              #
#'              # switch(polarity,
#'              #        "positive" = {
#'              #          sample <- sample_pos
#'              #        },
#'              #        "negative" = {
#'              #          sample <- sample_neg
#'              #        },
#'              #        "both" = {
#'              #          sample <- sample_pos %>% dplyr::bind_rows(sample_neg)
#'              #        })
#'              #
#'              # if (polarity == "positive" | polarity == "both"){
#'              #   result_annotation_pos <- result_annotation_pos %>%
#'              #     dplyr::mutate(peak_name = paste(peak_name, 'POS', sep = '_'))
#'              # }
#'              #
#'              # if (polarity == "negative" | polarity == "both"){
#'              #   result_annotation_neg <- result_annotation_neg %>%
#'              #     dplyr::mutate(peak_name = paste(peak_name, 'NEG', sep = '_'))
#'              # }
#'              #
#'              # # retrieve and merge result_annotation_pos and result_annotation_neg
#'              # switch(polarity,
#'              #        "positive" = {
#'              #          result_annotation <- result_annotation_pos
#'              #        },
#'              #        "negative" = {
#'              #          result_annotation <- result_annotation_neg
#'              #        },
#'              #        "both" = {
#'              #          result_annotation <- result_annotation_pos %>%
#'              #            dplyr::bind_rows(result_annotation_neg)
#'              #        })
#'
#'              # scale sample
#'              if (is_scale){
#'                switch(scale_method,
#'                       "auto" = {sample <- t(apply(sample, 1, function(x) {x/sd(x)}))},
#'                       "pareto" = {sample <- t(apply(sample, 1, function(x) {x/sqrt(sd(x))}))},
#'                       "z_score" = {sample <- t(apply(sample, 1, function(x) {(x-mean(x))/sd(x)}))}
#'                )
#'              }
#'
#'              cat("\n");cat("Transform metabolite information to quantitative pathway information. \n")
#'              cat('\n');cat("Select quantification metabolites... \n")
#'              # select quantification metabolites
#'              temp_pathway <- which(names(lib_pathway) == species) %>% lib_pathway[[.]]
#'              all_ids <- unlist(temp_pathway) %>% unname() %>% unique()
#'              all_ids_inchikey1 <- mapInChIKey1AndIds(items = all_ids,
#'                                                      mode = 'kegg_inchikey1',
#'                                                      ...)
#'              all_ids_inchikey1 <- all_ids_inchikey1[!is.na(all_ids_inchikey1)] %>% unique()
#'
#'              records_inchikey1 <- lapply(all_ids_inchikey1, function(x){
#'                temp_result_annotation <- result_annotation %>%
#'                  dplyr::filter(inchikey1 == x) %>%
#'                  dplyr::arrange(confidence_level) %>%
#'                  dplyr::filter(confidence_level == confidence_level[1])
#'
#'                temp_result_annotation
#'              })
#'
#'              names(records_inchikey1) <- all_ids_inchikey1
#'
#'              # remove the node which have no mapped peak, those node may be hidden nodes is network
#'              records_inchikey1 <- records_inchikey1[which(unlist(lapply(records_inchikey1, nrow)) != 0)]
#'
#'              # for each metabolite, select the peak which has the most intense intensity
#'              records_inchikey1 <- lapply(records_inchikey1, function(x){
#'                if(nrow(x) == 1) {return(x)}
#'                idx <- match(x$peak_name, rownames(sample))
#'                temp_int <- apply(sample[idx,], 1, mean)
#'                selected_peak <- which.max(temp_int) %>% x[.,]
#'              })
#'
#'              # remove id which are not in metabolomics.network
#'              network_inchikey1 <- mapInChIKey1AndIds(items = igraph::V(metabolic_network)$name,
#'                                                      # lib_cpd = lib_cpd,
#'                                                      mode = 'kegg_inchikey1', ...)
#'              records_inchikey1 <- (names(records_inchikey1) %in% network_inchikey1) %>%
#'                which() %>%
#'                records_inchikey1[.]
#'              records_inchikey1 <- records_inchikey1 %>% dplyr::bind_rows()
#'
#'              # output node quantitative information
#'              quanti_inchikey1 <- records_inchikey1$inchikey1
#'              quanti_peak_name <- records_inchikey1$peak_name
#'              quanti_id <- records_inchikey1$id
#'              quanti_cpd_name <- records_inchikey1$name
#'              quanti_confidence_level <- records_inchikey1$confidence_level
#'              quanti_id_synonyms <- mapInChIKey1AndIds(items = quanti_inchikey1,
#'                                                       # lib_cpd = lib_cpd,
#'                                                       mode = 'inchikey1_kegg',
#'                                                       ...)
#'              quanti_cpd_name_synonyms <- mapInChIKey1AndIds(items = quanti_inchikey1,
#'                                                             # lib_cpd = lib_cpd,
#'                                                             mode = 'inchikey1_kegg_name',
#'                                                             ...)
#'              quanti_sample_data <- match(records_inchikey1$peak_name, rownames(sample)) %>%
#'                sample[.,] %>%
#'                tibble::as_tibble()
#'
#'              metabolite_quantitative_data <- tibble::tibble(
#'                inchikey1 = quanti_inchikey1,
#'                peak_name = quanti_peak_name,
#'                id = quanti_id,
#'                cpd_name = quanti_cpd_name,
#'                id_synonyms = quanti_id_synonyms,
#'                cpd_name_synonyms = quanti_cpd_name_synonyms,
#'                confidence_level = quanti_confidence_level
#'              ) %>% dplyr::bind_cols(quanti_sample_data)
#'
#'              readr::write_csv(metabolite_quantitative_data,
#'                               file.path(path, "pathway_metabolite_quantitative_result.csv"))
#'
#'
#'              cat('\n');cat("Merge metabolites for each pathways... \n")
#'              # transform metabolite information to group information
#'              temp_pathway <- which(names(lib_pathway) == species) %>% lib_pathway[[.]]
#'              temp_metabolite_quant_data <- metabolite_quantitative_data %>% tidyr::separate_rows('id_synonyms')
#'              quanti_pathway_data <- lapply(temp_pathway, function(x){
#'                temp_idx <- match(x, temp_metabolite_quant_data$id_synonyms)
#'                temp_idx <- temp_idx[!is.na(temp_idx)]
#'                if(length(temp_idx) < 3) {return(NULL)}
#'                temp_peak_name <- temp_metabolite_quant_data$peak_name[temp_idx]
#'                temp_sample <- sample[match(temp_peak_name, rownames(sample)),,drop = FALSE]
#'              })
#'
#'              remove_idx <- which(unlist(lapply(quanti_pathway_data, is.null)))
#'              if (length(remove_idx > 0)){
#'                quanti_pathway_data <- quanti_pathway_data[-remove_idx]
#'              }
#'
#'              # get module quantitative information
#'              quanti_pathway_data <- lapply(quanti_pathway_data, function(x){
#'                switch(represent_method,
#'                       'mean' = {
#'                         result <- x %>%
#'                           apply(., 2, mean)
#'                       },
#'                       'sum' = {
#'                         result <- x %>%
#'                           apply(., 2, sum)
#'                       },
#'                       'median' = {
#'                         result <- x %>%
#'                           apply(., 2, median)
#'                       })
#'
#'                return(result)
#'              })
#'
#'              pathway_name <- names(quanti_pathway_data)
#'              quanti_pathway_data <- quanti_pathway_data %>% dplyr::bind_rows()
#'
#'              quanti_pathway_data <- quanti_pathway_data %>%
#'                dplyr::mutate(pathway_name = pathway_name) %>%
#'                dplyr::select(pathway_name, dplyr::everything())
#'
#'              readr::write_csv(quanti_pathway_data,
#'                               file.path(path,
#'                                         "pathway_quantitative_result.csv"))
#'
#'              cat('\n');cat("Perform univariate test for quantitative pathways. \n")
#'              # browser()
#'              quanti_pathway_data2 <- quanti_pathway_data %>%
#'                tidyr::separate(pathway_name, into = c('pathway_name', 'pathway_id'), sep = ';') %>%
#'                tidyr::separate(pathway_name, into = c('pathway_name', 'temp'), sep = ' - ') %>%
#'                dplyr::mutate(pathway_name = paste(pathway_name, pathway_id, sep = ';')) %>%
#'                dplyr::select(-c('temp', 'pathway_id')) %>%
#'                tibble::column_to_rownames(var = 'pathway_name')
#'
#'              stat_result_pathway <- doStatAnalysis(sample = quanti_pathway_data2,
#'                                                    sample_info = sample_info,
#'                                                    comp_group = comp_group,
#'                                                    uni_test = uni_test,
#'                                                    correct = correct_p,
#'                                                    p_cutoff = p_cutoff,
#'                                                    by_what = by_what)
#'
#'              readr::write_csv(stat_result_pathway,
#'                               file.path(path,
#'                                         "stat_pathway_quantitative_result.csv"))
#'
#'
#'            })
