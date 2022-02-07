################################################################################
# SpecAnnotationClass ----------------------------------------------------------

setClass(Class = "SpecAnnotationClass",
         representation(peak_info = "list",
                        annotation_result = "data.frame")
)


setMethod(f = "show",
          signature = "SpecAnnotationClass",
          definition = function(object) {
            cat("-----------Meta information------------\n")
            cat("Feature:", object@peak_info$name, "\n")
            cat("m/z:",object@peak_info$mz, "\n")
            cat("RT (s):",object@peak_info$rt, "\n")
            cat("CCS (A2):",object@peak_info$ccs, "\n")
            cat("Median peak area:",object@peak_info$intmed, '\n\n')

            cat("-----------Annotation result------------\n")

            if(nrow(object@annotation_result) == 0) {
              cat("No annotation.")
            }else{
              for(i in 1:nrow(object@annotation_result)){
                cat("\n")
                cat("Annotation:", i, '\n')
                cat("ID:", object@annotation_result$id[i], "\n")
                cat("Name:", object@annotation_result$name[i], "\n")
                cat("Formula:", object@annotation_result$formula[i], "\n")
                cat("SMILES:", object@annotation_result$smiles[i], "\n")
                cat("InChIKey:", object@annotation_result$inchikey[i], "\n")
                cat("InChIKey1:", object@annotation_result$inchikey1[i], "\n")
                cat("Adduct:", object@annotation_result$adduct[i], "\n")
                cat("Lib m/z:", object@annotation_result$mz_lib[i], "\n")
                cat("Lib RT:", object@annotation_result$rt_lib[i], "\n")
                cat("Lib CCS:", object@annotation_result$ccs_lib[i], "\n")
                cat("m/z error (ppm):", object@annotation_result$mz_error[i], "\n")
                cat("RT error (s):", object@annotation_result$rt_error[i], "\n")
                cat("RT score:", object@annotation_result$rt_score[i], "\n")
                cat("CCS error (%):", object@annotation_result$ccs_error[i], "\n")
                cat("CCS score:", object@annotation_result$ccs_score[i], "\n")
                cat("MS2 score (forward):", object@annotation_result$msms_score_forward[i], "\n")
                cat("MS2 score (reverse):", object@annotation_result$msms_score_reverse[i], "\n")
                cat("Number of matched fragments:", object@annotation_result$msms_matched_fra[i], "\n\n")
              }

            }
          }
)



#' @title convertMs1Data2SpecAnnotationClass
#' @description Convert Ms1DataInfor to AnnotationResult class
#' @author Zhiwei Zhou
#' @param ms1_info
#' @export

# test <- Ms1DataInfo2AnnotationResult(ms1_info = ms1_info)


setGeneric(name = 'convertMs1Data2SpecAnnotationClass',
           def = function(
             ms1_info
           ){
             cat('Convert to SpecAnnotationClass...\n')
             temp_annotation_result <- tibble::tibble(idx = numeric(),
                                                      id = character(),
                                                      name = character(),
                                                      formula = character(),
                                                      smiles = character(),
                                                      inchikey = character(),
                                                      inchikey1 = character(),
                                                      adduct = character(),
                                                      mz_lib = numeric(),
                                                      rt_lib = numeric(),
                                                      ccs_lib = numeric(),
                                                      mz_error = numeric(),
                                                      # mz_score=mz_score,
                                                      rt_error = numeric(),
                                                      rt_score = numeric(),
                                                      ccs_error = numeric(),
                                                      ccs_score = numeric(),
                                                      msms_score_forward = numeric(),
                                                      msms_score_reverse = numeric(),
                                                      msms_matched_frag = numeric())


             result <- pbapply::pblapply(seq_along(ms1_info$name), function(i){
               peak_info <- ms1_info[i,] %>% as.list()
               new(Class = "SpecAnnotationClass",
                   peak_info = peak_info,
                   annotation_result = temp_annotation_result)
             })

             names(result) <- ms1_info$name

             return(result)
           }
)




################################################################################
# PeakInfo ---------------------------------------------------------------------
#' @title PeakInfo
#' @description The PeakInfo class is the S4 class for Peak annotation information.
#' @slot name Peak name.
#' @slot mz Peak m/z.
#' @slot rt Peak rentention time (RT).
#' @slot ms2 MS/MS spectrum of peak.
#' @slot annotation Annotation information of peak. It is a list contains different annotation result.
#' @export

setClass("PeakInfo",
         representation(name = "character",
                        mz = "numeric",
                        rt = "numeric",
                        ccs = "numeric",
                        ms2 = "data.frame",
                        annotation = "list")

)



setMethod(f = "show",
          signature   = "PeakInfo",
          definition = function(object) {
            cat("Peak information\n")
            cat("--------------\n")
            cat("Peak name:", object@name)
            cat("\n")
            cat("mz:", object@mz)
            cat("\n")
            cat("rt:", object@rt)
            cat("\n")
            cat("ccs:", object@ccs)
            cat("\n")
            cat("MS/MS spectrum:", ifelse(nrow(object@ms2) > 0, "Yes\n", "No\n"))
            # cat(ifelse(nrow(object@ms2) > 0, "Yes\n", "No\n"))
            # format(object@ms2)
            cat('\n')
            cat("--------------\n")
            cat("Annotation information\n")
            cat("--------------\n")
            if(length(object@annotation) == 0){
              cat("No annotation.")
            }else{
              for(i in 1:length(object@annotation)){
                cat("\n")
                cat("Annotation:",i)
                cat("\n")
                cat("Type:", object@annotation[[i]]$type)
                cat("\n")
                cat("From:", object@annotation[[i]]$From)
                cat("\n")
                cat("From which peak:", object@annotation[[i]]$From.peak)
                cat("\n")
                cat("Step:", object@annotation[[i]]$step)
                cat("\n")
                cat("Annotation result (KEGG ID or EMRN ID):", object@annotation[[i]]$to)
                cat("\n")
                cat("Annotation level:", object@annotation[[i]]$level)
                cat("\n")
                cat("Is seed or not:", object@annotation[[i]]$as.seed)
                cat("\n")
                cat("In which round as seed:", object@annotation[[i]]$as.seed.round)
                cat("\n")
                cat("Isotpe:", object@annotation[[i]]$isotope)
                cat("\n")
                cat("Adduct:", object@annotation[[i]]$adduct)
                cat("\n")
                cat("Charge:", object@annotation[[i]]$charge)
                cat("\n")
                cat("Formula:", object@annotation[[i]]$formula)
                cat("\n")
                cat("m/z error (ppm):", object@annotation[[i]]$mz.error)
                cat("\n")
                cat("RT error (%/s):", object@annotation[[i]]$rt.error)
                cat("\n")
                cat("CCS error (%):", object@annotation[[i]]$ccs.error)
                cat("\n")
                cat("Intensity ratio error:", object@annotation[[i]]$int.error)
                cat("\n")
                cat("MS/MS similarity:", object@annotation[[i]]$ms2.sim)
                cat("\n")
                cat("Matched fragments:", object@annotation[[i]]$nfrag)
                cat("\n")
                cat("Score:", object@annotation[[i]]$score)
                cat("\n")
              }
            }

          }
)


setGeneric(
  name = "getInfo",
  def  = function(object,...) {
    standardGeneric("getInfo")
  }
)


setMethod(f = "getInfo",
          signature  = "PeakInfo",
          definition = function(object) {
            name <- object@name
            mz <- object@mz
            rt <- object@rt
            ccs <- object@ccs
            annotation <- object@annotation
            if (length(annotation) == 1){
              annotation <- unlist(annotation)
            } else {
              idx.max.score <- which.max(unlist(lapply(annotation, function(x) {x$score})))
              annotation <- unlist(annotation[[idx.max.score]])
              names(annotation)[19] <- "Score"
            }
            result <- c(name, mz, rt, ccs, annotation)
            names(result)[1:4] <- c("name", "mz", "rt", 'ccs')
            return(result)
          }
)



################################################################################
# showTags2 --------------------------------------------------------------------
#' @title showTags2
#' @description Show detail information of tags2 data.
#' @author Xiaotao Shen, Zhiwei Zhou
#' @param tags2 The tags2 data.
#' @param slot The inforamtion.
#' @return The information of tags2 data.
#' @export

setGeneric(name = "showTags2",
           def = function(tags2,
                          slot = c("name", "mz", "rt", "ccs",
                                   "annotation.len", "annotation.id",
                                   "annotation.from", "annotation.from.peak",
                                   "annotation.score", "annotation.level",
                                   "annotation.step", "seed.number",
                                   "annotation.type","annotation.formula",
                                   "as.seed")){
             slot <- match.arg(slot)
             if(slot == "name"){
               name <- unlist(lapply(tags2, function(x) {x@name}))
               return(name)
             }

             if(slot == "mz"){
               mz <- unlist(lapply(tags2, function(x) {x@mz}))
               return(mz)
             }

             if(slot == "rt"){
               rt <- unlist(lapply(tags2, function(x) {x@rt}))
               return(rt)
             }

             if(slot == "ccs"){
               ccs <- unlist(lapply(tags2, function(x) {x@ccs}))
               return(ccs)
             }


             annotation.len <- unlist(lapply(tags2, function(x) {length(x@annotation)}))

             if(slot == "annotation.len"){
               return(annotation.len)
             }

             index <- which(annotation.len > 0)
             if(length(index) == 0) {cat('No peaks have annotation.\n')
             }

             if(slot == "annotation.id"){
               id <- lapply(tags2[index], function(x) {
                 annotation <- x@annotation
                 lapply(annotation, function(x) {x$to})
               })
               names(id) <- index
               return(id)
             }

             if(slot == "annotation.type"){
               type <- lapply(tags2[index], function(x) {
                 annotation <- x@annotation
                 lapply(annotation, function(x) {x$type})
               })
               names(type) <- index
               return(type)
             }


             if(slot == "annotation.from"){
               from <- lapply(tags2[index], function(x) {
                 annotation <- x@annotation
                 lapply(annotation, function(x) {x$From})
               })
               names(from) <- index
               return(from)
             }

             if(slot == "annotation.formula"){
               formula <- lapply(tags2[index], function(x) {
                 annotation <- x@annotation
                 lapply(annotation, function(x) {x$Formula})
               })
               names(formula) <- index
               return(formula)
             }


             if(slot == "annotation.from.peak"){
               from.peak <- lapply(tags2[index], function(x) {
                 annotation <- x@annotation
                 lapply(annotation, function(x) {x$From.peak})
               })
               names(from.peak) <- index
               return(from.peak)
             }

             if(slot == "annotation.score"){
               score <- lapply(tags2[index], function(x) {
                 annotation <- x@annotation
                 lapply(annotation, function(x) {x$score})
               })
               names(score) <- index
               return(score)
             }

             if(slot == "annotation.level"){
               level <- lapply(tags2[index], function(x) {
                 annotation <- x@annotation
                 lapply(annotation, function(x) {x$level})
               })
               names(level) <- index
               return(level)
             }

             if(slot == "annotation.step"){
               step <- lapply(tags2[index], function(x) {
                 annotation <- x@annotation
                 lapply(annotation, function(x) {x$step})
               })
               names(step) <- index
               return(step)
             }

             if(slot == "seed.number"){
               seed.number <- lapply(tags2[index], function(x) {
                 annotation <- x@annotation
                 sum(unlist(lapply(annotation, function(x) {x$as.seed})))
               })
               seed.number <- unlist(seed.number)
               names(seed.number) <- index
               return(seed.number)
             }


             if(slot == "as.seed"){
               as.seed <- lapply(tags2[index], function(x) {
                 annotation <- x@annotation
                 unlist(lapply(annotation, function(x) {x$as.seed}))
               })


               seed.round <- lapply(tags2[index], function(x) {
                 annotation <- x@annotation
                 unlist(lapply(annotation, function(x) {x$as.seed.round}))
               })

               names(as.seed) <- index
               names(seed.round) <- index
               result <- list(as.seed, seed.round)
               names(result) <- c("seed", "round")
               return(result)
             }


           })



################################################################################
# filterPeak -------------------------------------------------------------------
setGeneric(
  name = "filterPeak",
  def  = function(object,...) {
    standardGeneric("filterPeak")
  }
)

# title filterPeak
# description Filter peakInfo data.
# author Xiaotao Shen
# \email{shenxt@@sioc.ac.cn}
# param object peakInfo data.
# param score.thr Score cutoff.
# return peakInfo data.


setMethod(f = "filterPeak",
          signature  = "PeakInfo",
          definition = function(object,
                                score_thr = 0) {
            annotation <- object@annotation
            if(length(annotation) != 0){
              score <- unlist(lapply(annotation, function(x) {x$score}))
              annotation <- annotation[order(score, decreasing = TRUE)]
              score <- unlist(lapply(annotation, function(x) {x$score}))
              if(length(which(score >= score_thr)) > 0){
                annotation <- annotation[which(score >= score_thr)]
                object@annotation <- annotation
              }else{
                object@annotation <- list()
              }

            }
            return(object)
          }
)


################################################################################
# MetDNA2Annotation ------------------------------------------------------------

setClass(Class = "MetDNA2Annotation",
         representation(peak_info = 'list',
                        initial_seed_annotation = 'data.frame',
                        recursive_annotation = 'data.frame',
                        credential_annotation = 'list'))



#' @title generateMetDNA2Annotation
#' @description Convert Ms1DataInfor to AnnotationResult class
#' @author Zhiwei Zhou
#' @param ms1_info
#' @export

# test <- generateMetDNA2Annotation(result_initial_annotation = result_initial_annotation)

setGeneric(name = 'generateMetDNA2Annotation',
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



setMethod(f = "show",
          signature = "MetDNA2Annotation",
          definition = function(object) {
            cat("----------------- Meta information ------------------\n")
            cat("Feature:", object@peak_info$name, "\n")
            cat("m/z:",object@peak_info$mz, "\n")
            cat("RT (s):",object@peak_info$rt, "\n")
            cat("CCS (A2):",object@peak_info$ccs, "\n")
            # cat("Median peak area:",object@peak_info$intmed, '\n\n')

            cat("----------- Initial seed annotation result ------------\n")

            if (nrow(object@initial_seed_annotation) == 0) {
              cat("No initial seed annotation.\n")
            } else {
              cat('# Initial seed annotation:', nrow(object@initial_seed_annotation), '\n')
            }

            cat("------------- Recursive annotation result --------------\n")

            if (nrow(object@recursive_annotation) == 0) {
              cat("No recursive annotation.\n")
            } else {
              cat('# Recursive annotation:', nrow(object@recursive_annotation), '\n')
            }

            cat("------------- Credential result --------------\n")

            if (nrow(object@credential_annotation$peak_group) == 0) {
              cat("No recursive annotation.\n")
            } else {
              cat('# Annotated peak group records:', nrow(object@credential_annotation$peak_group), '\n')
              cat('# Predicted formulas: ', nrow(object@credential_annotation$formula_prediction), '\n')
            }

          }
)
