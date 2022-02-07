#' @title mapProgress
#' @author Zhiwei Zhou
#' @param n
#' @export
#' @example
#' progress <- mapProgress(n = length(1:10))
#' temp <- purrr::map(seq_along(1:10), function(i){
#'   Sys.sleep(1)
#'   mapProgressPrint(progress = progress)
#'   return(NULL)
#' })

# progress <- mapProgress(n = 10)
setGeneric(name = 'mapProgress',
           def = function(n, ...){
             dplyr::progress_estimated(n = n,
                                       ...)
           })


#' @title mapProgressPrint
#' @author Zhiwei Zhou
#' @param n
#' @export
#' @example
#' progress <- mapProgress(n = length(1:10))
#' temp <- purrr::map(seq_along(1:10), function(i){
#'   Sys.sleep(1)
#'   mapProgressPrint(progress = progress)
#'   return(NULL)
#' })

# mapProgressPrint()
setGeneric(name = 'mapProgressPrint',
           def = function(progress){
             progress$tick()$print()

           })




setGeneric(name = 'getMzRange',
           def = function(mz,
                          ppm = 25,
                          mz_ppm_thr = 400){
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



#' @title SXTMTmatch
#' @description Match two data according to mz and RT.
#' @author Xiaotao Shen, Zhiwei Zhou
#' @param data1 First data for matching, first column must be mz and seconod column must be rt.
#' @param data2 Second data for matching, first column must be mz and seconod column must be rt.
#' @param mz.tol mz tol for ms1 and ms2 data matching. ppm
#' @param rt.tol RT tol for ms1 and ms2 data matching. s or %
#' @return Return a result which give the matching result of data1 and database.
# export

setGeneric(name = "SXTMTmatch",
           def = function(data1,
                          data2,
                          mz.tol = 25,
                          #rt.tol is relative
                          rt.tol = 10,
                          ccs.tol = NULL,
                          rt.error.type = c("relative", "abs")){
             rt.error.type <- match.arg(rt.error.type)

             if (nrow(data1) == 0 | nrow(data2) == 0) {
               result <- NULL
               return(result)
             }

             if (length(ccs.tol) == 0) {
               info1 <- data1[,c(1,2)]
               info1 <- apply(info1, 1, list)

               mz2 <- as.numeric(data2[, 1])
               rt2 <- as.numeric(data2[, 2])
             } else {
               info1 <- data1[,c(1,2,3)]
               info1 <- apply(info1, 1, list)

               mz2 <- as.numeric(data2[, 1])
               rt2 <- as.numeric(data2[, 2])
               ccs2 <- as.numeric(data2[, 3])
             }


             result <- pbapply::pblapply(info1, function(x) {
               temp.mz1 <- x[[1]][[1]]
               temp.rt1 <- x[[1]][[2]]
               mz.error <- abs(temp.mz1 - mz2) * 10 ^ 6 / temp.mz1
               if (rt.error.type == "relative") {
                 rt.error <- abs(temp.rt1 - rt2) * 100 / temp.rt1
               } else {
                 rt.error <- abs(temp.rt1 - rt2)
               }

               if (length(ccs.tol) == 0) {
                 j <- which(mz.error <= mz.tol & rt.error <= rt.tol)

                 if (length(j) == 0) {
                   matrix(NA, ncol = 7)
                 } else {
                   cbind(j, temp.mz1, mz2[j], mz.error[j], temp.rt1, rt2[j], rt.error[j])
                 }

               } else {
                 temp.ccs1 <- x[[1]][[3]]
                 ccs.error <- abs(temp.ccs1 - ccs2)*100/temp.ccs1

                 j <- which(mz.error <= mz.tol & rt.error <= rt.tol & ccs.error <= ccs.tol)

                 if (length(j) == 0) {
                   matrix(NA, ncol = 10)
                 } else {
                   cbind(j, temp.mz1, mz2[j], mz.error[j], temp.rt1, rt2[j], rt.error[j], temp.ccs1, ccs2[j], ccs.error[j])
                 }
               }

             })

             # add a column of number
             if (length(result) == 1) {
               result <- cbind(1,result[[1]])
             } else {
               result <- mapply(function(x,y){list(cbind(x,y))},
                                x <- 1:length(info1),
                                y = result)
               result <- do.call(rbind, result)
             }


             if (length(ccs.tol) == 0) {
               result <- matrix(result[which(!apply(result,1,function(x) any(is.na(x)))),], ncol = 8)
               if(nrow(result) == 0) return(NULL)
               colnames(result) <-
                 c("Index1",
                   "Index2",
                   "mz1",
                   "mz2",
                   "mz error",
                   "rt1",
                   "rt2",
                   "rt error")

               return(result)
             } else {
               result <- matrix(result[which(!apply(result,1,function(x) any(is.na(x)))),], ncol = 11)
               if(nrow(result) == 0) return(NULL)
               colnames(result) <-
                 c("Index1",
                   "Index2",
                   "mz1",
                   "mz2",
                   "mz error",
                   "rt1",
                   "rt2",
                   "rt error",
                   'ccs1',
                   'ccs2',
                   'ccs error')

               return(result)

             }

           })








################################################################################
# sumFormula -------------------------------------------------------------------

#' @title sumFormula
#' @description Combine metabolite and adduct as a new sum formula.
#' If there are no enough element to remove, return NA.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param formula The formula of metabolite.
#' @param adduct The adduct of metabolite.
#' @return  A sum formula.
#' @export

setGeneric(name = "sumFormula",
           def = function(formula = "C9H11NO2",
                          adduct = "M-H2O+H"){

             if(is.na(formula)) return(NA)
             if(is.na(adduct)) return(formula)
             if(adduct == "M+" | adduct == "M-" | adduct == 'M'){
               return(formula)
             }

             formula1 <- splitFormula(formula)
             adduct1 <- strsplit(x = adduct, split = "\\-|\\+")[[1]][-1]
             polymer <- as.numeric(gsub(pattern = "M", replacement = "",
                                        strsplit(x = adduct, split = "\\-|\\+")[[1]][1]))
             if (is.na(polymer)) polymer <- 1

             plusorminus <- strsplit(x = adduct, split = "")[[1]]
             plusorminus <- grep("\\+|\\-", plusorminus, value = TRUE)

             formula1$number <- formula1$number * polymer

             adduct1 <- mapply(function(x, y){
               temp <- splitFormula(x)
               temp$number <- temp$number * ifelse(y == "+", 1, -1)
               list(temp)
             },
             x = adduct1,
             y = plusorminus)

             adduct1 <- do.call(rbind, adduct1)

             formula <- rbind(formula1, adduct1)
             rownames(formula) <- NULL

             unique.element <- unique(formula$element.name)
             if(length(unique.element) == nrow(formula)){
               if(any(formula$number < 0)) {
                 return(NA)
               }else{
                 formula$number[formula$number==1] <- "W"
                 formula <- paste(paste(formula$element.name, formula$number, sep = ""), collapse = "")
                 formula <- strsplit(formula, split = "")[[1]]
                 formula[formula == "W"] <- ""
                 formula <- paste(formula, collapse = "")
                 return(formula)
               }
             }else{
               formula <- lapply(unique.element, function(x){
                 formula[formula$element.name == x,,drop = FALSE]
               })

               formula <- lapply(formula, function(x){
                 data.frame(unique(x$element.name), sum(x$number),
                            stringsAsFactors = FALSE)
               })

               formula <- do.call(rbind, formula)
               formula <- formula[formula[,2] != 0,]
               colnames(formula) <- c("element.name", "number")
               if(any(formula$number < 0)) {return(NA)}else{
                 formula$number[formula$number==1] <- "W"
                 formula <- paste(paste(formula$element.name, formula$number, sep = ""), collapse = "")
                 formula <- strsplit(formula, split = "")[[1]]
                 formula[formula == "W"] <- ""
                 formula <- paste(formula, collapse = "")
                 return(formula)
               }
             }
           })


#' @title splitFormula
#' @description Split a formula into element and number.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param formula The formula of metabolite.
#' @return  A splited formula.
#' @export

setGeneric(name = "splitFormula",
           def = function(formula = "C9H11NO2"){
             temp.formula <- strsplit(formula, split = "")[[1]]

             number <- NULL
             for(i in 1:length(temp.formula)){
               if(length(grep("[0-9]{1}", temp.formula[i])) == 0){break}
               number[i] <- temp.formula[i]
             }

             if(!is.null(number)) {
               number <- as.numeric(paste(number, collapse = ""))
             }else{
               number <- 1
             }
             ##first select the Na, Cl and so on element
             idx1 <- gregexpr("[A-Z][a-z][0-9]*", formula)[[1]]
             len1 <- attributes(idx1)$match.length
             ##no double element
             if(idx1[1] == -1) {
               double.formula <- matrix(NA, ncol = 2)
               formula1 <- formula
             }else{
               double.letter.element <- NULL
               double.number <- NULL
               remove.idx <- NULL
               for (i in 1:length(idx1)) {
                 double.letter.element[i] <- substr(formula, idx1[i], idx1[i] + len1[i] - 1)
                 if(nchar(double.letter.element[i]) == 2){
                   double.number[i] <- 1
                 }else{
                   double.number[i] <- as.numeric(substr(double.letter.element[i], 3, nchar(double.letter.element[i])))
                 }
                 double.letter.element[i] <- substr(double.letter.element[i], 1, 2)
                 remove.idx <- c(remove.idx, idx1[i] : (idx1[i] + len1[i] - 1))
               }

               double.formula <- data.frame(double.letter.element,
                                            double.number, stringsAsFactors = FALSE)
               formula1 <- strsplit(formula, split = "")[[1]]
               formula1 <- formula1[-remove.idx]
               formula1 <- paste(formula1, collapse = "")
             }

             ## no one element
             if(formula1 == ""){
               one.formula <- matrix(NA, ncol = 2)
             }else{
               idx2 <- gregexpr("[A-Z][0-9]*", formula1)[[1]]
               len2 <- attributes(idx2)$match.length
               one.letter.element <- NULL
               one.number <- NULL
               for (i in 1:length(idx2)) {
                 one.letter.element[i] <- substr(formula1, idx2[i], idx2[i] + len2[i] - 1)
                 if(nchar(one.letter.element[i]) == 1){
                   one.number[i] <- 1
                 }else{
                   one.number[i] <- as.numeric(substr(one.letter.element[i], 2, nchar(one.letter.element[i])))
                 }
                 one.letter.element[i] <- substr(one.letter.element[i], 1, 1)
               }
               one.formula <- data.frame(one.letter.element, one.number,
                                         stringsAsFactors = FALSE)
             }

             colnames(double.formula) <- colnames(one.formula) <- c("element.name","number")
             formula <- rbind(double.formula, one.formula)
             formula <- formula[!apply(formula, 1, function(x) any(is.na(x))),]

             formula <- formula[order(formula$element.name),]
             formula$number <- formula$number * number
             unique.element <- unique(formula$element.name)
             if(length(unique.element) == nrow(formula)){
               return(formula)
             }else{
               formula <- lapply(unique.element, function(x){
                 formula[formula$element.name == x,,drop = FALSE]
               })

               formula <- lapply(formula, function(x){
                 data.frame(unique(x$element.name), sum(x$number),
                            stringsAsFactors = FALSE)
               })

               formula <- do.call(rbind, formula)
               colnames(formula) <- c("element.name", "number")
               return(formula)
             }
           })


#' @title pasteElement
#' @description Paste formula and element.
#' Combine metabolite and adduct as a new sum formula.
#' If there are no enough element to remove, return NA.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param formula The formula of metabolite.
#' @param element The element.
#' @param mode Add or remove a module
#' @export
#' @return  A formula.


pasteElement <- function(formula = "C9H11NO2",
                         element = "H",
                         mode = c("plus", "minus")){

  mode <- match.arg(mode)
  formula <- splitFormula(formula = formula)
  element <- splitFormula(formula = element)


  ## mode = plus
  if(mode == "plus"){
    for (i in 1:nrow(element)){
      temp.name <- as.character(element[i,1])
      temp.number <- as.numeric(element[i,2])
      temp.idx <- match(temp.name, formula[,1])
      if(is.na(temp.idx)) {
        formula <- rbind(formula, element[i,])
      }else{
        formula[temp.idx, 2] <- formula[temp.idx, 2] + temp.number
      }
    }
  }else{
    for (i in 1:nrow(element)){
      temp.name <- as.character(element[i,1])
      temp.number <- as.numeric(element[i,2])
      temp.idx <- match(temp.name, formula[,1])
      if(is.na(temp.idx)) {
        # warning("Formula has no element in adduct!\n")
        return(NA)
      }else{
        formula[temp.idx,2] <- formula[temp.idx,2] - temp.number
        if(formula[temp.idx,2] < 0) {
          # warning("Formula has no enough element in adduct!\n")
          return(NA)
        }
      }
    }
  }

  ###return formula
  formula <- as.data.frame(formula)
  formula <- formula[formula[,2] != 0, , drop = FALSE]
  formula <- c(t(formula))
  formula <- gsub(pattern = " ", replacement = "", x = formula)
  formula <- formula[formula != "1"]
  formula <- paste(formula, collapse = "")
  return(formula)
}





#' @title checkElement
#' @description Check a formula can add one adduct or not.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param formula The formula of metabolite.
#' @param adduct The adduct.
#' @return  valid, return TRUE; invalid.
#' @export

# node.adduct <- adduct_table$name
# node.formula <- "CH2O"
# sapply(node.adduct, function(x) {cat(x, '\n'); checkElement(node.formula, x)})
# adduct <- 'M-'
# formula <- 'CH2O'

setGeneric(name = "checkElement",
           def = function(formula = "C9H11NO2",
                          adduct = "M-H2O+H"){
             formula1 <- splitFormula(formula)
             # fix checkElement for adducts 'M+' and 'M-', 'M'
             if (adduct %in% c('M+', 'M-', 'M')) {
               return(TRUE)
             }
             adduct1 <- strsplit(x = adduct, split = "\\-|\\+")[[1]][-1]
             plusorminus <- strsplit(x = adduct, split = "")[[1]]
             plusorminus <- grep("\\+|\\-", plusorminus, value = TRUE)
             if(all(plusorminus == "+")) return(TRUE)

             adduct1 <- mapply(function(x, y){
               temp <- splitFormula(x)
               temp$number <- temp$number * ifelse(y == "+", 1, -1)
               list(temp)
             },
             x = adduct1,
             y = plusorminus)

             adduct1 <- do.call(rbind, adduct1)

             formula <- rbind(formula1, adduct1)
             rownames(formula) <- NULL

             unique.element <- unique(formula$element.name)
             if(length(unique.element) == nrow(formula)){
               if(any(formula$number < 0)) {return(FALSE)}else{return(TRUE)}
             }else{
               formula <- lapply(unique.element, function(x){
                 formula[formula$element.name == x,,drop = FALSE]
               })

               formula <- lapply(formula, function(x){
                 data.frame(unique(x$element.name), sum(x$number),
                            stringsAsFactors = FALSE)
               })

               formula <- do.call(rbind, formula)
               colnames(formula) <- c("element.name", "number")
               if(any(formula$number < 0)) {return(FALSE)}else{return(TRUE)}
             }
           })




#' @title generateMSP
#' @author Zhiwei Zhou
#' @description Convert MS/MS spectra to MSP files
#' @param file_name Required. The file name of MSP. The suffix of file name must be ".msp".
#' @param cmp_name Required
#' @param precusormz Required
#' @param spec Required. It should be in a dataframe, 1st: "mz", 2nd: "intensity"
#' @param adduct Default: NULL
#' @param instrument_type Default: NULL
#' @param instrument Default: NULL
#' @param smiles Default: NULL
#' @param inchikey Default: NULL
#' @param inchikey1 Default: NULL
#' @param formula Default: NULL
#' @param polarity 'Positive' or 'Negative'
#' @param ce Default: NULL
#' @param rt Default: NULL
#' @param ccs Default: NULL
#' @param zhulib_id Default: NULL
#' @param kegg_id Default: NULL
#' @param hmdb_id Default: NULL
#' @param pubchem_id Default: NULL
#' @param links Default: ''
#' @param comment Default: ''
#' @example
#' setwd('F:/01 MetIMMS/00 data processing/190515 external validation msms data extraction/')
#'     GenerateMSP(file_name = 'zhumetlib_validation_pos_20v_190520.msp',
#'                 cmp_name = external_data_pos$compound_name[i],
#'                 precusormz = external_data_pos$mz[i],
#'                 adduct = external_data_pos$adducts[i],
#'                 instrument_type = 'LC-ESI-qTof',
#'                 instrument = 'LC-ESI-qTof',
#'                 smiles = external_data_pos$smiles[i],
#'                 inchikey = external_data_pos$inchikey[i],
#'                 inchikey1 = external_data_pos$inchikey1[i],
#'                 formula = external_data_pos$formula[i],
#'                 polarity = 'Positive',
#'                 ce = '20',
#'                 ccs = external_data_pos$CCS[i],
#'                 zhulib_id = external_data_pos$matched_zhulib_id[i],
#'                 pubchem_id = external_data_pos$pubchem_cid[i],
#'                 comment = paste(external_data_pos$id[i], external_data_pos$source[i], sep = ' '),
#'                 spec = temp_spec)



setGeneric(name = 'generateMSP',
           def = function(
             file_name = "./test.msp",
             cmp_name = NULL,
             precusormz = NULL,
             adduct=NULL,
             instrument_type=NULL,
             instrument=NULL,
             smiles=NULL,
             inchikey=NULL,
             inchikey1=NULL,
             formula=NULL,
             polarity=c('positive', 'negative'),
             ce=NULL,
             rt=NULL,
             ccs=NULL,
             zhulib_id=NULL,
             kegg_id=NULL,
             hmdb_id=NULL,
             pubchem_id=NULL,
             links='',
             comment='',
             spec=NULL
           ){

             if (!stringr::str_detect(file_name, '.msp')) {
               stop('The suffix of file name must be .msp\n')
             }

             if (is.null(cmp_name)) {
               stop('Please input cmp_name\n')
             }

             if (is.null(precusormz)) {
               stop('Please input precusormz\n')
             }

             if (is.null(spec)) {
               stop('Please input spec\n')
             }

             polarity = match.arg(polarity)


             if (ncol(spec)!=2) {
               stop('Please check the format of spectrum. It should be in a dataframe, 1st: "mz", 2nd: "intensity" \n')
             }

             if (!all(colnames(spec)==c('mz', 'intensity'))) {
               stop('Please check the format of spectrum. It should be in a dataframe, 1st: "mz", 2nd: "intensity" \n')
             }



             # write into msp
             file_result <- file(description = file_name, open = "a")
             cat('NAME: ', cmp_name, '\n', sep = '', file = file_result)
             cat('PRECURSORMZ: ', precusormz, '\n', sep = '', file = file_result)

             if (!is.null(adduct)) {
               cat('PRECURSORTYPE: ', adduct, '\n', sep = '', file = file_result)
             }

             if (!is.null(instrument_type)) {
               cat('INSTRUMENTTYPE: ', instrument_type, '\n', sep = '', file = file_result)
             }

             if (!is.null(instrument)) {
               cat('INSTRUMENT: ', instrument, '\n', sep = '', file = file_result)
             }

             if (!is.null(smiles)) {
               cat('SMILES: ', smiles, '\n', sep = '', file = file_result)
             }

             if (!is.null(inchikey)) {
               cat('INCHIKEY: ', inchikey, '\n', sep = '', file = file_result)
             }

             if (!is.null(inchikey1)) {
               cat('INCHIKEY1: ', inchikey1, '\n', sep = '', file = file_result)
             }

             if (!is.null(formula)) {
               cat('FORMULA: ', formula, '\n', sep = '', file = file_result)
             }

             if (!is.null(polarity)) {
               cat('IONMODE: ', polarity, '\n', sep = '', file = file_result)
             }

             if (!is.null(ce)) {
               cat('COLLISIONENERGY: ', ce, '\n', sep = '', file = file_result)
             }

             if (!is.null(rt)) {
               cat('RETENTIONTIME: ', rt, '\n', sep = '', file = file_result)
             }

             if (!is.null(ccs)) {
               cat('COLLISIONCROSSSECTION: ', ccs, '\n', sep = '', file = file_result)
             }

             if (!is.null(zhulib_id)) {
               cat('ZHULAB: ', zhulib_id, '\n', sep = '', file = file_result)
             }

             if (!is.null(kegg_id)) {
               cat('KEGG: ', kegg_id, '\n', sep = '', file = file_result)
             }

             if (!is.null(hmdb_id)) {
               cat('HMDB: ', hmdb_id, '\n', sep = '', file = file_result)
             }

             if (!is.null(pubchem_id)) {
               cat('PUBCHEM: ', pubchem_id, '\n', sep = '', file = file_result)
             }

             cat('Links: ', links, '\n', sep = '', file = file_result)
             cat('Comment: ', comment, "\n", sep = '', file = file_result)

             cat('Num Peaks: ',  nrow(spec),  '\n',  sep = '', file = file_result)

             for (i in 1:nrow(spec)) {
               cat(paste(as.numeric(round(spec[i,1], digits = 4)),
                         as.numeric(round(spec[i,2], digits = 2)),
                         collapse = ' '),
                   '\n', sep = '', file = file_result)
             }

             cat('\n', file = file_result)

             close(file_result)
           }
)


#' @title getPostfix
#' @author Zhiwei Zhou
#' @description get the postfix of file
#' @param x file_name
#' @export
#' @examples
#' getPostfix('abc.mgf')
# getPostfix('abc.mgf')

setGeneric(name = "getPostfix", def = function(x){
  unlist(lapply(strsplit(x = x, split = "\\."),  function(y) {
    if(length(y) == 1) return(NA)
    y[[length(y)]]}
  ))
})


#' @title splitInchiKey
#' @author Zhiwei Zhou
#' @description split inchikey to 3 parts
#' @return a data.frame with 4 columns
#' @export
#' @example
#' splitInchiKey(inchikey = 'VGONTNSXDCQUGY-RRKCRQDMSA-N')

setGeneric(name = 'splitInchiKey',
           def = function(
             inchikey = 'AEMRFAOFKBGASW-UHFFFAOYSA-N'
           ){
             result <- stringr::str_split_fixed(inchikey, pattern = '-', n = 3)
             result <- data.frame(result, stringsAsFactors = F)
             colnames(result) <- c('inchikey1', 'inchikey2', 'inchikey3')

             result <- data.frame(result,
                                  inchikey = inchikey,
                                  stringsAsFactors = F)

             return(result)
           }
)




#' @title 'retieveUpperPath'
#' @author Zhiwei Zhou
#' @param path
#' @examples
#' retieveUpperPath("/home/zhouzw/Data_processing/20210626_metdna2_merge_pos_neg/BOTH/")
#' @export

setGeneric(name = 'retieveUpperPath',
           function(
             path = "/home/zhouzw/Data_processing/20210626_metdna2_merge_pos_neg/BOTH/"
           ){
             path <- strsplit(path, split = "/")[[1]][-length(strsplit(path, split = "/")[[1]])] %>%
               paste(collapse = "/")

             return(path)

           })


#' @title removeFiles
#' @author Zhiwei Zhou
#' @param files
#' @param path
#' @export

setGeneric(name = 'removeFiles',
           def = function(files, path){
             if ((length(files) != length(path)) & (length(path) != 1)) {
               stop('The length of path should equal 1 or be same as length of files\n')
             }

             if (length(path) == 1) {
               purrr::walk(seq_along(files), function(i){
                 try(unlink(x = file.path(path, files[i]), recursive = TRUE), silent = TRUE)
               })
             }

             if (length(files) == length(path)) {
               purrr::walk(seq_along(files), function(i){
                 try(unlink(x = file.path(path[i], files[i]), recursive = TRUE), silent = TRUE)
               })
             }

})
