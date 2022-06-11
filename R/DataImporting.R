################################################################################
# checkQuality -----------------------------------------------------------------

#' @title checkQuality
#' @description Check ms1_file.
#' @author Xiaotao Shen, Zhiwei Zhou
#' @param ms1_file The name of ms1 ms1_file.
#' @param sample_info_file The name of sample information.
#' @param ms2_type The type of MS2 ms1_file.
#' @param path The work directory.
#' @return Check result.
#' @export
#'

# checkQuality(ms1_file = "data.csv",
#              sample_info_file = "sample.info.csv",
#              ms2_type = "mgf",
#              path = "/home/zhouzw/Data_processing/20210522_debug/10_nist_urine_neg_emrn0_without_mz_segment/")
#
# ms1_file <- "data.csv"
# sample_info_file <- 'sample_info_file.csv'
# ms2_type <- 'mgf'
# path <- '/home/zhouzw/Data_processing/20201102_metdna2_pos_development'

setGeneric(name = "checkQuality",
           def = function(ms1_file = "data.csv",
                          sample_info_file = "sample.info.csv",
                          ms2_type = c("mgf", "mzXML", "msp", 'cef'),
                          instrument = c("SciexTripleTOF", "AgilentQTOF", "BrukerQTOF", "ThermoOrbitrap",
                                         'ThermoExploris', "AgilentDTIMMS", "BrukerTIMS", 'WatersQTOF','WatersTWIMMS'),
                          path = "."){

             ms1_file_name <- ms1_file
             sample_info_file_name <- sample_info_file
             ms1_file_record <- NULL
             sample_info_file_record <- NULL
             ms2_file_record <- NULL

             ms2_type <- match.arg(ms2_type)
             cat("Read ms1_file.\n")
             cat("--------------------------------------------------------------\n")

             temp_error <- showError(
               ms1_data <- readr::read_csv(file.path(path, ms1_file), col_types = readr::cols(),
                                           progress = TRUE),
               error_info = paste("Error: There is no", ms1_file, "in", path, ". Please check it.\n"))
             if (class(temp_error) == "try-error") {
               cat(paste("Error: There is no", ms1_file, "in", path, ". Please check it.\n"),
                   file = file.path(path, "run.log.txt"), append = TRUE)
               stop("Error")
             }
             ms1_data <- as.data.frame(ms1_data)

             temp_error <- showError(
               sample_info_data <- readr::read_csv(file.path(path, sample_info_file), col_types = readr::cols(),
                                                   progress = TRUE),
               error_info = paste("Error: There is no", sample_info_file, "in", path, ". Please check it.\n"))
             if (class(temp_error) == "try-error") {
               cat(paste("Error: There is no", sample_info_file, "in", path, ". Please check it.\n"),
                   file = file.path(path, "run.log.txt"), append = TRUE)
               stop("Error")
             }
             sample_info_data <- as.data.frame(sample_info_data)

             # check the sample_info_file have - or not, if it has -, change it to .s
             if (length(grep(pattern = "-", sample_info_data$sample.name)) > 0){
               colnames(ms1_data) <- gsub(pattern = "-", replacement = ".", x = colnames(ms1_data))
               sample_info_data$sample.name <- gsub(pattern = "-", replacement = ".", x = sample_info_data$sample.name)
               readr::write_csv(x = ms1_data, path = file.path(path, ms1_file_name))
               readr::write_csv(x = sample_info_data, path = file.path(path, sample_info_file_name))
             }

             # check ms1_file, NA or space
             if(sum(is.na(ms1_data)) > 0){
               cat("Error: There are", sum(is.na(ms1_data)), "NAs in you ms1_file.\n")
               ms1_file_record <- c(ms1_file_record, "Error")
             }else{
               # cat("OK: There are no NAs in you ms1_file.\n")
               ms1_file_record <- c(ms1_file_record, "OK")
             }

             if(ifelse(is.na(sum(ms1_data == "") > 0), FALSE, sum(ms1_data == "") > 0)){
               cat("Error: There are", sum(ms1_data == ""), "spaces in you ms1_file.\n")
               ms1_file_record <- c(ms1_file_record, "Error")
             }else{
               # cat("OK: There are no spaces in you ms1_file.\n")
               ms1_file_record <- c(ms1_file_record, "OK")
             }

             # check non-alphanumeric character
             if (any(stringr::str_detect(ms1_data$name, '[^\\w]'))) {
               # cat('Error: non-alphanumeric character existed in name column in ms1_file')
               # ms1_file_record <- c(ms1_file_record, "Error")

               cat(paste('Error: non-alphanumeric character existed in name column in ms1_file'),
                   file = file.path(path, "run.log.txt"), append = TRUE)
               stop('Error: non-alphanumeric character existed in name column in ms1_file')
             }


             # impute NA and space as 0
             ms1_data1 <- ms1_data
             ms1_data1[is.na(ms1_data1)] <- 0
             ms1_data1[ms1_data1 == ""] <- 0


             temp_error <- showError(
               ms1_data1 <- ms1_data1[,setdiff(colnames(ms1_data1), c("mz", "rt", "name"))],
               error_info = "Error: The ms1_data don't contains 'name', 'mz' and 'rt', please check it.")
             if (class(temp_error) == "try-error") {
               cat("Error: The ms1_data don't contains 'name', 'mz' and 'rt', please check it.",
                   file = file.path(path, "run.log.txt"), append = TRUE)
               stop("Error")
             }

             temp_error <- showError(
               zero.per <- apply(ms1_data1, 1, function(x){
                 sum(x == 0)/ncol(ms1_data1)
               }),
               error_info = "Error: The ms1_data contains other non-numeric columns, please check it.")
             if(class(temp_error) == "try-error") {
               if (class(temp_error) == "try-error") {
                 cat("Error: The ms1_data contains other non-numeric columns, please check it.",
                     file = file.path(path, "run.log.txt"), append = TRUE)
                 stop("Error")
               }
             }


             if(sum(zero.per == 1) > 0){
               cat("Error: The ms1_data contain some feature with all zero peak-area, please check it.",
                   file = file.path(path, "run.log.txt"), append = TRUE)
               stop("Error: The ms1_data contain some feature with all zero peak-area, please check it.")
             }

             if(sum(zero.per > 0.5) > 0){
               cat("Warning: There are peaks with zero ratio > 50% in you ms1_file.\n")
               ms1_file_record <- c(ms1_file_record, "Warning")
             }

             # check ms1_data component
             ms1_data_col_name <- colnames(ms1_data)
             if (length(ms1_data_col_name) < 7){
               cat("Warning: There are less than 4 samples in your ms1_data, please check it.\n")
               ms1_file_record <- c(ms1_file_record, "Warning")
             } else {
               cat("There are", length(ms1_data_col_name) - 3, "samples in your ms1_data.\n")
               ms1_file_record <- c(ms1_file_record, "OK")
             }


             if (!(instrument %in% c("AgilentDTIMMS", "BrukerTIMS", "WatersTWIMMS"))) {
               if (ms1_data_col_name[1] != "name"){
                 cat("Error: The first column of ms1_file is not 'name'.\n")
                 ms1_file_record <- c(ms1_file_record, "Error")
               } else {
                 # cat("OK: The first column of ms1_file is name.\n")
                 ms1_file_record <- c(ms1_file_record, "OK")
               }

               if(ms1_data_col_name[2] != "mz"){
                 cat("Error: The second column of ms1_file is not 'mz'.\n")
                 ms1_file_record <- c(ms1_file_record, "Error")
               } else {
                 # cat("OK: The second column of ms1_file is mz.\n")
                 ms1_file_record <- c(ms1_file_record, "OK")
               }

               if(ms1_data_col_name[3] != "rt"){
                 cat("Error: The third column of ms1_file is not 'rt'.\n")
                 ms1_file_record <- c(ms1_file_record, "Error")
               } else {
                 # cat("OK: The third column of ms1_file is not rt.\n")
                 ms1_file_record <- c(ms1_file_record, "OK")
               }
             } else {
               if (ms1_data_col_name[1] != "mz"){
                 cat("Error: The first column of ms1_file is not 'mz'.\n")
                 ms1_file_record <- c(ms1_file_record, "Error")
               } else {
                 # cat("OK: The first column of ms1_file is name.\n")
                 ms1_file_record <- c(ms1_file_record, "OK")
               }

               if(ms1_data_col_name[2] != "rt"){
                 cat("Error: The second column of ms1_file is not 'rt'.\n")
                 ms1_file_record <- c(ms1_file_record, "Error")
               } else {
                 # cat("OK: The second column of ms1_file is mz.\n")
                 ms1_file_record <- c(ms1_file_record, "OK")
               }

               if(ms1_data_col_name[3] != "ccs"){
                 cat("Error: The third column of ms1_file is not 'ccs'.\n")
                 ms1_file_record <- c(ms1_file_record, "Error")
               } else {
                 # cat("OK: The third column of ms1_file is not rt.\n")
                 ms1_file_record <- c(ms1_file_record, "OK")
               }
             }

             # check ms1_file, RT
             if (any(ms1_data_col_name == "rt")){
               rt <- as.numeric(ms1_data$rt)

               if (!(instrument %in% c("AgilentDTIMMS", "BrukerTIMS", "WatersTWIMMS"))) {
                 if (max(rt) < 60){
                   cat("Warning: Please confirm that the rt is measured in seconds,\n")
                   ms1_file_record <- c(ms1_file_record, "Warning")
                 }
               } else {
                 if (max(rt) > 60){
                   cat("Warning: Please confirm that the rt is measured in minutes,\n")
                   ms1_file_record <- c(ms1_file_record, "Warning")
                 }
               }

               cat("RT range is", range(rt), ".\n")
             }

             cat("--------------------------------------------------------------\n")
             # check sample_info_data
             if (ncol(sample_info_data) > 2){
               cat("Error: The smple_info has more than two columns. Please check it.\n")
               sample_info_file_record <- c(sample_info_file_record, "Error")
             }

             if (colnames(sample_info_data)[1] != "sample.name"){
               cat("Error: The first column name of sample_info_file must be 'sample.name'. Please check it.\n")
               sample_info_file_record <- c(sample_info_file_record, "Error")
             } else {
               sample_info_file_record <- c(sample_info_file_record, "OK")
             }

             if (colnames(sample_info_data)[2] != "group"){
               cat("Error: The second column name of sample_info_file must be 'group'. Please check it.\n")
               sample_info_file_record <- c(sample_info_file_record, "Error")
             } else {
               sample_info_file_record <- c(sample_info_file_record, "OK")
             }

             if (sum(is.na(sample_info_data)) > 0){
               cat("Error: There are", sum(is.na(sample_info_data)), "NAs in you sample_info_file.\n")
               sample_info_file_record <- c(sample_info_file_record, "Error")
             } else {
               # cat("OK: There are no NAs in you sample_info_file.\n")
               sample_info_file_record <- c(sample_info_file_record, "OK")
             }

             if (ifelse(is.na(sum(sample_info_data == "") > 0), FALSE, sum(sample_info_data == "") > 0)){
               cat("Error: There are", sum(sample_info_file == ""), "spaces in you sample_info_file.\n")
               sample_info_file_record <- c(sample_info_file_record, "Error")
             } else {
               # cat("OK: There are no spaces in you sample_info_file.\n")
               sample_info_file_record <- c(sample_info_file_record, "OK")
             }


             group <-  table(as.character(sample_info_data[,2]))
             group1 <- matrix(group, nrow = 1)
             colnames(group1) <- names(group)
             rownames(group1) <- "Number"

             cat("Group information:\n")
             print(group1)

             if(any(group1[1,] < 2)){
               cat("Warning: ","Group", group1[1,which(group1[,1] < 3)],
                   ifelse(length(which(group1[,1] < 3)) > 1, "have", "has"),
                   "less than 3 samples, please check it.\n"
               )
               sample_info_file_record <- c(sample_info_file_record, "Warning")
             }else{
               sample_info_file_record <- c(sample_info_file_record, "OK")
             }


             if (!(instrument %in% c("AgilentDTIMMS", "BrukerTIMS", "WatersTWIMMS"))) {
               ms1_data_col_name <- setdiff(ms1_data_col_name, c("name", "mz", "rt"))
             } else {
               ms1_data_col_name <- setdiff(ms1_data_col_name, c("mz", "rt", "ccs"))
             }

             sample_name <- as.character(sample_info_data[,1])
             if (any(sort(ms1_data_col_name) != sort(sample_name))){
               cat("Error: The sample names in ms1_file and sample_info_file are not same.\n")
               sample_info_file_record <- c(sample_info_file_record, "Error")
             } else {
               sample_info_file_record <- c(sample_info_file_record, "OK")
             }

             # if the name of sample is begin with number, add "Sample" before it
             if (!(instrument %in% c("AgilentDTIMMS", "BrukerTIMS", "WatersTWIMMS"))) {
               ms1_data_col_name <- setdiff(ms1_data_col_name, c("name", "mz", "rt"))
             } else {
               ms1_data_col_name <- setdiff(ms1_data_col_name, c("mz", "rt", "ccs"))
             }
             sample_name <- as.character(sample_info_data[,1])

             if (any(substr(x = sample_name, start = 1, stop = 1) %in% 0:9)){
               cat("Warning: The names of samples should not start with a number.\n")
               sample_info_file_record <- c(sample_info_file_record, "Warning")

               idx1 <- which(substr(x = sample_name, start = 1, stop = 1) %in% 0:9)
               idx2 <- which(substr(x = ms1_data_col_name, start = 1, stop = 1) %in% 0:9)

               sample_name[idx1] <- paste("Sample", sample_name[idx1], sep = "")
               ms1_data_col_name[idx2] <- paste("Sample", ms1_data_col_name[idx2], sep = "")

               sample_info_data[,1] <- sample_name
               if (!(instrument %in% c("AgilentDTIMMS", "BrukerTIMS", "WatersTWIMMS"))) {
                 colnames(ms1_data)[-which(colnames(ms1_data) %in% c("name", "mz", "rt"))] <- ms1_data_col_name
               } else {
                 colnames(ms1_data)[-which(colnames(ms1_data) %in% c("mz", "rt", "ccs"))] <- ms1_data_col_name
               }
               readr::write_csv(ms1_data, file.path(path, "data.csv"))
               readr::write_csv(sample_info_data, file.path(path, "sample.info.csv"))

             } else {
               sample_info_file_record <- c(sample_info_file_record, "OK")
             }

             cat("--------------------------------------------------------------\n")
             # ms2.file
             # switch (ms2_type,
             #         "mgf" = {ms2_file <- grep("mgf", dir(path), value = TRUE)},
             #         "msp" = {ms2_file <- grep("msp", dir(path), value = TRUE)},
             #         "mzXML" = {ms2_file <- grep("mzXML", dir(path), value = TRUE)},
             #         "cef" = {ms2_file <- grep("cef", dir(path), value = TRUE)})

             switch(ms2_type,
                    'mgf' = {
                      ms2_file <- list.files(path = path,
                                             pattern = '.+\\.mgf$',
                                             recursive = TRUE,
                                             ignore.case = TRUE)
                      # ms2_file <- file.path(path, ms2_file)
                    },
                    'msp' = {
                      ms2_file <- list.files(path = path,
                                             pattern = '.+\\.msp$',
                                             recursive = TRUE,
                                             ignore.case = TRUE)
                      # ms2_file <- file.path(path, ms2_file)
                    },
                    'mzXML' = {
                      ms2_file <- list.files(path = path,
                                             pattern = '.+\\.mzXML$',
                                             recursive = TRUE,
                                             ignore.case = TRUE)
                      # ms2_file <- file.path(path, ms2_file)
                    },
                    'cef' = {
                      ms2_file <- list.files(path = path,
                                             pattern = '.+\\.cef$',
                                             recursive = TRUE,
                                             ignore.case = TRUE)
                      # ms2_file <- file.path(path, ms2_file)
                    })

             if (length(ms2_file) == 0){
               cat("Error: There are no ms2 ms1_file.\n")
               ms2_file_record <- c(ms2_file_record, "Error")
             } else {
               cat("There are", length(ms2_file), "ms2 ms1_file.\n")
               ms2_file_record <- c(ms2_file_record, "OK")
               cat("Total size of ms2 ms1_file is", sum(file.size(file.path(path, ms2_file))/(1024*1024)), "M.\n")
             }

             cat("--------------------------------------------------------------\n")
             cat("Summary:\n")
             stat <- lapply(list(ms1_file_record, sample_info_file_record, ms2_file_record),
                            function(x) {
                              return(c(ifelse(all(x!="Error"), "Valid", "Invalid"),
                                       sum(x == "OK"), sum(x == "Warning"), sum(x == "Error")))
                            })
             stat <- do.call(rbind, stat)
             colnames(stat) <- c("Check result","OK", "Warning", "Error")
             rownames(stat) <- c("ms1_data", "sample_info_data", "ms2_data")

             print(stat, quote = FALSE)
             cat("\n")

             cat("\n")
             cat("ms1_file:\n")
             if(all(ms1_file_record != "Error")){
               cat("ms1_file is valid.\n")
             }else{
               if (sum(ms1_file_record == "Warning") > 0){
                 cat("There", ifelse(sum(ms1_file_record == "Warning") > 1, "are", "is"),
                     sum(ms1_file_record == "Warning"),
                     ifelse(sum(ms1_file_record == "Warning") > 1, "Warnings", "Warning"),
                     "in your ms1_file. Please check it according to the information.\n")
               }

               if (sum(ms1_file_record == "Error") > 0){
                 cat("There", ifelse(sum(ms1_file_record == "Error") > 1, "are", "is"),
                     sum(ms1_file_record == "Error"),
                     ifelse(sum(ms1_file_record == "Error") > 1, "Errors", "Error"),
                     "in your ms1_file. Please check it according to the information.\n")
               }

             }



             cat("\n")
             cat("sample_info_file:\n")
             if (all(sample_info_file_record != "Error")){
               cat("sample_info_file is valid.\n")
             } else {
               if(sum(sample_info_file_record == "Warning") > 0){
                 cat("There", ifelse(sum(sample_info_file_record == "Warning") > 1, "are", "is"),
                     sum(sample_info_file_record == "Warning"),
                     ifelse(sum(sample_info_file_record == "Warning") > 1, "Warnings", "Warning"),
                     "in your sample_info_file. Please check it according to the information.\n")
               }

               if(sum(sample_info_file_record == "Error") > 0){
                 cat("There", ifelse(sum(sample_info_file_record == "Error") > 1, "are", "is"),
                     sum(sample_info_file_record == "Error"),
                     ifelse(sum(sample_info_file_record == "Error") > 1, "Errors", "Error"),
                     "in your sample_info_file. Please check it according to the information.\n")
               }

             }

             cat("\n")
             cat("ms2_file:\n")
             if(all(ms2_file_record != "Error")){
               cat("ms2_file are valid.\n")
             }else{
               if(sum(ms2_file_record == "Warning") > 0){
                 cat("There", ifelse(sum(ms2_file_record == "Warning") > 1, "are", "is"),
                     sum(ms2_file_record == "Warning"),
                     ifelse(sum(ms2_file_record == "Warning") > 1, "Warnings", "Warning"),
                     "in your ms2_file. Please check it according to the information.\n")
               }

               if(sum(ms2_file_record == "Error") > 0){
                 cat("There", ifelse(sum(ms2_file_record == "Error") > 1, "are", "is"),
                     sum(ms2_file_record == "Error"),
                     ifelse(sum(ms2_file_record == "Error") > 1, "Errors", "Error"),
                     "in your ms2.file. Please check it according to the information.\n")
               }

             }
             stat <- stat
           })



#   checkQualitySampleGroup -----------------------------------------------------
#' @title checkQualitySampleGroup
#' @author Zhiwei Zhou
#' @param ms1_file "data.csv"
#' @param sample_info_file "sample.info.csv"
#' @param path '.'

setGeneric(name = 'checkQualitySampleGroup',
           def = function(
             ms1_file = "data.csv",
             sample_info_file = "sample.info.csv",
             group,
             path = '.'
           ){
             old_sample_info_file_name <- sample_info_file
             old_ms1_file <- ms1_file
             temp_error <- showError(
               sample_info_file <- read.csv(file.path(path, sample_info_file), stringsAsFactors = FALSE),
               error_info = paste("Error: There is no", sample_info_file, "in", path, ". Please check it.\n"))

             if(class(temp_error) == "try-error") {
               cat(paste("Error: There is no", sample_info_file, "in", path, ". Please check it.\n"),
                   file = file.path(path, "run.log.txt"), append = TRUE)
               stop("Error")
             }

             if (missing(group)) {
               cat("You don't provide the group names.\n")
               cat("You don't provide the group names.\n", file = file.path(path, "run.log.txt"), append = TRUE)
               group <- sort(unique(as.character(sample_info_file[,2])))
               if (length(group) == 2) {
                 cat("The group:", group, ", will be used.\n")
                 cat("The group:", group, ", will be used.\n", file = file.path(path, "run.log.txt"), append = TRUE)
               } else {
                 # cat(paste("Error: There are", length(group), "groups (", group,
                 #           ") in your sample_info_file. Please provide the group names you want to process.\n"),
                 #     file = file.path(path, "run.log.txt"), append = TRUE)
                 # stop(paste("There are", length(group), "groups (", group,
                 #            ") in your sample_info_file. Please provide the group names you want to process.\n"))
                 temp_group <- group[1:2]
                 cat(paste("Warning: There are", length(group), "groups (", group,
                           ") in your sample_info_file. The first two groups ", paste0(temp_group[1], ' and ', temp_group[2]), " will be used here"),
                     file = file.path(path, "run.log.txt"), append = TRUE)
                 cat(paste("There are", length(group), "groups (", group,
                           ") in your sample_info_file. The first two groups ", paste0(temp_group[1], ' and ', temp_group[2]), " will be used here"))
               }
             } else {
               temp.idx <- which(!(group %in% as.character(sample_info_file[,2])))
               if(length(temp.idx) > 0){
                 cat(paste('Error: ',group[temp.idx], "is/are not in your sample_info_file. (", unique(sample_info_file[,2]), ")\n"),
                     file = file.path(path, "run.log.txt"), append = TRUE)
                 stop(paste(group[temp.idx], "is/are not in your sample_info_file. (", unique(sample_info_file[,2]), ")\n"))
               }
             }

             return(group)

             # cat('The quality check of sample information is done\n\n', file = file.path(path, "run.log.txt"), append = TRUE)

           })





#   showError ------------------------------------------------------------------
showError <- function(expr,
                      error_info = "error"){
  process.result <- try(expr,
                        silent = TRUE)
  if(class(process.result) == "try-error"){
    cat(error_info, "\n")
    process.result <- process.result
  }else{
    process.result <- "Right"
  }
}



# if(check.data){
#   switch(ms2.type,
#          "mgf" = {ms2.file <- grep("mgf", dir(path), value = TRUE)},
#          "msp" = {ms2.file <- grep("msp", dir(path), value = TRUE)},
#          "mzXML" = {ms2.file <- grep("mzXML", dir(path), value = TRUE)})
#
#
#   check.result <- try({checkData(data = ms1_file,
#                                  sample.info = old.sample.info.name,
#                                  ms2.type = ms2.type,
#                                  path = path)}, silent = TRUE)
#
#   if(class(check.result) == "try-error"){
#     cat(check.result[[1]], file = file.path(log.path, "run.log.txt"), append = TRUE)
#     stop(check.result[[1]])
#   }
#
#
#   if(any(as.numeric(check.result[,4]) > 0)){
#     cat('Error: Please check your data to make sure that they are valid.\n', file = file.path(log.path, "run.log.txt"), append = TRUE)
#     stop('Please check your data to make sure that they are valid.\n')
#   }
# }else{
#   cat('Skip data check.\n', file = file.path(log.path, "run.log.txt"), append = TRUE)
#   cat("Skip data check.\n")
# }


#




# checkPara --------------------------------------------------------------------

#' @title checkPara
#' @author Zhiwei Zhou
#' @description check parameters

setGeneric(name = 'checkPara',
           def = function(
             path,
             lib = c('zhumetlib_qtof', 'zhumetlib_orbitrap', 'fiehnHilicLib'),
             polarity = c("positive", "negative"),
             instrument = c("SciexTripleTOF", "AgilentQTOF", "BrukerQTOF", "ThermoOrbitrap",
                            'ThermoExploris', "AgilentDTIMMS", "BrukerTIMS", 'WatersQTOF','WatersTWIMMS'),
             column = c("hilic", "rp"),
             ce = c("30", "10", "20", "35,15", "40", "50",
                    'NCE10', 'NCE20', 'NCE30', 'NCE40', 'NCE50',
                    'SCE20_30_40%', "SNCE20_30_40%"),
             method_lc = c('Amide12min', 'Amide23min', 'Other', 'MetlinRP', 'RP12min'),
             is_rt_calibration = TRUE
           ){
             # # parameter check and decision
             # instrument <- match.arg(instrument)
             # polarity <- match.arg(polarity)
             # column <- match.arg(column)
             # ce <- match.arg(ce)
             # lib <- match.arg(lib)
             # method_lc <- match.arg(method_lc)
             # direction <- match.arg(direction)

             # check path of working directory
             if (stringr::str_detect(path, pattern = ' ')) {
               cat(paste("Error: please remove space in the working path. \n"),
                   file = file.path(path, "run.log.txt"), append = TRUE)
               stop("Error")
             }


             # check library
             if (instrument %in% c("SciexTripleTOF", "AgilentQTOF", "BrukerQTOF", "AgilentDTIMMS", "BrukerTIMS", "WatersTWIMMS", "WatersQTOF")) {
              if (!(lib == 'zhumetlib_qtof')) {
                stop('Please check the parameter (lib), zhumetlib_qtof is required for QTOF instrument\n')
              }
             } else {
               if (!(lib %in% c('zhumetlib_orbitrap', 'fiehnHilicLib'))) {
                 stop('Please check the parameter (lib), zhumetlib_orbitrap/fiehnHilicLib is required for Obitrap instrument\n')
               }
             }

             # check ce
             switch (lib,
                     'zhumetlib_qtof' = {
                       if (!(ce %in% c("30", "10", "20", "35,15", "40", "50"))) {
                         stop('Please check collision energy, the ce', ce, 'is not included in the library ', lib, '\n')
                       }
                     },
                     'fiehnHilicLib' = {
                       if (!(ce %in% c("SNCE20_30_40%"))) {
                         stop('Please check collision energy, the ce', ce, 'is not included in the library ', lib, '\n')
                       }
                     },
                     'zhumetlib_orbitrap' = {
                       if (!(ce %in% c("30", "10", "20",  "40", "50",
                                       'NCE10', 'NCE20', 'NCE30', 'NCE40', 'NCE50',
                                       'SCE20_30_40%', "SNCE20_30_40%"))) {
                         stop('Please check collision energy, the ce', ce, 'is not included in the library ', lib, '\n')
                       }
                     }
             )

             # check column
             switch (column,
                     'hilic' = {
                       if (method_lc %in% c('Other', 'MetlinRP') & is_rt_calibration) {
                         stop('Please check the parameter is_rt_calibration, the RT calibration only support Amide12min/Amide23min/RP12min')
                       }
                     },

                     'rp' = {
                       if (method_lc %in% c('Amide12min', 'Amide23min')) {
                         stop('Please check the parameter method_lc, reverse column should choose MetlinRP, RP12min or Other')
                       }

                       if (is_rt_calibration) {
                         if (method_lc == 'Other' & is_rt_calibration) {
                           stop('Please check the parameter is_rt_calibration, the RT calibration only support Amide12min and Amide23min')
                         }
                       }
                     }
             )

           })


################################################################################
# readMs2 ----------------------------------------------------------------------

#' @title readMs2
#' @author Zhiwei Zhou
#' @param ms2_file MS/MS file names
#' @param ms2_type 'mgf', 'cef', 'msp', 'MzXML'
#' @export

# path <- 'H:/00_projects/03_MetDNA2/00_data/20200922_metdna2_development/pos_201007_metdna2/'
# temp_files <- dir(path) %>% stringr::str_detect('\\.mgf') %>% dir(path)[.]
# ms2_data <- readMs2(ms2_file = file.path(path, temp_files),
#                     ms2_type = 'mgf')
#
# path <- '/home/zhouzw/Data_processing/20210315_tims_data_processing_metdna/'
# temp_files <- dir(path) %>% stringr::str_detect('\\.mgf') %>% dir(path)[.]
# test <- readMs2(ms2_file = file.path(path, temp_files),
#                 ms2_type = 'mgf',
#                 instrument = 'IMMS')

setGeneric(name = 'readMs2',
           def = function(
             ms2_file,
             ms2_type = c('mgf', 'cef', 'msp', 'mzXML'),
             ...
           ){
             ms2_type <- match.arg(ms2_type)

             cat('Read MS2...\n\n')

             switch(ms2_type,
                    'mgf' = {
                      pbapply::pboptions(type='timer', char='+')
                      ms2_data <- pbapply::pblapply(ms2_file, function(x){
                        readMGF(x, ...)
                      })
                    },

                    'cef' = {
                      pbapply::pboptions(type='timer', char='+')
                      ms2_data <- pbapply::pblapply(ms2_file, function(x){
                        readCEF(x, ...)
                      })
                    },

                    'msp' = {
                      pbapply::pboptions(type='timer', char='+')
                      ms2_data <- pbapply::pblapply(ms2_file, function(x){
                        readMSP(x, ...)
                      })
                    },

                    'mzXML' = {
                      pbapply::pboptions(type='timer', char='+')
                      ms2_data <- pbapply::pblapply(ms2_file, function(x){
                        readMzXML(x, ...)
                      })

                    }
             )

             ms2_data <- do.call(c, ms2_data)
           }
)




#   readMGF ----------------------------------------------------------------------

#' @title readMGF
#' @author Zhiwei Zhou
#' @param file the name of mgf file
#' @export

# path <- 'H:/00_projects/03_MetDNA2/00_data/20200922_metdna2_development/pos_201007_metdna2'
# temp_files <- dir(path) %>% stringr::str_detect('\\.mgf') %>% dir(path)[.]
# test <- readMGF(file.path(path, temp_files[1]))

# path <- '/home/zhouzw/Data_processing/20210315_tims_data_processing_metdna/'
# temp_files <- dir(path) %>% stringr::str_detect('\\.mgf') %>% dir(path)[.]
# test <- readMGF(file.path(path, temp_files[1]))

setGeneric('readMGF',
           def = function(file,
                          ...){

             mgf.data <- ListMGF(file)
             whether_has_ccs <- stringr::str_detect(mgf.data[[1]], pattern = '^CCS') %>% any()

             if (!whether_has_ccs) {
               # mgf.data <- ListMGF(file)
               # nl.spec <- grep('^\\d', mgf.data)
               nl.spec <- lapply(mgf.data, function(x) grep('^\\d', x)) # select spec index
               info.mz <- lapply(mgf.data, function(x) grep('^PEPMASS', x, value = T))
               info.rt <- lapply(mgf.data, function(x) grep('^RTINSECONDS', x, value = T))

               info.mz <- sapply(info.mz, function(x){
                 temp <- gsub(pattern = "\\w+=", "", x)
                 if (stringr::str_detect(temp, pattern = "\ ")) {
                   temp <- as.numeric(substr(temp, start = 1, stop = regexpr("\ ", temp)))
                 }

                 temp <- as.numeric(temp)
                 # temp <- round(as.numeric(temp), digits = 4)
               })

               info.rt <- unlist(info.rt)
               # info.rt <- round(as.numeric(gsub(pattern = "\\w+=", "", info.rt)), digits = 0)
               info.rt <- as.numeric(gsub(pattern = "\\w+=", "", info.rt))

               spec <- mapply(function(x, y){
                 list(do.call(rbind, strsplit(x[y], split = " ")))
               },
               x = mgf.data,
               y = nl.spec)

               spec <- lapply(seq(length(info.mz)), function(i){
                 # data.frame(mz=as.numeric(spec[1:10, i]), intensity=as.numeric(spec[11:20, i]), stringsAsFactors = F)
                 data.frame(mz=as.numeric(spec[[i]][, 1]), intensity=as.numeric(spec[[i]][, 2]), stringsAsFactors = F)
               })

               ms2 <- lapply(seq(length(info.mz)), function(i){
                 temp.info <- c(info.mz[i], info.rt[i])
                 names(temp.info) <- c("mz", "rt")
                 temp.spec <- spec[[i]]
                 temp.result <- list(info=temp.info, spec=temp.spec)
               })

               return(ms2)
             } else {
               # mgf.data <- ListMGF(file)
               # nl.spec <- grep('^\\d', mgf.data)
               nl.spec <- lapply(mgf.data, function(x) grep('^\\d', x)) # select spec index
               info.mz <- lapply(mgf.data, function(x) grep('^PEPMASS', x, value = T))
               info.rt <- lapply(mgf.data, function(x) grep('^RTINSECONDS', x, value = T))
               info.ccs <- lapply(mgf.data, function(x) grep('^CCS', x, value = T))

               info.mz <- sapply(info.mz, function(x){
                 temp <- gsub(pattern = "\\w+=", "", x)
                 if (stringr::str_detect(temp, pattern = "\ ")) {
                   temp <- as.numeric(substr(temp, start = 1, stop = regexpr("\ ", temp)))
                 }

                 temp <- as.numeric(temp)
                 # temp <- round(as.numeric(temp), digits = 4)
               })

               info.rt <- unlist(info.rt)
               # info.rt <- round(as.numeric(gsub(pattern = "\\w+=", "", info.rt)), digits = 0)
               info.rt <- as.numeric(gsub(pattern = "\\w+=", "", info.rt))

               info.ccs <- unlist(info.ccs)
               info.ccs <- as.numeric(gsub(pattern = "\\w+=", "", info.ccs))


               spec <- mapply(function(x, y){
                 list(do.call(rbind, strsplit(x[y], split = " ")))
               },
               x = mgf.data,
               y = nl.spec)

               spec <- lapply(seq(length(info.mz)), function(i){
                 # data.frame(mz=as.numeric(spec[1:10, i]), intensity=as.numeric(spec[11:20, i]), stringsAsFactors = F)
                 data.frame(mz=as.numeric(spec[[i]][, 1]), intensity=as.numeric(spec[[i]][, 2]), stringsAsFactors = F)
               })

               ms2 <- lapply(seq(length(info.mz)), function(i){
                 temp.info <- c(info.mz[i], info.rt[i], info.ccs[i])
                 names(temp.info) <- c("mz", "rt", 'ccs')
                 temp.spec <- spec[[i]]
                 temp.result <- list(info=temp.info, spec=temp.spec)
               })

               return(ms2)
             }

           })


#' @title ListMGF
#' @author Zhiwei Zhou
#' @param file the name of mgf file

setGeneric('ListMGF',
           def = function(file){
             mgf.data <- readLines(file)
             nl.rec.new <- 1
             idx.rec <- 1
             rec.list <- list()
             for(nl in 1:length(mgf.data))
             {
               if(mgf.data[nl]=="END IONS")
               {
                 rec.list[idx.rec] <- list(Compound = mgf.data[nl.rec.new : nl])
                 nl.rec.new <- nl + 1
                 idx.rec <- idx.rec + 1
               }
             }
             rec.list
           })



#   readMSP ----------------------------------------------------------------------
#' #' @title readMSP
#' #' @description  Read a MSP file and return a list of spectra for all feature with feature information
#' #' @param file path of the msp file
#' #' @export
#'
#' setGeneric('readMSP', function(file) {
#'   msp.data.list <- ListDB(file)
#'   nr.num.pk <- grep('Num Peaks', msp.data.list[[1]])
#'   info.spec <- lapply(msp.data.list, function(msp.data) {
#'     info.list <- strsplit(msp.data[1:nr.num.pk], split = ': ')
#'     info <- do.call(cbind, info.list)
#'     colnames(info) <- info[1, ]
#'     mz <- round(as.numeric(info[2,3]), digits = 4)
#'     rt <- round(as.numeric(strsplit(info[2,"Comment"], split = "_")[[1]][1])*60, digits = 0)
#'     name <- info[2,"Comment"]
#'
#'     # info <- data.frame(mz=mz, rt=rt)
#'     rownames(info) <- name
#'
#'     spec.list <- strsplit(msp.data[(nr.num.pk+1):length(msp.data)], ' |\\t')
#'     spec.list <- lapply(spec.list, as.numeric)
#'     spec <- do.call(rbind, spec.list)[, c(1, 2), drop = FALSE]
#'     colnames(spec) <- c('mz', 'intensity')
#'
#'     list('info' = info,
#'          'spec' = spec)
#'
#'   })
#'
#'   info.spec
#' })
#'
#' setGeneric('ListDB', function(file) {
#'   msp.data <- readLines(file)
#'   nl.db.new <- 1
#'   idx.db <- 1
#'   db.list <- list()
#'   len.data <- length(msp.data)
#'   for(nl in 1:len.data)
#'   {
#'     if(msp.data[nl]=="") {
#'       db.list[idx.db] <- list(Compound=msp.data[nl.db.new:(nl-1)])
#'       nl.db.new <- nl + 1
#'       idx.db <- idx.db + 1
#'     } else if (nl == len.data){
#'       db.list[idx.db] <- list(Compound=msp.data[nl.db.new:nl])
#'     }
#'   }
#'   db.list
#' })



#' @title readMSP
#' @description read MSP spectra files
#' @author Zhiwei Zhou
#' \email{zhouzw@@sioc.ac.cn}
#' @param file the file name
#' @param mode standard: extract name, mz and RT from the file, which fit for MSP data exported from QI software; all: extract all information in the file
#' @return
#' A list object. info: MS1 information, including ms1, rt etc.; spec: the MS/MS spectrum
#' @examples
#' test <- readMSP(file = 'F:/01 MetIMMS/00 data processing/190515 external validation msms data extraction/zhumetlib_validation_pos_20v_190520.msp', mode = 'all')

setGeneric('readMSP', function(file,
                               mode = c('all', 'standard'),
                               source = c('MetAnalyzer', 'MSDIAL', 'Other')) {
  # devtools::use_package('dplyr')

  mode <- match.arg(mode)
  source <- match.arg(source)
  msp.data.list <- ListDB(file)
  nr.num.pk <- grep('Num Peaks', stringr::str_to_title(msp.data.list[[1]]))
  info.spec <- lapply(msp.data.list, function(msp.data) {
    info.list <- strsplit(msp.data[1:nr.num.pk], split = ': ')
    info <- do.call(cbind, info.list)
    colnames(info) <- info[1, ]

    if (mode=='standard') {
      mz <- round(as.numeric(info[2,3]), digits = 4)
      rt <- round(as.numeric(strsplit(info[2, "Comment"], split = "_")[[1]][1])*60, digits = 0)
      name <- info[2, "Comment"]
      # info <- matrix(c(mz, rt), ncol = 2)
      info <- data.frame(mz=mz, rt=rt)
      rownames(info) <- name
      # colnames(info) <- c("mz", "rt")
    } else {
      info <- as.data.frame(tibble::as.tibble(info))
      info <- info[-1,,drop=F]

      if ('comment' %in% tolower(colnames(info))) {
        info$Comment <- stringr::str_replace(info$Comment, pattern = '\t', replacement = '')
      }

      if (info$NAME == 'Unknown') {
        source <- 'MSDIAL'
      }

      # change colnames for MSP file from MetAnalyzer
      switch (source,
              'MetAnalyzer' = {
                info <- info %>%
                  dplyr::rename(name = NAME,
                                mz = PRECURSORMZ,
                                npeaks = `Num Peaks`) %>%
                  dplyr::mutate(mz = as.numeric(mz))
              },
              'MSDIAL' = {
                info <- info %>%
                  dplyr::rename(name = NAME,
                                mz = PRECURSORMZ,
                                npeaks = `Num Peaks`) %>%
                  dplyr::mutate(mz = as.numeric(mz),
                                name = Comment)
              }
      )

      rownames(info) <- NULL
    }

    # if NULL spectra exported, return NULL (changed for MSDIAL)
    if (length(msp.data) <= nr.num.pk) {
      return(NULL)
    }

    spec.list <- strsplit(msp.data[(nr.num.pk+1):length(msp.data)], ' |\\t')
    spec.list <- lapply(spec.list, as.numeric)
    spec <- do.call(rbind, spec.list)[, c(1, 2), drop = FALSE]
    colnames(spec) <- c('mz', 'intensity')
    # spec <- list('spec' = spec)

    # list('info' = info[-1, , drop = FALSE],
    #      'spec' = spec)

    list('info' = info,
         'spec' = spec)

  })

  # remove NULL spectra
  idx_rm <- sapply(info.spec, function(x){length(x) == 0}) %>% which()
  if (length(idx_rm) > 0) {
    info.spec <- info.spec[-idx_rm]
  }

  return(info.spec)
})



setGeneric('ListDB', function(file) {
  msp.data <- readLines(file)
  nl.db.new <- 1
  idx.db <- 1
  db.list <- list()
  len.data <- length(msp.data)
  for(nl in 1:len.data)
  {
    if(msp.data[nl]=="") {
      db.list[idx.db] <- list(Compound=msp.data[nl.db.new:(nl-1)])
      nl.db.new <- nl + 1
      idx.db <- idx.db + 1
    } else if (nl == len.data){
      db.list[idx.db] <- list(Compound=msp.data[nl.db.new:nl])
    }
  }
  db.list
})


#   readCEF ----------------------------------------------------------------------
#' @title readCEF
#' @description  Read a CEF file and return a list of spectra for all feature
#'  with feature information
#' @param file path of the CEF file
#' @export

# file <- 'H:/00_projects/03_MetDNA2/00_data/20201012_metdna2_developmemt/CEF/mice_liver_pos/QC_01_MA-d3-c3_SR.cef'
# test <- readCEF(file)

setGeneric('readCEF', function(file) {
  cef.data.list <- ListCEF(file)
  # nr.num.pk: the number of
  # nr.num.pk <- grep('Num Peaks', cef.data.list[[1]])
  info.spec <- lapply(cef.data.list, function(cef.data) {

    # Extract basic information
    #  line 13, extract mz
    #  line 2, extract string rt, ccs "rt=\"\\d+.\\d+", example "rt=\"5.934"
    #    rt in seconds

    mz <- stringr::str_extract(string = cef.data[13], "\\d+.\\d+")
    mz <- round(as.numeric(mz), 4)

    rt <- stringr::str_extract(string = cef.data[2], "rt=\"\\d+.\\d+")
    ccs <- stringr::str_extract(string = cef.data[2], "ccs=\"\\d+.\\d+")

    rt <- round(as.numeric(stringr::str_extract(string = rt, pattern = "\\d+.\\d+"))*60, 0)
    ccs <- round(as.numeric(stringr::str_extract(string = ccs, pattern = "\\d+.\\d+")), 1)

    # info <- data.frame(mz=mz, rt=rt, ccs=ccs)
    info <- c(mz, rt, ccs)
    names(info) <- c('mz', 'rt', 'ccs')

    # extract MS/MS spectrum -------------------
    idx.spec.start <- stringr::str_locate(string = cef.data, pattern = "\\s+<MSPeaks>")
    idx.spec.end <- stringr::str_locate(string = cef.data, pattern = "\\s+</MSPeaks>")

    # search the label "<MSPeaks>"
    # 1st is the product ion, 2nd is the isopote

    idx.spec.start <- which(idx.spec.start[,1]>0)[1]+1
    idx.spec.end <- which(idx.spec.end[,1]>0)[1]-1

    spec.list.mz<- stringr::str_extract(cef.data[idx.spec.start:idx.spec.end],
                                        pattern = "x=\"\\d+.\\d+")
    spec.list.int <- stringr::str_extract(cef.data[idx.spec.start:idx.spec.end],
                                          pattern = "y=\"\\d+.\\d+")

    spec.list.mz <- round(as.numeric(stringr::str_extract(spec.list.mz, "\\d+.\\d+")),4)
    spec.list.int <- round(as.numeric(stringr::str_extract(spec.list.int, "\\d+.\\d+")),0)


    spec <- data.frame(mz=spec.list.mz, intensity=spec.list.int)
    spec <- as.matrix(spec)

    list('info' = info,
         'spec' = spec)

  })

  return(info.spec)
})



#' @title ListCEF
#' @author Zhiwei Zhou
#' @param file the name of CEF file

setGeneric('ListCEF',
           def = function(file){
             # cef.data <- readLines(file)
             cef.data <- readr::read_lines(file)

             idx.start <- stringr::str_locate(string = cef.data, "<Compound mppid=\"")
             idx.end <- stringr::str_locate(string = cef.data, "</Compound>")

             idx.start <- which(as.numeric(idx.start[,1])>0)
             idx.end <- which(as.numeric(idx.end[,1])>0)

             rec.list <- mapply(function(x, y){
               result <- cef.data[x:y]
               return(result)
             },
             x = idx.start,
             y = idx.end)

             return(rec.list)
           })


#   readMzXML ----------------------------------------------------------------------

#' @title readMzXML
#' @author Zhiwei Zhou, Xiaotao Shen
#' @param file
#' @export
# file <- 'H:/00_projects/03_MetDNA2/00_data/20201012_metdna2_developmemt/mzXML/POS/Sample2-W30 POS06-POS-W30 POS06.mzXML'
# test <- readMzXML(file)

setGeneric('readMzXML',
           def = function(file){
             mzxml.data <- mzR::openMSfile(file)
             mzxml.info <- mzR::header(mzxml.data)
             mzxml.peak <- mzR::peaks(mzxml.data)

             ms2.idx <- which(mzxml.info$msLevel == 2)
             ms2.info <- mzxml.info[ms2.idx, c("precursorMZ", "retentionTime")]
             ms2.info <- apply(ms2.info, 1, list)
             ms2.spec <- mzxml.peak[ms2.idx]

             result.ms2 <- mapply(function(x, y){
               # temp_info <- data.frame(mz = x[[1]][1],
               #                         rt = x[[1]][2],
               #                         stringsAsFactors = FALSE)

               temp_info <- as.numeric(c(x[[1]][1], x[[1]][2]))
               names(temp_info) <- c("mz", "rt")

               temp_spec <- y
               colnames(temp_spec) <- c("mz", "intensity")

               result <- list('info' = temp_info,
                              'spec' = temp_spec)

               return(result)
             },
             x = ms2.info,
             y = ms2.spec,
             SIMPLIFY = FALSE)

             return(result.ms2)

           })


