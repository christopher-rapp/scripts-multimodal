#
#
#
#
# # readPSD_NAS
#
#
# rawDATA <- read.delim("/Users/christopherrapp/Downloads/Ebas_250609_1955/US9050R.20220322000000.20230901184052.smps.particle_number_size_distribution.aerosol.73h.1h.US09L_TSI_3080L_SPL.US09L_smps04.lev2.nas",
#                       header = F)
#
# library(stringr)
# library(data.table)
# library(logr)
# library(dplyr)
# library(tidyr)
# library(ggplot2)
# library(patchwork)
# library(lubridate)
#
# line.header = as.numeric(str_extract(rawDATA[1, 1], "(\\d+)\\s\\d+", group = T))
#
# metaDATA <- as.vector(rawDATA[1:(as.numeric(line.header) - 1), 1])
#
# # Create a dataframe
# {
#   data <- rawDATA[as.numeric(line.header):nrow(rawDATA), ]
#   data <- read.table(text = data, header = T)
# }
#
# header.c <- rawDATA[line.header, ]
#
# var1 <- str_which(metaDATA, "end_time")
# var2 <- str_which(metaDATA, "numflag")
#
# variables.c <- rawDATA[var1:var2, ]
#
#
# date.range <- metaDATA[7]
# time.format <- metaDATA[9]
#
# date1 <- str_extract(date.range, "\\d{4}\\s{1}\\d{2}\\s{1}\\d{2}")
# date2 <- str_extract(date.range, "\\s{1}(\\d{4}\\s{1}\\d{2}\\s{1}\\d{2})", group = T)
#
# date1 <- as.POSIXct(date1, format = "%Y %m %d", tz = "UTC")
# date2 <- as.POSIXct(date2, format = "%Y %m %d", tz = "UTC")
#
# if (length(variables.c) != ncol(data)){
#   variables.c <- append(paste0("start_time of measurement, ", time.format), variables.c)
# }
#
# keep.c <- str_which(variables.c, "mean")
#
# # This keeps the two time columns
# keep.c <- append(1:2, keep.c)
#
# # This keeps the flag variable
# keep.c <- append(keep.c, length(variables.c))
#
# data <- data[, keep.c]
#
#
# tmp.c <- variables.c[keep.c]
#
#
# bin.ix <- str_which(tmp.c, "D=\\d+")
# binned.c <- tmp.c[bin.ix]
#
# binned.c <- as.numeric(str_extract(binned.c, "(\\d*\\.\\d*)"))
#
# tmp.c[bin.ix] <- binned.c
#
# if (any(str_which(tmp.c, "hPa|pressure|Pressure|pres."))){
#
#   tmp.c[str_which(tmp.c, "hPa|pressure|Pressure|pres.")] <- "pressure (hPa)"
# }
#
#
# if (any(str_which(tmp.c, "temperature|Temperature|Temp."))){
#
#   tmp.c[str_which(tmp.c, "temperature|Temperature|Temp.")] <- "temperature"
#
#   if (mean(data[, str_which(tmp.c, "temperature|Temperature|Temp.")]) > 50){
#     tmp.c[str_which(tmp.c, "temperature|Temperature|Temp.")] <- "temperature (K)"
#   } else {
#     tmp.c[str_which(tmp.c, "temperature|Temperature|Temp.")] <- "temperature (C)"
#   }
# }
#
# tmp.c[1] <- "starttime"
# tmp.c[2] <- "endtime"
#
# setnames(data, new = tmp.c)
#
#
# as.POSIXct(data$starttime)
#
# paste0(date1, " ", data$starttime)
#
# data$starttime <- as.POSIXct(date1) + as.difftime(data$starttime, units = "days")
# data$endtime <- as.POSIXct(date1) + as.difftime(data$endtime, units = "days")
#
# test <- multimodal.fitting(data, log.path = "/Users/christopherrapp/Documents", frequency = 5, 30, 10, 1, 1000, 0.05, F)
#
# passlist <- purrr::map(test, 1)
#
# pass.ls <- test[which(passlist == T)]
#
# plotlist <- purrr::map(pass.ls, 2)
#
# plotlist[[8]]
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
