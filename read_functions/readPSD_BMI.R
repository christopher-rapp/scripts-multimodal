readPSD_BMI <- function(import.path, tz){

  # -------------------------------------------------------------------------- #
  ##### SECTION: File Wrangling #####
  #'

  {
    # List all level0 csv files
    files.bmi <- list.files(path = import.path,
                            recursive = T,
                            full.names = T)

    # Detect files of specific type
    files.MONO_DATA <- files.bmi[grepl(paste("MONO_DATA", collapse = '|'), ignore.case = T, files.bmi)]
    files.SEMS_CONC <- files.bmi[grepl(paste("SEMS_CONC", collapse = '|'), ignore.case = T, files.bmi)]
    files.SEMS_AUX <- files.bmi[grepl(paste("SEMS_AUX", collapse = '|'), ignore.case = T, files.bmi)]
    files.SEMS_RAW <- files.bmi[grepl(paste("SEMS_RAW", collapse = '|'), ignore.case = T, files.bmi)]
    files.SEMS_VOLTS <- files.bmi[grepl(paste("SEMS_VOLTS", collapse = '|'), ignore.case = T, files.bmi)]
  }

  # ------------------------------------------------------------------------ #
  ##### SECTION: Create Export Structure #####
  #'

  data.ls <- lapply(files.SEMS_CONC, function(x) {

    # Use data.table's fread to read in data
    # Fastest method in reading CSV's and least memory consuming
    # fill MUST equal FALSE with SPIN files
    tmp.df <- data.table::fread(paste0(x),
                                na.strings = c("", "NA", "NaN"),
                                skip = "#StartDate",
                                fill = TRUE,
                                strip.white = TRUE,
                                stringsAsFactors = FALSE)

    if (length(str_which(colnames(tmp.df), "Bin\\d{1,3}_")) != 0){
      setnames(tmp.df, new = str_replace(colnames(tmp.df), "Bin\\d{1,3}_", ""))
    } else {
      tmp.df <- NULL
    }

    filename.c <- str_extract(x, "(?<=CONC_)\\d{6}\\w{1}\\d{6}")

    if (nrow(tmp.df) >= 1){

      tmp.df <- tmp.df %>%
        mutate("Units" = "dNdlogDp")

      if (tz == "UTC") {

        # Create time objects
        tmp.df$`UTC Time` = as.POSIXct(paste0(tmp.df$`#StartDate`, ' ', tmp.df$`StartTime`),
                                       format = '%y%m%d %H:%M:%S', tz = tz)

      } else {

        # Create time objects
        tmp.df$`Local Time` = as.POSIXct(paste0(tmp.df$`#StartDate`, ' ', tmp.df$`StartTime`),
                                         format = '%y%m%d %H:%M:%S', tz = tz)

        # Create time objects
        tmp.df$`UTC Time` = with_tz(tmp.df$`Local Time`, tzone = "UTC")
      }

      # Retrieve median diameters and associated columns
      bins.ix = str_which(colnames(tmp.df), "\\d+\\.\\d+")
      bins.nm = as.numeric(colnames(tmp.df)[bins.ix])

      # Calculate the size range within the file
      tmp.df <- tmp.df %>%
        mutate(`Size Range` = bins.nm[length(bins.nm)] - bins.nm[1])
    }

    # Identify column classes of dataframe for usage later
    col.classes <- lapply(tmp.df, class)

    # Find columns containing the POSIX class
    time.index = c(grep("POSIXt", col.classes))

    tmp.df <- tmp.df %>%
      select(all_of(time.index), everything())
  })

  # Remove empty data.frames
  data.ls <- purrr::discard(data.ls, ~nrow(.) == 0)

  return(data.ls)
}
