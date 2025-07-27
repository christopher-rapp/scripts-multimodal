readPSD_TSI <- function(import.path, tz){

  # -------------------------------------------------------------------------- #
  ##### SECTION: File Wrangling #####
  #'

  # List all level0 csv files
  files.tsi <- list.files(path = import.path,
                          recursive = T,
                          full.names = T)

  data.ls <- lapply(files.tsi, function(x) {

    if (any(str_detect(readLines(paste0(x)), "Laser, Flow,")) == T){

      print(paste0("WARNING: ", "Unallowed delimiter detected in status flag for ", x, ". Removing commas and retrying..."))

      tmp.catch <- paste0(str_replace(readLines(paste0(x)), "Laser, Flow,", "Laser Flow"), collapse = "\n")

      tmp.df = fread(tmp.catch, header = TRUE, stringsAsFactors = FALSE,
                     fill = TRUE, skip = 'Sample #', check.names = FALSE)

      print("Done")
    } else {

      tmp.df <- fread(paste0(x), header = TRUE, stringsAsFactors = FALSE,
                      fill = TRUE, skip = 'Sample #', check.names = FALSE)
    }

    issues = NULL
    issues = str_which(tmp.df$`td(s)`, "flow")
    issues = append(issues, str_which(tmp.df$`td(s)`, "Flow"))

    if (length(issues) > 0){

      # This is used to correct the embedded delimiter issue with
      # Status Flag "CPC Laser, Flow, or Temp out of range causes this problem
      for (i in issues) {

        # Index two columns ahead and shift values to the left
        start = which(c(colnames(tmp.df)) == "td(s)") + 2
        end = ncol(tmp.df)

        for (j in start:end){

          tmp.df[i, (j - 2)] <- tmp.df[i, j, with = FALSE]
        }
      }

      rm(i, j, start, end)
    }

    # Extract string of column names
    tmp.nm = colnames(tmp.df)

    # Change any non-valid locale column headers
    if (any(!validUTF8(tmp.nm)) == TRUE){

      # Replace non-locale characters
      tmp.nm = str_replace(tmp.nm, '<', '')
      tmp.nm = str_replace(tmp.nm, 'cm�', 'cc')
      tmp.nm = str_replace(tmp.nm, '≥', '')
      tmp.nm = str_replace(tmp.nm, '\\?', '')
      tmp.nm = str_replace(tmp.nm, '\\. ', '.')
      tmp.nm = str_replace(tmp.nm, ' \\(', '(')

      # Rename columns without erroneous column characters
      setnames(tmp.df, tmp.nm)
    }

    tmp.df <- tmp.df %>%
      mutate("Units" = "dNdlogDp")

    if (tz == "UTC") {

      # Create time objects
      tmp.df$`UTC Time` = as.POSIXct(paste0(tmp.df$Date, ' ', tmp.df$`Start Time`),
                                     format = '%m/%d/%y %H:%M:%S', tz = tz)

    } else {

      # Create time objects
      tmp.df$`Local Time` = as.POSIXct(paste0(tmp.df$Date, ' ', tmp.df$`Start Time`),
                                       format = '%m/%d/%y %H:%M:%S', tz = tz)

      # Create time objects
      tmp.df$`UTC Time` = with_tz(tmp.df$`Local Time`, tzone = "UTC")
    }

    # Retrieve median diameters and associated columns
    bins.ix = str_which(colnames(tmp.df), "\\d+\\.\\d+")
    bins.nm = as.numeric(colnames(tmp.df)[bins.ix])

    # Calculate the size range within the file
    tmp.df <- tmp.df %>%
      mutate(`Size Range` = bins.nm[length(bins.nm)] - bins.nm[1])

    # Identify column classes of dataframe for usage later
    col.classes <- lapply(tmp.df, class)

    # Find columns containing the POSIX class
    time.index = c(grep("POSIXt", col.classes))

    tmp.df <- tmp.df %>%
      select(all_of(time.index), everything())

    return(tmp.df)
  })

  # Remove empty data.frames
  data.ls <- purrr::discard(data.ls, ~nrow(.) == 0)
}
