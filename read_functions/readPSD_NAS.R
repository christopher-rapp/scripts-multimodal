readPSD_NAS <- function(filepath){

  # List files if multiple
  files <- list.files(filepath, full.names = T)

  if (length(files) == 0){
    files.nas = filepath
  }

  if (length(files) >= 1){
    files.nas <- files
  }

  data.nas <- lapply(files.nas, function(x){

    # Read file
    tmp.nas <- read.delim(x, header = F)

    # Identify where the variable headers are located
    header.ix = as.numeric(str_extract(tmp.nas[1, 1], "(\\d+)\\s\\d+", group = T))

    # Extract
    header.c <- tmp.nas[header.ix, ]

    # Retrieve metadata
    meta.data <- as.vector(tmp.nas[1:(as.numeric(header.ix) - 1), 1])

    # Create a dataframe
    data <- tmp.nas[as.numeric(header.ix):nrow(tmp.nas), ]
    data <- read.table(text = data, header = T)

    # Identify standard variable labels
    var1 <- str_which(meta.data, "end_time")
    var2 <- str_which(meta.data, "numflag")

    # Variable names
    variables.c <- tmp.nas[var1:var2, ]

    # Time variables
    # Based on standard NASA-AMES format
    {
      date.range <- meta.data[7]
      time.format <- meta.data[9]

      # Regular expressions to catch dates
      date1 <- str_extract(date.range, "\\d{4}\\s{1}\\d{2}\\s{1}\\d{2}")
      date2 <- str_extract(date.range, "\\s{1}(\\d{4}\\s{1}\\d{2}\\s{1}\\d{2})", group = T)

      # Format to UTC
      date1 <- as.POSIXct(date1, format = "%Y %m %d", tz = "UTC")
      date2 <- as.POSIXct(date2, format = "%Y %m %d", tz = "UTC")
    }

    if (length(variables.c) != ncol(data)){
      variables.c <- append(paste0("start_time of measurement, ", time.format), variables.c)
    }

    # For PSD's retrieve the mean
    {
      keep.c <- str_which(variables.c, "mean")

      # This keeps the two time columns
      keep.c <- append(1:2, keep.c)

      # This keeps the flag variable
      keep.c <- append(keep.c, length(variables.c))

      data <- data[, keep.c]

      # Which variables are labeled with mean
      tmp.c <- variables.c[keep.c]
    }

    # Binned data
    {
      # Use regex
      bin.ix <- str_which(tmp.c, "D=\\d+")
      binned.c <- tmp.c[bin.ix]

      # Extract bins
      binned.c <- as.numeric(str_extract(binned.c, "(\\d*\\.\\d*)"))

      # Setnames as the bins
      tmp.c[bin.ix] <- binned.c
    }

    # Temperature and pressure if present
    {
      if (any(str_which(tmp.c, "hPa|pressure|Pressure|pres."))){

        tmp.c[str_which(tmp.c, "hPa|pressure|Pressure|pres.")] <- "pressure (hPa)"
      }


      if (any(str_which(tmp.c, "temperature|Temperature|Temp."))){

        tmp.c[str_which(tmp.c, "temperature|Temperature|Temp.")] <- "temperature"

        if (mean(data[, str_which(tmp.c, "temperature|Temperature|Temp.")]) > 50){
          tmp.c[str_which(tmp.c, "temperature|Temperature|Temp.")] <- "temperature (K)"
        } else {
          tmp.c[str_which(tmp.c, "temperature|Temperature|Temp.")] <- "temperature (C)"
        }
      }
    }

    # Rename time variables
    tmp.c[1] <- "starttime"
    tmp.c[2] <- "endtime"

    # Rename dataframe
    setnames(data, new = tmp.c)

    # Format time to POSIX
    {
      data$starttime <- as.POSIXct(date1) + as.difftime(data$starttime, units = "days")
      data$endtime <- as.POSIXct(date1) + as.difftime(data$endtime, units = "days")
    }

    return(data)
  })

  return(data.nas)
}