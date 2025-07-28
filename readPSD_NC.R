
readPSD_NC <- function(filepath){

  # List files if multiple
  files <- list.files(filepath, full.names = T)

  if (length(files) == 0){
    files.nc = filepath
  }

  if (length(files) > 1){
    files.nc <- files
  }

  data.nc <- lapply(files.nc, function(x){

    tmp.nc <- nc_open(x)

    # Time variable
    {
      time <- ncvar_get(tmp.nc, "time")
      time.unit <- ncatt_get(tmp.nc, "time", "units")

      # Use CFtime package to retrieve timestamp
      time <- CFtime::CFtime(time.unit$value, calendar = "standard", time) # convert time to CFtime class

      # Convert to POSIX format
      time <- lubridate::as_datetime(time$as_timestamp())
    }

    # Concentration values
    {
      conc <- t(ncvar_get(tmp.nc, 'dN_dlogDp'))
      conc.qc <- t(ncvar_get(tmp.nc, 'qc_dN_dlogDp'))
    }

    # Mobility diameters
    {
      bins <- ncvar_get(tmp.nc, 'diameter_mobility_bounds')
      bins = round(rowMeans(t(bins)), 2)
    }

    # Create dataframe
    {
      tmp.df <- as.data.frame(conc)

      # Setnames as diameter midpoints
      colnames(tmp.df) <- round(bins, 2)

      # Apply QC
      tmp.df <- tmp.df * as.data.frame(!conc.qc)

      # Remove empty columns
      tmp.df <- purrr::discard(tmp.df, ~all(is.na(.)))

      # Add time column to dataframe
      tmp.df <- cbind("Time" = time, tmp.df)
    }
  })

  return(data.nc)
}




