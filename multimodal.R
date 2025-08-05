#' Title: Multimodal lognormal particle size distribution fitting algorithm
#'
#' @author Christopher Rapp - first and corresponding author, Fred Brechtel, and Daniel Cziczo
#'
#' @description
#' A function that fits a multimodal lognormal particle size distribution or MM-LPSD to
#' a measured dataset using an iterative subtractive Levenberg–Marquardt NLS
#' algorithm. The curve fitting is based solely on the theoretical lognormal PSD
#' assumed for both ambient and laboratory aerosol.
#'
#' @param data a data.frame with a column of time data formatted to POSIX CT
#' and binned concentrations with median Dp as name. Does not need to be in nm, but limits need to be consistent
#' @param log.path
#' @param frequency
#' @param max.iterations
#' @param max.modes
#' @param lower.limit
#' @param upper.limit
#' @param accuracy
#' @param verbose
#'
#'
#' @returns A 6-level list containing:
#' 1. Success flag either TRUE or FALSE for easy filtering, if true # of modes
#' 2. Interpolated data between lower.limit and upper.limit for each mode;
#' 3. Mode parameters and performance report;
#' 4. Plot list for visual check
#'
#' @import logr
#'
#' @export
#'
#' @examples

# max.iterations = 30
# max.modes = 10
# lower.limit = 10
# upper.limit = 1500
# accuracy = 0.05
# verbose = T
# log.path = '/Users/christopherrapp/Documents'
#
# data = rawSEMS.df


multimodal.fitting <- function(data, log.path, labeling, frequency, max.iterations, max.modes, lower.limit, upper.limit, NMRSE.threshold, FVU.threshold, verbose){

  # Default arguments ----------------------------------------------------------

  if(missing(labeling)) labeling <- T
  if(missing(max.iterations)) max.iterations <- 20
  if(missing(max.modes)) max.modes <- 5
  if(missing(lower.limit)) lower.limit <- 10
  if(missing(upper.limit)) upper.limit <- 1000
  if(missing(NMRSE.threshold)) NMRSE.threshold <- 0.05
  if(missing(FVU.threshold)) FVU.threshold <- 20
  if(missing(verbose)) verbose <- FALSE
  if(missing(log.path)) log.path <- tempdir()

  {
    {
      # Identify column classes of dataframe for usage later
      col.classes <- lapply(data, class)

      # Convert to data table for column indexing
      data <- as.data.table(data)

      # Find columns containing the POSIX class
      time.index = c(grep("POSIXt", col.classes))[1]

      if (length(time.index) >= 1){

        timestamp = unlist(data[1, time.index, with = F])

        filename <- format(as.POSIXct(timestamp, origin = "1970-01-01"), "%Y%m%d%H%M%S")
      } else {
        filename <- "NA"
      }

      log.path <- file.path(paste0(log.path, "/multimodal", filename, "_",  format(Sys.time(), "%Y%m%d%H%M%S"), ".log"))
      print(paste0("Log Path: ", log.path))

      LOG <- logr::log_open(log.path, show_notes = F)
    }

    # Return objects -------------------------------------------------------------

    flag.control <- TRUE
    export.gg <- list()
    export.df <- list()
    export.ft <- list()
    export.lm <- list()
    export.pf <- list()

    # PRE-PROCESSING -----------------------------------------------------------

    # Identify column classes of dataframe for usage later
    col.classes <- lapply(data, class)

    {
      # Temporal data ----------------------------------------------------------

      {
        # Find columns containing the POSIX class
        time.index = c(grep("POSIXt", col.classes))

        if (length(time.index) >= 1){

          # See if the time zone is UTC
          if (any(which(str_detect(colnames(data), pattern = c("UTC|GMT"))))){

            timezone <- "UTC"

            # Select UTC instance
            time.index <- which(str_detect(colnames(data), pattern = c("UTC|GMT")))
          } else{

            timezone <- "Unknown"

            # Select first instance
            time.index <- time.index[1]
          }

          tmp.print <- paste0("Current Dataset Time: ",
                              lubridate::as_datetime(as.numeric(first(data[, time.index, with = F]))),
                              " ", timezone)

          if (verbose){

            # Print to console
            print(tmp.print)

            # Print to log
            log_print(tmp.print, console = FALSE)
          } else {
            # Print to log
            log_print(tmp.print, console = FALSE)
          }

        } else {

          tmp.print <- "The data set does not contain any detectable POSIXct POSIXt time objects, please format your data according to the README documentation"

          # Print to console
          print(tmp.print)

          # Print to log
          log_print(tmp.print, console = FALSE)
        }
      }

      # Sampling frequency detection (mins)
      {
        # Calculate the time difference between samples in minutes
        tmp.diff <- difftime(data[[time.index]][-1], data[[time.index]][-nrow(data)], units = "mins")
        time.interval <- round(mean(tmp.diff, na.rm = T), 1)

        tmp.print <- paste0("Dataset sampling frequency is ", time.interval, " min")

        if (verbose){
          # Print to console
          print(tmp.print)

          # Print to log
          log_print(tmp.print, console = FALSE)
        } else {
          # Print to log
          log_print(tmp.print, console = FALSE)
        }
      }

      # INTERVAL AVERAGING -----------------------------------------------------

      {
        if (any(c(is.na(frequency), is.null(frequency), is.nan(frequency)))){
          Group = lubridate::as_datetime(as.numeric(first(data[, time.index, with = F])))
        } else {
          Group = round_date(data[[time.index]], paste(frequency, " mins"))
        }

        # Averaging frequency for use in multimodal analysis
        tmp.df <- data %>%
          mutate(`Group` = Group, .after = everything())

        # Split by group
        data.ls <- split(tmp.df, tmp.df$Group)

        # Remove empty lists to prevent errors downstream
        data.ls <- purrr::discard(data.ls, ~nrow(.) == 0)
      }
    }

    # PROCESSING ---------------------------------------------------------------

    export.list <- lapply(data.ls, function(z){

      # Used for logging
      date.time.c <- lubridate::as_datetime(as.numeric(first(z[, time.index, with = F])))
      date.time.end.c <- lubridate::as_datetime(as.numeric(last(z[, time.index, with = F])))

      # Number of observations per group
      n.obs <- nrow(z)

      # Size bins --------------------------------------------------------------
      {
        # Retrieve binned diameters and associated columns
        bins.ix = str_which(colnames(z), "\\d+(\\.\\d+)?$")
        bins.nm = as.numeric(colnames(z)[bins.ix])

        # Subset all non binned data from dataframe
        tmp.df <- z %>%
          select(!all_of(bins.ix))

        # Calculate log difference of diameters
        {
          # Stagger the bins i.e trim end points for each to subtract
          # Trims first bin as it typically has the highest density of measurements
          tmp1 <- as.numeric(bins.nm[-c(length(bins.nm))])/1000
          tmp2 <- as.numeric(bins.nm[-1])/1000

          # Find log difference between bin diameters
          dlogDp = log10(tmp2) - log10(tmp1)
        }

        # Create a vector to select columns
        # Due to needing to calculate dlogDp the last column is dropped
        select.c <- bins.ix[1:(length(bins.ix) - 1)]
      }

      # Number Density Matrix --------------------------------------------------
      {
        # Subset data and create a numeric matrix
        dNdlogDp.mat <- z[, select.c, with = F]
        dNdlogDp.mat <- sapply(dNdlogDp.mat, as.numeric)

        # Calculate number density by multipling by the log difference to convert
        # dNdlogDp to dN
        dN.mat = t(t(dNdlogDp.mat)*dlogDp)

        # Remove empty rows
        dN.mat <- dN.mat[complete.cases(dN.mat), ]

        # This is to prevent single modes of organic particles to throw an error
        if (is.null(nrow(dNdlogDp.mat))){
          dNdlogDp <- dNdlogDp.mat
        } else {
          dNdlogDp <- colMeans(dNdlogDp.mat, na.rm = T)
        }

        # This is to prevent single modes of organic particles to throw an error
        if (is.null(nrow(dN.mat))){
          dN <- dN.mat
        } else {
          dN <- colMeans(dN.mat, na.rm = T)
        }

        # Test if sum of all particle density is lower than 100 n/cc
        if (sum(dN.mat) < 100){

          flag.control <- FALSE

          tmp.print <- paste0(date.time.c, ": Total concentration is lower than 100 n/cc, model fitting is not possible for this dataset")

          if (verbose){
            # Print to console
            print(tmp.print)

            # Print to log
            log_print(tmp.print, console = FALSE)
          } else {
            # Print to log
            log_print(tmp.print, console = FALSE)
          }

          {
            stats.c <- c("Pearson Correlation", "RMSE", "NRMSE", "dN RMSE", "dN NRMSE", "Students T Test", "Chi-Squared")
            stats.nm <- rep(NA, 7)

            export.pf <- t(as.data.frame(stats.nm))
            colnames(export.pf) <- stats.c
            rownames(export.pf) <- NULL
          }

          export.df <- data.frame(
            "Dp" = as.numeric(names(dN)),
            "Predicted dNdlogDp" = NA,
            "Predicted dN" = NA,
            "Actual dNdlogDp" = dNdlogDp,
            "Actual dN" = dN,
            "Residual dNdlogDp" = NA,
            "Residual dN" = NA,
            "Ratio" = NA,
            check.names = F
          )

          export.ls <- list("pass" = flag.control, "plot" = NULL, "data" = export.df, "predict" = NULL, "fits" = NULL, "evaluation" = export.pf)
          return(export.ls)
        }
      }

      # MODEL SETUP ------------------------------------------------------------
      {
        # Create temporary dataset used to initialize the fitting algorithm
        tmp.data = data.frame(Dp = as.numeric(names(dN)),
                              dlogDp = dlogDp,
                              dN = dN,
                              dNdlogDp = dNdlogDp)

        peak.fitting <- list()
        model.fitting <- list()
        export.lm <- list()
        slopes <- list()

        # Initial parameters to initialize while loop
        FVU = 100
        FVU.i = 100
        i = 1
        m = 1
      }

      # Peak Identified Levenberg-Marquardt NLS Algorithm ----------------------
      {
        # Start loop
        while (FVU > 1 & i < max.iterations & length(purrr::compact(model.fitting)) < max.modes){

          {
            # First iteration
            if (i == 1){
              tmp.data$residual = tmp.data$dNdlogDp
            }

            # Find all peaks
            {
              # Peak identification
              peaks.df <- as.data.frame(pracma::findpeaks(tmp.data$residual,
                                                          minpeakdistance = 5,
                                                          sortstr = T
              ))

              # Use indices to match corresponding diameters
              peaks.df$Max <- peaks.df$V1
              peaks.df$Mode <- tmp.data$Dp[peaks.df$V2]
              peaks.df$Lower <- tmp.data$Dp[peaks.df$V3]
              peaks.df$Upper <- tmp.data$Dp[peaks.df$V4]
              peaks.df$Width = peaks.df$V4 - peaks.df$V3

              # Stop code from continuing if there are no ramps detected which means the dataset will be empty and continue to throw errors
              if (m > nrow(peaks.df)) {

                tmp.print <- paste0(date.time.c, ": Terminating loop, no remaining peaks")

                if (verbose){
                  # Print to console
                  print(tmp.print)

                  # Print to log
                  log_print(tmp.print, console = FALSE)
                } else {
                  # Print to log
                  log_print(tmp.print, console = FALSE)
                }

                break
              }

              # Parameters used to select data
              peaks.df <- peaks.df[m, ]
              range.c <- peaks.df$V3:peaks.df$V4
            }

            # Peak statistics
            {
              # Total Number Concentration
              {
                # Find both dNdlogDp and dlogDp to convert to dN
                tmp.dNdlogDp = tmp.data$residual[range.c]
                tmp.dlogDp = tmp.data$dlogDp[range.c]

                # Find the original diameters corresponding to the peak range
                tmp.Dp = tmp.data$Dp[range.c]

                # Add column for total number concentration (N)
                peaks.df$N = sum(tmp.dNdlogDp*tmp.dlogDp)

                if (sum(tmp.dNdlogDp*tmp.dlogDp) > 0){

                  # Add column for total number concentration (N)
                  peaks.df$N = sum(tmp.dNdlogDp*tmp.dlogDp)
                } else {

                  tmp.print <- paste0(date.time.c, ": Current Loop Iteration: ", i, ", Negative Concentration!")

                  if (verbose){
                    # Print to console
                    print(tmp.print)

                    # Print to log
                    log_print(tmp.print, console = FALSE)
                  } else {
                    # Print to log
                    log_print(tmp.print, console = FALSE)
                  }

                  m = m + 1
                  i = i + 1

                  next
                }
              }

              # GSD calculation
              {
                # Prevent negative values from appearing in GSD calculation
                # This also prevents curves that are predicting negative concentrations
                tmp.dNdlogDp[which(tmp.dNdlogDp < 0)] <- 0
                tmp.dlogDp[which(tmp.dlogDp < 0)] <- 0

                # Calculate GSD
                # Warnings are suppressed for function output
                GSD = suppressWarnings({
                  10^(sqrt(sum(tmp.dNdlogDp*((log10(tmp.Dp) - log10(peaks.df$Mode))^2))/(peaks.df$N - 1)))
                })

                # If GSD fails, increase iteration counter and try the second peak
                if (is.nan(GSD)){

                  tmp.print <- paste0(date.time.c, ": Current Loop Iteration: ", i, ", GSD Error!")

                  if (verbose){
                    # Print to console
                    print(tmp.print)

                    # Print to log
                    log_print(tmp.print, console = FALSE)
                  } else {
                    # Print to log
                    log_print(tmp.print, console = FALSE)
                  }

                  m = m + 1
                  i = i + 1

                  next
                } else {

                  peaks.df$GSD = GSD
                }
              }
            }

            # Temporary dataframe for LM-NLS model
            # This modifies the "view" NLS has when trying to curve fit
            tmp <- data.frame(x = tmp.data$Dp[range.c], y = tmp.data$residual[range.c])

            # LM_NLS -----------------------------------------------------------------
            # Try fitting an NLS model to the residuals (or original data if first iteration)
            {
              NLS.break <<- FALSE

              # Try fitting an NLS model to the residuals (or original data if first iteration)
              tryCatch(expr = {

                NLS.MODEL <- suppressWarnings(minpack.lm::nlsLM(y ~ dNdlogDp.PSD(x, N, GSD, Dpg), data = tmp,
                                                                start = list(N = peaks.df$N, GSD = peaks.df$GSD, Dpg = peaks.df$Mode),
                                                                trace = F))

              },
              error = function(e){
                NLS.break <<- TRUE
              }
              )

              # Stop code from continuing if there are no ramps detected which means the dataset will be empty and continue to throw errors
              if (NLS.break) {

                tmp.print <- paste0(date.time.c, ": Terminating loop ", i, ", model fitting failure")

                if (verbose){
                  # Print to console
                  print(tmp.print)

                  # Print to log
                  log_print(tmp.print, console = FALSE)
                } else {
                  # Print to log
                  log_print(tmp.print, console = FALSE)
                }

                break
              }

              # Predict
              tmp.predict = data.frame(x = seq(lower.limit, upper.limit, by = 0.01))
              tmp.predict[[paste0("Mode ", i)]] <- predict(NLS.MODEL, newdata = tmp.predict)
            }

            tmp.xL <- round(tmp.predict$x, 2)
            tmp.xR <- round(tmp.data$Dp, 2)

            index.match <- which(tmp.xL %in% tmp.xR)

            predict.diameter = tmp.predict[index.match, 1]
            predict.dNdlogDp = tmp.predict[index.match, 2]

            if (length(predict.diameter) != nrow(tmp.data)){

              tmp.print <- paste0(date.time.c, ": Error, please modify lower and upper limits to accommadate data set")

              # Print to console
              print(tmp.print)

              # Print to log
              log_print(tmp.print, console = FALSE)

              break
            }

            # Evaluation parameters
            p.N <- summary(NLS.MODEL)$parameters[, 4][1]
            p.GSD <- summary(NLS.MODEL)$parameters[, 4][2]
            p.Dpg <- summary(NLS.MODEL)$parameters[, 4][3]

            peak.height.ratio = max(tmp.predict[, 2])/peaks.df$Max

            if (peak.height.ratio > 1.5){

              tmp.data$residual[which(tmp.data$Dp < peaks.df$Upper & tmp.data$Dp > peaks.df$Lower)] <- 0

              tmp.print <- paste0(date.time.c, ": Current Loop Iteration: ", i, ", Peak estimation is 150% higher than data, resetting to zero and retrying")

              if (verbose){
                # Print to console
                print(tmp.print)

                # Print to log
                log_print(tmp.print, console = FALSE)
              } else {
                # Print to log
                log_print(tmp.print, console = FALSE)
              }
            } else if (isTRUE(any(p.N >= 0.05, p.GSD >= 0.05, p.Dpg >= 0.05)) == T){

              tmp.print <- paste0(date.time.c, ": Current Loop Iteration: ", i, ", LM-NLS significance greater than 0.05")

              if (verbose){
                # Print to console
                print(tmp.print)

                # Print to log
                log_print(tmp.print, console = FALSE)
              } else {
                # Print to log
                log_print(tmp.print, console = FALSE)
              }

              break
            } else {

              tmp.diff <- tmp.data$residual - predict.dNdlogDp

              tmp.data$residual <- tmp.diff

              # Add model values to output list
              model.fitting[[i]] <- tmp.predict

              # Add peak values to output list
              peak.fitting[[i]] <- peaks.df

              # Add model values to output list
              export.lm[[i]] <- NLS.MODEL

              FVU.i = round((var(tmp.data$residual)/var(tmp.data$dNdlogDp))*100, 2)
            }

            if (FVU.i < FVU){

              FVU <- FVU.i
            }

            tmp.print <- paste0(date.time.c, ": Current Loop Iteration: ", i, ", Remaining Variance: ", FVU, "%",
                                ", Number of Modes: ", nrow(rbindlist(purrr::compact(peak.fitting))))

            # Move to next peak
            i <- i + 1

            if (verbose){

              # Print to console
              print(tmp.print)

              # Print to log
              log_print(tmp.print, console = F)
            } else {
              # Print to log
              log_print(tmp.print, console = FALSE)
            }
          }
        }
      }

      # Post-Processing --------------------------------------------------------
      {
        # Output list of identified peak data from loop
        {
          # Bind list output
          tmp.ls <- purrr::compact(export.lm)

          if (length(tmp.ls) == 0){

            export.ls <- list("pass" = flag.control, "plot" = NULL, "data" = export.df, "predict" = NULL, "fits" = NULL, "evaluation" = export.pf)
            return(export.ls)
          }

          tmp.df <- data.frame(matrix(nrow = length(tmp.ls), ncol = 11))
          for (s in 1:length(tmp.ls)){

            tmp.df[s, 1] <- paste0("Mode ", s)
            tmp.df[s, 2] <- c(coef(tmp.ls[[s]]))[1]
            tmp.df[s, 3] <- c(coef(tmp.ls[[s]]))[2]
            tmp.df[s, 4] <- c(coef(tmp.ls[[s]]))[3]
            tmp.df[s, 5] <- BIC(tmp.ls[[s]])
            tmp.df[s, 6] <- deviance(tmp.ls[[s]])
            tmp.df[s, 7] <- sum((tmp.ls[[s]]$m$getEnv()$y - mean(tmp.ls[[s]]$m$getEnv()$y))^2)
            tmp.df[s, 8] <- 1 - tmp.df[s, 6]/tmp.df[s, 7]
            tmp.df[s, 9] <- signif(summary(tmp.ls[[s]])$parameters[, 4][1], 1)
            tmp.df[s, 10] <- signif(summary(tmp.ls[[s]])$parameters[, 4][2], 2)
            tmp.df[s, 11] <- signif(summary(tmp.ls[[s]])$parameters[, 4][3], 3)
          }

          setnames(tmp.df, new = c("Mode Label", "N", "GSD", "Dpg", "BIC", "RSS", "TSS", "R2", "N T pval", "GSD T pval", "Dpg T pval"))

          peaks.df <- rbindlist(purrr::compact(peak.fitting))

          # Drop first 4 columns of data
          peaks.df <- peaks.df %>%
            select(Max:Width)

          #
          export.lm <- cbind(tmp.df, peaks.df) %>%
            relocate(Max:Width, .after = Dpg)
        }

        # Output list of model fitting from loop
        {
          model.fitting <- purrr::compact(model.fitting)
          output <- do.call(cbind, model.fitting)
          output <- output[append(1, str_which(colnames(output), 'Mode'))]
          output <- round(output, 2)

          setnames(output, old = colnames(output)[str_which(colnames(output), pattern = "Mode \\d{1,2}")],
                   new = paste0("Mode ", seq_along((str_which(colnames(output), pattern = "Mode \\d{1,2}")))))

          setnames(output, old = "x", new = "Dp")

          if (ncol(output) > 2){
            output$dNdlogDp <- rowSums(output[, 2:ncol(output)])
          } else {
            output$dNdlogDp <- output[, 2]
          }

          export.ft <- output
        }
      }

      # Statistics -------------------------------------------------------------
      {

        # Calculate various statistical measures

        # Matching predictions to measurements
        {
          # Find the matched diameters so measurements to predicted binned diameters is correct
          match.ix <- which(output$Dp %in% tmp.data$Dp)

          # Diameters
          Dp <- output$Dp[match.ix]

          # Predicted Lognormal Concentration
          predicted.dNdlogDp <- output$dNdlogDp[match.ix]

          # Predicted Concentration
          predicted.dN <- predicted.dNdlogDp*dlogDp

          # Actual lognormal concentration
          actual.dNdlogDp <- tmp.data$dNdlogDp

          # Actual concentration
          actual.dN <- tmp.data$dN
        }

        # Pearson Correlation
        {
          stats.R2 = round(cor(predicted.dN, actual.dN), 4)
        }

        # RMSE
        {
          stats.RMSE = round(Metrics::rmse(actual.dNdlogDp, predicted.dNdlogDp), 2)
          dN.RMSE = round(Metrics::rmse(actual.dN, predicted.dN), 2)
        }

        # Max min normalized RMSE
        {
          NRMSE = stats.RMSE/(max(actual.dNdlogDp) - min(actual.dNdlogDp))
          dN.NRMSE = dN.RMSE/(max(actual.dN) - min(actual.dN))
        }

        # Significance testing
        # NOTE: This is not a recommended evaluation method as it has demonstrated to
        # pass or fail even the best fittings
        # USE AT YOUR OWN RISK
        {
          stats.STTEST <- t.test(actual.dNdlogDp, predicted.dNdlogDp, paired = T)
          stats.STTEST <- round(stats.STTEST$p.value, 4)

          stats.CHI <- suppressWarnings(chisq.test(actual.dNdlogDp, predicted.dNdlogDp))
          stats.CHI <- round(stats.CHI$p.value, 4)
        }

        export.df <- data.frame(
          "Dp" = Dp,
          "Predicted dNdlogDp" = predicted.dNdlogDp,
          "Predicted dN" = predicted.dN,
          "Actual dNdlogDp" = actual.dNdlogDp,
          "Actual dN" = actual.dN,
          "Residual dNdlogDp" = actual.dNdlogDp - predicted.dNdlogDp,
          "Residual dN" = actual.dN - predicted.dN,
          "Ratio" = predicted.dNdlogDp/actual.dNdlogDp,
          check.names = F
        )

        stats.c <- c("Pearson Correlation", "RMSE", "NRMSE", "dN RMSE", "dN NRMSE", "Students T Test", "Chi-Squared")
        stats.nm <- c(stats.R2, stats.RMSE, NRMSE, dN.RMSE, dN.NRMSE, stats.STTEST, stats.CHI)

        export.pf <- t(as.data.frame(stats.nm))
        colnames(export.pf) <- stats.c
        rownames(export.pf) <- NULL

        rm(stats.nm, stats.c)

        if (verbose){
          print(paste0("Concentration RMSE: ", stats.RMSE, " n/cc"))
        }

        if (FVU > FVU.threshold){
          flag.control <- FALSE
        }

        if (NRMSE > NMRSE.threshold){
          flag.control <- FALSE
        }
      }

      # Plotting ---------------------------------------------------------------
      {
        plot.df <- output %>%
          pivot_longer(
            cols = !c("Dp"),
            names_to = "Mode",
            values_to = "Concentration"
          )

        setnames(plot.df, old = "Dp", new = "Diameter")

        # Replace label with total for plotting
        plot.df[plot.df == "dNdlogDp"] <- "Total"

        plot.df <- plot.df %>%
          mutate(`Mode` = relevel(factor(`Mode`, ordered = F), ref = "Total"))

        # Use log distributed values rather than numerical to improve performance
        # Needs several orders of magnitude less data to produce smooth curves
        x = 10^seq(log10(lower.limit), log10(upper.limit), length.out = 1000)

        x.limits = 10^((log10(lower.limit)):(log10(upper.limit)))
        y.limits = c(-1*max(c(export.df$`Predicted dNdlogDp`, export.df$`Actual dNdlogDp`))/10,
                     max(c(export.df$`Predicted dNdlogDp`, export.df$`Actual dNdlogDp`))/10)

        y.breaks = unique(c(pretty(export.df$`Predicted dNdlogDp`), pretty(export.df$`Actual dNdlogDp`)))

        x.label.positions = round(log10(x.limits))

        breaks <- NULL
        labels <- NULL
        for (i in x.label.positions){
          breaks <- append(breaks, i)
          labels <- append(labels, as.character(10^(i)))
        }

        cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

        cbPalette <- colorRampPalette(cbPalette)(max.modes)

        cm.palette <- c("black", cbPalette)

        plot.df <- plot.df %>%
          filter(`Diameter` <= x.limits[length(x.limits)] & `Diameter` >= x.limits[1])

        error <- output %>%
          mutate(Min = if_else(dNdlogDp - stats.RMSE < 0, 0, dNdlogDp - stats.RMSE)) %>%
          mutate(Max = dNdlogDp + stats.RMSE)

        labels.df <- export.lm %>%
          select(`Mode Label`, `Dpg`, `N`, `GSD`, `Max`)

        labels.df$`Mode Label` <- str_replace(labels.df$`Mode Label`, "Mode ", "")

        # Plot title
        if (n.obs == 1 & any(c(is.na(frequency), is.null(frequency), is.nan(frequency)))){
          plot.title = paste0(date.time.c)
        } else {
          plot.title = paste0(date.time.c, " - ", date.time.end.c)
        }

        text.labels <- plot.df %>%
          group_by(Mode) %>%
          arrange(desc(Concentration)) %>%
          slice(1) %>%
          ungroup()

        text.labels <- text.labels %>%
          filter(Mode != "Total") %>%
          mutate(label = paste0("(", str_replace(Mode, "Mode ", ""), ")"))

        if (labeling == FALSE){
          text.labels <- text.labels %>%
            mutate(label = "")
        }

        # Plotting each panel seperately then using patchwork to join
        {
          top.gg <- ggplot(data = tmp.data, aes(x = Dp, y = dNdlogDp)) +
            geom_point(shape = 1) +
            geom_line(data = plot.df, aes(x = Diameter, y = Concentration, color = Mode), linewidth = 0.5) +
            geom_ribbon(data = error, aes(x = Dp, ymin = Min, ymax = Max), inherit.aes = F, fill = "gray75", alpha = 0.5) +
            scale_x_log10(breaks = 10^breaks,
                          labels = labels,
                          limits = c(min(x.limits), max(x.limits))) +
            scale_y_continuous(breaks = y.breaks, limits = range(y.breaks)) +
            ggrepel::geom_text_repel(data = text.labels, aes(x = `Diameter`, y = `Concentration`, label = label),
                                     nudge_x = 0.1, nudge_y = y.limits[2]*0.1, box.padding = 1) +
            labs(title = plot.title,
                 subtitle = paste0("Pass: ", flag.control, ", Concentration RMSE: ", stats.RMSE, " n/cc", ", NRMSE: ", round(NRMSE, 2), ", (n = ", n.obs, ")")) +
            ylab(expression("dN/dlog"[10]*"D"[p]*"  ["*"cm"^-3*"]")) +
            xlab(expression("D"[p]*"  [nm]")) +
            scale_color_manual(values = cm.palette) +
            theme(
              plot.title = element_text(hjust = 0.5, size = 16),
              plot.subtitle = element_text(hjust = 0.5, size = 12),
              panel.background = element_rect(fill = "white"),
              panel.grid.major = element_line(colour = "grey80", linewidth = 0.1),
              panel.grid.minor = element_line(colour = "grey80", linewidth = 0.01),
              panel.border = element_rect(colour = "black", fill = NA),
              axis.title.y = element_text(angle = 90, vjust = 2),
              plot.margin = unit(c(1, 1, 0.1, 1), "cm"),
              legend.background = element_blank(),
              legend.box.background = element_blank(),
              legend.key = element_blank(),
              legend.title = element_blank(),
              legend.position = "right"
            ) + guides(x = guide_axis_logticks(long = 3)) +
            coord_cartesian(clip = "off")

          mid.gg <- ggplot(export.df, aes(x = Dp, y = `Residual dNdlogDp`)) +
            geom_point(shape = 1) +
            geom_segment(aes(x = Dp, xend = Dp, y = 0, yend = `Residual dNdlogDp`)) +
            scale_x_log10(breaks = 10^breaks,
                          labels = labels,
                          limits = c(min(x.limits), max(x.limits))) +
            scale_y_continuous(limits = range(export.df$`Residual dNdlogDp`)) +
            labs(title = paste0("Fraction Variance Unexplained: ", FVU)) +
            ylab(expression("Residual (dN/dlog"[10]*"D"[p]*")")) +
            xlab(expression("D"[p])) +
            scale_color_manual(values = cm.palette) +
            theme(
              plot.title = element_text(hjust = 0.5, size = 12),
              plot.subtitle = element_text(hjust = 0.5, size = 12),
              panel.background = element_rect(fill = "white"),
              panel.grid.major = element_line(colour = "grey80", linewidth = 0.1),
              panel.grid.minor = element_line(colour = "grey80", linewidth = 0.01),
              panel.border = element_rect(colour = "black", fill = NA),
              axis.title.y = element_text(angle = 90, vjust = 2),
              plot.margin = unit(c(1, 1, 0.1, 1), "cm"),
              legend.background = element_blank(),
              legend.box.background = element_blank(),
              legend.key = element_blank(),
              legend.title = element_blank(),
              legend.position = "right"
            ) + guides(x = guide_axis_logticks(long = 3))

          bot.gg <- ggplot(export.df, aes(x = `Predicted dN`, y = `Actual dN`)) +
            geom_point(shape = 1) +
            geom_abline(slope = 1) +
            scale_y_continuous(n.breaks = 5) +
            xlab(expression("Predicted Concentration "~ "n cm"^-3)) +
            ylab(expression("Actual Concentration "~ "n cm"^-3)) +
            scale_color_manual(values = cm.palette) +
            labs(title = paste0("Pearson Correlation: ", stats.R2)) +
            theme(
              plot.title = element_text(hjust = 0.5, size = 12),
              plot.subtitle = element_text(hjust = 0.5, size = 12),
              plot.caption = element_text(hjust = 0, vjust = -5),
              panel.background = element_rect(fill = "white"),
              panel.grid.major = element_line(colour = "grey80", linewidth = 0.1),
              panel.grid.minor = element_line(colour = "grey80", linewidth = 0.01),
              panel.border = element_rect(colour = "black", fill = NA),
              axis.title.x = element_text(angle = 0, vjust = -1),
              axis.title.y = element_text(angle = 90, vjust = 2),
              plot.margin = unit(c(1, 1, 1, 1), "cm"),
              legend.background = element_blank(),
              legend.box.background = element_blank(),
              legend.key = element_blank(),
              legend.title = element_blank(),
              legend.position = "right"
              )
        }

        # Merge plots
        export.gg <- top.gg / mid.gg / bot.gg
      }

      # END  -----------------------------------------------------------------------
      export.ls <- list("pass" = flag.control, "plot" = export.gg, "data" = export.df, "predict" = export.ft, "fits" = export.lm, "evaluation" = export.pf)
    })

    # Close log
    log_close()
  }

  return(export.list)
}




#' Lognormal particle size distribution
#' Uses natural log as formally defined then converts to base 10 logarithm
#'
#' @param dx
#' @param N
#' @param GSD
#' @param Dpg
#'
#' @returns Numeric vector
#'
dNdlogDp.PSD <- function(dx, N, GSD, Dpg){

  A = N/(((2*pi)^(1/2))*log(GSD))
  B = -1*(log(dx)-log(Dpg))^2
  C = 2*(log(GSD)^2)

  result = 2.303*A*exp(B/C)

  return(result)
}
