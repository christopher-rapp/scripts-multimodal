#' Multi-modal lognormal particle size distribution algorithm
#'
#' @param data.ls
#' @param max.iterations
#'
#' @returns
#' @export
#'
#' @examples
multimodal.fitting <- function(data.ls, max.iterations){

  if(missing(max.iterations)) max.iterations <- 10

  f.dNdlogDp.curve = function(dx, N, GSD, Dpg){

    A = N/(((2*pi)^(1/2))*dx*log(GSD))
    B = -1*(log(dx)-log(Dpg))^2
    C = 2*(log(GSD)^2)

    result = 2.303*dx*A*exp(B/C)

    return(result)
  }

  export.gg <- list()
  export.ls <- list()
  for (z in 1:length(data.ls)){

    {
      # Retrieve data from function input
      {
        data = data.ls[[z]]

        time.index = str_which(colnames(data), paste0(c("UTC", "GMT"), collapse = '|'))

        if (time.index >= 1){
          print(paste0("Current Dataset: ", as.POSIXct(as.numeric(first(data[, ..time.index])), tz = "UTC")))
        } else {
          print("Please check if a UTC time is available in your datasets")
        }

        # Find numeric bins
        bins.ix <- suppressWarnings(which(!is.na(as.numeric(colnames(data)))))
        bins.nm <- as.numeric(colnames(data)[bins.ix])

        # Subset all non binned data from dataframe
        tmp.df <- data %>%
          select(!all_of(bins.ix))

        # Calculate log difference of diameters
        {
          # Stagger the bins i.e trim end points for each to subtract
          tmp1 <- as.numeric(bins.nm[-c(length(bins.nm))])/1000
          tmp2 <- as.numeric(bins.nm[-1])/1000

          # Find log difference between bin diameters
          dlogDp = log10(tmp2) - log10(tmp1)
        }

        # Create a vector to select columns
        # Due to needing to calculate dlogDp the last column is dropped
        select.c <- bins.ix[1:(length(bins.ix) - 1)]

        # Subset data and create a numeric matrix
        dNdlogDp.mat <- data[, ..select.c]
        dNdlogDp.mat <- sapply(dNdlogDp.mat, as.numeric)

        {
          # Remove rows with erroneous values from the temporary numeric matrix
          remove.c <- which(rowSums(dNdlogDp.mat, na.rm = T) > 1e10)

          # Filter out extremely high erroneous values
          dNdlogDp.mat[remove.c, ] <- NA
        }

        # Calculate number distribution by multipling by the log difference to convert dNdlogDp to dN
        dN.mat = sweep(dNdlogDp.mat, 2, dlogDp, FUN = "*")

        # Remove empty rows
        dN.mat <- dN.mat[complete.cases(dN.mat), ]
      }

      # PRE-MODEL FILTERING
      {
        if (sum(dN.mat) < 100){
          next
        }
      }

      # Peak Identification
      # Levenberg-Marquardt Nonlinear Least-Squares Algorithm
      {
        # MODEL SETUP
        {
          peak.fitting <- list()
          model.fitting <- list()

          # Create temporary dataset used to initialize the fitting algorithm
          tmp.data = data.frame(diameter = as.numeric(names(colMeans(dN.mat, na.rm = T))), dlogDp = dlogDp, dN = colMeans(dN.mat, na.rm = T), dNdlogDp = colMeans(dNdlogDp.mat, na.rm = T))

          # Initial parameters to initialize loop
          FUV = 100
          max.iterations = 10
          i = 1
        }

        # MODEL RUN
        while (FUV > 1 & i < max.iterations){

          {
            if (i == 1){
              tmp.data$residual = tmp.data$dNdlogDp
            }

            # Find all peaks
            {
              # Peak identification
              peaks.df <- as.data.frame(pracma::findpeaks(tmp.data$residual, minpeakdistance = 5, sortstr = T, npeaks = 1))

              # Use indices to match to corresponding diameters
              peaks.df$Mode <- tmp.data$diameter[peaks.df$V2]
              peaks.df$Lower <- tmp.data$diameter[peaks.df$V3]
              peaks.df$Upper <- tmp.data$diameter[peaks.df$V4]
            }

            # Peak statistics
            {
              # Total Number Concentration
              {
                # Find both dNdlogDp and dlogDp to convert to dN
                tmp.dNdlogDp = tmp.data$residual[peaks.df$V3:peaks.df$V4]
                tmp.dlogDp = tmp.data$dlogDp[peaks.df$V3:peaks.df$V4]

                peaks.df$Max = max(tmp.data$dNdlogDp[peaks.df$V3:peaks.df$V4])

                # Add column for total number concentration (N)
                peaks.df$N = sum(tmp.dNdlogDp*tmp.dlogDp)

                # Prevent negative values from appearing in GSD calculation
                tmp.dNdlogDp[which(tmp.dNdlogDp < 0)] <- 0
                tmp.dlogDp[which(tmp.dlogDp < 0)] <- 0
              }

              # GSD calculation
              {
                tmp.Dp = tmp.data$diameter[peaks.df$V3:peaks.df$V4]

                GSD = suppressWarnings({
                  10^(sqrt(sum(tmp.dNdlogDp*((log10(tmp.Dp) - log10(peaks.df$Mode))^2))/(peaks.df$N - 1)))
                })

                peaks.df$GSD = GSD
              }
            }

            # Temporary dataframe for LM-NLS model
            # This modifies the "view" NLS has when trying to curve fit
            tmp <- data.frame(x = tmp.data$diameter[peaks.df$V3:peaks.df$V4], y = tmp.data$residual[peaks.df$V3:peaks.df$V4])

            # Try fitting an NLS model to the residuals (or original data if first iteration)
            {
              NLS.break <<- FALSE

              # Try fitting an NLS model to the residuals (or original data if first iteration)
              tryCatch(expr = {
                NLS.MODEL <- suppressWarnings(minpack.lm::nlsLM(y ~ f.dNdlogDp.curve(x, N, GSD, Dpg), data = tmp, start = list(N = peaks.df$N, GSD = peaks.df$GSD, Dpg = peaks.df$Mode)))
              },
              error = function(e){NLS.break <<- TRUE}
              )

              # Stop code from continuing if there are no ramps detected which means the dataset will be empty and continue to throw errors
              if (NLS.break) {break("Terminating loop, model fitting failure")}

              # Predict
              tmp.predict = data.frame(x = seq(0.01, 1500, by = 0.01))
              tmp.predict[[paste0("Mode ", i)]] <- predict(NLS.MODEL, newdata = tmp.predict)
            }
          }

          tmp.xL <- round(tmp.predict$x, 2)
          tmp.xR <- round(tmp.data$diameter, 2)

          index.match <- which(as.character(tmp.xL) %in% as.character(tmp.xR))

          predict.diameter = tmp.predict[index.match, 1]
          predict.dNdlogDp = tmp.predict[index.match, 2]

          peak.height.ratio = max(tmp.predict[, 2])/peaks.df$Max

          if (peak.height.ratio > 1.5){

            tmp.data$residual[which(tmp.data$diameter < peaks.df$Upper & tmp.data$diameter > peaks.df$Lower)] <- 0

          } else {

            tmp.data$residual <- tmp.data$residual - predict.dNdlogDp

            # Add model values to output list
            model.fitting[[i]] <- tmp.predict

            # Add peak values to output list
            peak.fitting[[i]] <- peaks.df[, -c(1:4)]
          }

          FUV = round((var(tmp.data$residual)/var(tmp.data$dNdlogDp))*100, 2)

          print(paste0("Current Iteration: ", i, ", Remaining Variance: ", FUV, "%"))
          i <- i + 1
        }
      }

      # Output list of peak data from loop
      {
        peak.fitting <- purrr::compact(peak.fitting)
        peaks.df <- rbindlist(peak.fitting)
        peaks.df <- peaks.df %>%
          select(`Lower`, `Mode`, everything()) %>%
          arrange(desc(`Mode`))
      }

      # Output list of model fitting from loop
      {
        model.fitting <- purrr::compact(model.fitting)
        output <- do.call(cbind, model.fitting)
        output <- output[append(1, str_which(colnames(output), 'Mode'))]
        output <- round(output, 2)

        setnames(output, old = colnames(output)[str_which(colnames(output), pattern = "Mode \\d{1,2}")],
                 new = paste0("Mode ", seq_along((str_which(colnames(output), pattern = "Mode \\d{1,2}")))))

        if (ncol(output) > 2){
          output$Total <- rowSums(output[, 2:ncol(output)])
        } else {
          output$Total <- output[, 2]
        }

      }

      diameter <- output$x[which(as.character(output$x) %in% as.character(tmp.data$diameter))]
      predicted <- output$Total[which(as.character(output$x) %in% as.character(tmp.data$diameter))]
      actual <- tmp.data$dNdlogDp

      performance.df <- data.frame(predicted, actual) %>%
        mutate(`ratio` = predicted/actual) %>%
        mutate(`difference` = actual - predicted) %>%
        mutate(`diameter` = diameter)

      R2 = round(cor(performance.df$predicted, performance.df$actual), 4)

      # Plotting
      {
        plot.df <- output %>%
          pivot_longer(
            cols = !c("x"),
          )

        setnames(plot.df, new = c("Diameter", "Mode", "Concentration"))

        # Use log distributed values rather than numerical to improve performance
        # Needs several orders of magnitude less data to produce smooth curves
        x = 10^seq(log10(1), log10(1000), length.out = 1000)

        x.limits = 10^((1):(3))

        x.label.positions = log10(x.limits)

        breaks <- NULL
        labels <- NULL
        for (i in x.label.positions){
          breaks <- append(breaks, i)
          labels <- append(labels, as.character(10^(i)))
        }

        cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

        cm.palette <- append(cbPalette[1:(ncol(output) - 2)], "black")

        top.gg <- ggplot(plot.df, aes(x = Diameter, y = Concentration, color = Mode)) +
          geom_point(data = tmp.data, aes(x = diameter, y = dNdlogDp), inherit.aes = F, shape = 1) +
          geom_line(linewidth = 0.5) +
          scale_x_log10(breaks = 10^breaks,
                        labels = labels,
                        limits = c(min(x.limits), max(x.limits))) +
          scale_y_continuous(n.breaks = 5) +
          ylab(expression("dN/dlog"[10]*"D"[p])) +
          scale_color_manual(values = cm.palette) +
          theme(
            plot.title = element_text(hjust = 0.5, size = 16),
            plot.subtitle = element_text(hjust = 0.5, size = 12),
            panel.background = element_rect(fill = "white"),
            panel.grid.major = element_line(colour = "grey80", linewidth = 0.1),
            panel.grid.minor = element_line(colour = "grey80", linewidth = 0.01),
            panel.border = element_rect(colour = "black", fill = NA),
            axis.title.x = element_blank(),
            axis.title.y = element_text(angle = 90, vjust = 2),
            plot.margin = unit(c(1, 1, 0.1, 1), "cm"),
            legend.background = element_blank(),
            legend.box.background = element_blank(),
            legend.key = element_blank(),
            legend.title = element_blank(),
            legend.position = "right"
          ) + guides(x = guide_axis_logticks(long = 3))

        mid.gg <- ggplot(performance.df, aes(x = diameter, y = difference)) +
          geom_point(shape = 1) +
          geom_segment(aes(x = diameter, xend = diameter, y = 0, yend = difference)) +
          scale_x_log10(breaks = 10^breaks,
                        labels = labels,
                        limits = c(min(x.limits), max(x.limits))) +
          scale_y_continuous(n.breaks = 10) +
          labs(title = paste0("Fraction Unexplained Variance: ", FUV)) +
          ylab("Difference") +
          scale_color_manual(values = cm.palette) +
          theme(
            plot.title = element_text(hjust = 0.5, size = 16),
            plot.subtitle = element_text(hjust = 0.5, size = 12),
            panel.background = element_rect(fill = "white"),
            panel.grid.major = element_line(colour = "grey80", linewidth = 0.1),
            panel.grid.minor = element_line(colour = "grey80", linewidth = 0.01),
            panel.border = element_rect(colour = "black", fill = NA),
            axis.title.x = element_blank(),
            axis.title.y = element_text(angle = 90, vjust = 2),
            plot.margin = unit(c(1, 1, 0.1, 1), "cm"),
            legend.background = element_blank(),
            legend.box.background = element_blank(),
            legend.key = element_blank(),
            legend.title = element_blank(),
            legend.position = "right"
          ) + guides(x = guide_axis_logticks(long = 3))

        bot.gg <- ggplot(performance.df, aes(x = predicted, y = actual)) +
          geom_point(shape = 1) +
          geom_abline(slope = 1) +
          scale_y_continuous(n.breaks = 5) +
          xlab(expression("Predicted Concentration "~ "n cm"^-3)) +
          ylab(expression("Actual Concentration "~ "n cm"^-3)) +
          scale_color_manual(values = cm.palette) +
          labs(title = paste0("Pearson Correlation: ", R2)) +
          theme(
            plot.title = element_text(hjust = 0.5, size = 16),
            plot.subtitle = element_text(hjust = 0.5, size = 12),
            panel.background = element_rect(fill = "white"),
            panel.grid.major = element_line(colour = "grey80", linewidth = 0.1),
            panel.grid.minor = element_line(colour = "grey80", linewidth = 0.01),
            panel.border = element_rect(colour = "black", fill = NA),
            axis.title.x = element_text(angle = 0, vjust = -1),
            axis.title.y = element_text(angle = 90, vjust = 2),
            plot.margin = unit(c(1, 1, 0.1, 1), "cm"),
            legend.background = element_blank(),
            legend.box.background = element_blank(),
            legend.key = element_blank(),
            legend.title = element_blank(),
            legend.position = "right"
          )

        export.gg[[z]] <- top.gg / mid.gg / bot.gg

        export.ls[[z]] <- list(export.gg, output, performance.df)
      }
    }
  }

  return(list(export.gg, export.ls))
}
