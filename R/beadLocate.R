#' Detect changes in bead location during a cruise in order to build a filter parameter table. 
#' 
#' @param cruise SeaFlow cruise directory name
#' @param inst Instrument serial number
#' @param save_path Path to save plots and *.filter-table.tsv files 
#' @param k Rolling median window.  Default = 7.
#' @param q_out Quantile % used to define a jump in value.  Defines the lower and upper tails. Default = 0.25.
#' @param fsc_spread is the maximum range in fsc of the identified bead peak. Default = 3000, which is quite small but seems to work.
#' @param show_plots Logical to show plots or not and save to save_path.  Default = TRUE.
#' @return bead_locs
#' @usage beadLocate(cruise, save_path, k = 7, show_plots = FALSE)

beadLocate <- function(cruise, inst, save_path, k = 7, q_out = 0.25, fsc_spread = 3000, show_plots = TRUE){
  library(popcycle)
  library(tidyverse)
  library(arrow)
  library(zoo)
  library(grid)
  library(gridExtra)
  
  path <- paste0('/Users/annettehynes/Library/CloudStorage/GoogleDrive-ahynes@uw.edu/Shared drives/SeaFlow-VCT/snakemake/')

  file <- paste0(path, cruise, "/beadfinder/", cruise, ".summary.round2.parquet") # beadfinder output
  bd <- arrow::read_parquet(file)
  bd <- bd[!is.na(bd$peak_pe_location), ] # Don't consider times when no peak is found
  bd <- bd[!is.na(bd$fsc_small_2Q), ]
  bd <- bd[!is.na(bd$D1_2Q), ]
  bd <- bd[!is.na(bd$D2_2Q), ]
  bd <- bd[(bd$fsc_small_range <= fsc_spread), ]
  
  ### rolling median:  FSC
  roll <- zoo::rollmedian(bd$peak_fsc_small_med, k = k)
  roll_date <-  zoo::rollmedian(bd$date, k = k)  # is this the best way to get the rolling date value for plotting?
  
  # Changepoint analysis.
  dy <- diff(roll)
  thresholds <- quantile(dy, c(q_out/100, 1 - q_out/100), na.rm=TRUE)
  jumps <- which(dy < thresholds[1] | dy > thresholds[2]) + 1
  
  ### rolling median:  D1 & D2/
  roll1 <- zoo::rollmedian(bd$D1_2Q, k = k)
  roll2 <- zoo::rollmedian(bd$D2_2Q, k = k)
  
  # Changepoint analysis.
  dy1 <- diff(roll1)
  thresholds <- quantile(dy1, c(q_out/100, 1 - q_out/100), na.rm=TRUE)
  jumps1 <- which(dy1 < thresholds[1] | dy1 > thresholds[2]) + 1
  
  dy2 <- diff(roll2)
  thresholds <- quantile(dy2, c(q_out/100, 1 - q_out/100), na.rm=TRUE)
  jumps2 <- which(dy2 < thresholds[1] | dy2 > thresholds[2]) + 1
  
  all_jumps <- sort(unique(c(jumps, jumps1, jumps2)))
  dj <- diff(all_jumps)
  if (min(dj) <= 2) all_jumps <- all_jumps[-(which(dj <= 2) + 1)] # Choose the smaller of two indices that are very close together
  
  ### Plotting ###
  
  g1 <- ggplot(bd) +
    geom_point(aes(x = date, y = peak_fsc_small_med), pch = 16, color = "red", alpha = 0.25) +
    geom_vline(xintercept = roll_date[all_jumps], color = "black") +
    geom_vline(xintercept = roll_date[jumps], color = "red3", linewidth = 1) +
    theme_bw(base_size = 18) +
    labs(x = "", y = "FSC")
  
  g2 <- ggplot(bd) +
    geom_point(aes(x = date, y = D1_2Q), pch = 16, color = "red", alpha = 0.25) +
    geom_vline(xintercept = roll_date[all_jumps], color = "black") +
    geom_vline(xintercept = roll_date[jumps1], color = "red3", linewidth = 1) +
    theme_bw(base_size = 18) +
    labs(x = "", y = "D1")
  
  g3 <- ggplot(bd) +
    geom_point(aes(x = date, y = D2_2Q), pch = 16, color = "red", alpha = 0.25) +
    geom_vline(xintercept = roll_date[all_jumps], color = "black") +
    geom_vline(xintercept = roll_date[jumps2], color = "red3", linewidth = 1) +
    theme_bw(base_size = 18) +
    labs(x = "", y = "D2")
  
  # Find medians in break points and build filter parameters. Beadfinder finds a much tighter peak 
  # than what is found manually with inflection_point(), so instead of using the 1st and 3rd quartiles 
  # to define the 2.5 and 97.5 quantile, I used the IQR above and below the median.
  
  limits <- c(1, all_jumps, length(roll))
  x <- rep(NA, length(all_jumps)+1)
  xend <- rep(NA, length(all_jumps)+1)
  y.hat <- rep(NA, length(all_jumps)+1)
  fsc_iqr <-  rep(NA, length(all_jumps)+1)
  y1.hat <- rep(NA, length(all_jumps)+1)
  D1_iqr <-  rep(NA, length(all_jumps)+1)
  y2.hat <- rep(NA, length(all_jumps)+1)
  D2_iqr <-  rep(NA, length(all_jumps)+1)
  
  filter_table <- NULL
  filter_plan <- data.frame(start_date = as.POSIXct(c(bd$date[1], roll_date[all_jumps])), filter_id = NA)
  
  for (i in 1:(length(all_jumps)+1)) {
    j0 <- limits[i]
    j1 <- limits[i+1]-1
    x[i] <- roll_date[j0]
    xend[i] <- roll_date[j1]
    y.hat[i] <- median(bd$peak_fsc_small_med[j0:j1])
    y1.hat[i] <- median(bd$D1_2Q[j0:j1])
    y2.hat[i] <- median(bd$D2_2Q[j0:j1])
    fsc_iqr[i] <- median(bd$fsc_small_IQR[j0:j1])
    D1_iqr[i] <- median(bd$D1_IQR[j0:j1])
    D2_iqr[i] <- median(bd$D2_IQR[j0:j1])
    
    ip <- data.frame(quantile = c(2.5, 50, 97.5))   # inflection point
    ip$fsc <- c(y.hat[i] - fsc_iqr[i], y.hat[i], y.hat[i] + fsc_iqr[i])
    ip$d1 <- c(y1.hat[i] + D1_iqr[i], y1.hat[i], y1.hat[i] - D1_iqr[i])
    ip$d2 <- c(y2.hat[i] + D2_iqr[i], y2.hat[i], y2.hat[i] - D2_iqr[i])
    
    ip[ip < 0] <- 0   # When IQR is larger than the median, can get negative results. 

    filter_params <- popcycle::create_filter_params(inst, fsc = ip$fsc, d1 = ip$d1, d2 = ip$d2, min_d1 = -5000, min_d2 = -5000, width = 15000)
    filter_id <- uuid::UUIDgenerate()
    date.stamp <- to_date_str(lubridate::now("UTC"))
    
    # From popcycle::save_filter_params
    df <- data.frame()
    for (quantile in filter_params$quantile) {
      p <- filter_params[filter_params$quantile == quantile, ]
      if (nrow(p) > 1) {
        stop("Duplicate quantile rows found in parameters passed to save_filter_params()")
      }
      df <- rbind(df, cbind(id = filter_id, date = date.stamp, 
                            quantile = quantile, beads_fsc_small = p$beads.fsc.small, 
                            beads_D1 = p$beads.D1, beads_D2 = p$beads.D2, width = p$width, 
                            notch_small_D1 = p$notch.small.D1, notch_small_D2 = p$notch.small.D2, 
                            notch_large_D1 = p$notch.large.D1, notch_large_D2 = p$notch.large.D2, 
                            offset_small_D1 = p$offset.small.D1, offset_small_D2 = p$offset.small.D2, 
                            offset_large_D1 = p$offset.large.D1, offset_large_D2 = p$offset.large.D2))
    }
    
    filter_table <- rbind(filter_table, df)
    filter_plan$filter_id[i] <- filter_id 
  } # end break point loop
  
  # Save filter table and filter plan table
  
  file_ft <- paste0(save_path, cruise, ".filter_params.filter.tsv")
  write.table(filter_table, file = file_ft, row.names = FALSE, sep = "\t")
  
  file_fplan <- paste0(save_path, cruise, ".filter_params.filter_plan.tsv")
  write.table(filter_plan, file = file_fplan, row.names = FALSE, sep = "\t")
  
  segments <- data.frame(x = as.POSIXct(x), xend = as.POSIXct(xend), fsc = y.hat, D1 = y1.hat, D2 = y2.hat)
  g1 <- g1 + geom_segment(data = segments, aes(x = x, y = y.hat, xend = xend, yend = y.hat), color = "black", linewidth = 1)
  g2 <- g2 + geom_segment(data = segments, aes(x = x, y = y1.hat, xend = xend, yend = y1.hat), color = "black", linewidth = 1)
  g3 <- g3 + geom_segment(data = segments, aes(x = x, y = y2.hat, xend = xend, yend = y2.hat), color = "black", linewidth = 1)
  
  if (show_plots){
    if(!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
  
    figname <- paste0(save_path, cruise, "_roll_med_bead_locations.png")
    png(figname, width = 1000, height = 1000)
      gridExtra::grid.arrange(g1, g2, g3, nrow = 3)
    dev.off()
  }
  
  
  # Output 
  bead_locs <- data.frame(start_date = c(bd$date[1], roll_date[all_jumps]), 
        fsc = y.hat, D1 = y1.hat, D2 = y2.hat, fsc_IQR = fsc_iqr, D1_IQR = D1_iqr, D2_IQR = D2_iqr)
  
}
