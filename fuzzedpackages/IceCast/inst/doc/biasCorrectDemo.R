## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE) 

## ----libraries, message = FALSE, warning = FALSE-------------------------
library("IceCast")
library("fields")

## ----quickstart approach, eval = FALSE-----------------------------------
#  ##Not run##
#  quick_run(obs_NCDF = "/obs.nc", pred_NCDF = "/pred.nc", pred_years = 2008,
#           start_year = 1993, month = 2, output_file = "/outputFile.nc", level = 15,
#           dat_type_obs = "bootstrap")

## ----bootstrap read-in, eval = FALSE-------------------------------------
#  ##Not run##
#  raw_data <- read_bootstrap("bt_198301_n07_v02_n.bin")

## ----how obs data obtained, eval = FALSE---------------------------------
#  ##Not run##
#  observed <- read_monthly_BS(start_year = 1993,end_year = 2007,
#                              file_folder = "myFilePath/", version = 2)
#  obs_sept <- observed[, 9, , ] #Use September data only

## ----plot Sep 2007 obs, fig.height = 5, fig.width = 5, fig.align = "center"----
image.plot(obsSep2006_2007[length(2006:2007),,], 
           main = "Observed Sea Ice Concentration \n September 2007", 
           xaxt = "n", yaxt = "n")

## ----plot Sep 2007 pred, fig.height = 5, fig.width = 5, fig.align = "center"----
image.plot(sipSep2006_2007[2,,], 
           main = "Predicted Sea Ice Concentration \n September 2007 (2.5-month lead time)",
           xaxt = "n", yaxt = "n")

## ----get land, fig.height = 5, fig.width = 5, fig.align = "center"-------
plot(land, col = "grey", main = "Land")

## ----plot misc polygons, fig.height = 5, fig.width = 5, fig.align = "center"----
plot(land, col = "grey", main = "Seas of the Arctic")
plot(bg_water, col = "black", add = T)
plot(all_regions, col = "blue", add = T)
legend("bottom", fill = c("blue", "black"), cex = 0.75,
       legend = c("Regions", "Outside Regions"))

## ----plot regions, fig.height = 5, fig.width = 5, fig.align = "center"----
colors <- c("darkblue", "green", "blue", "red", "orange", "yellow", "purple", 
            "pink",   "lightgreen", "brown", "tan", "darkgreen", "hotpink", 
            "navy", "beige", "darkblue", "green", "blue", "red", "orange",
            "yellow")
plot(land, col = "grey", main = "Mapping Lines & Regions")
nReg <- length(reg_info$regions)
for (i in 1:nReg) {
  plot(reg_info$regions[[i]], add = T, lwd = 1.5)
}
for (i in 2:nReg) {
  plot(reg_info$start_lines[[i]], col = colors[i], add = T, lwd = 2)
}

## ----plot central Arctic, fig.height = 5, fig.width = 5, fig.align = "center"----
#Find angle of mapping line (and color code)
nLines <- length(reg_info$lines[[1]])
ang <- rep(NA, nLines)
for (i in 1:nLines) {
  temp <- reg_info$lines[[1]][[i]]@lines[[1]]@Lines[[1]]@coords
  nTemp <- nrow(temp)
  ang[i] <- atan2(temp[nTemp, 2] - temp[1, 2], temp[nTemp, 1] - temp[1, 1])
} 
bp <- seq(-pi, pi, length.out = 65)
angCol <- rainbow(length(ang))

#plot region and lines 
plot(reg_info$regions[[1]], main = "Central Arctic Boundary Lines")
plot(land, add = T, col = "grey")
for (s in 1:length(reg_info$lines[[1]])) {
  plot(reg_info$lines[[1]][[s]], col = angCol[s], add = T)
}


## ---- mapAnch example, fig.align = "center", fig.height = 7, fig.width = 7, fig.align = "center"----
obs <- get_region(dat = obsSep2006_2007[length(2006:2007), ,],
                 dat_type = "bootstrap", level = 15)
obs_map <- get_map(ice = obs, plotting = TRUE, reg_info,
                 main = "Observed Mapping \n September 2007")

## ---- createMapping Feb example, fig.align = "center", fig.width = 8, fig.height = 4----
par(mfrow = c(2, 2), oma = rep(0, 4), mar = c(1, 1, 2, 1))
discrep_demo1 <- create_mapping(start_year = 2007, end_year = 2007,
                                obs_start_year = 2006, pred_start_year = 2006,
                                observed = obsSep2006_2007, predicted = sipSep2006_2007[],
                                reg_info,  month = 9, level = 15,
                                dat_type_obs = "bootstrap", dat_type_pred = "simple",
                                plotting = TRUE)

## ----look at a mapped list, fig.height = 4, fig.width = 4----------------
 head(discrep_demo1$pred_list[[1]][1,,])

## ----remove discrep demo, include = FALSE--------------------------------
rm(discrep_demo1)

## ----finding all maps, eval = FALSE--------------------------------------
#  ##Not run##
#  discrep <- create_mapping(start_year = 1993, end_year = 2007,
#                            obs_start_year = 1993, pred_start_year = 1993,
#                            observed = observed[,month,,],predicted = sip[,month,,],
#                            reg_info, month,level = 15, dat_type_obs = "bootstrap",
#                            dat_type_pred = "simple", plotting = TRUE)

## ----bias correct month--------------------------------------------------
adj <- contour_shift(maps = discrep, predicted = sipSep2008, bc_year = 2008,
                     pred_start_year = 2008, reg_info,
                     level = NA, dat_type_pred = "simple")

## ----get obs and raw-----------------------------------------------------
obs <- get_region(dat = obsSep2008, dat_type = 'bootstrap', level = 15)
un_adj <- get_region(dat = sipSep2008, dat_type = 'gfdl', level = 15)

## ----plot results, fig.align = "center", fig.height = 4.5, fig.width = 4.5----
plot(land, col = "grey", border = F,
     main = "Predicted vs. Bias−Corrected Contours \n September 2008 (2.5-Month Lead Time)")
plot(obs, col = "lightblue", add = T, border = F)
plot(un_adj, border = "red", add = T)
plot(adj, border = "navy", add = T)

## ----find error Areas----------------------------------------------------
over_est_un_adj <- gDifference(obs, un_adj)
under_est_un_adj <- gDifference(un_adj, obs)
over_est_adj <- gDifference(obs, adj)
under_est_adj <- gDifference(adj, obs)

## ----plot over- and under-estimated regions, fig.align = "center", fig.height = 4, fig.width = 4----
par(mfrow = c(1, 2), oma = rep(0, 4), mar = rep(0, 4))
#Unadjusted
plot(land, col = "grey", border = FALSE, main = "Error Regions:\n Unadjusted")
plot(obs, col = "lightblue", border = F, add = T)
plot(over_est_un_adj, col = "green", border = F, add = T)
plot(under_est_un_adj, col = "yellow", border = F, add = T)
plot(un_adj, add = T, border = "red")

#bias-corrected
plot(land, col = "grey", border = FALSE, main = "Error Regions:\n Bias-corrected")
plot(obs, col = "lightblue", border = F, add = T)
plot(over_est_adj, col = "green", border = F, add = T)
plot(under_est_adj, col = "yellow", border = F, add = T)
plot(adj, add = T, border = "navy")

## ----calculate area difference-------------------------------------------
un_adj_IIEE <- get_area(over_est_un_adj) + get_area(under_est_un_adj)
adj_IIEE <- get_area(over_est_adj) + get_area(under_est_adj)
IIEE_red <- (un_adj_IIEE - adj_IIEE)/1e5 #in 10^5 km
IIEE_red

## ----percent error reduction---------------------------------------------
per_red <- 100*(un_adj_IIEE - adj_IIEE)/un_adj_IIEE
per_red

