## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----load package, message = FALSE, warning = FALSE----------------------
library("IceCast")
library("viridis")
library("fields")

## ----time info-----------------------------------------------------------
month <- 9
train_start_year <- 1993
train_end_year <- 2007
level <- 15
train_bc_start_year <- 1993

## ----mcmc info-----------------------------------------------------------
n_iter <- 100
burn_in <- 20

## ----obs data info, eval = FALSE-----------------------------------------
#  ##Not run##
#  obs_file_path <- 'path_to_observation_data'
#  bs_version <- 3.1
#  dat_type_obs <- "bootstrap"
#  obs <- read_monthly_BS(start_year = train_bc_start_year, end_year = train_end_year,
#                         version = bs_version, file_folder = obs_file_path)

## ----create mapping, eval = FALSE----------------------------------------
#  ## Not run ##
#  obs_maps <- create_mapping(start_year = train_bc_start_year, end_year = train_end_year,
#                             obs_start_year = train_bc_start_year, pred_start_year = NULL,
#                             observed = obs[,month,,], predicted = NULL,
#                             reg_info, month, level, dat_type_obs,
#                             dat_type_pred = NULL, obs_only = TRUE)

## ----compute y_obs-------------------------------------------------------
y_obs <- y_obs(maps = obs_maps, reg_info)

## ----identify regions----------------------------------------------------
temp <- to_fit(y_obs, reg_info)
regs_to_fit <- temp$regs_to_fit
full <- temp$full

## ----run mcmc------------------------------------------------------------
res <- list()
for (r in regs_to_fit)  {
    start_time <- proc.time()
	  res[[r]] <- fit_cont_pars(r, n_iter, y_obs, reg_info)
    end_time <- proc.time()
    elapse_time <- end_time - start_time
    print(sprintf("MCMC for region %i finished, elapsed time %f", r, elapse_time[3]))
}

## ----estimate parameters-------------------------------------------------
pars <- list()
for (r in regs_to_fit) {
	pars[[r]] <- calc_pars(res[[r]], burn_in, w = res[[r]]$w)
}

## ----store gen info------------------------------------------------------
gen_info <- list("regs_to_fit" = regs_to_fit, "full" = full, "pars" = pars,
                 "obs_maps" = obs_maps,"train_start_year" = train_start_year,
                 "train_end_year" = train_end_year)

## ------------------------------------------------------------------------
forecast_year <- gen_info$train_end_year + 1

## ----stat_bin forecast---------------------------------------------------
indiv_stat_bin <- list()
for (r in gen_info$regs_to_fit) {
  indiv_stat_bin[[r]] <- gen_cont(r, pars_r = gen_info$pars[[r]], reg_info,
                                  stat_only = TRUE, mean_only = TRUE)
  print(sprintf("stat_bin  region %i contours generated", r))
}

## ------------------------------------------------------------------------
stat_bin <- merge_conts(conts = indiv_stat_bin, full = gen_info$full)
stat_bin <- conv_to_grid(stat_bin[[1]])

## ----colors--------------------------------------------------------------
n_color <- 8
colors <- viridis(n_color)

## ----plot stat_bin forecast, fig.height = 5, fig.width  = 5--------------
stat_bin_vis <- stat_bin #convert na's to a number for visualization (only!)
stat_bin_vis[is.na(stat_bin)] <- 1 + 1/n_color
image(stat_bin_vis, col = c(colors, "grey"), xaxt = "n", yaxt = "n",
      main = sprintf("Statistical Binary Forecast \n Month: %i, Year: %i",
                     month, forecast_year))
legend("topright", fill = c(colors[1], colors[length(colors)], "grey"),
       legend = c("no sea ice", "sea ice", "land"))

## ----gen stat probs------------------------------------------------------
n_gen <- 5
indiv_stat_prob <- list()
for (r in gen_info$regs_to_fit) {
  indiv_stat_prob[[r]] <- gen_cont(r, pars_r = gen_info$pars[[r]], reg_info,
                                   n_gen = n_gen, stat_only = TRUE)
  print(sprintf("stat_prob region %i contours generated", r))
}

## ----merge stat probs----------------------------------------------------
stat_prob <- merge_conts(conts = indiv_stat_prob, full = gen_info$full)
stat_prob <- prob_map(merged = stat_prob)

## ----plot stat_prob forecast, fig.height = 5, fig.width  = 5-------------
par(oma = c(0, 0, 0, 4))
stat_prob_vis <- stat_prob #convert na's to a number for visualization (only!)
stat_prob_vis[is.na(stat_prob)] <- 1 + 1/n_color
image(stat_prob_vis, col = c(colors, "grey"), xaxt = "n", yaxt = "n",
     main = sprintf("Statistical Probabilistic Forecast \n Month: %i, Year: %i",
                     month, forecast_year))
legend("topright", fill = c("grey"), legend = c("land"))
par(oma = c(0, 0, 0, 1))
image.plot(stat_prob, col = colors, legend.only = TRUE)

## ------------------------------------------------------------------------
lag <- 2
init_month <- get_init_month(month, lag)

## ------------------------------------------------------------------------
ecmwf_start_year <- 1993
dat_type_pred <- "simple"

## ----compute mappings for prediction, eval = FALSE-----------------------
#  ## Not run##
#  pred_maps <- create_mapping(start_year = train_bc_start_year,
#                              end_year = forecast_year - 1,
#                              obs_start_year = NULL,
#                              pred_start_year = ecmwf_start_year,
#                              observed = NULL, predicted = ecmwf_bin,
#                              reg_info, month, level, dat_type_obs = NULL,
#                              dat_type_pred = "simple", pred_only = TRUE)

## ----combine maps--------------------------------------------------------
maps <- pred_maps
maps$obs_list <- gen_info$obs_maps$obs_list

## ----contour shifting----------------------------------------------------
ppe_bin_poly <- contour_shift(maps,
                                 predicted = ecmwf_bin[length(ecmwf_start_year:forecast_year),,],
                                 bc_year = forecast_year, pred_start_year = ecmwf_start_year,
                                 reg_info, level, dat_type_pred)
ppe_bin <- conv_to_grid(ppe_bin_poly)

## ----plot ppe bin, fig.height = 5, fig.width  = 5------------------------
ppe_bin_vis <- ppe_bin #convert na's to a number for visualization (only!)
ppe_bin_vis[is.na(ppe_bin)] <- 1 + 1/n_color
image(ppe_bin_vis, col = c(colors, "grey"), xaxt = "n", yaxt = "n",
      main = sprintf("Post-Processed Ensemble Binary Forecast
                    \n Month: %i, Year: %i, Initialized Month: %i",
                      month, forecast_year, init_month))
legend("topright", fill = c(colors[1], colors[length(colors)], "grey"),
       legend = c("no sea ice", "sea ice", "land"))

## ----map ppe binary forecast---------------------------------------------
map_curr <- get_map(ice = ppe_bin_poly, reg_info)

## ----hybrd prob contours-------------------------------------------------
indiv_ppe_prob <- list()
for (r in gen_info$regs_to_fit) {
  indiv_ppe_prob[[r]] <- gen_cont(r, pars_r = gen_info$pars[[r]], reg_info,
                                     n_gen, map_pred_r = map_curr[[r]])
  print(sprintf("ppe_prob region %i contours generated", r))
}

## ----merge ppe prob contours---------------------------------------------
ppe_prob <- merge_conts(conts = indiv_ppe_prob, full = gen_info$full)
ppe_prob <- prob_map(merged = ppe_prob)

## ----plot ppe prob, fig.height = 5, fig.width  = 6-----------------------
ppe_prob_vis <- ppe_prob #convert na's to a number for visualization (only!)
ppe_prob_vis[is.na(ppe_prob)] <- 1 + 1/n_color
par(oma = c(0, 0, 0, 4))
image(ppe_prob_vis, col = c(colors, "grey"), xaxt = "n", yaxt = "n",
      main = sprintf("Post-Processed Ensemble Probabilistic Forecast
                  \n Month: %i, Year: %i, Initialized Month: %i",
                   month, forecast_year, init_month))
legend("topright", fill = c("grey"), legend = c("land"))
par(oma = c(0, 0, 0, 1))
image.plot(ppe_prob, col = colors, legend.only = TRUE)

