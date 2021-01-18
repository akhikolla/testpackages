#' Generate the contours for a particular region given the model prediction
#' @title Generate contours
#' @param r integer indicating the number of the region in which the contours
#'          should be generated
#' @param pars_r List of parameter information for region r. The list should contain
#'            two elements, \code{muEst} and \code{sigmaEst}, which give
#'            estimates for the \eqn{\mu} and \eqn{\Sigma} parameters used in
#'            generating contours. Typically obtained from the \code{calc_pars}
#'            function
#' @param reg_info  a \code{reg_info} list (see documentation for \code{reg_info})
#' @param n_gen integer specifying the number of contours to be generated, must
#'              be at least 2
#' @param map_pred_r output of \code{get_map} function applied to
#'                   \code{SpatialPolygons} object corresponding to an intial
#'                   forecast (typically a bias-corrected dynamic ensemble forecast)
#' @param stat_only boolean indicating that forecast is purely statistical
#'                  (no dynamic ensemble model forecast considered)
#' @param mean_only boolean indicating that only the mean contour will be
#'                  computed rather than distribution
#' @param eff_zero how close a generated vector needs to be to zero to be counted as a
#' zero, defaults to 12.5
#' @param stat_only_trend boolean indicating if a trend adjustment should
#' be applied when \code{stat_only = TRUE}. Defaults to true
#' @importFrom rgeos gIntersection gDifference
#' @importFrom MASS mvrnorm
#' @importFrom sp SpatialPolygons
#' @examples
#' \dontrun{
#' #statistical binary, region 1
#' stat_bin_1 <- gen_cont(r = 1, pars_r = pars_1, reg_info,
#'                        stat_only = TRUE, mean_only = TRUE)
#'
#' #statistical probabilistic, region 1, 2 generated contours
#' stat_prob_1 <- gen_cont(r = 1, pars_r = pars_1, reg_info,
#'                         n_gen = 2, stat_only = TRUE)
#'
#'  #hybrid probabilistic, region 1, 2 generated contours
#'  hybrid_prob_1 <- gen_cont(r = 1, pars_r = pars_1, reg_info,
#'                            n_gen = 2, map_pred_r = map_curr_1)
#'  }
gen_cont <- function(r, pars_r, reg_info, n_gen = NULL, map_pred_r = NULL,
                     stat_only = FALSE, mean_only = FALSE, eff_zero = 12.5,
                     stat_only_trend = TRUE) {
  #must generate at least 2 contours
  if (!mean_only) {
    stopifnot(n_gen > 1)
  }

  #should not be providing dynamic forecast prediction if making a mean-only
  #stat-only forecast
  if (stat_only) {
    stopifnot(is.null(map_pred_r))
  }
  if (mean_only) {
    stopifnot(is.null(map_pred_r))
  }

  #Relevant values
  angs_r <- reg_info$angs[[r]]
  loop_r <- reg_info$loop[[r]]

  #get line lengths from prediction
  if (stat_only) {
    y_pred <- pars_r$mu_est
  } else {
    y_pred <- sqrt((map_pred_r[,5])^2 + (map_pred_r[,6])^2)
  }

  #generate x and y's
  if (mean_only) {
    y_gen <- matrix(ncol = 1, data = y_pred)
    n_gen <- 1
  } else {
    x_gen = mvrnorm(n_gen, y_pred, pars_r$sigma_est)
    y_gen <- apply(x_gen, 1, function(s){censor(s, reg_info$dist[[r]])})
    y_gen[y_gen <= eff_zero] <- 0
  }

  #make contours out of x's and y's
  conts_r <- list()
  for (k in 1:n_gen) {
    new_pts <- reg_info$start_coords[[r]] + cbind(y_gen[,k]*cos(angs_r), y_gen[,k]*sin(angs_r))
    conts_r[[k]] <- make_polygons(r, my_end = new_pts, poly_name = "new", loop_r)
  }
  return(conts_r)
}

#' Take in unbounded line lengths (x-values) and truncate them based on provided
#' bounds to create bounded lengths (y-values)
#' @title Truncate simulated line lengths based on a list of bounds
#' @param x a vector of generated line lengths
#' @param bounds a vector which gives the lengths of the end and breakpoints
#'               for each x value
#' @return vector of new line lengths
censor <- function(x, bounds) {
  stopifnot(length(x) == length(bounds))
  n_lines <- length(x)
  y <- x
  for (k in 1:n_lines) {
    m <- bounds[[k]] #bounds associated with line k
    n_bound <- length(m)
    if (x[k] < m[1]) {#lower than min length, set to zero
      y[k] <- 0
    } else if (x[k] >= m[n_bound]) {#higher than max, set to b
      y[k] <- m[n_bound]
    } else if (n_bound > 2) {
      #cycle over non-endpoint bounds, assign to closest one
      for (i in seq(2, n_bound - 1, by = 2)) {
        if ((x[k] >= m[i]) & (x[k] <= m[i] + (m[i + 1] - m[i])/2)) {
          y[k] <- m[i]
        } else if ((x[k] <= m[i + 1]) & (x[k] >= m[i] + (m[i + 1] - m[i])/2)) {
          y[k] <- m[i + 1]
        }
      }
    }
  }
  return(y)
}

#' Merge generated contours for all regions together
#' @title Merge contours
#' @param conts list of contours organized as a list of regions by a list of
#' years by a list of samples
#' @param full \code{SpatialPolygons} object for area to be included in all
#'             generated contours
#' @return Returns a list of contours organized as a list of years by a list of
#'  samples
#' @importFrom raster aggregate
#' @importFrom maptools spRbind
#' @export
merge_conts <- function(conts, full) {
  first <- min(which(sapply(conts, function(x){!is.null(x)})))
  n_gen <- length(conts[[first]])
  merged <- list()
  reg_fit <- which(sapply(conts, function(x){!is.null(x)}))
  for (i in 1:n_gen) {
    temp <- NULL
    for (r in reg_fit) {
      if (i <= length(conts[[r]])) {#if last polygon is empty, no value for index i
        if (!is.null(conts[[r]][[i]])) { #ignore NULL (empty) polygons
          if (!is.null(temp) && !is.null(conts[[r]][[i]])) { #add new polygon
            temp <- aggregate(spRbind(temp, conts[[r]][[i]]))
          } else { #first polygon
            temp <- conts[[r]][[i]]
          }
        }
        temp@polygons[[1]]@ID <- "genCont"
      }
    }
    if (is(full)[1] == "SpatialPolygons") {
      temp <- aggregate(spRbind(temp, full))
    } else {
      temp <- aggregate(temp)
    }
    merged[[i]] <- temp
  }
  return(merged)
}


#' Takes in list of polygon objects from merged function and produces a
#' map of probabilities
#' @title Get probabilities on a grid from contours
#' @param merged list of contours organized as a list of years by a list of
#'               samples
#' @param nX dimension in the x (defaults to value for Northern Polar
#'           stereographic grid: 304)
#' @param nY dimension in the y (defaults to value for Northern Polar
#'           stereographic grid: 448)
#' @return array of dimension number of years by longitude by latitude that
#' gives the proportion of contours in which the grid box is ice-covered
#' @export
#' @examples
#' \dontrun{ probs <- prob_map(merged) }
#'
prob_map <- function(merged, nX = 304, nY = 448) {
  n_gen <- length(merged)
  indiv <- array(dim = c(n_gen, nX, nY))
  for (i in 1:n_gen) {
    indiv[i,,] <- conv_to_grid(merged[[i]])
  }
  map <- apply(indiv, 2:3, mean)
  return(map)
}

#' Trend Adjustment For \code{mu}
#' @param obs_list partial output of \code{get_map function},
#' \code{maps$obs_list[[r]]}, where r is the region of interest
#' @param forecast_year year to be forecast
#' @param train_start_year first year in training period
#' @param train_end_year last year in training period
#' @return vector of the length of the number of lines in the mapping that
#' represent by what factor each estimated \code{mu} should be adjusted
ts_adj_mu <- function(obs_list, forecast_year, train_start_year, train_end_year) {
  mid_year <- (train_start_year + train_end_year)/2

  #indices of data to use
  year_ind <- train_start_year:train_end_year
  train_ind <- 1:(length(year_ind))

  #Check for and remove any indices for years with no observed data
  na_only <- apply(obs_list, 1, function(x){all(is.na(x))})
  missing <- which(na_only)
  if (length(missing) > 0) {
    train_ind <- train_ind[-missing]
  }

  #corresponding years to use
  if (length(missing) > 0) {
    year_ind <- year_ind[-missing]
  }

  n_lines <- dim(obs_list)[2]
  adj_prop <- rep(NA, n_lines)
  ##Compute adjustment for all lines
  for (s in 1:n_lines) {
    #calculate length of lines
    x_temp_obs <- obs_list[train_ind, s, 5]
    y_temp_obs <- obs_list[train_ind, s, 6]
    obs_length <- sqrt(x_temp_obs^2 + y_temp_obs^2)

    #use robust regression and compute estimated length at forecast year
    #divided by estimate at mid-point year of training period
    if (length(unique(obs_length)) > 1) {
      lm_obs <- suppressWarnings(rlm(obs_length ~ year_ind))
      length_pred <- predict(lm_obs,
                              newdata = data.frame(year_ind = c(mid_year, forecast_year)))
      length_pred[length_pred < 0] <- 0
      adj_prop[s] <- length_pred[2]/length_pred[1]

     } else {
      adj_prop[s] <- 1
    }
  }
  return(adj_prop)
}



