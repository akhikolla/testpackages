#' Apply contour-shifting to bias correct a predicted contour
#' using existing mappings.
#' @title Apply contour-shifting to bias correct
#' @param maps object obtained from the \code{create_mapping} function
#'             (see details)
#' @param predicted array of predicted values of dimension year x month x
#'                  longitude x latitude
#' @param bc_year year to be bias-corrected
#' @param pred_start_year year prediction array starts in
#' @param reg_info  a \code{reg_info} list (see documentation for \code{reg_info})
#' @param level concentration level for which to build contour
#' @param dat_type_pred string indicating the format of the prediction: either
#'                    "gfdl" or "simple" (see details)
#' @param my_land_mat binary matrix specifying land locations
#' @param my_land \code{SpatialPolygons} corresponding to the land
#' @param n_train_years number of years prior to the current year used in fitting
#'                    the bias correction
#' @return \code{SpatialPolygons} object of the adjusted region
#' @export
#' @details The \code{map} parameter is a list of length four that has the form
#'          of a list obtained from running the \code{create_mapping} function.
#'          The values \code{start_year} and \code{endYear} give the
#'          first and last year that were mapped. The variables
#'          \code{obs_list} and \code{pred_list} are lists of arrays with one
#'          3-dimensional array for each region. The first dimension of each array
#'          coresponds to the year, the second dimension corresponds to the lines
#'          on which the region is being mapped, and the third dimension
#'          corresponds to the variables of interest. The first and second
#'          dimension are indexed sequentially. The variables for the third dimension
#'          are for the fixed points' x-coordinates, the fixed points' y-coordinates,
#'          the mapped points' x-coordinates, the mapped points' y-coordinates,
#'          the length of the mapping vectorsin the x-direction,
#'          the length of the vectors in the y-direction, and the angles of the
#'          mapping vectors.

#'          The predicted data array, \code{predicted}, should be a single array
#'          of dimension: years x longitude (304) x latitude (448). If
#'          \code{dat_type_pred = ``simple"}, the values in the array should
#'          indicate whether each grid box is categorized to contain ice
#'          (1: ice-covered, 0: no ice, NA: land). If
#'          \code{dat_type_pred =``gfdl"} the
#'          values in the \code{predicted} array correspond to the raw ice
#'          concentrations values predicted (including indicators for missing
#'          data, land etc.) formatted as in the CM2.5 Forecast-oriented
#'          Low-Ocean Resolution (FLOR) model produced by the National Oceanic
#'          and Atmospheric Administration’s Geophysical Fluid Dynamics
#'          Laboratory and converted to a Polar Stereographic grid
#'          (Vecchi et al. 2014; Msadek et al. 2014). Weights for converting to
#'          a polar stereograhic grid were obtained from the spherical
#'          coordinate remapping and interpolation package (SCRIP) (Jones 1997).
#'
#'
#' @references
#'
#' Jones, P.W. "A user’s guide for SCRIP: A spherical coordinate remapping and interpolation package."
#' Los Alamos National Laboratory, Los Alamos, NM (1997).
#'
#' Msadek, R., et al.
#' \href{http://onlinelibrary.wiley.com/doi/10.1002/2014GL060799/full}{"Importance of initial conditions in seasonal predictions of Arctic sea ice extent."}
#'  Geophysical Research Letters 41.14 (2014): 5208-5215.
#'
#' Vecchi, Gabriel A., et al.
#' \href{http://journals.ametsoc.org/doi/abs/10.1175/JCLI-D-14-00158.1}{"On the seasonal forecasting of regional tropical}
#' cyclone activity." Journal of Climate 27.21 (2014): 7994-8016.
#'
#' @importFrom raster aggregate
#' @importFrom sp disaggregate Polygon
#' @importFrom rgeos gArea
#' @importFrom MASS rlm
#' @importFrom stats predict
#' @examples
#' \dontrun{
#' adj <- contour_shift(maps = discrep, predicted = emFeb2012, bc_year = 2012,
#'                      pred_start_year = 1980, reg_info, level = 15,
#'                      dat_type_pred = "gfdl")
#' plot(land, col = "grey", border = FALSE)
#' plot(adj, add = TRUE, col = "blue")
#' }
contour_shift <- function(maps, predicted, bc_year, pred_start_year, reg_info,
                          level, dat_type_pred, my_land_mat = land_mat,
                          my_land = land, n_train_years = NULL) {

  ##Read-in and format prediction
  raw <- get_region(dat = predicted, dat_type = dat_type_pred, level)

  #indices of data to use
  if (is.null(n_train_years)) {
    train_ind <- (1:(bc_year - maps$start_year))
  } else {
    train_ind <- ((bc_year - maps$start_year - n_train_years  + 1):(bc_year - maps$start_year))
  }

  #Check for and remove any indices for years with no observed data
  na_only <- sapply(maps$obs_list, function(x){apply(x, 1, function(y){all(is.na(y))})})
  missing <- which(apply(na_only, 1, function(x){all(x)}))
  if (length(missing) > 0) {
    train_ind <- train_ind[-missing]
  }

  #corresponding years to use
  if (is.null(n_train_years)) {
    year_ind <- maps$start_year:(bc_year - 1)
  } else {
    year_ind <- (bc_year -n_train_years):(bc_year - 1)
  }
  if (length(missing) > 0) {
    year_ind <- year_ind[-missing]
  }

  n_reg <- length(reg_info$regions)

  ##find mappings for raw prediction
  raw_map <- get_map(ice = raw, reg_info)

  end <- list()
  ##Adjust each region
  for (r in 1:n_reg) {
    angs_r <- reg_info$angs[[r]]
    n_lines_r <- length(reg_info$lines[[r]])
    start_coords_r <- reg_info$start_coords[[r]]
    ##adjust Central Arctic region
    end_r <- matrix(ncol = 2, nrow = n_lines_r)
    for (s in 1:n_lines_r) {
      #calculate length of lines
      x_temp_raw <- maps$pred_list[[r]][train_ind, s, 5]
      y_temp_raw <- maps$pred_list[[r]][train_ind, s, 6]
      raw_length <- sqrt(x_temp_raw^2 + y_temp_raw^2)
      x_temp_obs <- maps$obs_list[[r]][train_ind, s, 5]
      y_temp_obs <- maps$obs_list[[r]][train_ind, s, 6]
      obs_length <- sqrt(x_temp_obs^2 + y_temp_obs^2)

      #use robust regression to get adjustment
      if (length(unique(raw_length)) > 1) {
        lm_raw <- suppressWarnings(rlm(raw_length ~ year_ind))
        raw_pred <- predict(lm_raw, newdata = data.frame(year_ind = bc_year))
      } else {
        raw_pred <- raw_length[1]
      }
      if (length(unique(obs_length)) > 1) {
        lm_obs <- suppressWarnings(rlm(obs_length ~ year_ind))
        obs_pred <- predict(lm_obs, newdata = data.frame(year_ind = bc_year))
      } else {
        obs_pred <- obs_length[1]
      }
      #calculate new length of mapping vectors
      x_raw <- raw_map[[r]][s, 5]
      y_raw <- raw_map[[r]][s, 6]
      raw_dist <- sqrt(x_raw^2 + y_raw^2)
      adj_dist <- raw_dist + (obs_pred - raw_pred)

      #new x and y coordinates
      end_r[s, 1] <- start_coords_r[s, 1] + adj_dist*cos(angs_r[s])
      end_r[s, 2] <- start_coords_r[s, 2] + adj_dist*sin(angs_r[s])
    }
    end[[r]] <- end_r
  }

  ##make polygon for central Arctic region
  adj <- make_polygons(r = 1, my_end = end[[1]], poly_name = "new1",
                       loop_r = reg_info$loop[1])

  ##make adjusted polygons for all other regions
  for (r in 2:n_reg) {
    new <- make_polygons(r, my_end = end[[r]], poly_name = sprintf("new%i", r),
                         loop_r  = reg_info$loop[r])
    if (!is.null(new)) {
      adj <- aggregate(spRbind(adj, new))
      adj@polygons[[1]]@ID <- "new"
    }
  }

  #Add back ice islands
  raw <- disaggregate(raw)
  for (r in 1:n_reg) {
    ice_curr <- keep_poly(gIntersection(raw, reg_info$regions[[r]]))
    if (!is.null(ice_curr)) {
      ice_curr <- disaggregate(ice_curr)
      isl <- ice_curr[!gIntersects(ice_curr, reg_info$start_lines[[r]], byid = T)]
      if (!is.null(isl)) {
        spRbind(adj, isl)
      }
    }
  }

  return(adj)
}
