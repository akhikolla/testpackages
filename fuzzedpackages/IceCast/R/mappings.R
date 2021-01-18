#' Find the mapping vectors for one observation or prediction.
#' @title Map one observation or prediction
#' @param ice \code{SpatialPolygons} object corresponding to the region of ice
#' @param reg_info  \code{reg_info} list (see \code{reg_info} documentation)
#' @param plotting boolean indicating if map should be plotted
#' @param main string specifying the name for the plot
#' @param my_land \code{SpatialPolygons} object corresponding to the land
#' @return List of the length of the number of regions. Each item in the list is
#'         a matrix. Each row of each matrix corresponds to a point in the
#'         region's line. The six columns give the fixed point's x-coordinate,
#'         the fixed point's y-coordinate, the mapped point's x-coordinate, the
#'         mapped point's y-coordinate, the length of the mapping fvectors in the
#'         x-direction, and the length of the vectors in the y-direction.
#' @importFrom sp SpatialPoints disaggregate
#' @importFrom rgeos gIntersection gIntersects gLineMerge
#' @importFrom maptools spRbind
#' @importFrom methods as
#' @importFrom graphics lines points
#' @importFrom sp plot
#' @importFrom utils combn
#' @examples
#' \dontrun{
#' obs <- get_region(dat = obsFeb19811982[1,,], dat_type = "bootstrap",
#'                level = 15)
#' obs_map <- get_map(ice = obs, plotting = TRUE, reg_info,
#'                  main = "Observed Mapping \n February 1985")
#' }
#' @export
get_map <- function(ice, reg_info, plotting = FALSE, main = "", my_land = land) {

  ##cycle through regions to get their mappings
  first <- TRUE
  n_reg <- length(reg_info$regions)
  map_list <- list()
  for (r in 1:n_reg) {
    ice_curr <- keep_poly(gIntersection(ice, reg_info$regions[[r]]))
    n_lines_r <- length(reg_info$lines[[r]])
    map <- matrix(nrow = n_lines_r, ncol = 6)
    colnames(map) <- c("fromX", "fromY", "toX", "toY", "deltaX","deltaY")
    map[, c("fromX", "fromY")] <- reg_info$start_coords[[r]]
    if (!is.null(ice_curr)) {
      ice_curr <- disaggregate(keep_poly(ice_curr))
      if (r != 1) {
        inter_test <- gIntersects(ice_curr, reg_info$start_lines[[r]], byid = T)
        if (any(inter_test)) {
          ice_curr <- ice_curr[inter_test]
        } else {
          ice_curr <- NULL
        }
      } else {
        ice_curr <- ice_curr[which.max(gArea(ice_curr, byid = T))]
      }
      if (!is.null(ice_curr)) {
        for (l in 1:n_lines_r) {
          inter_pts <- keep_line(gIntersection(ice_curr, reg_info$lines[[r]][[l]]))
          if (!is.null(inter_pts)) {
            inter_pts <- get_coords(inter_pts)
            ep <- which.max(apply(inter_pts, 1,
                                  function(x){get_dist(x, reg_info$start_coords[[r]][l,])}))
            map[l, c("toX", "toY")] <- inter_pts[ep,]
          } else {
            map[l,c("toX", "toY")] <- map[l, c("fromX", "fromY")]
          }
        }
      } else {
        map[,c("toX", "toY")] <- map[, c("fromX", "fromY")]
      }
    } else {
      map[,c("toX", "toY")] <- map[, c("fromX", "fromY")]
    }

    map[, "deltaX"] <- map[, "toX"] - map[, "fromX"]
    map[, "deltaY"] <- map[, "toY"] - map[, "fromY"]

    ##optional plotting of mapping
    if (first & plotting) {
      plot(my_land, col = "grey", main = main)
      first <- FALSE
    }

    #plotting every nth point on the fixed line and mapping lines
    if (plotting) {
      #divide up (-pi, pi) into even break points
      if (!is.null(ice_curr)) {
        plot(ice_curr, lwd = 2, add = T)
      }
      points(map[,c("fromX", "fromY")], col = "black",
             pch = 20, cex = .25)
      #Note: Appropriate-sized arrows can be used in place of lines with the
      #Arrows function in the  "shape" package
      sapply((1:n_lines_r), function(s){
        lines(x = c(map[s,"fromX"], map[s,"toX"]),
              y = c(map[s,"fromY"], map[s,"toY"]), lwd = .7)})
    }
    map_list[[r]] <- map
  }
  return(map_list)
}

#' The function evenly spaces the number of points that are on one line,
#' \code{pred_l}, on a different line, \code{obs_l}
#' @title Space points along a line
#' @param pred_l predicted line (n1 x 2 matrix of coordinates)
#' @param obs_l predicted line (n2 x 2 matrix of coordinates)
#' @param plotting boolean indicating whether maps should be plotted
#' @return n x 2 matrix of evenly-spaced coordinates
#' @export
#' @examples
#' line_space <- int_line(predLEx, obsLEx, plotting = TRUE)
int_line <- function(pred_l, obs_l, plotting = FALSE) {
  #number of points in prediction and observation
  n_pred <- nrow(pred_l); n_obs <- nrow(obs_l)
  #matrix giving where new points should be mapped to
  map <- matrix(nrow = n_pred, ncol = 2, data = NA)
  colnames(map) <- c("x", "y")

  ##flip indexing of one line if needed
  if (get_dist(obs_l[1, ], pred_l[1, ]) >
      get_dist(obs_l[1, ], pred_l[nrow(pred_l),])) {
    obs_l <- obs_l[nrow(obs_l):1, ]
  }

  ##identify where points on prediction should map onto observation
  #same number of observation as prediction points, just map 1-to-1
  if (n_pred == n_obs) {
    map <- obs_l
  } else  {
    if ((n_obs < n_pred) || (n_pred != 2)) { #typical case
      #assign first and last points to stay matching
      map[1, ] <- pred_l[1,]; map[n_pred,] <- pred_l[n_pred,]
      #proportion of sections of the line used per mapping
      step <- (n_obs - 1)/(n_pred - 1)
      used <- step #keep track of what proportion of the line has been used
      for (i in 2:(n_pred - 1)) {
        #indices of what points the dividing line for each section is between
        s <- floor(used); e  <- ceiling(used)
        #assign mappings,indexing of matrices is off by one from indices of line
        map[i, 1] <- obs_l[s + 1, 1] +
          (used - s)*(obs_l[e + 1, 1] - obs_l[s + 1, 1])
        map[i, 2] <- obs_l[s + 1, 2] +
          (used - s)*(obs_l[e + 1, 2] - obs_l[s + 1, 2])
        used <- used + step #update portion of line used
      }
    } else { #have to be careful when the two
      step <- (n_obs - 1)/3 #proportion of sections of the line used per mapping
      used <- step #set initial proportion used
      for (i in 1:n_pred) {
        #indices of what points the dividing line for each section is between
        s <- floor(used); e  <- ceiling(used)
        #assign mappings,indexing of matrices is off by one from indices of line
        map[i, 1] <- obs_l[s + 1, 1] +
          (used - s)*(obs_l[e + 1, 1] - obs_l[s + 1, 1])
        map[i, 2] <- obs_l[s + 1, 2] +
          (used - s)*(obs_l[e + 1, 2] - obs_l[s + 1, 2])
        used <- used + step #update portion of line used
      }
    }
  }

  ##optional plotting of maps
  if (plotting) {
    plot(pred_l, col = "blue",
         xlim = c(min(c(obs_l[,1], pred_l[,1])), max(c(obs_l[,1], pred_l[,1]))),
         ylim = c(min(c(obs_l[,2], pred_l[,2])), max(c(obs_l[,2], pred_l[,2]))))
    points(pred_l, type = "l", col = "blue")
    points(obs_l, col = "black")
    points(obs_l, col = "black", type = "l")
    points(map, col = "red")
  }
  return(map)
}




#' Finds all the mappings for a set of observations and predictions often over
#' multiple years
#' @title Map a set of observations and predictions
#' @param start_year first year to be mapped
#' @param end_year last year to be mapped
#' @param obs_start_year year in which observation array starts
#' @param pred_start_year year in which prediction array starts
#' @param observed array of observed values of dimension year x longitude x
#'                 latitude
#' @param predicted array of predicted values of dimension year x longitude x
#'                  latitude
#' @param reg_info a \code{reg_info} list (see documentation for \code{reg_info})
#' @param month month under consideration
#' @param level concentration level for which to build contour
#' @param dat_type_obs string of either "bootstrap" or "simple" indicating the
#'                   file type of the observation (see details)
#' @param dat_type_pred string of either "gfdl" or "simple" indicating the file
#'                   type of the prediction (see details)
#' @param plotting boolean indicatng whether maps should be plotted (defaults to
#'                 false)
#' @param obs_only indicator to run mapping only for observations
#' @param pred_only indicator to run mapping only for predictions
#' @param nX dimension in the x (defaults to value for Northern Polar
#'           stereographic grid: 304)
#' @param nY dimension in the y (defaults to value for Northern Polar
#'           stereographic grid: 448)
#' @param xmn min x value (defaults to value for Northern Polar stereographic
#'            grid: -3850)
#' @param xmx max x value (defaults to value for Northern Polar
#'            stereographic grid: 3750)
#' @param ymn min y value (defaults to value for Northern Polar
#'            stereographic grid: -5350)
#' @param ymx max y value (defaults to value for Northern Polar
#'            stereographic grid: 5850)
#' @return \code{map} object (see details)
#' @details The object \code{maps} is obtained from running the
#'          \code{create_mapping} function. It is a list of four objects. The
#'          first two items in the list, \code{start_year} and \code{end_year},
#'          give the first and last year that were mapped. The second two items,
#'          \code{obs_list} and \code{pred_list}, are lists of arrays with one
#'          3-dimensional array for each region. The first dimension is for the
#'          year. The other two dimensions are for the fixed points'
#'          y-coordinates, the mapped points' x-coordinates, the mapped points'
#'          y-coordinates, the length of the mapping vectors in the x-direction,
#'          the length of the vectors in the y-direction, and the angles of the
#'          mapping vectors.
#'
#'
#'         For \code{dat_type_obs = "simple"} and \code{dat_type_pred = "simple"}
#'         the values in the \code{observed} and \code{predicted} arrays are
#'         indicators of whether the grid box contains ice (1: ice-covered,
#'         0: no ice, NA: land). If \code{datTypePred = "gfdl"} or
#'         \code{dat_type_obs = "bootstrap"}, the values in the \code{observed}
#'         and \code{predicted} arrays correspond to the raw ice concentrations
#'         values observed or predicted (including indicators for missing data,
#'         land etc.). If \code{datTypePred = "gfdl"}, the predictions are
#'         formatted as in the CM2.5 Forecast-oriented Low-Ocean Resolution
#'         (FLOR) model produced by the National Oceanic and Atmospheric
#'         Administration’s Geophysical Fluid Dynamics Laboratory and converted
#'         to a Polar Stereographic grid (Vecchi et al. 2014; Msadek et al. 2014).
#'         If \code{datTypeObs = "bootstrap"} the array values are assumed to be
#'         from the monthly sea ice concentration  obtained from the National
#'         Aeronautics and Space Administration (NASA) satellites Nimbus-7 SMMR
#'         and DMSP SSM/I-SSMIS and processed by the bootstrap algorithm.
#'         Weights for converting to a polar stereograhic grid were obtained
#'         from the spherical coordinate remapping and interpolation package
#'         (SCRIP) (Jones 1997).
#'
#' @references Comiso, J., 2017: Bootstrap sea ice concentrations
#'             from Nimbus-7 SMMR and DMSP SSM/I-SSMIS. version 3.
#'             Boulder, Colorado USA: NASA National Snow and Ice Data Center
#'             Distributed Active Archive Center.
#'             doi: \url{https://doi.org/10.5067/7Q8HCCWS4I0R}
#'
#'            CM2.5 Forecast-oriented Low-Ocean Resolution (FLOR) model:
#'            Vecchi, Gabriel A., et al.
#'            \href{http://journals.ametsoc.org/doi/abs/10.1175/JCLI-D-14-00158.1}{"On the seasonal forecasting of regional tropical}
#'            cyclone activity." Journal of Climate 27.21 (2014): 7994-8016.
#'
#'            Msadek, R., et al.
#'           \href{http://onlinelibrary.wiley.com/doi/10.1002/2014GL060799/full}{"Importance of initial conditions in seasonal predictions of
#'            Arctic sea ice extent."} Geophysical Research Letters
#'            41.14 (2014): 5208-5215.
#'
#'           National Center for Atmospheric Research, 2017: Earth system grid
#'           at NCAR. \url{https://www.earthsystemgrid.org/home.html}.
#'
#'           Jones, P.W. "A user’s guide for SCRIP: A spherical coordinate
#'           remapping and interpolation package." Los Alamos National
#'           Laboratory, Los Alamos, NM (1997).
#'
#' @importFrom sp disaggregate
#' @importFrom rgeos gArea gIntersection
#' @importFrom methods as
#' @importFrom graphics lines
#' @export
#' @examples \dontrun{
#' create_mapping(start_year = 1981, end_year = 1981, obs_start_year = 1981,
#'                pred_start_year = 1980, observed = obsFeb19811982,
#'                predicted = emFeb19811982, reg_info = reg_info, month = 2,
#'                level = 15, dat_type_obs = "bootstrap", dat_type_pred = "gfdl",
#'                plotting = TRUE) }
create_mapping <- function(start_year, end_year, obs_start_year, pred_start_year,
                           observed, predicted, reg_info, month, level,
                           dat_type_obs, dat_type_pred, plotting = FALSE,
                           obs_only = FALSE, pred_only = FALSE,
                           nX = 304, nY = 448, xmn = -3850,
                           xmx = 3750, ymn = -5350, ymx = 5850) {
  month_lab <- c("January", "February", "March", "April", "May", "June",
                "July", "August", "September", "October", "November",
                "December")

  n_reg <- length(reg_info$regions)
  n_years <- length(start_year:end_year)
  n_lines <- sapply(reg_info$lines, length)

  #storage matrices
  pred_list <- obs_list <- list()
  for (r in 1:n_reg) {
    pred_list[[r]] <- obs_list[[r]]  <- array(dim = c(n_years, n_lines[r], 6))
  }

  #years of interest
  years <- start_year:end_year
  if (!pred_only) {
    missing <- which(apply(observed, 1, function(x){all(is.na(x))}))
    if (length(missing) > 0) {
      years <- years[-missing]
    }
  }

  ##store mappings for all years
  for (i in years) {
    #read-in and format observation and prediction
    if (!pred_only) {
      obs <- get_region(dat = observed[i - obs_start_year + 1, ,],
                        dat_type = dat_type_obs, level = level)
    }
    if (!obs_only) {
      raw <- get_region(dat = predicted[i - pred_start_year + 1, ,],
                        dat_type = dat_type_pred, level = level)
    }

    #map regions
    if (!pred_only) {
      obs_map <- get_map(ice = obs, plotting, reg_info = reg_info,
                       main = sprintf("Observed Mapping:\n %s %i",
                                      month_lab[month], i))
    }
    if (!obs_only) {
      raw_map <- get_map(ice = raw, plotting, reg_info = reg_info,
                       main = sprintf("Predicted Mapping:\n %s %i",
                                      month_lab[month], i))
    }

    #store results
    index <- i - start_year + 1
    for (r in 1:n_reg) {
      if (!pred_only) {
        obs_list[[r]][index,,] <- obs_map[[r]]
      }
      if (!obs_only) {
        pred_list[[r]][index,,] <- raw_map[[r]]
      }
    }

    print(sprintf("mapping complete for year %i", i))
  }

  if (obs_only) {
    pred_list <- NULL
  }
  if (pred_only) {
    obs_list <- NULL
  }

  return(list("start_year" = start_year, "end_year" = end_year,
              "pred_list" = pred_list, "obs_list" = obs_list))
}
