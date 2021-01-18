#' Convert \code{SpatialPolygons} object to a binary grid. Grid boxes whose
#' centers are part ofthe \code{SpatialPolygons} are given value 1 and all other
#'  grid boxes are given value 0. Land values are set to NA.
#' @title Convert \code{SpatialPolygons} object to a grid
#' @param x \code{SpatialPolygons} object
#' @param my_land_mat binary matrix specifying land locations
#' @export
#' @importFrom sp SpatialPolygons Polygons
#' @importFrom raster rasterize as.matrix
conv_to_grid <- function(x, my_land_mat = land_mat) {
  x <- aggregate(x)
  n_poly <- length(x@polygons[[1]]@Polygons)

  poly <- NULL
  for (i in 1:n_poly) {
    temp <- SpatialPolygons(list(Polygons(list(Polygon(
                           x@polygons[[1]]@Polygons[[i]]@coords)),
                           ID =  sprintf("temp%i", i))))
    temp <- untwist(temp, poly_name = sprintf("temp%i", i))
    if (!is.null(poly) & !is.null(temp)) {
      poly <- spRbind(poly, temp)
    } else if (is.null(poly) & !is.null(temp)){
      poly <- temp
    }
    if (!is.null(poly)) {
      poly <- aggregate(poly) #makes into one polygon, so binary
    }
  }
  #raster has columns and rows flipped from typical orientation
  rast <- raster(nrows = 448, ncols = 304, xmn = -3850, xmx = 3750,
                 ymn = -5350, ymx = 5850)
  rast <- rasterize(poly, rast, fun = max, background = 0)
  rast <- as.matrix(rast)
  rast <- t(rast)[,448:1]#fix weird orientation of raster
  rast[which(my_land_mat == 1, arr.ind = T)] <- NA
  return(rast)
}

#' Reads in netCDF files of observations and predictions, performs bias
#' correction, and exports a new netCDF file with bias-corrected predictions
#' @title Simple evaluation of contour-shifting
#' @param obs_NCDF filepath for observed data array (see details for info about
#'                array structure)
#' @param pred_NCDF filepath for predicted data array (see details for info about
#'                 array structure)
#' @param pred_years vectors of years for which to make prediction
#' @param start_year first year to use when learning model
#' @param month month of prediction
#' @param output_file filepath for where bias-corrected netCDF file should be
#'                   stored
#' @param dat_type_obs string of either "bootstrap" or "simple" indicating the
#'                   file type of the observation (see details for info about
#'                   array structure)
#' @param level concentration level for which to build contour
#' @param n_train_years number of prior years used in training bias correction
#' @importFrom ncdf4 nc_open ncvar_get nc_close nc_create ncdim_def
#'             ncvar_def ncvar_put
#' @details The predicted data array, \code{pred_NCDF}, should be a netCDF file
#'          with a single array of dimension: years x longitude (304) x
#'          latitude (448). The variable should be named \code{ice_ind}. The
#'          values in the array should indicate whether each grid box is
#'          categorized to contain ice (1: ice-covered, 0: no ice, NA: land).
#'          The observed data array, \code{obs_NCDF}, should be a netCDF file
#'          with a single array of dimension: years x longitude (304) x
#'          latitude (448). The observed data array, \code{obs_NCDF}, can be
#'          formatted the same as \code{pred_NCDF} if
#'          \code{dat_type_obs = "simple"}. Alternatively, if
#'          \code{dat_type_obs = "bootstrap"} the array values can be ice
#'          concentration values obtained from the National Aeronautics and Space
#'          Administration (NASA) satellites Nimbus-7 SMMR and DMSP SSM/I-SSMIS
#'          and processed by the bootstrap algorithm. Data should be retained in
#'          the same format as given by bootstrap (including indicators for
#'          missing data, land etc.). The variable should be named "conc".
#'
#' @references Comiso, J., 2017: Bootstrap sea ice concentrations
#'             from Nimbus-7 SMMR and DMSP SSM/I-SSMIS. version 3.
#'             {Boulder, Colorado USA: NASA National Snow and Ice Data Center
#'             Distributed Active Archive Center}
#'
#' @return netCDF file of dimension years by longitude (304) by latitude (448)
#'         with indicators for where ice is predicted after bias correction.
#'         (1: ice-covered, 0: not ice, NA: land). Grid boxes will be
#'         categorized as ice if their centers are ice covered (within R the
#'         bias-corrected contours are not restricted to align to a grid).
#' @examples
#' \dontrun{
#' quick_run(obs_NCDF = "/obs.nc", pred_NCDF = "/pred.nc",
#'          pred_years = c(2001:2013), start_year = 1980, month = 2,
#'          output_file = "/outputFile.nc", level = 15, dat_type_obs = "simple")
#' }
#' @export
quick_run <- function(obs_NCDF, pred_NCDF, pred_years, start_year, month,
                      output_file, level, dat_type_obs = "bootstrap",
                      n_train_years = NULL) {
  #extract input dimensions
  obs <- nc_open(obs_NCDF)
  if (dat_type_obs == "bootstrap") {
    obsMat <- ncvar_get(obs, "conc")
  } else {
    obsMat <- ncvar_get(obs, "iceInd")
  }
  obs_start_year <- obs$dim$year$vals[1]
  pred <- nc_open(pred_NCDF)
  pred_mat <- ncvar_get(pred, "iceInd")
  pred_start_year <- pred$dim$year$vals[1]

  #prepare for output
  year_dim <- ncdim_def("year", "years", pred_years)
  lon_dim <- ncdim_def("lon", "longitude", 1:304)
  lat_dim <- ncdim_def("lat", "latitude", 1:448)
  n_pred_years <- length(pred_years)
  output <- array(dim = c(n_pred_years, 304, 448))

  #run mappings for all years
  discrep <- create_mapping(start_year, end_year = max(pred_years) - 1,
                           obs_start_year,pred_start_year,
                           observed = obsMat[,month,,],
                           predicted = pred_mat[,month,,],
                           reg_info, month, level, dat_type_obs,
                           dat_type_pred = "simple")

  #Bias correct predictions
  bg_water <- conv_to_grid(bg_water)
  for (k in 1:length(pred_years)) {
    adj <- contour_shift(maps = discrep,
                         predicted = pred_mat[length(pred_start_year:pred_years[k]),
                                            month,,],
                         bc_year = pred_years[k], pred_start_year,
                         reg_info, level, dat_type_pred = "simple",
                         n_train_years)
    adj <- conv_to_grid(adj)
    adj[bg_water == 1] <- 2
    output[k, ,] <- adj
    print(sprintf("Correction of year %i completed", pred_years[k]))
  }

  #define variables
  ice_ind_def <- ncvar_def("iceInd", "indicator", list(year_dim, lon_dim, lat_dim),
                       longname = "Indicator of if pixel is ice covered
                       (0: not ice, 1: ice, NA: land, 2: Outside region")

  #create netCDF file and arrays
  nc_out <- nc_create(output_file, ice_ind_def)

  #put variables
  ncvar_put(nc_out, ice_ind_def, output)

  #close file (writing file to disk)
  nc_close(nc_out)
}


