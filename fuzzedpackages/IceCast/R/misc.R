#Load data files that are needed globally (avoids CRAN error about global
#binding)
globalVariables(c("all_regions", "bg_water", "land", "land_mat", "reg_info"))

#' @useDynLib IceCast
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
NULL

#' Finds the euclidean distance between two points (ignoring projection)
#' @title Find euclidean distance
#' @param p1 vector giving the x and y coordinate pair for the first point
#' @param p2 vector giving the x and y coordinate pair for the second point
#' @examples get_dist(c(1, 2), c(3, 4))
#' @return distance value
get_dist <- function(p1, p2) {
  sqrt(((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2))
}

#' Remove holes from a \code{SpatialPolygons} object.  Note that this function
#' differs from the function \code{findHoles} in that
#' it only removes holes contained within the polygon itself, not gaps between
#' the polygon and region boundaries
#' @title Remove holes in a polygon
#' @param my_poly \code{SpatialPolygon} object
#' @param poly_name character string to name polygon (defaults to "notSpecified")
#' @return \code{SpatialPolygon} object with holes removed
#' @export
#' @importFrom sp Polygon Polygons SpatialPolygons
#' @importFrom raster aggregate
#' @importFrom maptools spRbind
#' @importFrom methods is
#' @examples
#' with_holes <- bg_water[2]
#' plot(with_holes, col = "blue", main = "Polygon with Holes")
#' no_holes <- rm_holes(with_holes)
#' plot(no_holes, col = "blue", main = "Holes removed")
rm_holes<- function(my_poly, poly_name = "notSpecified") {
  ##check my_poly input
  if (!(is(my_poly)[1] == "SpatialPolygons")) {
    stop("myPoly is not a polygon object")
  }

  ##remove holes
  my_poly <- suppressWarnings(aggregate(my_poly)) #give whole polygon one ID
  #identify holes
  holes <- sapply(my_poly@polygons[[1]]@Polygons, function(x){x@hole})
  for (i in which(!holes)) {
    #current part of the polygon
    coords <- my_poly@polygons[[1]]@Polygons[i][[1]]@coords
    temp <- SpatialPolygons(list(Polygons(list(Polygon(coords, hole = FALSE)),
                                          sprintf("temp%i", i))))
    if (i == 1) { #make final version of the polygon if first iteration
      new_poly <-  temp
    } else { #add to existing final version of the polygon
      new_poly <- spRbind(new_poly, temp)
    }
  }
  new_poly <- aggregate(new_poly)
  new_poly@polygons[[1]]@ID <- poly_name

  return(new_poly)
}

#' Keep only \code{SpatialPolygons} from a spatial object.
#' @title Keep only spatial polygons
#' @param my_poly \code{SpatialCollections}, \code{SpatialPolygons},
#'                \code{SpatialPoints}, or \code{SpatialLines} object
#' @return \code{SpatialPolygons} object
#' @importFrom methods is
#' @export
#' @examples
#' par(mfrow = c(1, 2))
#' plot(spatialCollEx, col = "blue", main = "Spatial Collections Object")
#' poly_only <- keep_poly(spatialCollEx)
#' plot(poly_only, col = "blue", main = "Spatial Polygon Only")
keep_poly <- function(my_poly) {
  if (is.null(my_poly)) {
    return(NULL)
  } else if (is(my_poly)[1] == "SpatialCollections") {
    my_poly <- my_poly@polyobj
  } else if (is(my_poly)[1] == "SpatialPoints") {
    return(NULL)
  } else if (is(my_poly)[1] == "SpatialLines") {
    return(NULL)
  } else if (is(my_poly)[1] == "SpatialPolygons") {
    return(my_poly)
  }else {
    stop("Incorrect object type")
  }
  return(my_poly)
}

#' Keep only \code{SpatialLines} from a spatial object.
#' @title Keep only spatial lines
#' @param my_poly \code{SpatialCollections}, \code{SpatialPolygons},
#'                \code{SpatialPoints}, or \code{SpatialLines} object
#' @return \code{SpatialPolygons} object
#' @importFrom methods is
#' @export
#' @examples
#' par(mfrow = c(1, 2))
#' plot(spatialCollEx, col = "blue", main = "Spatial Collections Object")
#' line_only <- keep_line(spatialCollEx)
#' plot(line_only, col = "blue", main = "Spatial Line Only")
keep_line <- function(my_poly) {
  if (is.null(my_poly)) {
    return(NULL)
  } else if (is(my_poly)[1] == "SpatialCollections") {
    my_poly <- my_poly@lineobj
  } else if (is(my_poly)[1] == "SpatialPoints") {
    return(NULL)
  } else if (is(my_poly)[1] == "SpatialLines") {
    return(my_poly)
  } else {
    stop("Incorrect object type")
  }
  return(my_poly)
}

#' Get coordinates from a spatial object of lines and points. There is no
#' ordering of points returned. Note: This differs from \code{extract_coords} in
#' that the ordering of the points is NOT considered.
#' @title Extract coordinates from a spatial object of lines and points
#' @param my_points spatial object of type \code{SpatialCollections},
#'                 \code{SpatialPoints}, or \code{SpatialLines}
#' @return n x 2 matrix of coordinates
#' @importFrom rgeos gLineMerge
#' @importFrom sp disaggregate
#' @importFrom methods as is
#' @examples
#' #Load sample line
#' ex_line <- as(rm_holes(bg_water[2]), "SpatialLines")
#' get_coords(ex_line)
get_coords <- function(my_points) {
  ##split apart lines and points
  line <- points <- NULL
  first <- TRUE
  if (is(my_points)[1] == "SpatialCollections") {
    line <- my_points@lineobj
    if (!is.null(line)) {
      line <- disaggregate(gLineMerge(my_points@lineobj))
    }
    points <- my_points@pointobj
  } else if (is(my_points)[1] == "SpatialPoints") {
    points <- my_points
  } else if (is(my_points)[1] == "SpatialLines") {
    line <- disaggregate(my_points)
  }

  ##add coordinates from line objects
  if (!is.null(line)) {
    n1 <- length(line@lines)
    for (i in 1:n1) {
      temp <- line@lines[[i]]
      n2 <- length(temp)
      for (j in 1:n2) {
        if (first) {
          coords <- temp@Lines[[j]]@coords
          first <- FALSE
        } else {
          coords <- rbind(coords, temp@Lines[[j]]@coords)
        }
      }
    }
  }

  ##add individual points
  if (!is.null(points)) {
    if (first) {
      coords <- points@coords
      first <- FALSE
    } else {
      coords <- rbind(coords, points@coords)
    }
  }

  return(coords)
}

#' Function to find to which matrix indices coordinates correspond (on a 304 x
#'  448 grid)
#' @title Find indices in matrix
#' @param coords coordinates of interest
#' @param xmn min x (defaults to value for Northern Polar stereographic grid:
#'           -3850)
#' @param ymn min y (defaults to value for Northern Polar stereographic grid:
#'            -5350)
#' @return n x 2 matrix of coordinates on a 304 x 448 grid
#' @examples
#' dat <- matrix(nrow = 2, ncol = 2, data = c(-2000, 0, 300, 1000))
#' get_ind(dat)
get_ind <- function(coords, xmn = -3850, ymn = -5350) {
  cbind(round((coords[,1] - xmn)/25 + 1, 0),
        round((coords[,2] - ymn)/25 + 1, 0))
}

#' Function to extract coordinates from a \code{SpatialLines} object.
#' If there are breaks in the line, this function connects the closest points to
#' create one line. Note: This differs from the function \code{getCoords} in
#' that the ordering of the points is considered.
#' @title Function to extract coordinates.
#' @param x \code{SpatialLines} or \code{SpatialPolygons} object
#' @return n x 2 matrix of coordinates
#' @importFrom methods as is slot
#' @examples
#' coords <- extract_coords(reg_info$regions[[3]])
#' par(mfrow = c(1, 2))
#' plot(reg_info$regions[[3]], main = "Polygon Object")
#' plot(coords, type = "p", main = "Coordinates", pch = 20)
extract_coords <- function(x){
  ##convert spatial polygon to spatial line object if needed
  if (is(x)[1] == "SpatialPolygons") {
    x <- as(x, "SpatialLines")
  }

  ##extract coordinates into a list of lines where each line is composed of a
  #list of segments
  res <- lapply(slot(x, "lines"), function(x)
                        lapply(slot(x, "Lines"), function(y) slot(y, "coords")))
  ##go through each line segment and extract the coordinates in order
  n_lines <- length(res)
  lines <- list()
  for (i in 1:n_lines) {
    n_segs <- length(res[[i]])
    used <- rep(FALSE, n_segs)
    #list of first coordinates in each segment
    first <- lapply(res[[i]], function(x){x[1, ]})
    #list of last coordinates in each segment
    last <- lapply(res[[i]], function(x){x[nrow(x),]})
    coords <- res[[i]][[1]] #start coordinate matrix with first line
    used[1] <- TRUE
    while (any(!used)) {
      #find segments that start or end as close as possible to the beginning or
      #end of the existing coordinate matrix
      match_first_dist <- unlist(lapply(last,
                                      function(x){get_dist(x, coords[1,])}))
      match_last_dist <-unlist(lapply(first,
                              function(x){get_dist(x, coords[nrow(coords),])}))
      match_first_dist[used] <- Inf; match_last_dist[used] <- Inf
      test <- which.min(c(min(match_first_dist), min(match_last_dist)))
      if (test == 1) { #add segment to the start of the coordinate matrix
        add_top <- which.min(match_first_dist)
        temp <- res[[i]][[add_top]]
        coords <- rbind(temp[1:(nrow(temp) - 1), ], coords)
        used[add_top] <- TRUE
      } else { #assign coordinate to the end of the coordinate matrix
        add_bottom <- which.min(match_last_dist)
        temp <- res[[i]][[add_bottom]]
        coords <- rbind(coords, temp[2:nrow(temp), ])
        used[add_bottom] <- TRUE
      }
    }
    lines[[i]] <- coords
  }

  #if there is more than one line, put lines together into one line
  if (n_lines > 1) {
    used <- rep(FALSE, n_lines)
    #list of first coordinates in each segment
    first <- lapply(lines, function(x){x[1, ]})
    #list of last coordinates in each segment
    last <- lapply(lines, function(x){x[nrow(x),]})
    coords <- lines[[1]] #start coordinate matrix with first line
    used[1] <- TRUE
    while (any(!used)) {
      #find lines that start or end as close as possible to the beginning or
      #end of the existing coordinate matrix
      match_first_dist <- unlist(lapply(last,
                                      function(x){get_dist(x, coords[1,])}))
      match_last_dist <- unlist(lapply(first,
                               function(x){get_dist(x, coords[nrow(coords),])}))
      match_first_dist[used] <- Inf; match_last_dist[used] <- Inf
      test <- which.min(c(min(match_first_dist), min(match_last_dist)))
      if (test == 1) { #add line to the start of the coordinate matrix
        add_top <- which.min(match_first_dist)
        temp <- lines[[add_top]]
        coords <- rbind(temp[1:(nrow(temp) - 1), ], coords)
        used[add_top] <- TRUE
      } else { #add line to the bottom of the coordinate matrix
        add_bottom <- which.min(match_last_dist)
        temp <- lines[[add_bottom]]
        coords <- rbind(coords, temp[2:nrow(temp), ])
        used[add_bottom] <- TRUE
      }
    }
  }

  return(coords)
}

#' Takes in a matrix and returns a \code{SpatialPolygon} object
#' representing regions fitting some criteria. Typically these regions are
#' either where the sea ice concentration is above a certain level or
#' where there is land.
#' @title Get polygons corresponding to regions
#' @param dat matrix of one of the allowed data types ("gfdl", "bootstrap", or
#'            "simple)  (see details)
#' @param dat_type string indicating the format of the data: either "gfdl",
#'                "bootstrap", or "simple" (see details)
#' @param level concentration level of interest
#' @param my_land_mat binary matrix specifying land locations
#' @param my_all_regions \code{SpatialPolygons} object specifying region that will
#'                      be considered
#' @param use_all boolean, if true indicates to use the full area (overrides
#'                \code{land_mat})
#' @param land_ind boolean, if true indicates that the region of interest is the
#'                land
#' @param xmn min x dimension (defaults to value for polar stereographic grid: -3850)
#' @param xmx max x dimension (defaults to value for polar stereographic grid: 3750)
#' @param ymn min y dimension (defaults to value for polar stereographic grid: -5350)
#' @param ymx max y dimension (defaults to value for polar stereographic grix: 5850)
#' @details For \code{datType = "simple"}  the values in the \code{dat} matrix
#'          are indicators of whether the grid box contains ice (1: ice-covered,
#'          0: no ice, NA: land). If \code{datType = "gfdl"} or
#'          \code{datType  = "bootstrap"}, the values in the matrix correspond
#'          to the raw ice concentrations values observed or predicted
#'          (including indicators for missing data, land etc.). If
#'          \code{datType = "gfdl"}, the predictions are formatted as in the
#'          CM2.5 Forecast-oriented Low-Ocean Resolution (FLOR) model produced
#'          by the National Oceanic and Atmospheric Administration’s Geophysical
#'          Fluid Dynamics Laboratory converted to a Polar Stereographic grid
#'          (Vecchi et al. 2014; Msadek et al. 2014). If
#'            \code{datType = "bootstrap"} the array values are formatted the same
#'          as the ice concentration values obtained from the  National
#'          Aeronautics and Space Administration (NASA) satellites Nimbus-7
#'          SMMR and DMSP SSM/I-SSMIS and processed by the bootstrap algorithm.
#' @references Bootstrap sea ice concentration:
#'             Comiso, J., 2017: Bootstrap sea ice concentrations
#'             from Nimbus-7 SMMR and DMSP SSM/I-SSMIS. version 3.
#'             {Boulder, Colorado USA: NASA National Snow and Ice Data Center
#'             Distributed Active Archive Center}
#'
#'            CM2.5 Forecast-oriented Low-Ocean Resolution (FLOR) model:Vecchi,
#'            Gabriel A., et al.
#'            \href{http://journals.ametsoc.org/doi/abs/10.1175/JCLI-D-14-00158.1}{"On the seasonal forecasting of regional tropical cyclone activity."}
#'            Journal of Climate 27.21 (2014): 7994-8016.
#'
#'            Msadek, R., et al.
#'            \href{http://onlinelibrary.wiley.com/doi/10.1002/2014GL060799/full}{"Importance of initial conditions in seasonal predictions of Arctic sea ice extent."}
#'            Geophysical Research Letters
#'            41.14 (2014): 5208-5215.
#' @return region of interest as a \code{SpatialPolygons} object
#' @importFrom raster raster rasterToPolygons aggregate
#' @importFrom rgeos gIntersection
#' @importFrom sp disaggregate
#' @importFrom maptools unionSpatialPolygons
#' @importFrom methods is
#' @export
#' @examples
#' \dontrun{
#' obs_example <-  get_region(dat = obsFeb2012, dat_type = "bootstrap", level = 15)
#' plot(land, col = 'grey', border = FALSE)
#' plot(obs_example, col = "lightblue", add = TRUE)
#' }
get_region <- function(dat,  dat_type, level = NULL, my_land_mat = land_mat,
                       my_all_regions = all_regions, use_all = FALSE,
                       land_ind = FALSE, xmn = -3850, xmx = 3750,
                       ymn = -5350, ymx = 5850) {
  #check input
  if (!(dat_type == "bootstrap" || dat_type == "gfdl" || dat_type == "simple")) {
    stop("datType must be 'bootstrap', 'gfdl', or 'simple'")
  }

  if (dat_type != "simple") {
    ##temporary matrix to store region before converting to a polygon
    temp <- matrix(nrow = nrow(dat), ncol = ncol(dat), data = NA)
  } else {
    temp <- dat
  }

  ##find region (if not land)
  if (dat_type == "gfdl" & land_ind == FALSE) { #ice region, GFDL data
    temp[dat >= level/100] <- 1 #above ice concentration threshold
    temp[dat < level/100] <- 0 #below ice concentration threshold
    temp[is.na(dat)] <- 0 #land is not part of contiguous region
    #ice region, bootstrap data
  } else if (dat_type == "bootstrap" & land_ind == FALSE) {
     #convert specified land in land_mat to land indicator for gfdl (120)
    dat[my_land_mat == 1] <- 12 #land is not ice
    temp[dat >= level & dat <= 100] <- 1 #above ice concentration threshold
    temp[dat == 110] <- 1 #assume satelite hole is ice
    temp[dat == 120] <- 0 #land is not part of the main contiguous area
    temp[dat < level] <- 0 #below ice concentration threshold
  ##Find region (if land)
  } else if (dat_type == "gfdl" & land_ind == TRUE) {
    temp[dat < 1] <- 0 #set ice to zero
    temp[is.na(dat)] <- 1 #set land to zero
  } else if (dat_type == "bootstrap" & land_ind == TRUE) {
    dat[my_land_mat == 1] <- 120 #convert to land boundaries
    temp[dat <= 100] <- 0 #ice and ocean are not land
    temp[dat == 110] <- 0 #assume satelite hole is not on land
    temp[dat == 120] <- 1 #land
  }

  ##convert to polygon object
  temp2 <- t(temp[,ncol(temp):1]) #orientation used by raster
  poly <- raster(temp2, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx) #make raster
  #make set of grid boxes corresponding to polygon
  poly <- rasterToPolygons(poly, fun = function(x){x > 0})
  if (land_ind == TRUE) { #polygon ID
    ID <- rep("land", length(poly))
  } else {
    ID <- rep("ice", length(poly))
  }
  #aggregate into one polygon with boundaries
  polyUnion <- unionSpatialPolygons(poly, ID)
  #seperate into each polygon with boundaries
  indiv <- disaggregate(polyUnion, hole = FALSE)

  ##finalize object
  if (!land_ind & !use_all) {
    indiv <- gIntersection(indiv, my_all_regions)
  }
  if (is(indiv)[1] == "SpatialCollections") {
    indiv <- indiv@polyobj
  }

  return(indiv)
}

#' Determine initialization month based on month being forecast and lag.
#' Considers lags up to 11 months in advance.
#' @title Get initialization month
#' @param month forecast month (integer from 1 to 12 corresponing to month of year)
#' @param lag months in advance prediction is being made (integer from 1 to 11).
#' @return integer corresponding to the initialization month
#' @details Note that this calculation assumes that the prediction for a month
#'          is its on first day. This differs from the labeling used in Director
#'          et al. (2017) which rounds up to the nearest full month.
#' @export
#' @examples
#' init_month <- get_init_month(month = 10, lag = 4)
#' init_month
get_init_month <- function(month, lag) {
  ##check "month" input
  if (!(month%%1== 0) || month < 0 || month > 12) {
    stop("month is not an integer from 0 to 12")
  }
  ##check "lag" input
  if (!(lag%%1  == 0) || lag < 0 || lag > 11) {
    stop("lag is not an integer from 0 to 11")
  }

  ##calculate initialization month
  init_month <- month - lag
  if (init_month < 1) {
    init_month <- 12 - (-init_month)
  }

  return(init_month)
}

