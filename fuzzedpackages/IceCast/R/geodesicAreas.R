#' Caclulate the geodesic areas of \code{SpatialPolygons} object on the Nothern
#' Hemisphere Polar Stereographic grid projection
#' @title Calculate geodesic area
#' @param poly \code{SpatialPolygons} object for which to calculate area
#' @param byid boolean indicating whether areas should be calculated for each
#'             polygon individually or for the whole object together
#' @details Area calculations are for the Polar stereographic grid with major
#'          axis of 6378273m and ellipsoid flattening of 1/298.2794111.
#' @references Information on Polar Stereographic North projection:
#'             \url{https://nsidc.org/data/polar-stereo/ps_grids.html}
#' @return Area of polygon (or vector of areas if \code{byid} is set to \code{TRUE})
#' @importFrom sp SpatialPolygons
#' @importFrom rgeos gArea
#' @importFrom methods is
#' @export
#' @examples
#' get_area(reg_info$regions[[1]])
#' get_area(land, byid = TRUE)
get_area <- function(poly, byid = FALSE) {

  #if myPoly is null, return zero area
  if (is.null(poly)) {
    return(0)
  }

  #check for appropriate data
  if (!(is(poly)[1] == "SpatialPolygons")) {
    stop("myPoly is not a polygon object")
  }

  #check byid input
  if (!is.logical(byid)) {
    stop("byid must be a logical value (TRUE or FALSE)")
  }

  if (requireNamespace("geosphere", quietly = TRUE)) {
    n_ID <- length(poly@polygons) #number of polygons ID's
    a_ID <- rep(NA, n_ID) #vector to store area of each polygons
    for (i in 1:n_ID) {
      #Get coordinates of each polygon and record whether it's a hole
      holes <- sapply(poly[i]@polygons[[1]]@Polygons,
                      function(x){x@hole})
      if (length(holes) > 1) {
        coords <- lapply(poly[i]@polygons[[1]]@Polygons,
                         function(x){x@coords})
      } else {
        coords <- list(poly[i]@polygons[[1]]@Polygons[[1]]@coords)
      }
      n_sec <- length(holes)

      ##calculate area using areaPolygon
      a_ID_curr <- 0 #area for the particular myPoly
      for (j in 1:n_sec) {
        #get geodetic latitude and longitude
        temp2 <- mapxy(coords[[j]][,1], coords[[j]][,2])
        #calc area
        a_curr <- geosphere::areaPolygon(cbind(temp2$alon, temp2$alat),
                                        6378273, f = 1/298.2794111)/(1000^2)
        if (holes[j]) { #substract off area if its a hole
          a_ID_curr <- a_ID_curr - a_curr
        } else { #otherwise add areas
          a_ID_curr <- a_ID_curr + a_curr
        }
      }
      a_ID[i] <- a_ID_curr
    }

    #return results
    if (byid) {
      return(a_ID)
    } else { #summed
      return(sum(a_ID))
    }
  } else {
    print("WARNING: Geosphere Package is not installed. Planar areas,
          not geodetic areas were calculated. For geodetic areas,
          install Geosphere package and re-run get_area")
    if (byid) {
      return(gArea(poly, byid = TRUE))
    } else {
      return(gArea(poly))
    }
  }
}


#' Get corresponding latitude and longitude values for coordinates on a Polar
#' Stereographic North projection grid
#' @title Get geodetic latitudes and longitudes
#' @param X Polar Stereographic X Coordinate (km)
#' @param Y Polar Stereographic Y Coordinate (km)
#' @param sgn indicator for Northern hemisphere (defaults to 1)
#' @param slat standard latitude (defaults to 70)
#' @param re Earth's radius (defaults to 6378.273)
#' @param e2 eccentricity squared (defaults to 0.006693883)
#' @param degrees boolean indicating whether result should be
#'                returned in degrees or radians
#' @return list with elements \code{coords$aLat}, the geodetic latitude
#'         (degrees, +90 to -90), and \code{coords$aLon}, the geodetic longitude
#'         (degrees, -180 to 180)
#' @references The equations for this calculation are from Snyder, J. P., 1982,
#'             Map Projections Used by the U.S. Geological Survey, Geological
#'             Survey Bulletin 1532, U.S. Government Printing Office.
#'             See JPL Technical Memorandum 3349-85-101 for further details.
#' @examples
#' new <- mapxy(100, 300)
#' new$alat
#' new$alon
mapxy <- function(X, Y, sgn = 1, slat = 70,  re = 6378.273, e2 = 0.006693883,
                  degrees = TRUE){
  pi = 3.141592654
  E <-  sqrt(e2)
  SL <-  slat*pi/180.
  delta <- 0

  rho <- sqrt(X^2 + Y^2)
  CM <- cos(SL)/sqrt(1 - e2*(sin(SL)^2))
  T <- tan((pi/4)-(SL/(2)))/((1 - E*sin(SL))/(1 + E*sin(SL)))^(E/2)

  if (abs(slat - 90) < 1e-5) {
    T <- rho*sqrt((1 + E)^(1 + E)*(1 - E)^(1 - E))/2/re
  } else {
    T <- rho*T/(re*CM)
  }

  chi <- (pi/2) - 2*atan(T)
  a_lat <- (chi + ((e2/2) + (5*e2^2/24) + (e2^3/12))*sin(2*chi) +
             ((7*e2^2/48)+(29*e2^3/240))*sin(4*chi) +
             (7*e2^3/120)*sin(6*chi))
  a_lat <- sgn*a_lat;
  a_long <- atan2(sgn*X, -sgn*Y)
  a_long <- sgn*a_long

  if (any(rho <= 0.1)) {
    temp <- which(rho <= 0.1)
    a_lat[temp] <- (pi/2)*sgn
    a_long[temp] <- 0
  }

  #Calculated longitude is rotated pi/4 clockwise from NSIDC polar stereographic
  #grid
  a_long <- a_long - pi/4;
  a_long[a_long < -pi] <- 2*pi + a_long[a_long < -pi]

  if (degrees == TRUE) {
    a_lat <- (180/pi)*a_lat
    a_long <- (180/pi)*a_long - delta
  }

  coords <- list("alat" = a_lat, "alon" = a_long)

  return(coords)
}

