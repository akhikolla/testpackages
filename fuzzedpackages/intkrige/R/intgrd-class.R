### CREATE CLASSES ###
#=============================================================================
# This is the "intgrd" class, which extends the "SpatialPixelsDataFrame"
# class. Many of the generic methods defined here were
# adapted from code in the sp package.
# Source code for the sp package can be obtained from
#
# https://github.com/edzer/sp
#
# We add a slot for intervals that can be called by plot,
# print, and summary functions.
NULL

# Necessary package imports and inclusions
#' @import sp
#' @import gstat
#' @include generics.R
NULL
#=============================================================================

### CREATE METHODS FOR GENERICS ###
#=============================================================================
# This method ensures that the "$" operator looks through the data,
# coordinate, and interval slots when looking for a subset.
# Adapted from:
# - https://github.com/edzer/sp/blob/master/R/SpatialPoints-methods.R

#' @name $
#' @rdname extract-methods
#' @aliases $,intgrd-method
setMethod("$", "intgrd",
          function(x, name) {
            if (name %in% coordnames(x))
              return(x@coords[,name])
            if (name %in% colnames(x@interval))
              return(x@interval[, name])
            if (!("data" %in% slotNames(x)))
              stop("no $ method for object without attributes")
            x@data[[name]]
          }
)

# This method stops users from editing the coordinates or intervals
# using the standard "$" operator.
# Adapted from:
# - https://github.com/edzer/sp/blob/master/R/Spatial-methods.R

#' @name $<-
#' @rdname extract-methods
#' @aliases $<-,intgrd-method
setMethod("$<-", "intgrd",
          function(x, name, value) {
            if (name %in% coordnames(x))
              stop(paste(name,
                         "is a coordinate name, please choose another name"))
            # Addition for coordinate slot.
            if (name %in% colnames(x@interval))
              stop(paste(name, "is currently assigned to the interval slot,
                         please choose another name"))
            if (!("data" %in% slotNames(x))) {
              df = list(value); names(df) = name
              return(addAttrToGeom(x, data.frame(df), match.ID = FALSE))
              # stop("no $<- method for object without attributes")
            }
            #if (is.list(value))
            #	warning("assigning list or data.frame to attribute vector")
            x@data[[name]] = value
            x
          }
)

# Method to properly subset and intsp object. The native functions fail to
# subset the interval slot.
#' @rdname extract-methods
#' @aliases [,intgrd-method
setMethod("[", signature(x = "intgrd",
                         i = "ANY",
                         j = "missing",
                         drop = "missing"),
          function(x, i, j, ..., drop) {
            if (!("data" %in% slotNames(x)))
              stop("no [ method for object without attributes")

            # Subset the contents of each slot.
            x@coords <- x@coords[i, ]
            x@data <- x@data[i, ]
            x@interval <- x@interval[i, ]
            x
          }
)

#' @name interval
#' @rdname interval-methods
#' @aliases interval,intgrd-method
setMethod("interval", "intgrd", function(x) {
  return(x@interval)
})

#' @name interval<-
#' @rdname interval-methods-assign
#' @aliases interval<-,intgrd-method
setMethod("interval<-", "intgrd",
          function(x, value) {

            # A null or na value will cause the intgrd object to revert to
            # its parent class.
            if(is.null(value) | all(is.na(value))){
              # If the interval slot is not empty, return these values
              # to the data frame.
              if (nrow(x@interval) > 0) {
                x@data <- cbind(x@data, as.data.frame(x@interval))
              }

              x <- as(x, "SpatialPixelsDataFrame")

              return(x)
            }else if(inherits(value, "character")){
              # If a character string is provided, then search the data frame
              # for the values to place in the interval slot.
              if(length(value) != 2){
                stop("exactly two variable names are needed
                     to define an interval")
              }

              # If the interval slot is not empty, return these values
              # to the data frame.
              if (nrow(x@interval) > 0) {
                x@data <- cbind(x@data, as.data.frame(x@interval))
              }

              # Ensure that only numeric variables are placed in the interval
              if (!inherits(x@data[, value[1]], "numeric") |
                  !inherits(x@data[, value[2]], "numeric")) {
                stop("non-numeric input detected")
              }

              # Fill the interval slot with the requested columns.
              x@interval <-
                matrix(c(x@data[, value[1]], x@data[, value[2]]), ncol = 2)

              # Remove the variables that are now in the interval slot from the
              # data frame.
              x@data[, value[1]] <- NULL
              x@data[, value[2]] <- NULL

              # Preserve the column names in the interval slot.
              colnames(x@interval) <- c(value[1], value[2])

              # If a matrix is input, allow the interval values to be replaced.
            }else if(inherits(value, "matrix")){
              # Ensure that the number of rows in the replacement
              # matches the number of rows expected.
              if(nrow(value) != nrow(coordinates(x@coords)) |
                 ncol(value) != ncol(coordinates(x@coords))){
                stop("matrix dimensions must match the slot dimensions")
              }

              # Replace the values in the interval slot.
              x@interval <- value

            }else{
              stop("\'value\' must be a vector of column names or a matrix")
            }

            # Check to ensure that the lower enpoints are defined first
            if(any(x@interval[, 1] > x@interval[, 2])){
              stop("detected at least one lower endpoint
                            greater than upper endpoint")
            }

            return(x)
          })

#' @name interval<-
#' @rdname interval-methods-assign
#' @aliases interval<-,SpatialPixelsDataFrame-method
setMethod("interval<-", "SpatialPixelsDataFrame",
          function(x, value) {

            if(is.null(value) | all(is.na(value))){
              warning("No interval provided, interval not specified.")
              return(x)
            }

            x <- as(x, "intgrd")

            if(inherits(value, "character")){
              # Ensure that only numeric variables are placed in the interval
              if (!inherits(x@data[, value[1]], "numeric") |
                  !inherits(x@data[, value[2]], "numeric")) {
                stop("non-numeric input detected")
              }

              # Fill the interval slot with the requested columns.
              x@interval <-
                matrix(c(x@data[, value[1]], x@data[, value[2]]), ncol = 2)

              # Remove the variables that are now in the interval slot from the
              # data frame.
              x@data[, value[1]] <- NULL
              x@data[, value[2]] <- NULL

              # Preserve the column names in the interval slot.
              colnames(x@interval) <- c(value[1], value[2])

            }else if(inherits(value, "matrix")){
              # Ensure that the number of rows in the replacement
              # matches the number of rows expected.
              if(nrow(value) != nrow(coordinates(x@coords)) |
                 ncol(value) != ncol(coordinates(x@coords))){
                stop("matrix dimensions must match the slot dimensions")
              }

              # Replace the values in the interval slot.
              x@interval <- value

            }else{
              stop("\'value\' must be a vector of column names or a matrix")
            }

            # Check to ensure that the lower enpoints are defined first
            if(any(x@interval[, 1] > x@interval[, 2])){
              stop("detected at least one lower endpoint
                            greater than upper endpoint")
            }

            return(x)
          })

# This method redefines how to print the intgrd object to the screen.
# Adapted from the sp package
#' Print the contents of an intgrd object
#'
#' This function extends print.sp by including a display and summary of the
#' interval slot for the object.
#'
#' @param x An object of class \code{intgrd}.
#' @param ... Additional arguments to \code{\link{print}}.
#' @param digits Determines how numeric values are printed to the screen
#'   (default from \code{sp} package).
#' @return Prints object to the screen, identical to
#'   \code{\link[sp]{SpatialPoints-class}}, as well as summary statistics for the
#'   interval slot.
#' @method print intgrd
#' @export
print.intgrd = function(x, ..., digits = getOption("digits")) {
  cc = substring(paste(as.data.frame(
    t(signif(coordinates(x), digits)))),2,999)
  int = paste("[", as.data.frame(
    t(signif(interval(x), digits))), "]", sep = "")
  int <- gsub(int, pattern = "[()c]", replacement = "")

  # drop = false ensures the resul remains a data frame and
  # cannot be converted to a vector.
  df = data.frame("coordinates" = cc, "interval" = int, x@data)

  colnames(df) <- c("coordinates", "interval", colnames(x@data))
  print(df, ..., digits = digits)
}
#' Extension of the show function for intgrd objects
#' @param object and object of class intgrd
#' @method show intgrd
setMethod("show", "intgrd", function(object) print(object))

# This method redefines how to print the head the intgrd object to the screen.
# Adapted from SP
#' Print the head of an intgrd object.
#'
#' This function extends print.sp by including a display and summary of the
#' interval slot for the object.
#'
#' @param x An object of class \code{intgrd}.
#' @param n Number of rows to print to the screen.
#' @param ... Additional arguments to \code{print}.
#' @param digits Determines how numeric values are printed to the screen
#'   (default from \code{sp} package).
#' @return Prints a subset of the object observations to the screen,
#'   identical to \code{\link[sp]{SpatialPoints-class}}, as well as
#'   summary statistics for the interval slot.
#' @method head intgrd
#' @export
head.intgrd = function(x, n = 6, ..., digits = getOption("digits")){
  cc = substring(paste(as.data.frame(
    t(signif(sp::coordinates(x), digits)))),2,999)
  int = paste("[", as.data.frame(
    t(signif(intkrige::interval(x), digits))), "]", sep = "")
  int <- gsub(int, pattern = "[()c]", replacement = "")
  df = data.frame("coordinates" = cc, "interval" = int, x@data)

  colnames(df) <- c("coordinates", "interval", colnames(x@data))
  print(df[1:n, ], ..., digits = digits)
}
setMethod("head", "intgrd", function(x, ...) head.intgrd(x, ...))

# This method redefines how to print the tail of the intgrd object to the screen.
# Adapted from SP
#' Print the tail of an intgrd object.
#'
#' This function extends print.sp by including a display and summary of the
#' interval slot for the object.
#'
#' @param x An object of class \code{intgrd}.
#' @param n The number of rows to print to the screen.
#' @param ... Additional arguments to \code{tail}.
#' @param digits Determines how numbers are displayed to the screen
#'   (default taken from package \code{sp}).
#' @return Prints a subset of the object observations to the screen,
#'   identical to \code{\link[sp]{SpatialPoints-class}}, as well as
#'   summary statistics for the interval slot.
#' @method tail intgrd
#' @export
tail.intgrd = function(x, n = 6, ..., digits = getOption("digits")){
  cc = substring(paste(as.data.frame(
    t(signif(sp::coordinates(x), digits)))),2,999)
  int = paste("[", as.data.frame(
    t(signif(intkrige::interval(x), digits))), "]", sep = "")
  int <- gsub(int, pattern = "[()c]", replacement = "")
  df = data.frame("coordinates" = cc, "interval" = int, x@data)

  colnames(df) <- c("coordinates", "interval", colnames(x@data))
  print(df[(nrow(df)-(n-1)):nrow(df), ], ..., digits = digits)
}
setMethod("tail", "intgrd", function(x, ...) tail.intgrd(x, ...))


# This method adapts the summary.spatial to include a covariance matrix
# for the interval-center and radii in the output.
#' Summarize the contents of an \code{intgrd} object,
#' including special summaries for the interval slot.
#'
#' @param object An object of class \code{intgrd}.
#' @param ... Additional arguments to \code{\link{summary}}.
#' @return Prints a summary of the object observations to the screen,
#'   identical to \code{\link[sp]{SpatialPoints-class}}, as well as
#'   summary statistics for the interval slot.
#' @method summary intgrd
#'
#' @export
summary.intgrd = function(object, ...) {
  obj = list()
  obj[["class"]] = class(object)
  obj[["bbox"]] = bbox(object)
  obj[["is.projected"]] = is.projected(object)
  obj[["proj4string"]] = object@proj4string@projargs
  # Add the covariance matrix for center/radius interaction.
  obj[["vcov"]] <- stats::cov(interval(object))
  obj[["itvl"]] <- summary(interval(object))
  obj[["grid"]] = gridparameters(object)
  if ("data" %in% slotNames(object) & ncol(object@data) > 0)
    obj[["data"]] = summary(object@data)
  class(obj) = "summary.intgrd"
  obj
}
setMethod("summary", "intgrd", summary.intgrd)

# This method defines how intgrd objects are printed to the screen.
#' Print the object summary to the screen.
#'
#' @param x An object of class \code{intgrd}.
#' @param ... Additional arguments to \code{\link{print}}.
#' @return Prints a subset of the object observations to the screen,
#'   identical to \code{\link[sp]{SpatialPoints-class}}, as well as
#'   summary statistics for the interval slot.
#'
#' @export
print.summary.intgrd = function(x, ...) {
  cat(paste("Object of class ", x[["class"]], "\n", sep = ""))
  cat("Coordinates:\n")
  print(x[["bbox"]], ...)
  cat(paste("Is projected:", x[["is.projected"]], "\n"))
  #    cat(paste("proj4string : [", x[["proj4string"]], "]\n", sep=""))
  pst <- paste(strwrap(x[["proj4string"]]), collapse="\n")
  if (nchar(pst) < 40) cat(paste("proj4string : [", pst, "]\n", sep=""))
  else cat(paste("proj4string :\n[", pst, "]\n", sep=""))
  if (!is.null(x$npoints)) {
    cat("Number of points: ")
    cat(x$npoints)
    cat("\n")
  }
  if (!is.null(x$n.polygons)) {
    cat("Number of polygons: ")
    cat(x$n.polygons)
    cat("\n")
  }
  if (!is.null(x$grid)) {
    cat("Grid attributes:\n")
    print(x$grid, ...)
  }
  if (!is.null(x[["vcov"]])) {
    cat("center-radius covariance matrix:\n")
    print(x[["vcov"]], ...)
  }
  if (!is.null(x[["itvl"]])) {
    cat("enpoints:\n")
    print(x[["itvl"]], ...)
  }
  if (!is.null(x$data)) {
    cat("Data attributes:\n")
    print(x$data, ...)
  }
  invisible(x)
} # No "setMethod" as this directly calls intgrd.


#' @name intvariogram
#' @rdname intvariogram-methods
#' @aliases intvariogram,intgrd-method
setMethod("intvariogram", "intgrd",
          function(x, formulas = list(center ~ 1, radius ~ 1), ...){

            # First ensure that the center and
            # radius are included in the proper formulas
            # (When strings are converted to strings, mathematical operators
            # act as a string split. Because we simply need "center" and
            # "radius" to appear somewhere in the text, we use the "any"
            # function to return one logical)
            check1 <- any(regexpr(pattern = "center",
                                  text = formulas[[1]][[2]]) > 0)
            check2 <- any(regexpr(pattern = "radius",
                                  text = formulas[[2]][[2]]) > 0)

            if(!check1){
              stop("Formula one must contain \'center\'
                            in the dependent variable slot.")
            }
            if(!check2){
              stop("Formula two must contain \'radius\'
                            in the dependent variable slot.")
            }

            x$center <- (interval(x)[, 1] + interval(x)[, 2]) / 2
            x$radius <- (interval(x)[, 2] - interval(x)[, 1]) / 2

            g1 <- gstat::gstat(NULL, "center",
                               formula = formulas[[1]], data = x, ...)
            g2 <- gstat::gstat(g1, "radius",
                               formula = formulas[[2]], data = x, ...)

            gv <- gstat::variogram(g2, ...)

            # Have the intvariogram class inherit the
            # original variogram class
            class(gv) <- c("intvariogram", class(gv))

            return(gv)
          })

#' @name as.data.frame
#' @rdname interval.as.data.frame-methods
#' @aliases as.data.frame,intgrd-method
setMethod("as.data.frame", "intgrd", function(x){
  interval(x) <- NULL
  return(as.data.frame(x))
})


#' Create an interval plot for spatial grid.
#'
#' Calls \code{\link[sp]{spplot}} to plot the locations, centers, and
#' radii of an \code{intgrd} object in a single figure.
#'
#' @param x An object of class \code{intgrd}.
#' @param beside Tf true, center and radius plotted side by side
#'  if false, center and radius are plotted in a single figure with the center
#'  plotted using color and the radius plotted using circles circumscribed
#'  within each grid cell.
#' @param circleCol If beside=TRUE, the color of the circles
#'  that will be circumscribed within each grid cell.
#' @param minRad The minimum value of the radius in the circles
#'  drawn to represent the interval radii. Must be a number between 0 and 1
#'  where approaching 0 results in a point being drawn in the center of the grid,
#'  while approaching 1 results in every circle being circumscribed in their
#'  respective grid cell (which is not very interesting).
#' @param ... Additional arguments to \code{\link[sp]{spplot}}.
#' @method plot intgrd,missing
setMethod("plot",
          signature = c("intgrd", "missing"),
          function(x, beside = TRUE, circleCol = "black",
                   minRad = 0.25, ...){


            x$center <- (interval(x)[, 1] + interval(x)[, 2]) / 2
            x$radius <- (interval(x)[, 2] - interval(x)[, 1]) / 2

            interval(x) <- NULL

            if(beside){

              pc <- sp::spplot(x, zcol = "center",
                               main = "center", ...)
              pr <- sp::spplot(x, zcol = "radius",
                               main = "radius", ...)
              return(gridExtra::grid.arrange(pc, pr, ncol = 2))

            }else{

              if(minRad < 0 | minRad > 1){
                stop("minRad must be an argument between 0 and 1")
              }

              # Determine max radius and include in plot title
              maxrad <- max(x$radius, na.rm = TRUE)
              minrad <- min(x$radius, na.rm = TRUE)

              # Now extract the coordinates and
              # resolution of the grid.
              tempcoords <- sp::coordinates(x)
              x2 <- as(x, "RasterLayer")
              res <- raster::res(x2)

              # Determine the proportion of the distance between max and min
              radScale <- (x$radius - minrad)/(maxrad - minrad)

              # Create circular objects with radii proportional
              # to the value of the radius.
              rads <- 1:360*pi/180
              ll <- nrow(tempcoords)
              templines <- vector("list", ll)
              for(i in 1:ll){
                lon <- tempcoords[i, 1] +
                  (radScale[i] + (1-radScale[i])*minRad)*(res[1]/2)*cos(rads)
                lat <- tempcoords[i, 2] +
                  (radScale[i] + (1-radScale[i])*minRad)*(res[2]/2)*sin(rads)

                templines[[i]] <- Lines(Line(cbind(lon, lat)), ID = as.character(i))
              }

              # Place radius legend in each corner
              radLeg <- vector("list", 8)
              count = 1
              for(i in 1:4){
                for(j in 1:2){
                  if(i == 1){
                    xanchor <- x@bbox[1, 1] - res[1]/2
                    yanchor <- x@bbox[2, 1] - res[2]/2
                  }else if(i == 2){
                    xanchor <- x@bbox[1, 1] - res[1]/2
                    yanchor <- x@bbox[2, 2] + res[2]/2
                  }else if(i == 3){
                    xanchor <- x@bbox[1, 2] + res[1]/2
                    yanchor <- x@bbox[2, 2] + res[2]/2
                  }else{
                    xanchor <- x@bbox[1, 2] + res[1]/2
                    yanchor <- x@bbox[2, 1] - res[2]/2
                  }

                  if(j ==1){
                    xleg <- (res[1]/2)*cos(rads) + xanchor
                    yleg <- (res[2]/2)*sin(rads) + yanchor
                  }else{
                    xleg <- (res[1]/2)*minRad*cos(rads) + xanchor
                    yleg <- (res[2]/2)*minRad*sin(rads) + yanchor
                  }

                  radLeg[[count]] <- Lines(Line(cbind(xleg, yleg)),
                                           ID = count)
                  count <- count + 1
                }
              }


              radLeg <- SpatialLines(radLeg)
              tester <- SpatialLines(templines)

              # Extract the resolution of the grid
              pc <- sp::spplot(x, zcol = "center",
                               main = paste("Min radius: ", round(minrad, 3), "; ",
                                            "Max radius: ", round(maxrad, 3), sep = ""),
                               xlim = c(x@bbox[1, 1] - res[1],
                                        x@bbox[1, 2] + res[1]),
                               ylim = c(x@bbox[2, 1] - res[2],
                                        x@bbox[2, 2] + res[2]),
                               ...)
              p2 <- pc +
                latticeExtra::layer(sp.lines(tester, col = circleCol),
                                    data = list(tester = tester,
                                                circleCol = circleCol))+
                latticeExtra::layer(sp.lines(radLeg, col = circleCol),
                                    data = list(radLeg = radLeg,
                                                circleCol = circleCol))

              return(p2)
            }

          }
)
