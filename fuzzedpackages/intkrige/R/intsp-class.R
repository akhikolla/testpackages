### CREATE CLASSES ###
#=============================================================================
# This is the "intsp" class, which extends the "SpatialPointsDataFrame"
# class in the sp package. Many of the generic methods defined here were
# adapted from code in the sp package.
# Source code for the sp package can be obtained from
#
# https://github.com/edzer/sp
#
# We add a slot for intervals that can be called by plot,
# print, and summary functions.
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
#' @aliases $,intsp-method
setMethod("$", "intsp",
          function(x, name) {
            if (name %in% sp::coordnames(x))
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
#' @aliases $<-,intsp-method
setMethod("$<-", "intsp",
          function(x, name, value) {
            if (name %in% sp::coordnames(x))
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
#' @aliases [,intsp-method
setMethod("[", signature(x = "intsp",
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

# This method stops users from editing the coordinates or intervals
# using the standard "$" operator.
# Adapted from:
# - https://github.com/edzer/sp/blob/master/R/Spatial-methods.R
#' @name $<-
#' @rdname extract-methods
#' @aliases $<-,intsp-method
setMethod("$<-", "intsp",
          function(x, name, value) {
            if (name %in% sp::coordnames(x))
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

# This method returns the interval slot.
#' @name interval
#' @rdname interval-methods
#' @aliases interval,intsp-method
setMethod("interval", "intsp", function(x) {
  return(x@interval)
})

#' @name interval<-
#' @rdname interval-methods-assign
#' @aliases interval<-,intsp-method
setMethod("interval<-", "intsp",
          function(x, value) {

            # A null or na value will cause the intsp object to revert to
            # its parent class.
            if(is.null(value) | all(is.na(value))){
              # If the interval slot is not empty, return these values
              # to the data frame.
              if (nrow(x@interval) > 0) {
                x@data <- cbind(x@data, as.data.frame(x@interval))
              }

              x <- as(x, "SpatialPointsDataFrame")

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
#' @aliases interval<-,SpatialPointsDataFrame-method
setMethod("interval<-", "SpatialPointsDataFrame",
          function(x, value) {

            if(is.null(value) | all(is.na(value))){
              warning("No interval provided, interval not specified.")
              return(x)
            }

            x <- as(x, "intsp")

            if(inherits(value, "character")){
              # Ensure that only numeric variables are placed in the interval
              if (!inherits(x@data[, value[1]], "numeric") |
                  !inherits(x@data[, value[2]], "numeric") ) {
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

# This method redefines how to print the intsp object to the screen.
# Adapted from print generic in the sp package.
#' Print the contents of an \code{intsp} object
#'
#' This function extends printing methods in the \code{sp} package
#'   by including a display and summary of the interval slot for
#'   the object.
#'
#' @param x An object of class \code{intsp}.
#' @param ... Additional arguments to \code{\link{print}}.
#' @param digits Determines how numbers are displayed on the screen.
#'   Default option taken from \code{sp} package.
#' @return Prints object to the screen, identical to
#'   \code{\link[sp]{SpatialPoints-class}}, as well as summary statistics for the
#'   interval slot.
#' @method print intsp
#' @export
print.intsp = function(x, digits = getOption("digits"), ...) {
  cc = substring(paste(as.data.frame(
    t(signif(sp::coordinates(x), digits)))),2,999)
  int = paste("[", as.data.frame(
    t(signif(intkrige::interval(x), digits))), "]", sep = "")
  int <- gsub(int, pattern = "[()c]", replacement = "")

  # drop = false ensures the resul remains a data frame and
  # cannot be converted to a vector.
  df = data.frame("coordinates" = cc, "interval" = int, x@data)

  colnames(df) <- c("coordinates", "interval", colnames(x@data))
  print(df, ..., digits = digits)
}
#' Extension of the show function for intgrd objects
#' @param object An object of class \code{intsp}.
#' @method show intsp
setMethod("show", "intsp", function(object) print.intsp(object))

# This method redefines how to print the head the intsp object to the screen.
# Adapted from sp package source code.
# This method redefines how to print the heaf intsp object to the screen.
# Adapted from sp package source code.
#' Print the head of an \code{intsp} object.
#'
#' This function extends printing methods in the \code{sp} package
#'   by including a display and summary of the head of the
#'   interval slot for the object.
#'
#' @param x An object of class \code{intsp}.
#' @param n The number of rows to print to the screen.
#' @param ... Additional arguments to \code{\link{print}}.
#' @param digits Determines how values are printed to the screen (default
#'   taken from sp package).
#' @return Prints a subset of the object to the screen, identical to
#'   \code{\link[sp]{SpatialPoints-class}}, as well as summary statistics for the
#'   interval slot.
#' @method head intsp
#' @export
head.intsp <- function(x, n = 6, ..., digits = getOption("digits")) {
  cc = substring(paste(as.data.frame(
    t(signif(sp::coordinates(x), digits)))),2,999)
  int = paste("[", as.data.frame(
    t(signif(intkrige::interval(x), digits))), "]", sep = "")
  int <- gsub(int, pattern = "[()c]", replacement = "")
  df = data.frame("coordinates" = cc, "interval" = int, x@data)

  colnames(df) <- c("coordinates", "interval", colnames(x@data))
  print(df[1:n, ], digits = digits, ...)
}
setMethod("head", "intsp", function(x, ...) head.intsp(x, ...))

# This method redefines how to print the tail of the intsp object to the screen.
# Adapted from sp package source code
#' Print the tail of an intsp object.
#'
#' This function extends print.sp by including a display and summary of the
#' interval slot for the object.
#'
#' @param x An object of class \code{intsp}.
#' @param n The number of rows to print to the screen.
#' @param ... Additional arguments to \code{\link{tail}}.
#' @param digits Determines how numbers are displayed (default taken from
#'   sp package).
#' @return Prints a subset of the object to the screen, identical to
#'   \code{\link[sp]{SpatialPoints-class}}, as well as summary statistics for the
#'   interval slot.
#' @method tail intsp
#' @export
tail.intsp = function(x, n = 6, ..., digits = getOption("digits")) {
  cc = substring(paste(as.data.frame(
    t(signif(sp::coordinates(x), digits)))),2,999)
  int = paste("[", as.data.frame(
    t(signif(intkrige::interval(x), digits))), "]", sep = "")
  int <- gsub(int, pattern = "[()c]", replacement = "")
  df = data.frame("coordinates" = cc, "interval" = int, x@data)

  colnames(df) <- c("coordinates", "interval", colnames(x@data))
  print(df[(nrow(df)-(n-1)):nrow(df), ], ..., digits = digits)
}
setMethod("tail", "intsp", function(x, ...) tail.intsp(x, ...))


# This method adapts the summary.spatial to include a covariance matrix
# for the interval-center and radii in the output.
#' Summarize the contents of an \code{intsp} object,
#'   including special summaries for the interval slot.
#'
#' @param object An object of class \code{intsp}.
#' @param ... Additional arguments to \code{\link{summary}}.
#' @return Prints a series of summaries to the screen, identical to
#'   \code{\link[sp]{SpatialPoints-class}}, as well as summary statistics for the
#'   interval slot.
#' @method summary intsp
#' @export
summary.intsp = function(object, ...) {
  obj = list()
  obj[["class"]] = class(object)
  obj[["bbox"]] = sp::bbox(object)
  obj[["is.projected"]] = sp::is.projected(object)
  obj[["proj4string"]] = object@proj4string@projargs
  # Add the covariance matrix for center/radius interaction.
  obj[["vcov"]] <- stats::cov(intkrige::interval(object))
  obj[["itvl"]] <- summary(intkrige::interval(object))
  obj[["npoints"]] = nrow(object@coords)
  if ("data" %in% slotNames(object) & ncol(object@data) > 0)
    obj[["data"]] = summary(object@data)
  class(obj) = "summary.intsp"
  obj
}
setMethod("summary", "intsp", summary.intsp)

# This method defines how intsp objects are printed to the screen.
#' Print the object summary to the screen.
#'
#' @param x An object an object of class \code{intsp}.
#' @param ... Additional arguments to \code{\link{print}}.
#' @return Prints a series of summaries to the screen, identical to
#'   \code{\link[sp]{SpatialPoints-class}}, as well as summary statistics for the
#'   interval slot.
#' @export
print.summary.intsp = function(x, ...) {
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
} # No "setMethod" as this directly calls intsp.

#' @name intvariogram
#' @rdname intvariogram-methods
#' @aliases intvariogram,intsp-method
setMethod("intvariogram", "intsp",
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

            x$center <- (intkrige::interval(x)[, 1] + intkrige::interval(x)[, 2]) / 2
            x$radius <- (intkrige::interval(x)[, 2] - intkrige::interval(x)[, 1]) / 2

            g1 <- gstat::gstat(NULL, "center", formula = formulas[[1]], data = x)
            g2 <- gstat::gstat(g1, "radius", formula = formulas[[2]], data = x)

            gv <- gstat::variogram(g2, ...)

            # Have the intvariogram class inherit the original variogram class
            class(gv) <- c("intvariogram", class(gv))

            return(gv)
          })

#' @name as.data.frame
#' @rdname interval.as.data.frame-methods
#' @aliases as.data.frame,intsp-method
setMethod("as.data.frame", "intsp", function(x){
  intkrige::interval(x) <- NULL
  return(as.data.frame(x))
})

#' Create an interval plot for spatial points.
#'
#' Calls \code{\link[sp]{spplot}} to plot the locations, centers, and
#'   radii of an interval-valued spatial data frame in a single figure.
#' @param x An object of class \code{intsp}.
#' @param locationsOnly If TRUE, simply plots geographic
#'   locations.
#' @param legend.positions The positions of the center and radius legend
#'   relative to the plotting window.
#' @param cuts The number of ranges of values to print in
#'   the center and radius legend respectively.
#' @param radSize A vector of length 2 indicating the range
#'   of point sizes to plot to visualize radii magnitudes.
#' @param pch The shape of the points
#'   (see \code{\link{plot}}).
#' @param alpha The transparency of the points.
#' @param ... Additional arguments to \code{\link[sp]{spplot}}.
#' @method plot \code{intsp},missing
setMethod("plot", signature = c("intsp", "missing"),
          function(x, locationsOnly = FALSE,
                   legend.positions = c("left", "right"),
                   cuts = c(5, 5), radSize = c(0.1, 3),
                   pch = 16, alpha = 0.5,
                   ...){

            if(locationsOnly){
              intkrige::interval(x) <- NULL
              test <- plot(x, ...)
              return(test)
            }

            if(length(legend.positions) != 2 |
               !inherits(legend.positions, "character")){
              stop("two character legend positions must be provided")
            }
            if(length(cuts) != 2 | !inherits(cuts, "numeric")){
              stop("two numeric cut arguments must be provided")
            }
            if(length(radSize) != 2 | !inherits(radSize, "numeric")){
              stop("two numeric size arguments must be provided")
            }
            if(radSize[2] < radSize[1]){
              stop("backwards range detected for radSize")
            }

            # First, compute the center and radius for the interval slot.
            x$center <- (x@interval[, 1] + x@interval[, 2]) / 2
            x$radius <- (x@interval[, 2] - x@interval[, 1]) / 2

            # remove the interval slot and convert back
            # to spatialpointsdataframe
            intkrige::interval(x) <- NULL

            # Determine the range of radius values.
            radRange <- c(min(x$radius), max(x$radius))
            radCuts <- cut(x$radius, cuts[2])

            # Map the range of radius values to the range of cex values provided
            # by the user.
            x$weights <- (x$radius - radRange[1]) / (radRange[2] - radRange[1])
            x$cex <- x$weights*radSize[1] + (1-x$weights)*radSize[2]

            key.rad <- list(title="radius",
                            points=list(pch=pch,
                                        cex = seq(radSize[1],
                                                  radSize[2],
                                                  length = cuts[2]),
                                        col="black"),
                            text=list(levels(radCuts)),
                            cex.title=1.25)

            if(legend.positions[2] == "right"){
              legend = list(right = list(fun = lattice::draw.key(key.rad)))
            }else if(legend.positions[2] == "left"){
              legend = list(left = list(fun = lattice::draw.key(key.rad)))
            }else if(legend.positions[2] == "bottom"){
              legend = list(bottom = list(fun = lattice::draw.key(key.rad)))
            }else if(legend.positions[2] == "top"){
              legend = list(top = list(fun = lattice::draw.key(key.rad)))
            }else{
              stop("invalid radius legend position supplied")
            }

            test <- sp::spplot(x, zcol = "pointDL",
                               cex = x$cex,
                               alpha = alpha, pch = pch,
                               legend = legend,
                               key.space = legend.positions[1],
                               auto.key = list(title = "center"),
                               cuts = cuts[1], ...)

            # Change the size of the points in the center
            # legend so they don't overplot
            # http://r-sig-geo.2731867.n2.nabble.com/spplot-size-of-plotting-symbol-in-legend-td2762610.html
            test$legend[[legend.positions[1]]]$args$key$points$cex <-
              rep(1.5, length(test$legend[[legend.positions[1]]]$args$key$points$cex))


            return(test)
          })
#=============================================================================
