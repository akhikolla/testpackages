# R package viscomplexr - phase portraits of functions in the
# complex number plane
# Copyright (C) 2020 Peter Biber
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>



# ------------------------------------------------------------------------------
# Package viscomplexr
# Main R file
# ------------------------------------------------------------------------------

#' @importFrom grDevices  as.raster hsv
#' @importFrom graphics   par plot rasterImage abline polypath text points
#' @importFrom stats      runif
#' @importFrom parallel   detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach    foreach
#' @importFrom foreach    getDoParWorkers
#' @importFrom foreach    registerDoSEQ
#' @importFrom foreach    %dopar%
#' @importFrom Rdpack     reprompt
#' @importFrom grDevices  col2rgb rgb

# in order to avoid package build warning for the i iterator
# in the foreach loops
utils::globalVariables("i")

# For including Rcpp
#' @useDynLib viscomplexr, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

# -------------------------------------------------------------------------------
# Pointer emulation after
# https://www.stat.berkeley.edu/~paciorek/computingTips/
#   Pointers_passing_reference_.html

newPointer <- function(inputValue) {
    object        <- new.env(parent = emptyenv())
    object$value  <- inputValue
    class(object) <- "pointer"
    return(object)
} # newPointer

# -------------------------------------------------------------------------------
# Function phaseColhsv

# Calculates a hsv color array based on an array of complex numbers, both arrays
# handed over as pointers.

# This function only takes into account the arguments of the complex numbers.
# It is called from phasePortrait when pType = "p".

# The checks for NaN values (obtained at definition gaps other than Inf,
# i.e. result NaN) are important in order to ensure proper execution of hsv.
# The color for NaN values is hsvNaN.
# Default is a neutral grey (h = 0, s = 0, v = 0.5).

phaseColhsv <- function(pCompArr, pHsvCol, stdSaturation = 0.8,
                        hsvNaN = c(0, 0, 0.5)) {

  names(hsvNaN) <- c("h", "s", "v")
  dims    <- dim(pCompArr$value)
  h       <- Arg(pCompArr$value)
  h       <- ifelse(is.nan(h), hsvNaN["h"], ifelse(h < 0, h + 2*pi, h) / (2*pi))
  v       <- ifelse(is.nan(pCompArr$value), hsvNaN["v"], 1)
  s       <- ifelse(is.nan(pCompArr$value), hsvNaN["s"], stdSaturation)

  pHsvCol$value <- array(hsv(h = h, v = v, s = s), dims)

  return(pHsvCol)

} # phaseColhsv

# ------------------------------------------------------------------------------
# Function phaseAngColhsv

# Works like phaseColhsv above, but additionally provides phase contour
# lines and shading. Called from phasePortrait when pType = "pa".

phaseAngColhsv <- function(pCompArr, pHsvCol, lambda = 7, pi2Div = 24,
                           argOffset = 0, stdSaturation = 0.8,
                           darkestShade = 0.1,
                           hsvNaN = c(0, 0, 0.5)) {

  names(hsvNaN) <- c("h", "s", "v")
  dims    <- dim(pCompArr$value)
  argmt   <- Arg(pCompArr$value)
  h       <- ifelse(is.nan(argmt), hsvNaN["h"],
                    ifelse(argmt < 0, argmt + 2*pi, argmt)/ (2*pi))
  v       <- ifelse(is.nan(pCompArr$value), hsvNaN["v"],
                    darkestShade + (1 - darkestShade) *
                      (((argmt - argOffset) / (2*pi / pi2Div)) %% 1)^(1/lambda))
  s       <- ifelse(is.nan(pCompArr$value), hsvNaN["s"], stdSaturation)

  pHsvCol$value <- array(hsv(h = h, v = v, s = s), dims)

  return(pHsvCol)

} # phaseAngColhsv

# ------------------------------------------------------------------------------
# Function phaseModColhsv

# Works like phaseColhsv above, but additionally provides contour
# lines and shading for the modulus. Called from phasePortrait
# when pType = "pm".

# The checks for NaN values (obtained at definition gaps other than Inf,
# i.e. result NaN) are important in order to ensure proper execution of hsv.
# The color for Nan values is hsvNaN.
# Default is a neutral grey (h = 0, s = 0, v = 0.5).
# Complex numbers with infinite modulus have a valid argument, so infinite
# modulus gets no shading at all (v = 1).

phaseModColhsv <- function(pCompArr, pHsvCol, lambda = 7, logBase = 2,
                           stdSaturation = 0.8, darkestShade = 0.1,
                           hsvNaN = c(0, 0, 0.5)) {

  names(hsvNaN) <- c("h", "s", "v")
  dims    <- dim(pCompArr$value)
  h       <- Arg(pCompArr$value)
  h       <- ifelse(is.nan(h), hsvNaN["h"], ifelse(h < 0, h + 2*pi, h)/ (2*pi))
  v       <- Mod(pCompArr$value)
  v       <- ifelse(is.nan(pCompArr$value), hsvNaN["v"],
                    ifelse(is.infinite(v), 1,
                           ifelse(v == 0, darkestShade,
                                  darkestShade + (1 - darkestShade) *
                                    (log(v, logBase) %% 1)^(1/lambda))))
  s       <- ifelse(is.nan(pCompArr$value), hsvNaN["s"], stdSaturation)

  pHsvCol$value <- array(hsv(h = h, v = v, s = s), dims)

  return(pHsvCol)

} # phaseModColhsv

# ------------------------------------------------------------------------------
# Function phaseModAngColhsv

# Works like phaseColhsv above, but additionally provides contour lines and
# shading for the modulus _and_ for the argument. Called from phasePortrait
# when pType = "pma".

# The checks for NaN values (obtained at definition gaps other than Inf,
# i.e. result NaN) are important in order to ensure proper execution of hsv.
# The color for NaN values is hsvNaN.
# Default is a neutral grey (h = 0, s = 0, v = 0.5).
# Complex numbers with infinite modulus have a valid argument, so infinite
# modulus gets no shading (v = 1) for modulus. But shading for the argument
# can still occur.

phaseModAngColhsv <- function(pCompArr, pHsvCol, lambda = 7, gamma = 9/10,
                              logBase = 2, pi2Div = 24, argOffset = 0,
                              stdSaturation = 0.8, darkestShade = 0.1,
                              hsvNaN = c(0, 0, 0.5)) {

  names(hsvNaN) <- c("h", "s", "v")
  dims    <- dim(pCompArr$value)
  argmt   <- Arg(pCompArr$value)
  h       <- ifelse(is.nan(argmt), hsvNaN["h"],
                    ifelse(argmt < 0, argmt + 2*pi, argmt)/(2 * pi))

  vMod    <- Mod(pCompArr$value)
  vMod    <- ifelse(is.nan(pCompArr$value), hsvNaN["v"],
                    ifelse(is.infinite(vMod), 1,
                           ifelse(vMod == 0, 0,
                                  (log(vMod, logBase) %% 1)^(1/lambda))))

  vAng    <- ifelse(is.nan(pCompArr$value), hsvNaN["v"],
                    (((argmt - argOffset) / (2 * pi / pi2Div)) %% 1)^(1/lambda))

  v       <- ifelse(is.nan(pCompArr$value), hsvNaN["v"],
                    darkestShade + (1 - darkestShade) *
                      (gamma * vMod * vAng  +
                         (1 - gamma) * (1 - (1 - vMod) * (1 - vAng))))

  s       <- ifelse(is.nan(pCompArr$value), hsvNaN["s"], stdSaturation)

  pHsvCol$value <- array(hsv(h = h, v = v, s = s), dims)

  return(pHsvCol)

} # phaseModAngColhsv

# ------------------------------------------------------------------------------
# Function buildArray

# Build an array of complex numbers which represents the complex plane.
# Each array cell represents a pixel and contains the complex number z
# associated with this pixel. These are the z-values to be transformed
# later by the function of interest.

# Computing times for 10x10 inch x 150x150 pixel = 2250000 have been
# ok so far. Thus, we compose the array as a set of vertically arranged
# blocks of about that size. Each block is saved as a temporary file in the
# tempDir folder. The function gives back all information required to locate,
# identify and process the temporary files in the correct order.

# The blocks are built and saved in a parallel loop.

buildArray <- function(widthPx, heightPx, xlim, ylim,
                       blockSizePx = 2250000,
                       tempDir,
                       verbose) {

  # How many blocks to build?
  linesPerBlock  <- blockSizePx %/% widthPx
  nBlocks        <- heightPx    %/% linesPerBlock
  linesRest      <- heightPx    %%  linesPerBlock
  nBlocks        <- nBlocks + 1 * (linesRest > 0)
  if(verbose) cat("\n.making", nBlocks, "blocks ... ")

  # First and last line number covered by each block
  iBlocks        <- c(1:nBlocks)
  lower          <- 1 + (iBlocks - 1) * linesPerBlock
  upper          <- iBlocks * linesPerBlock
  if(linesRest != 0) upper[nBlocks] <- lower[nBlocks] - 1 + linesRest
  uplow          <- cbind(lower, upper)

  # Random Code for all files to be saved during this run
  rndCode        <- substr(format(round(runif(1), 10), nsmall = 10), 3, 12)

  # Define z-Values of the Pixels
  xPxValVec <- (xlim[2] - xlim[1])/(widthPx  - 1) *
                  c(0:(widthPx  - 1)) + xlim[1]
  yPxValVec <- (ylim[2] - ylim[1])/(heightPx - 1) *
                  c((heightPx - 1):0) + ylim[1]

  # Find the xlim and ylim values for each Block, combine all meta
  # information about the blocks into one data.frame (metaZ).
  fileNames <- paste(formatC(lower, width = trunc(log10(lower[nBlocks])) + 1,
                       flag = "0"), "zmat", rndCode, ".RData", sep = "")

  # Define data.frame with meta information
  metaZ <- cbind(data.frame(fileNames = fileNames),
                 uplow,
                 xlim1 = xPxValVec[1],                 ylim1 = yPxValVec[upper],
                 xlim2 = xPxValVec[length(xPxValVec)], ylim2 = yPxValVec[lower])

  # Check for temporary directory, if it is not the tempdir() of the current
  # R session, create it if it does not exist
  if(tempDir != tempdir())
    if(!dir.exists(tempDir)) { dir.create(tempDir, recursive = TRUE) }

  # Parallel loop
  if(verbose) cat("parallel loop starting ... ")
  foreach(i = iBlocks) %dopar% {
    yPxValVecBlock <- rep(yPxValVec[c(metaZ[i,"lower"]:metaZ[i,"upper"])],
                          widthPx)
    xPxValVecBlock <- rep(xPxValVec,
                          each = metaZ[i,"upper"] - metaZ[i,"lower"] + 1)
    compVec <- complex(real = xPxValVecBlock, imaginary = yPxValVecBlock)
    compArr <- array(compVec, dim = c(metaZ[i,"upper"] - metaZ[i,"lower"] + 1,
                                      widthPx))
    save(compArr, file = paste(tempDir, metaZ[i,"fileNames"], sep = "/"))
    rm(list = c("compArr", "compVec", "xPxValVecBlock", "yPxValVecBlock"))
  } # foreach i
  if(verbose) cat("done.")

  # Arrange meta-information to be returned
  metaBack <- list(tempDir = tempDir, rndCode = rndCode, metaZ = metaZ)

  return(metaBack)

} # function buildArray

# -----------------------------------------------------------------------------
# rbindArraysbyPointer

# Takes a list of pointers to arrays of the same width (and same type), rbinds
# the arrays in order and gives back the pointer to the combined array
# while deleting all others. For reasons of convenience this one return value
# is a one-element list.

rbindArraysbyPointer <- function(pArray) {

  nn <- length(pArray)
  if(nn > 1) {
    for(i in (2:nn)) {
      # Always append the next one
      pArray[[1]]$value <- rbind(pArray[[1]]$value, pArray[[2]]$value)
      # Remove just appended element from list. Next element becomes #2.
      pArray[[2]]$value <- NULL
      pArray[[2]]       <- NULL
    } # for i
  } # if length(nn > 1)

  return(pArray[[1]]) # Useful though looking strange

} # function rbindArraysbyPointer

# ------------------------------------------------------------------------------
# verticalSplitIndex

# Returns the row indexes for splitting an array with nnRows as required for
# parallel processing with nnCores cores.
# Used in the functions complexArrayPlot and phasePortrait.

verticalSplitIndex <- function(nnRows, nnCores) {

  # nProcss: Number of processes to be parallelized
  nProcss        <- min(nnRows, nnCores)
  nRowsPerChunk  <- nnRows %/% nProcss
  nRest          <- nnRows %%  nProcss
  iProcss        <- c(1:nProcss)
  iPerChunk      <- c(1:nRowsPerChunk)
  lower          <- 1 + nRowsPerChunk * (iProcss - 1)
  upper          <- nRowsPerChunk * (iProcss)
  upper[nProcss] <- upper[nProcss] + nRest
  uplow          <- asplit(cbind(lower, upper), MARGIN = 1)

  return(uplow)

} # function verticalSplitIndex

# ------------------------------------------------------------------------------
# Function complexArrayPlot

# Displays an array of complex numbers in an existing plot.
# In order to do so,the temporary files that together form the array are
# read from disk one by one, but each one is processed in a parallel loop.
# The resulting array of hsv color values is finally plotted as
# a raster image.

complexArrayPlot <- function(zMetaInfrm, xlim, ylim,
                             pType = "pma", invertFlip = FALSE,
                             lambda = 7, gamma = 9/10, pi2Div = 9,
                             logBase = exp(2*pi/pi2Div),
                             argOffset = 0,
                             stdSaturation = 0.8,
                             darkestShade = 0.1,
                             hsvNaN = c(0, 0, 0.5),
                             asp = 1,
                             xlab = "", ylab = "",
                             verbose,
                             ...) {

  # Set up plot
  plot(NULL, xlim = xlim, ylim = ylim, asp = asp, xlab = xlab, ylab = ylab, ...)

  # Define call to color transformation function depending user's
  # choice of pType
  colCmd <- switch(pType,
                   "p"   = "phaseColhsv(pListCompArr[[i]],
                                        pHsvCol,
                                        stdSaturation = stdSaturation,
                                        hsvNaN = hsvNaN)",

                   "pm"  = "phaseModColhsv(pListCompArr[[i]],
                                           pHsvCol,
                                           lambda = lambda,
                                           logBase = logBase,
                                           stdSaturation = stdSaturation,
                                           darkestShade = darkestShade,
                                           hsvNaN = hsvNaN)",

                   "pa"  = "phaseAngColhsv(pListCompArr[[i]],
                                           pHsvCol,
                                           lambda = lambda,
                                           pi2Div = pi2Div,
                                           argOffset = argOffset,
                                           stdSaturation = stdSaturation,
                                           darkestShade = darkestShade,
                                           hsvNaN = hsvNaN)",

                   "pma" = "phaseModAngColhsv(pListCompArr[[i]],
                                              pHsvCol,
                                              lambda = lambda,
                                              gamma = gamma,
                                              pi2Div = pi2Div,
                                              logBase = logBase,
                                              argOffset = argOffset,
                                              stdSaturation = stdSaturation,
                                              darkestShade = darkestShade,
                                              hsvNaN = hsvNaN)"
  ) # switch


  # Obtain the names of the files to load and process
  zMetaInfrm$metaZ$wFileNames <- paste(zMetaInfrm$tempDir,
                                       zMetaInfrm$metaZ$wFileNames, sep = "/")

  # Run the color transformation function over each file
  pHsvCol <- lapply(c(1:nrow(zMetaInfrm$metaZ)),
                    function(i, zMetaInfrm, colCmd) {

      if(verbose) cat("\n.transforming block", i, "... ")
      # load a block (will soon become a list of pointers, hence the name)
      pListCompArr  <- get(load(zMetaInfrm$metaZ[i,]$wFileNames))
      # split it for parallel processing
      nCores   <- getDoParWorkers()
      uplow    <- verticalSplitIndex(nrow(pListCompArr), nCores)

      # here's the actual splitting, pListCompArr becomes a list of pointers
      pListCompArr  <- lapply(uplow, FUN = function(uplow, pListCompArr) {
        nwPtr <- newPointer(pListCompArr[c(uplow[1]:uplow[2]),])
        # if the split result has only one line, it will automatically become a
        # vector, which is undesired, because functions coming later require it
        # as a two-dimensional array. This is made sure here.
        if(uplow[1] == uplow[2]) {
          dim(nwPtr$value) <- c(1, length(nwPtr$value))
        }
        return(nwPtr)
      }, pListCompArr = pListCompArr)

      # Parallel loop transforming the chunks into a color raster each;
      # giving back a list of pointers to the rasters
      if(verbose) cat("parallel loop starting ... ")
      pHsvCol <- foreach(i = c(1:length(pListCompArr)),
                         .export  = c("phaseColhsv",
                                      "phaseModColhsv",
                                      "phaseAngColhsv",
                                      "phaseModAngColhsv",
                                      "logBase",
                                      "lambda",
                                      "gamma",
                                      "pi2Div",
                                      "stdSaturation",
                                      "darkestShade",
                                      "hsvNaN",
                                      "newPointer",
                                      "argOffset"),
                         .combine = c) %dopar% {

        pHsvCol  <- newPointer(NULL)
        eval(parse(text = colCmd))      # Does not require a return value,
                                        # changes color array via pointer
        pListCompArr[[i]]$value <- NULL # Reduced here, but removed after
                                        # the foreach loop
        return(pHsvCol)
      } # foreach
      if(verbose) cat("done.")

      # Remove the original list of array pointers
      rm(pListCompArr)

      # Combine the color arrays in the value of the first pointer.
      # Free the others (rbindArraysbyPointer).
      #   Enforce (one-element-) list in case there is only one value
      #   (i.e. if foreach loop was executed sequentially, one core only)
      if(length(pHsvCol) == 1) pHsvCol <- list(pHsvCol)
      pHsvCol <- rbindArraysbyPointer(pHsvCol)

      return(pHsvCol)
    }, # function in lapply
    zMetaInfrm = zMetaInfrm, colCmd = colCmd
  ) # lapply

  # Now combine all blocks into the big raster ...
  if(verbose) cat("\nCombine color rasters ... ")
  pHsvCol <- rbindArraysbyPointer(pHsvCol)
  if(verbose) cat("done.\n")

  # ... and plot it
  if(verbose) cat("Plotting raster image ... ")
  rasterImage(as.raster(pHsvCol$value), xlim[1], ylim[1], xlim[2], ylim[2])
  if(verbose) cat("done.\n")

  pHsvCol$value <- NULL
  rm(pHsvCol)

  return(NULL)

} # complexArrayPlot

# -----------------------------------------------------------------------------
# Function makeFunctionFromInput

# Transform an input FUN and moreArgs (a named list), both inputs to
# phasePortrait, into an executable function.

makeFunctionFromInput <- function(FUN, moreArgs = NULL) {

  # If match.fun() detects a function, give it back. If not, return NULL
  testFun <- tryCatch(match.fun(FUN), error = function(err) NULL)

  # If this does not work, maybe we have a useful character string
  if(is.null(testFun) & mode(FUN) == "character") {
    if(!is.null(moreArgs)) {
      moreArgString <- paste(",", names(moreArgs), collapse = "")
    }
    else {
      moreArgString <- ""
    }
    exprText <- paste("function(z", moreArgString, ") ", FUN, sep = "")
    testFun  <- eval(parse(text = exprText))
  } # if character

  # Test the function if the above resulted in something else
  # than NULL
  if(!is.null(testFun)) {
    # Arbitrary number for testing the function
    testNum <- complex(real = runif(1), imaginary = runif(1))
    if(!is.null(moreArgs)) testArgs      <- c(testNum, moreArgs)
    else                   testArgs      <- list(testNum)
    # if NULL comes back from this call, the function does not work
    testOut <- tryCatch(do.call(testFun, testArgs),
                        error = function(err) NULL)
    if(is.null(testOut)) testFun <- NULL
  } # if(!is.null(compFun))

  return(testFun)

} # makeFunctionFromInput


# -----------------------------------------------------------------------------
#' Create phase portraits of complex functions
#'
#' \code{phasePortrait} makes phase portraits of functions in the complex number
#' plane. It uses a technique often (but not quite correctly) called
#' \emph{domain coloring} (\url{https://en.wikipedia.org/wiki/Domain_coloring}).
#' While many varieties of this technique exist, this book relates closely to
#' the standards proposed by E. Wegert in his book \emph{Visual Complex
#' Functions} \insertCite{wegert_visualcpx_2012}{viscomplexr}. In a nutshell,
#' the argument (\code{\link{Arg}}) of any complex function value is displayed
#' as a color from the chromatic circle. The fundamental colors red, green, and
#' blue relate to the arguments (angles) of 0, 2/3pi, and 4/3pi (with smooth
#' color transitions in between), respectively. Options for displaying the
#' modulus (\code{\link{Mod}}) of the complex values and additional reference
#' lines for the argument are available. This function is designed for being
#' used inside the framework of R base graphics. It makes use of parallel
#' computing, and depending on the desired resolution it may create extensive
#' sets of large temporary files (see Details and Examples).
#'
#' This function is intended to be used inside the framework of R base graphics.
#' It plots into the active open graphics device where it will display the phase
#' plot of a user defined function as a raster image. If no graphics device is
#' open when called, the function will plot into the default graphics device.
#' This principle allows to utilize the full functionality of R base graphics.
#' All graphics parameters (\code{\link{par}}) can be freely set and the
#' function \code{phasePortrait} accepts all parameters that can be passed to
#' the \code{\link{plot.default}} function. This allows all kinds of plots -
#' from scientific representations with annotated axes and auxiliary lines,
#' notation, etc. to poster-like artistic pictures.\cr
#'
#' \describe{
#'   \item{Mode of operation}{After being called, \code{phasePortrait} gets the
#'   size in inch of the plot region of the graphics device it is plotting into.
#'   With the parameter \code{res} which is the desired plot resolution in dpi,
#'   the horizontal and vertical number of pixels is known. As \code{xlim} and
#'   \code{ylim} are provided by the user, each pixel can be attributed a
#'   complex number z from the complex plane. In that way a two-dimensional
#'   array is built, where each cell represents a point of the complex plane,
#'   containing the corresponding complex number z. This array is set up in
#'   horizontal strips (i.e. split along the imaginary axis), each strip
#'   containing approximately \code{blockSizePx} pixels. In a parallel computing
#'   loop, the strips are constructed, saved as temporary files and immediately
#'   deleted from the RAM in order to avoid memory overflow. After that, the
#'   strips are sequentially loaded and subdivided into a number of chunks that
#'   corresponds to the number of registered parallel workers (parameter
#'   \code{nCores}). By parallely processing each chunk, the function
#'   \code{f(z)} defined by the user in the argument \code{FUN} is applied to
#'   each cell of the strip. This results in an array of function values that
#'   has exactly the same size as the original strip. The new array is saved as
#'   a temporary file, the RAM is cleared, and the next strip is loaded. This
#'   continues until all strips are processed. In a similar way, all strips
#'   containing the function values are loaded sequentially, and in a parallel
#'   process the complex values are translated into colors which are stored in a
#'   raster object. While the strips are deleted from the RAM after processing,
#'   the color values obtained from each new strip are appended to the color
#'   raster. After all strips are processed, the raster is plotted into the plot
#'   region of the graphics device. If not explicitly defined otherwise by the
#'   user, all temporary files are deleted after that.
#'   }
#'   \item{Temporary file system}{By default, the above-mentioned temporary
#'   files are deleted after use. This will not happen, if the parameter
#'   \code{deleteTempFiles} is set to \code{FALSE} or if \code{phasePortrait}
#'   does not terminate properly. In both cases, you will find the files in the
#'   directory specified with the parameter \code{tempDir}. These files are
#'   \code{.RData} files, each one contains a two-dimensional array of complex
#'   numbers. The file names follow a strict convention, see the following
#'   examples:\cr\cr
#'   \code{0001zmat2238046385.RData}\cr
#'   \code{0001wmat2238046385.RData}\cr\cr
#'   Both names begin with '0001', indicating that the array's top line is the
#'   first line of the whole clipping of the complex number plane where the
#'   phase portrait relates to. The array which follows below can e.g. begin
#'   with a number like '0470', indicating that its first line is line number
#'   470 of the whole clipping. The number of digits for these line numbers is
#'   not fixed. It is determined by the greatest number required. Numbers with
#'   less digits are zero-padded. The second part of the file name is either
#'   \code{zmat} or \code{wmat}. The former indicates an array whose cells
#'   contain untransformed numbers of the complex number plane. The latter
#'   contains the values obtained from applying the function of interest to the
#'   first array. Thus, cells at the same position in both arrays exactly relate
#'   to each other. The third part of the file names is a ten-digit integer.
#'   This is a random number which all temporary files stemming from the same
#'   call of \code{phasePortrait} have in common. This guarantees that no
#'   temporary files will be confounded by the function, even if undeleted
#'   temporary files from previous runs are still present.
#'   }
#'   \item{HSV color model}{For color-coding the argument of a complex number,
#'   \code{phasePortrait} uses the \code{\link{hsv}} (hue, saturation, value)
#'   color model. Hereby, the argument is mapped to a position on the chromatic
#'   circle, where the fundamental colors red, green, and blue relate to the
#'   arguments (angles) of 0, 2/3pi, and 4/3pi, respectively. This affects only
#'   the hue component of the color model. The value component is used for
#'   shading modulus and/or argument zones. The saturation component for all
#'   colors can be defined with the parameter \code{stdSaturation}.
#'   }
#'   \item{Zone definitions and shading}{In addition to displaying colors for
#'   the arguments of complex numbers, zones for the modulus and/or the argument
#'   are shaded for \code{pType} other than "p". The modulus zones are defined
#'   in a way that each zone covers moduli whose logarithms to the base
#'   \code{logBase} have the same integer part. Thus, from the lower edge of one
#'   modulus zone to its upper edge, the modulus multiplies with the value of
#'   \code{logBase}. The shading of a modulus zone depends on the fractional
#'   parts \code{x} of the above-mentioned logarithms, which cover the interval
#'   \code{[0, 1[}.
#'   This translates into the value component \code{v} of the \code{\link{hsv}}
#'   color model as follows:\cr\cr
#'   \code{v = darkestShade + (1 - darkestShade) * x^(1/lambda)}\cr\cr
#'   where \code{darkestShade} and \code{lambda} are parameters that can be
#'   defined by the user. Modifying the parameters \code{lambda} and
#'   \code{darkestShade} is useful for adjusting contrasts in the phase
#'   portraits. The argument zone definition is somewhat simpler: Each zone
#'   covers an angle domain of \code{2*pi / pi2Div}, the "zero reference" for
#'   all zones being \code{argOffset}. The angle domain of one zone is linearly
#'   mapped to a value \code{x} from the range \code{[0, 1[}.
#'   The value component of the color to be displayed is calculated as a
#'   function of \code{x} with the same equation as shown above. In case the
#'   user has chosen \code{pType = "pma"}, x-values \code{xMod} and \code{xArg}
#'   are calculated separately for the modulus and the argument, respectively.
#'   They are transformed into preliminary v-values as follows:\cr\cr
#'   \code{vMod = xMod^(1/lambda)} and {vArg = xArg^(1/lambda)}\cr\cr
#'   From these, the final v value results as\cr\cr
#'   \code{v = darkestShade + (1-darkestShade) * (gamma * vMod * vArg +
#'   (1-gamma) * (1 - (1-vMod) * (1-vArg)))}\cr\cr
#'   The parameter \code{gamma} (range \code{[0, 1]}) determines they way how
#'   vMod and vArg are combined. The closer \code{gamma} is to one, the more
#'   the smaller of both values will dominate the outcome and vice versa.
#'   }
#'   \item{Defining more complicated functions as strings with
#'   \code{\link{vapply}}}{You might want to write and use functions which
#'   require more code than a single statement like \code{(z-3)^2+1i*z}. As
#'   mentioned in the description of the parameter \code{FUN}, we recommend to
#'   define such functions as separate objects and hand them over as such. There
#'   might be, however, cases, where it is more convenient, to define a function
#'   as a single long string, and pass this string to \code{FUN}.
#'   In order to make this work, \code{\link{vapply}} should be be used for
#'   wrapping the actual code of the function. This is probably not the use of
#'   \code{\link{vapply}} intended by its developers, but it works nicely and
#'   performs well. The character string has to have the following structure
#'   "vapply(z, function(z, \emph{other arguments if required}) \{\emph{define
#'   function code in here}\}, \emph{define other arguments here}, FUN.VALUE =
#'   complex(1))". See examples.
#'   }
#' }
#'
#'
#' @param FUN The function to be visualized. There are two possibilities to
#'   provide it, a quoted character string, or a function object.
#'   \describe{
#'   \item{Quoted character string}{ A function can be provided as a quoted
#'   character string containing an expression R can interpret as a function of
#'   a complex number z. Examples: "sin(z)", "(z^2 - 1i)/(tan(z))", "1/4*z^2 -
#'   10*z/(z^4+4)". Names of functions known in your R session can be used in a
#'   standalone way, without mentioning z, e.g. "sin", "tan", "sqrt". Obviously,
#'   this also works for functions you defined yourself, e.g.
#'   "myIncredibleFunction" would be valid if you coded a function with this
#'   name before. This is especially useful for functions which require
#'   additional parameters beside the complex number they are supposed to
#'   calculate with. Such arguments can be provided via the parameter
#'   \code{moreArgs}. One-liner expressions provided as strings are also
#'   compatible with \code{moreArgs} (see examples).
#'
#'   While it is not the way we recommend for most purposes, you can even define
#'   more complicated functions of your own as character strings. In this case,
#'   you need to use \code{\link{vapply}} as a wrapper for your actual function
#'   (see Details, and Examples). Such constructions allow to provide additional
#'   input variables as a part of the character string by using the
#'   \code{\link{vapply}}-mechanism (see Details and Examples). The helper
#'   function \code{\link{vector2String}}) can be useful for that matter.
#'   However, the parameter \code{moreArgs} is not applicable in this context.
#'   Probably, the most useful application of the function-as-string concept is
#'   when the user defined function, possibly including values for additional
#'   arguments, is to be pasted together at runtime.}
#'
#'   \item{Function object}{ It is also possible
#'   to directly provide function objects to \code{FUN}. This can be any
#'   function known to R in the current session. Simply put, for functions like
#'   sin, tan, cos, and sqrt you do not even have to quote their names when
#'   passing them to \code{phasePortrait}. Same applies to functions you defined
#'   yourself. It is also possible to hand over an anonymous function to
#'   \code{FUN} when calling \code{phasePortrait}. In all these cases, the
#'   parameter \code{moreArgs} can be used for providing additional arguments to
#'   \code{FUN}. In general, providing a function as an object, and using
#'   \code{moreArgs} in case additional arguments are required, is what we
#'   recommend for user-defined functions.}
#'   }
#'
#'   When executing \code{phasePortrait}, \code{FUN} is first evaluated with
#'   \code{\link{match.fun}}. If this is not successful, an attempt to interpret
#'   \code{FUN} as an expression will be made. If this fails,
#'   \code{phasePortrait} terminates with an error.
#'
#' @param moreArgs A named list of other arguments to FUN. The names must match
#'   the names of the arguments in FUN's definition.
#'
#' @param xlim The x limits (x1, x2) of the plot as a two-element numeric
#'   vector. Follows exactly the same definition as in
#'   \code{\link{plot.default}}. Here, \code{xlim} has to be interpreted as the
#'   plot limits on the real axis.
#'
#' @param ylim The y limits of the plot (y1, y2) to be used in the same way as
#'   \code{xlim}. Evidently, \code{ylim} indicates the plot limits on the
#'   imaginary axis.
#'
#' @param invertFlip If \code{TRUE}, the function is mapped over a z plane,
#'   which has been transformed to \code{1/z * exp(1i*pi)}. This is the
#'   projection required to plot the north Riemann hemisphere in the way
#'   proposed by \insertCite{wegert_visualcpx_2012;textual}{viscomplexr}, p. 41.
#'   Defaults to \code{FALSE}. If this option is chosen, the numbers at the
#'   axis ticks have another meaning than in the normal case. Along the real
#'   axis, they represent the real part of \code{1/z}, and along the imaginary
#'   axis, they represent the imaginary part of \code{1/z}. Thus, if you want
#'   annotation, you should choose appropriate axis labels like \code{xlab =
#'   Re(1/z)}, and \code{ylab = Im(1/z)}.
#'
#' @param res Desired resolution of the plot in dots per inch (dpi). Default is
#'   150 dpi. All other things being equal, \code{res} has a strong influence on
#'   computing times (double \code{res} means fourfold number of pixels to
#'   compute). A good approach could be to make a plot with low resolution (e.g.
#'   the default 150 dpi) first, adjust whatever required, and plot into a
#'   graphics file with high resolution after that.
#'
#' @param blockSizePx Number of pixels and corresponding complex values to be
#'   processed at the same time (see Details). Default is 2250000. This value
#'   gave good performance on older systems as well as on a high-end gaming
#'   machine, but some tweaking for your individual system might even improve
#'   things.
#'
#' @param tempDir \code{NULL} or a character string, specifying the name of the
#'   directory where the temporary files written by \code{phasePortrait} are
#'   stored. Default is \code{NULL}, which makes \code{phasePortrait} use the
#'   current R session's temporary directory. Note that if you specify another
#'   directory, it will be created if it does not exist already. Even though the
#'   temporary files are deleted after completing a phase portrait (unless the
#'   user specifies \code{deleteTempFiles = FALSE}, see below), the directory
#'   will remain alive even if has been created by \code{phasePortrait}.
#'
#' @param nCores Number of processor cores to be used in the parallel computing
#'   tasks. Defaults to the maximum number of cores available. Any number
#'   between 1 (serial computation) and the maximum number of cores available as
#'   indicated by \code{parallel::detectCores()} is accepted.
#'
#' @param pType One of the four options for plotting, "p", "pa", "pm", and "pma"
#'   as a character string. Defaults to "pma". Option "p" produces a mere phase
#'   plot, i.e. contains only colors for the complex numbers' arguments, but no
#'   reference lines at all. the option "pa" introduces shading zones that
#'   emphasize the arguments. These zones each cover an angle defined by
#'   \code{2*pi/pi2Div}, where p2Div is another parameter of this function (see
#'   there). These zones are shaded darkest at the lowest angle (counter
#'   clockwise). Option "pm" displays the modulus by indicating zones, where the
#'   moduli at the higher edge of each zone are in a constant ratio with the
#'   moduli at the lower edge of the zone. Default is a ratio of almost exactly
#'   2 (see parameter \code{logBase}) for details. At the lower edge, color
#'   saturation is lowest and highest at the higher edge (see parameters
#'   \code{darkestShade}, and \code{stdSaturation}). Option "pma" (default)
#'   includes both shading schemes.
#'
#' @param pi2Div Angle distance for the argument reference zones added for
#'   \code{pType = "pma"} or \code{pType = "pa"}. The value has to be given as
#'   an integer (reasonably) fraction of 2*pi (i.e. 360 degrees). The default is
#'   9; thus, reference zones are delineated by default in distances of 2*pi/9,
#'   i.e. (40 degrees), starting with 0, i.e. the color red if not defined
#'   otherwise with the parameter \code{argOffset}. In contrast to the borders
#'   delimiting the modulus zones, the borders of the reference zones for the
#'   argument always follow the same color (by definition).
#'
#' @param logBase Modulus ratio between the edges of the modulus reference zones
#'   in \code{pType} \code{"pm"} and \code{"pma"}. As recommended by
#'   \insertCite{wegert_visualcpx_2012;textual}{viscomplexr}, the default
#'   setting is \code{logBase = exp(2*pi/pi2Div)}. This relation between the
#'   parameters \code{logBase} and \code{pi2Div} ensures an analogue scaling of
#'   the modulus and argument reference zones (see Details). Conveniently, for
#'   the default \code{pi2Div = 9}, we obtain \code{logBase == 2.0099...},
#'   which is very close to 2. Thus, the modulus at the higher edge of a given
#'   zone is almost exactly two times the value at the lower edge.
#'
#' @param argOffset The (complex number) argument in radians counterclockwise,
#'   at which the argument reference zones are fixed. Default is 0, i.e. all
#'   argument reference zones align to the center of the red area.
#'
#' @param darkestShade Darkest possible shading of modulus and angle reference
#'   zones for \code{pType} \code{"pm"} and \code{"pma"}. It corresponds to the
#'   value "v" in the \code{\link{hsv}} color model. \code{darkestShade = 0}
#'   means no brightness at all, i.e. black, while \code{darkestShade = 1}
#'   indicates maximum brightness. Defaults to 0.1, i.e. very dark, but hue
#'   still discernible.
#'
#' @param lambda Parameter steering the shading interpolation between the higher
#'   and the lower edges of the the modulus and argument reference zones in
#'   \code{pType} \code{"pm"} and \code{"pm"}. Should be > 0, default and
#'   reference is \code{lambda = 7}. Values < 7 increase the contrast at the
#'   zone borders, values > 7 weaken the contrast.
#'
#' @param gamma Parameter for adjusting the combined shading of modulus and
#'   argument reference zones in \code{pType} \code{"pma"}. Should be in the
#'   interval \code{[0, 1]}. Default is 0.9. The higher the value, the more the
#'   smaller of both shading values will dominate the outcome and vice versa.
#'
#' @param stdSaturation Saturation value for unshaded hues which applies to the
#'   whole plot in \code{pType} \code{"p"} and to the (almost) unshaded zones in
#'   \code{pType} \code{"pm"} and \code{"p"}. This corresponds to the value "s"
#'   in the \code{\link{hsv}} color model. Must be between 0 and 1, where 1
#'   indicates full saturation and 0 indicates a neutral grey. Defaults to 0.8.
#'
#' @param hsvNaN \code{\link{hsv}} coded color for being used in areas where the
#'   function to be plotted is not defined. Must be given as a numeric vector
#'   with containing the values h, s, and v in this order. Defaults to
#'   \code{c(0, 0, 0.5)} which is a neutral grey.
#'
#' @param asp Aspect ratio y/x as defined in \code{\link{plot.window}}. Default
#'   is 1, ensuring an accurate representation of distances between points on
#'   the screen.
#'
#' @param deleteTempFiles If TRUE (default), all temporary files are deleted
#'   after the plot is completed. Set it on FALSE only, if you know exactly what
#'   you are doing - the temporary files can occupy large amounts of hard disk
#'   space (see details).
#'
#' @param noScreenDevice Suppresses any graphical output if TRUE. This is only
#'   intended for test purposes and makes probably only sense together with
#'   \code{deleteTempFiles == FALSE}. For dimensioning purposes,
#'   \code{phasePortrait} will use a 1 x 1 inch pseudo graphics device in this
#'   case. The default for this parameter is \code{FALSE}, and you should change
#'   it only if you really know what you are doing.
#'
#' @param autoDereg if TRUE, automatically register sequential backend after the
#'   phase portrait is completed. Default is FALSE, because registering a
#'   parallel backend can be time consuming. Thus, if you want make several
#'   phase portraits in succession, you should set \code{autoDereg} only for the
#'   last one, or simply type \code{foreach::registerDoSEQ} after you are done.
#'   In any case, you don't want to have an unused parallel backend lying about.
#'
#' @param verbose if TRUE (default), \code{phasePortrait} will continuously
#'   write progress messages to the console. This is convenient for normal
#'   purposes, as calculating larger phase portraits in higher resolution may
#'   take several minutes. The setting \code{verbose = FALSE}, will suppress any
#'   output to the console.
#'
#' @param ... All parameters accepted by the \code{\link{plot.default}}
#'   function.
#'
#'
#' @references
#'   \insertAllCited{}
#'
#'
#' @examples
#' # Map the complex plane on itself
#'
#' # x11(width = 8, height = 8)   # Screen device commented out
#'                                # due to CRAN test requirements.
#'                                # Use it when trying this example
#' phasePortrait("z", xlim = c(-2, 2), ylim = c(-2, 2),
#'               xlab = "real", ylab = "imaginary",
#'               verbose = FALSE, # Suppress progress messages
#'               nCores = 2)      # Max. two cores allowed on CRAN
#'                                # not a limit for your own use
#' \dontshow{
#' # R CMD check: make sure any open connections are closed afterward
#' foreach::registerDoSEQ()
#' doParallel::stopImplicitCluster()
#' }
#'
#'
#' # A rational function
#' \donttest{
#' # x11(width = 10, height = 8) # Screen device commented out
#'                               # due to CRAN test requirements.
#'                               # Use it when trying this example
#' phasePortrait("(2-z)^2*(-1i+z)^3*(4-3i-z)/((2+2i+z)^4)",
#'               xlim = c(-8, 8), ylim = c(-6.3, 4.3),
#'               xlab = "real", ylab = "imaginary",
#'               nCores = 2)     # Max. two cores allowed on CRAN
#'                               # not a limit for your own use
#'   \dontshow{
#'   # R CMD check: make sure any open connections are closed afterward
#'   foreach::registerDoSEQ()
#'   doParallel::stopImplicitCluster()
#'   }
#' }
#'
#'
#' # Different pType options by example of the sine function.
#' # Note the different equivalent definitions of the sine
#' # function in the calls to phasePortrait
#' \donttest{
#' # x11(width = 9, height = 9) # Screen device commented out
#'                              # due to CRAN test requirements.
#'                              # Use it when trying this example
#' op <- par(mfrow = c(2, 2), mar = c(2.1, 2.1, 2.1, 2.1))
#' phasePortrait("sin(z)", xlim = c(-pi, pi), ylim = c(-pi, pi),
#'               pType = "p",   main = "pType = 'p'",   axes = FALSE,
#'               nCores = 2) # Max. two cores on CRAN, not a limit for your use
#'   \dontshow{
#'   # R CMD check: make sure any open connections are closed afterward
#'   foreach::registerDoSEQ()
#'   doParallel::stopImplicitCluster()
#'   }
#' phasePortrait("sin(z)", xlim = c(-pi, pi), ylim = c(-pi, pi),
#'               pType = "pm",  main = "pType = 'pm'",  axes = FALSE,
#'               nCores = 2)
#'   \dontshow{
#'   # R CMD check: make sure any open connections are closed afterward
#'   foreach::registerDoSEQ()
#'   doParallel::stopImplicitCluster()
#'   }
#' phasePortrait("sin",    xlim = c(-pi, pi), ylim = c(-pi, pi),
#'               pType = "pa",  main = "pType = 'pa'",  axes = FALSE,
#'               nCores = 2)
#'   \dontshow{
#'   # R CMD check: make sure any open connections are closed afterward
#'   foreach::registerDoSEQ()
#'   doParallel::stopImplicitCluster()
#'   }
#' phasePortrait(sin,      xlim = c(-pi, pi), ylim = c(-pi, pi),
#'               pType = "pma", main = "pType = 'pma'", axes = FALSE,
#'               nCores = 2)
#'   \dontshow{
#'   # R CMD check: make sure any open connections are closed afterward
#'   foreach::registerDoSEQ()
#'   doParallel::stopImplicitCluster()
#'   }
#' par(op)}
#'
#'
#' # I called this one 'nuclear fusion'
#'
#' # x11(width = 16/9*8, height = 8) # Screen device commented out
#'                                   # due to CRAN test requirements.
#'                                   # Use it when trying this example
#' \donttest{
#' op <- par(mar = c(0, 0, 0, 0), omi = c(0.2, 0.2, 0.2, 0.2), bg = "black")
#' phasePortrait("cos((z + 1/z)/(1i/2 * (z-1)^10))",
#'               xlim = 16/9*c(-2, 2), ylim = c(-2, 2),
#'               axes = FALSE, xaxs = "i", yaxs = "i",
#'               nCores = 2) # Max. two cores allowed on CRAN
#'                           # not a limit for your own use
#'   \dontshow{
#'   # R CMD check: make sure any open connections are closed afterward
#'   foreach::registerDoSEQ()
#'   doParallel::stopImplicitCluster()
#'   }
#' par(op)}
#'
#'
#' # Passing function objects to phasePortrait:
#' # Two mathematical celebrities - Riemann's zeta function
#' # and the gamma function, both from the package pracma.
#' # R's built-in gamma is not useful, as it does not work
#' # with complex input values.
#' \donttest{
#' if(requireNamespace("pracma", quietly = TRUE)) {
#' # x11(width = 16, height = 8) # Screen device commented out
#'                               # due to CRAN test requirements.
#'                               # Use it when trying this example
#' op <- par(mfrow = c(1, 2))
#' phasePortrait(pracma::zeta,  xlim = c(-35, 15), ylim = c(-25, 25),
#'               xlab = "real", ylab = "imaginary",
#'               main = expression(zeta(z)), cex.main = 2,
#'               nCores = 2) # Max. two cores on CRAN, not a limit for your use
#'   \dontshow{
#'   # R CMD check: make sure any open connections are closed afterward
#'   foreach::registerDoSEQ()
#'   doParallel::stopImplicitCluster()
#'   }
#' phasePortrait(pracma::gammaz, xlim = c(-10, 10), ylim = c(-10, 10),
#'               xlab = "real", ylab = "imaginary",
#'               main = expression(Gamma(z)), cex.main = 2,
#'               nCores = 2) # Max. two cores allowed on CRAN
#'                           # not a limit for your own use
#'   \dontshow{
#'   # R CMD check: make sure any open connections are closed afterward
#'   foreach::registerDoSEQ()
#'   doParallel::stopImplicitCluster()
#'   }
#' }
#' }
#'
#'
#' # Using vapply for defining a whole function as a string.
#' # This is a Blaschke product with a sequence a of twenty numbers.
#' # See the documentation of the function vector2String for a more
#' # convenient space-saving definition of a.
#' # But note that a C++ version of the Blaschke product is available
#' # in this package (function blaschkeProd()).
#' \donttest{
#' # x11(width = 10, height = 8) # Screen device commented out
#'                               # due to CRAN test requirements.
#'                               # Use it when trying this example
#' phasePortrait("vapply(z, function(z, a) {
#'                 fct <- ifelse(abs(a) != 0,
#'                   abs(a)/a * (a-z)/(1-Conj(a)*z), z)
#'                 return(prod(fct))
#'               },
#'               a = c(0.12152611+0.06171533i,  0.53730315+0.32797530i,
#'                     0.35269601-0.53259644i, -0.57862039+0.33328986i,
#'                    -0.94623221+0.06869166i, -0.02392968-0.21993132i,
#'                     0.04060671+0.05644165i,  0.15534449-0.14559097i,
#'                     0.32884452-0.19524764i,  0.58631745+0.05218419i,
#'                     0.02562213+0.36822933i, -0.80418478+0.58621875i,
#'                    -0.15296208-0.94175193i, -0.02942663+0.38039250i,
#'                    -0.35184130-0.24438324i, -0.09048155+0.18131963i,
#'                     0.63791697+0.47284679i,  0.25651928-0.46341192i,
#'                     0.04353117-0.73472528i, -0.04606189+0.76068461i),
#'               FUN.VALUE = complex(1))",
#'               pType = "p",
#'               xlim = c(-4, 2), ylim = c(-2, 2),
#'               xlab = "real", ylab = "imaginary",
#'               nCores = 2) # Max. two cores allowed on CRAN
#'                           # not a limit for your own use
#'   \dontshow{
#'   # R CMD check: make sure any open connections are closed afterward
#'   foreach::registerDoSEQ()
#'   doParallel::stopImplicitCluster()
#'   }
#' }
#'
#'
#' # Much more elegant: Define the function outside.
#' # Here comes a Blaschke product with 200 random points.
#' \donttest{
#' # define function for calculating blaschke products, even
#' # possible as a one-liner
#' blaschke <- function(z, a) {
#'   return(prod(ifelse(abs(a) != 0, abs(a)/a * (a-z)/(1-Conj(a)*z), z)))
#' }
#' # define 200 random numbers inside the unit circle
#' n <- 200
#' a <- complex(modulus = runif(n), argument = runif(n)*2*pi)
#' # Plot it
#' # x11(width = 10, height = 8) # Screen device commented out
#'                               # due to CRAN test requirements.
#'                               # Use it when trying this example
#' phasePortrait(blaschke,
#'   moreArgs = list(a = a),
#'   pType = "p",
#'   xlim = c(-2.5, 2.5), ylim = c(-1.7, 1.7),
#'   xlab = "real", ylab = "imaginary",
#'   nCores = 2) # Max. two cores allowed on CRAN
#'               # not a limit for your own use
#'   \dontshow{
#'   # R CMD check: make sure any open connections are closed afterward
#'   foreach::registerDoSEQ()
#'   doParallel::stopImplicitCluster()
#'   }
#' }
#'
#'
#' # A hybrid solution: A one-liner expression given as a character string
#' # can be provided additional arguments with moreArgs
#' \donttest{
#' n <- 73
#' a <- complex(modulus = runif(n), argument = runif(n)*2*pi)
#' # x11(width = 10, height = 8) # Screen device commented out
#'                               # due to CRAN test requirements.
#'                               # Use it when trying this example
#' phasePortrait("prod(ifelse(abs(a) != 0,
#'   abs(a)/a * (a-z)/(1-Conj(a)*z), z))",
#'   moreArgs = list(a = a),
#'   pType = "p",
#'   xlim = c(-2.5, 2.5), ylim = c(-1.7, 1.7),
#'   xlab = "real", ylab = "imaginary",
#'   nCores = 1) # Max. two cores allowed on CRAN
#'               # not a limit for your own use
#'   \dontshow{
#'   # R CMD check: make sure any open connections are closed afterward
#'   foreach::registerDoSEQ()
#'   doParallel::stopImplicitCluster()
#'   }
#' }
#'
#'
#' # Note the difference in performance when using the C++ defined
#' # function blaschkeProd() provided in this package
#' \donttest{
#' n <- 73
#' a <- complex(modulus = runif(n), argument = runif(n)*2*pi)
#' # Plot it
#' # x11(width = 10, height = 8) # Screen device commented out
#'                               # due to CRAN test requirements.
#'                               # Use it when trying this example
#' phasePortrait(blaschkeProd,
#'   moreArgs = list(a = a),
#'   pType = "p",
#'   xlim = c(-2.5, 2.5), ylim = c(-1.7, 1.7),
#'   xlab = "real", ylab = "imaginary",
#'   nCores = 1) # Max. two cores allowed on CRAN
#'               # not a limit for your own use
#'   \dontshow{
#'   # R CMD check: make sure any open connections are closed afterward
#'   foreach::registerDoSEQ()
#'   doParallel::stopImplicitCluster()
#'   }
#' }
#'
#'
#' # Interesting reunion with Benoit Mandelbrot.
#' # The function mandelbrot() is part of this package (defined
#' # in C++ for performance)
#' \donttest{
#' # x11(width = 11.7, height = 9/16*11.7) # Screen device commented out
#'                                         # due to CRAN test requirements.
#'                                         # Use it when trying this example
#' op <- par(mar = c(0, 0, 0, 0), bg = "black")
#' phasePortrait(mandelbrot,
#'               moreArgs = list(itDepth = 100),
#'               xlim = c(-0.847, -0.403), ylim = c(0.25, 0.50),
#'               axes = TRUE, pType = "pma",
#'               hsvNaN = c(0, 0, 0), xaxs = "i", yaxs = "i",
#'               nCores = 1) # Max. two cores allowed on CRAN
#'                           # not a limit for your own use
#' par(op)
#'   \dontshow{
#'   # R CMD check: make sure any open connections are closed afterward
#'   foreach::registerDoSEQ()
#'   doParallel::stopImplicitCluster()
#'   }
#' }
#'
#'
#' # Here comes a Julia set.
#' # The function juliaNormal() is part of this package (defined
#' # in C++ for performance)
#' \donttest{
#' # x11(width = 11.7, height = 9/16*11.7) # Screen device commented out
#'                                         # due to CRAN test requirements.
#'                                         # Use it when trying this example
#' op <- par(mar = c(0, 0, 0, 0), bg = "black")
#' phasePortrait(juliaNormal,
#'   moreArgs = list(c = -0.09 - 0.649i, R_esc = 2),
#'   xlim = c(-2, 2),
#'   ylim = c(-1.3, 1.3),
#'   hsvNaN = c(0, 0, 0),
#'   nCores = 1) # Max. two cores allowed on CRAN
#'               # not a limit for your own use
#' par(op)
#'   \dontshow{
#'   # R CMD check: make sure any open connections are closed afterward
#'   foreach::registerDoSEQ()
#'   doParallel::stopImplicitCluster()
#'   }
#' }
#'
#'
#' @export

phasePortrait <- function(FUN, moreArgs = NULL, xlim, ylim,
                          invertFlip = FALSE,
                          res = 150,
                          blockSizePx = 2250000,
                          tempDir = NULL,
                          nCores = parallel::detectCores(),
                          pType = "pma",
                          pi2Div = 9,
                          logBase = exp(2*pi/pi2Div),
                          argOffset = 0,
                          darkestShade = 0.1,
                          lambda = 7, gamma = 0.9,
                          stdSaturation = 0.8,
                          hsvNaN = c(0, 0, 0.5),
                          asp = 1,
                          deleteTempFiles = TRUE,
                          noScreenDevice = FALSE,
                          autoDereg = FALSE,
                          verbose = TRUE,
                          ...) {

  # Bring the user's function definition in workable form
  compFun <- makeFunctionFromInput(FUN, moreArgs)
  if(is.null(compFun)) stop("\nFUN cannot be interpreted.")

  # Calculate pixel array size from plot region size in inch and the plot
  # range for the function given with xlim and ylim

  ## par("pin"): plot region size in inch; first is horizontal
  ## if noScreenDevice, region is set to 1 x 1 inch
  if(!noScreenDevice) regionPi  <- par("pin")
  else                regionPi  <- c(1, 1)

  xRange    <- abs(xlim[2] - xlim[1])
  yRange    <- abs(ylim[2] - ylim[1])

  yxRangeRatio <- yRange      / xRange
  yxPinchRatio <- regionPi[1] / regionPi[2]

  if(yxRangeRatio < yxPinchRatio) { # height is limiting
    heightPx <- res * regionPi[2]
    widthPx  <- res * regionPi[2] / yxRangeRatio
  } # if
  else { # width is limiting
    widthPx  <- res * regionPi[1]
    heightPx <- res * regionPi[1] * yxRangeRatio
  } #else

  widthPx  <- round(widthPx)
  heightPx <- round(heightPx)

  # In case of invertFlip == TRUE swap xlim
  if(invertFlip) {
    xlim  <- c(xlim[2], xlim[1])
  } # if invertFlip

  # Register parallel Cluster if required or change number of workers
  nWorkers   <- getDoParWorkers() # number registered
  availCores <- detectCores()     # number available
  nCores     <- min(max(nCores, 1), availCores) # register at least 1 :)
                                                # and not more than available
  if(nCores != 1) {
    if(nWorkers != nCores) {
      if(verbose) cat("\nRegistering parallel workers ... ")
      registerDoSEQ()    # Unregister parallel for the sake of safety before
      registerDoParallel(cores = nCores) # register with new number of cores
      if(verbose) cat(nCores, "parallel workers registered ...")
    }
    else {
      if(verbose)
        cat("\n", nCores, " parallel workers previously registered ...",
            sep = "")
    }
  }
  # Only one core desired
  else {
    registerDoSEQ()
    if(verbose)
      cat("\nnCores set to 1.",
          "Parallel loops will be executed sequentially ...")
  }

  # Make pixelwise array of z-Values (input values to function)
  if(verbose) cat("\nBuilding z plane array ...")
  if(is.null(tempDir)) tempDir <- tempdir()
  zMetaInfrm <- buildArray(widthPx, heightPx, xlim, ylim, blockSizePx, tempDir,
                           verbose)

  # This is where it really happens
  if(verbose) cat("\nEvaluation loop starting ... ")
  zMetaInfrm$metaZ$wFileNames <- vapply(c(1:nrow(zMetaInfrm$metaZ)),
                                        function(i, zMetaInfrm, compFun,
                                                 moreArgs) {

       if(verbose) cat("\n.processing block", i, "... ")
       fileName       <- paste(zMetaInfrm$tempDir,
                               zMetaInfrm$metaZ[i,]$fileName, sep = "/")
       z              <- get(load(fileName))

       # Split z vertically (by rows) into nCores chunks to be processed
       # in parallel
       # - here's some pre-work
       uplow <- verticalSplitIndex(dim(z)[1], nCores)

       # - here's the actual splitting, z becomes a list
       z <- lapply(uplow, FUN = function(uplow, z) {
          return(z[c(uplow[1]:uplow[2]),])
       }, z = z)

       # Construct function call to be evaluated inside the parallel loop
       if(is.null(moreArgs)) {
         vCall <- "vapply(z[[i]], compFun, FUN.VALUE = complex(1))"
       }
       else {
         vCall <- paste("vapply(z[[i]], compFun, FUN.VALUE = complex(1),",
                        paste(names(moreArgs), "=", moreArgs, collapse = ","),
                        ")")
       }
       vCall <- parse(text = vCall)

       # Run the evaluation parallel on each core and put it together again
       if(verbose) cat("parallel loop starting ... ")

       w <- foreach(i = c(1:length(z)), .combine = rbind,
                    .export = c("invertFlip", "compFun")) %dopar% {
         # Make sure dimensions are correct, because
         # one-line arrays can become vectors mysteriously ...
         if(length(dim(z[[i]])) < 2) dims <- c(1, length(z[[i]]))
         else                        dims <- dim(z[[i]])
         if(invertFlip) z[[i]]            <- Conj(1 / z[[i]])
         array(eval(vCall), dim = dims)
       } # foreach i
       if(verbose) cat("done.")

       rm(z) # discard z array

       wFileName <- paste(formatC(
         zMetaInfrm$metaZ[i,]$lower,
         width =
           trunc(log10(zMetaInfrm$metaZ$lower[nrow(zMetaInfrm$metaZ)])) + 1,
         flag = "0"
         ), # formatC
         "wmat", zMetaInfrm$rndCode, ".RData", sep = "")

       save(w, file = paste(zMetaInfrm$tempDir, wFileName, sep = "/"))
       rm(w)

       return(wFileName)
    }, # function FUN

    FUN.VALUE = character(1),
    zMetaInfrm = zMetaInfrm, compFun = compFun, moreArgs = moreArgs
  ) # vapply

  # Transform into color values and plot it
  if(!noScreenDevice) {
    if(verbose) cat("\nTransforming function values into colors ...")
    complexArrayPlot(zMetaInfrm, xlim, ylim, pType, invertFlip,
                     lambda, gamma, pi2Div, logBase,
                     argOffset, stdSaturation, darkestShade, hsvNaN,
                     verbose = verbose, ...)
  } # if(!noScreenDevice)
  else if(verbose) cat("\nNo plot is made (explicit wish of the user) ...")

  # Delete all temporary files ... or not
  if(deleteTempFiles) {
    if(verbose) cat("Deleting temporary files ... ")
    filesToDelete <- paste(zMetaInfrm$tempDir,
                           c(as.character(zMetaInfrm$metaZ$fileNames),
                             as.character(zMetaInfrm$metaZ$wFileNames)),
                           sep = "/")
    unlink(filesToDelete)
    if(verbose) cat("done.\n")
  } else {
    if(verbose)
      cat("\nTemporary files are NOT deleted (explicit wish of the user).\n")
  } # else (temp files ore not deleted)

  # If a parallel backend has been registered, keep or register a sequential
  # backend dependent on user settings
  if(nCores > 1) {
    if(!autoDereg) {
      if(verbose) cat("\nParallel backend with", nCores,
          "cores remains registered for convenience.")
      if(verbose)
        cat("\nCan be de-registered manually with",
            "'foreach::registerDoSEQ()'.\n")
    }
    else {
      foreach::registerDoSEQ()
      if(verbose) cat("\nSequential backend registered again.\n")
    }
  } # if nCores > 1

  invisible(TRUE) # For test purposes
} # phasePortrait

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

