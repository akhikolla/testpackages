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



# -----------------------------------------------------------------------------
# function phaseModColBw

# Calculates a hex two-color array based on an array of complex numbers, both
# arrays handed over as pointers.

# The function takes into account only the modulus values. The modulus of a
# complex number is attributed to a zone based on the input parameter logBase,
# and either assigned the first or the second value of the input variable
# bwCols. Only in cases where the modulus cannot be determined (NaNs or Inf),
# the third color in bwCols is used.

# In more detail, for an input number's modulus, the logarithm with base
# logBase is calculated and cut down to the next lower integer value. If this
# is an even number, the first color of bwCols is taken. In case of an odd
# number, the second color is used.

phaseModColBw <- function(pCompArr,
                          pBwCol,
                          logBase = exp(2*pi/18),
                          bwCols = c("black", "gray95", "gray")) {

  hexCols <- sapply(bwCols,
                    function(bwc) rgb(t(col2rgb(bwc)), maxColorValue = 255))

  dims    <- dim(pCompArr$value)

  intMod  <- floor(log(Mod(pCompArr$value), logBase))

  intIdx  <- intMod %% 2 + 1
  intIdx  <- ifelse(is.nan(intIdx), 3, intIdx)

  pBwCol$value <- array(hexCols[intIdx], dims)

  return(pBwCol)

} # phaseModColBw

# -----------------------------------------------------------------------------
# function phaseAngColBw

# Calculates a hex two-color array based on an array of complex numbers, both
# arrays handed over as pointers.

# The function takes into account only the argument values. The argument of a
# complex number is attributed to a zone based on the input parameters pi2Div
# and argOffset. Then, it is either assigned to the first or the second value
# of the input variable bwCols. Only in cases where the argeument cannot be
# determined (NaNs) the third color in bwCols is used.

# In more detail, the full angle (2*pi) is divided into p2Div zones, which are
# numbered from 0 to pi2Div - 1 with increasing angle. Even and odd zone numbers
# are attributed the first and the second color in bwCols, respectively.
# Usually, the input parameter pi2Div should be an even number in order to avoid
# the first and the last zone having the same color.

phaseAngColBw <- function(pCompArr,
                          pBwCol,
                          pi2Div = 18,
                          argOffset = 0,
                          bwCols = c("black", "gray95", "gray")) {

  hexCols <- sapply(bwCols,
                    function(bwc) rgb(t(col2rgb(bwc)), maxColorValue = 255))

  dims    <- dim(pCompArr$value)
  argmt   <- Arg(pCompArr$value)
  intArg  <- floor(ifelse(argmt - argOffset < 0, argmt + 2*pi, argmt) /
                     (2 * pi / pi2Div))

  intIdx  <- intArg %% 2 + 1
  intIdx  <- ifelse(is.nan(intIdx), 3, intIdx)

  pBwCol$value <- array(hexCols[intIdx], dims)

  return(pBwCol)

} # phaseAngColBw

# -----------------------------------------------------------------------------
# function phaseModAngColBw

# Calculates a hex two-color array based on an array of complex numbers, both
# arrays handed over as pointers.

# The function takes into account the modulus and the argument values and
# colors the resulting grid in a chessboard-like alternation using the first
# and the second color in the input variable bwCols. Only in cases where the
# modulus or the argument cannot be determined (NaNs or Inf), the third color
# in bwCols is used.

# In more detail, for an input number's modulus, the logarithm with base
# logBase is calculated and cut down to the next lower integer value. For the
# argument, the full angle (2*pi) is divided into p2Div zones, which are
# numbered from 0 to pi2Div - 1 with increasing angle. The sum of both integers
# is taken, and if it is an even or an odd number, the first or the second
# color from bwCols is used, respectively.

phaseModAngColBw <- function(pCompArr,
                             pBwCol,
                             pi2Div  = 18,
                             logBase = exp(2*pi/pi2Div),
                             argOffset = 0,
                             bwCols = c("black", "gray95", "gray")) {

  hexCols <- sapply(bwCols,
                    function(bwc) rgb(t(col2rgb(bwc)), maxColorValue = 255))

  dims    <- dim(pCompArr$value)
  argmt   <- Arg(pCompArr$value)
  intArg  <- floor(ifelse(argmt - argOffset < 0, argmt + 2*pi, argmt) /
                     (2 * pi / pi2Div))
  intMod  <- floor(log(Mod(pCompArr$value), logBase))

  intIdx  <- (intArg + intMod) %% 2 + 1
  intIdx  <- ifelse(is.nan(intIdx), 3, intIdx)

  pBwCol$value <- array(hexCols[intIdx], dims)

  return(pBwCol)

} # phaseModAngColBw

# -----------------------------------------------------------------------------
# Function complexArrayPlotBw

# Very much like the function complexArrayPlot, but tailored for two-color
# plots to be created by calling phasePortraitBw.

# Displays an array of complex numbers in an existing plot.
# In order to do so,the temporary files that together form the array are
# read from disk one by one, but each one is processed in a parallel loop.
# The resulting array of hex color values is finally plotted as
# a raster image.

complexArrayPlotBw <- function(zMetaInfrm,
                               xlim,
                               ylim,
                               bwType = "ma",
                               invertFlip = FALSE,
                               pi2Div = 18,
                               logBase = exp(2*pi/pi2Div),
                               argOffset = 0,
                               bwCols = c("black", "grey95", "grey"),
                               asp = 1,
                               xlab = "", ylab = "",
                               verbose,
                               ...) {

  # Set up plot
  plot(NULL, xlim = xlim, ylim = ylim, asp = asp, xlab = xlab, ylab = ylab, ...)

  # Define call to color transformation function depending user's
  # choice of bwType
  colCmd <- switch(bwType,
                   "m"  = "phaseModColBw(pListCompArr[[i]],
                                         pBwCol,
                                         logBase = logBase,
                                         bwCols = bwCols)",

                   "a"  = "phaseAngColBw(pListCompArr[[i]],
                                         pBwCol,
                                         pi2Div = pi2Div,
                                         argOffset = argOffset,
                                         bwCols = bwCols)",

                   "ma" = "phaseModAngColBw(pListCompArr[[i]],
                                            pBwCol,
                                            pi2Div = pi2Div,
                                            logBase = logBase,
                                            argOffset = argOffset,
                                            bwCols = bwCols)"
  ) # switch


  # Obtain the names of the files to load and process
  zMetaInfrm$metaZ$wFileNames <- paste(zMetaInfrm$tempDir,
                                       zMetaInfrm$metaZ$wFileNames, sep = "/")

  # Run the color transformation function over each file
  pBwCol <- lapply(c(1:nrow(zMetaInfrm$metaZ)),
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
      pBwCol <- foreach(i = c(1:length(pListCompArr)),
                        .export  = c("phaseModColBw",
                                     "phaseAngColBw",
                                     "phaseModAngColBw",
                                     "bwCols",
                                     "logBase",
                                     "pi2Div",
                                     "newPointer",
                                     "argOffset"),
                        .combine = c) %dopar% {

        pBwCol  <- newPointer(NULL)
        eval(parse(text = colCmd))      # Does not require a return value,
        # changes color array via pointer
        pListCompArr[[i]]$value <- NULL # Reduced here, but removed after
                                        # the foreach loop
        return(pBwCol)
      } # foreach
      if(verbose) cat("done.")

      # Remove the original list of array pointers
      rm(pListCompArr)

      # Combine the color arrays in the value of the first pointer.
      # Free the others (rbindArraysbyPointer).
      #   Enforce (one-element-) list in case there is only one value
      #   (i.e. if foreach loop was executed sequentially, one core only)
      if(length(pBwCol) == 1) pBwCol <- list(pBwCol)
      pBwCol <- rbindArraysbyPointer(pBwCol)

      return(pBwCol)
    }, # function in lapply
    zMetaInfrm = zMetaInfrm, colCmd = colCmd
  ) # lapply

  # Now combine all blocks into the big raster ...
  if(verbose) cat("\nCombine color rasters ... ")
  pBwCol <- rbindArraysbyPointer(pBwCol)
  if(verbose) cat("done.\n")

  # ... and plot it
  if(verbose) cat("Plotting raster image ... ")
  rasterImage(as.raster(pBwCol$value), xlim[1], ylim[1], xlim[2], ylim[2])
  if(verbose) cat("done.\n")

  pBwCol$value <- NULL
  rm(pBwCol)

  return(NULL)

} # complexArrayPlotBw

# -----------------------------------------------------------------------------
#' Create two-color phase portraits of complex functions
#'
#' \code{phasePortraitBw} allows for creating two-color phase portraits of
#' complex functions based on a polar chessboard grid (cf.
#' \insertCite{wegert_visualcpx_2012;textual}{viscomplexr}, p. 35). Compared to
#' the full phase portraits that can be made with \code{\link{phasePortrait}},
#' two-color portraits omit information. Especially in combination with full
#' phase portraits they can be, however, very helpful tools for interpretation.
#' Besides, two-color phase portraits have a special aesthetic appeal which is
#' worth exploring for itself. In its parameters and its mode of operation,
#' \code{phasePortraitBw} is very similar to \code{\link{phasePortrait}}.
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
#' notation, etc. to poster-like artistic pictures. The general mode of operation,
#' including the usage of parallel processing is exactly the same as with
#' \code{\link{phasePortrait}}, see details section there.
#'
#'
#'
#' @param FUN The function to be visualized. There are two possibilities to
#'   provide it, a quoted character string, or a function object. The quoted
#'   character string must contain an expression that can be interpreted by R as
#'   a function of a complex number \code{z} (like e.g. "sin(z)", "(z^2 -
#'   1i)/(tan(z))", "1/4*z^2 - 10*z/(z^4+4)"). See the documentation of
#'   \code{\link{phasePortrait}} for a complete presentation of all options.
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
#'   Defaults to \code{FALSE}. If this option is chosen, the numbers at the axis
#'   ticks have another meaning than in the normal case. Along the real axis,
#'   they represent the real part of \code{1/z}, and along the imaginary axis,
#'   they represent the imaginary part of \code{1/z}. Thus, if you want
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
#' @param bwType One of the three options for plotting, "m", "a", and "ma", to
#'   be provided as a character string. Defaults to "ma". This parameter has a
#'   comparable role to the parameter \code{pType} in
#'   \code{\link{phasePortrait}}. Option "m" produces a plot that colors modulus
#'   zones only. In more detail, for each input number's modulus, the logarithm
#'   with base \code{logBase} (see below) is calculated and cut down to the next
#'   lower integer value. If this is an even number, the first color given in
#'   \code{bwCols} (see below) is taken. In case of an odd number, the second
#'   color is used. Option "a" produces a plot that exclusively colors argument
#'   (phase angle) zones. To that end, the full angle (2*pi) is divided into
#'   \code{p2Div} (see below) zones, which are numbered from 0 to pi2Div - 1
#'   with increasing angle. Such an integer number is attributed to the complex
#'   number of interest according to the zone it falls into. Even and odd zone
#'   numbers are mapped to the first and the second color in \code{bwCols},
#'   respectively. For normal purposes, the input parameter \code{pi2Div} should
#'   be an even number in order to avoid the first and the last zone having the
#'   same color. With option "ma", a chessboard-like alternation of colors is
#'   displayed over the tiles formed by the intersecting modulus and argument
#'   zones (both determined separately as with the options "m" and "a").
#'
#' @param pi2Div Angle distance for the argument reference zones added for
#'   \code{pType = "pma"} or \code{pType = "pa"}. The value has to be given as
#'   an integer (reasonably) fraction of 2*pi (i.e. 360 degrees). Unlike with
#'   \code{\link{phasePortrait}}, the default is 18; thus, reference zones are
#'   delineated by default in distances of 2*pi/18, i.e. (20 degrees), starting
#'   with 0 if not defined otherwise with the parameter \code{argOffset}. While
#'   the default of \code{pi2Div} is 9 with \code{\link{phasePortrait}} for good
#'   reasons (see there), setting \code{pi2Div} to an odd number is usually not
#'   a good choice with two-color phase portraits, because the first and the
#'   last phase angle zone would get the same color. However, as \code{pi2Div}
#'   here defaults to double the value as with \code{\link{phasePortrait}}, both
#'   plot types can be nicely compared even when using their specific defaults
#'   of \code{pi2Div}.
#'
#' @param logBase Modulus ratio between the edges of the modulus zones in
#'   \code{bwType} \code{"m"} and \code{"ma"}. As recommended by
#'   \insertCite{wegert_visualcpx_2012;textual}{viscomplexr}, the default
#'   setting is \code{logBase = exp(2*pi/pi2Div)}. This relation between the
#'   parameters \code{logBase} and \code{pi2Div} ensures an analogue scaling of
#'   the modulus and argument reference zones (see Details section in the
#'   documentation of \code{\link{phasePortrait}}). Conveniently, for the
#'   default \code{pi2Div = 18}, we obtain \code{logBase == 1.4177...}, which is
#'   very close to the square root of 2. Thus, when crossing two modulus zones,
#'   the modulus at the higher edge of the second zone is almost exactly two
#'   times the value at the lower edge of the first zone.
#'
#' @param argOffset The (complex number) argument in radians counterclockwise,
#'   at which the argument (phase angle) reference zones are fixed, i.e. the
#'   lower angle of the first zone. Default is 0.
#'
#' @param bwCols Color definition for the plot provided as a character vector of
#'   length 3. Each element of the vector must be either a color name R
#'   recognizes, or a hexadecimal color string like e.g. "#00FF11". The first
#'   and the second color make the appearance of two-color phase portraits (see
#'   \code{bwType} above for details), while the third color is reserved for
#'   special cases, where the input value cannot sufficiently evaluated (NaNs,
#'   partly Inf). Defaults to c("black", "gray95", "gray"), which leads to an
#'   alternation of black and very light gray zones or tiles, and uses a neutral
#'   gray in special cases.
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
#'   \code{phasePortraitBw} will use a 1 x 1 inch pseudo graphics device in this
#'   case. The default for this parameter is \code{FALSE}, and you should change
#'   it only if you really know what you are doing.
#'
#' @param autoDereg if TRUE, automatically register sequential backend after the
#'   plot is completed. Default is FALSE, because registering a parallel backend
#'   can be time consuming. Thus, if you want make several phase portraits in
#'   succession, you should set \code{autoDereg} only for the last one, or
#'   simply type \code{foreach::registerDoSEQ} after you are done. In any case,
#'   you don't want to have an unused parallel backend lying about.
#'
#' @param verbose if TRUE (default), \code{phasePortraitBw} will continuously
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
#' @export
#'
#'
#' @examples
#' # Map the complex plane on itself
#'
#' # x11(width = 8, height = 8)     # Screen device commented out
#'                                  # due to CRAN test requirements.
#'                                  # Use it when trying this example
#' phasePortraitBw("z", xlim = c(-2, 2), ylim = c(-2, 2),
#'                 xlab = "real", ylab = "imaginary",
#'                 verbose = FALSE, # Suppress progress messages
#'                 nCores = 2)      # Max. two cores allowed on CRAN
#'                                  # not a limit for your own use
#' \dontshow{
#' # R CMD check: make sure any open connections are closed afterward
#' foreach::registerDoSEQ()
#' doParallel::stopImplicitCluster()
#' }
#'
#'
#'
#' # Sinus with default colors and default bwType ("ma")
#' \donttest{
#' # x11(width = 8, height = 8)       # Screen device commented out
#'                                    # due to CRAN test requirements.
#'                                    # Use it when trying this example
#' phasePortraitBw("sin(z)",
#'                 xlim = c(-pi, pi),
#'                 ylim = c(-pi, pi),
#'                 verbose = FALSE,
#'                 nCores = 2)        # Max. two cores allowed on CRAN
#'                                    # not a limit for your own use
#'   \dontshow{
#'   # R CMD check: make sure any open connections are closed afterward
#'   foreach::registerDoSEQ()
#'   doParallel::stopImplicitCluster()
#'   }
#' }
#'
#'
#'
#' # Sinus with custom colors and bwType "a"
#' \donttest{
#' # x11(width = 8, height = 8)       # Screen device commented out
#'                                    # due to CRAN test requirements.
#'                                    # Use it when trying this example
#' phasePortraitBw("sin(z)",
#'                 xlim = c(-pi, pi),
#'                 ylim = c(-pi, pi),
#'                 bwType = "a",
#'                 bwCols = c("darkgreen", "green", "gray"),
#'                 verbose = FALSE,
#'                 nCores = 2)        # Max. two cores allowed on CRAN
#'                                    # not a limit for your own use
#'   \dontshow{
#'   # R CMD check: make sure any open connections are closed afterward
#'   foreach::registerDoSEQ()
#'   doParallel::stopImplicitCluster()
#'   }
#' }
#'
#'
#'
#' # Sinus with custom colors and bwType "m"
#' \donttest{
#' # x11(width = 8, height = 8)       # Screen device commented out
#'                                    # due to CRAN test requirements.
#'                                    # Use it when trying this example
#' phasePortraitBw("sin(z)",
#'                 xlim = c(-pi, pi),
#'                 ylim = c(-pi, pi),
#'                 bwType = "m",
#'                 bwCols = c("darkblue", "skyblue", "gray"),
#'                 verbose = FALSE,
#'                 nCores = 2)        # Max. two cores allowed on CRAN
#'                                    # not a limit for your own use
#'   \dontshow{
#'   # R CMD check: make sure any open connections are closed afterward
#'   foreach::registerDoSEQ()
#'   doParallel::stopImplicitCluster()
#'   }
#' }
#'
#'
#'
#' # Map the complex plane on itself, show all bwType options
#' \donttest{
#' # x11(width = 8, height = 8)       # Screen device commented out
#'                                    # due to CRAN test requirements.
#'                                    # Use it when trying this example
#' op <- par(mfrow = c(2, 2), mar = c(4.1, 4.1, 1.1, 1.1))
#' for(bwType in c("ma", "a", "m")) {
#'   phasePortraitBw("z", xlim = c(-2, 2), ylim = c(-2, 2),
#'                   bwType = bwType,
#'                   xlab = "real", ylab = "imaginary",
#'                   verbose = FALSE, # Suppress progress messages
#'                   nCores = 2)      # Max. two cores allowed on CRAN
#'                                    # not a limit for your own use
#' }
#' # Add normal phase portrait for comparison
#' phasePortrait("z", xlim = c(-2, 2), ylim = c(-2, 2),
#'               xlab = "real", ylab = "imaginary",
#'               verbose = FALSE,
#'               pi2Div = 18,         # Use same angular division as default
#'                                    # in phasePortraitBw
#'               nCores = 2)
#' par(op)
#'
#'   \dontshow{
#'   # R CMD check: make sure any open connections are closed afterward
#'   foreach::registerDoSEQ()
#'   doParallel::stopImplicitCluster()
#'   }
#' }
#'
#'
#' # A rational function, show all bwType options
#' \donttest{
#' # x11(width = 8, height = 8)       # Screen device commented out
#'                                    # due to CRAN test requirements.
#'                                    # Use it when trying this example
#' funString <- "(z + 1.4i - 1.4)^2/(z^3 + 2)"
#' op <- par(mfrow = c(2, 2), mar = c(4.1, 4.1, 1.1, 1.1))
#' for(bwType in c("ma", "a", "m")) {
#'   phasePortraitBw(funString, xlim = c(-2, 2), ylim = c(-2, 2),
#'                   bwType = bwType,
#'                   xlab = "real", ylab = "imaginary",
#'                   verbose = FALSE, # Suppress progress messages
#'                   nCores = 2)      # Max. two cores allowed on CRAN
#'                                    # not a limit for your own use
#' }
#' # Add normal phase portrait for comparison
#' phasePortrait(funString, xlim = c(-2, 2), ylim = c(-2, 2),
#'               xlab = "real", ylab = "imaginary",
#'               verbose = FALSE,
#'               pi2Div = 18,         # Use same angular division as default
#'                                    # in phasePortraitBw
#'               nCores = 2)
#' par(op)
#'
#'   \dontshow{
#'   # R CMD check: make sure any open connections are closed afterward
#'   foreach::registerDoSEQ()
#'   doParallel::stopImplicitCluster()
#'   }
#' }
#'
#'
phasePortraitBw <- function(FUN, moreArgs = NULL, xlim, ylim,
                            invertFlip = FALSE,
                            res = 150,
                            blockSizePx = 2250000,
                            tempDir = NULL,
                            nCores = parallel::detectCores(),
                            bwType = "ma",
                            pi2Div = 18,
                            logBase = exp(2*pi/pi2Div),
                            argOffset = 0,
                            bwCols = c("black", "gray95", "gray"),
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
    function(i, zMetaInfrm, compFun, moreArgs) {

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
    complexArrayPlotBw(zMetaInfrm, xlim, ylim, bwType, invertFlip,
                       pi2Div, logBase, argOffset, bwCols,
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
} # phasePortraitBw

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------


