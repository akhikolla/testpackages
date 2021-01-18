#' savePlot.plotSummaryAlphaPart.R
#'
#' Save plot objects on the disk for permanent storage. Function \code{\link[grDevices]{savePlot}} from the \pkg{grDevices} package
#' works for current page on graphical device. This is an attempt to make this
#' function generic so that one can define \code{savePlot} methods for particular
#' needs.
#'
#' @seealso \code{\link[grDevices]{savePlot}} help page on the default \code{savePlot}
#' method in the \pkg{grDevices} package; \code{\link[AlphaPart]{savePlot.plotSummaryAlphaPart}}
#' help page on the method for the objects of \code{plotSummaryAlphaPart} class; and
#' \code{\link[AlphaPart]{plot.summaryAlphaPart}} for ploting results of summaryAlphaPart object.
#'
#' @param x Object on which to chose savePLot method.
#' @param filename Character, filename to save to.
#' @param type  Character, file/device type.
#' @param device Device, the device to save from.
#' @param pre.hook Function, call some code before calling print method for plots (see examples).
#' @param traitsAsDir Logical, should plots be saved within trait folders; the construction is
#' \code{file.path(dirname(file), trait, basename(file))}.
#' folders are created if they do not exist.
#' @param ... Arguments passed to \code{type} specific methods, say
#' \code{width} and \code{height} for \code{type="pdf"} etc.
#'
#' @example /inst/examples/examples_savePlot.plotsummaryAlphaPart.R
#'
#' @return Beside the side effect of saving plots to disk, filenames are printed on
#' screen during the process and at the end invisibly returned.
#'
#' @useDynLib AlphaPart, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom grDevices dev.cur
#' @importFrom grDevices dev.off


#' @export
savePlot  <- function ( ...)    {
    UseMethod("savePlot")
}


#' @describeIn savePlot Default \code{savePlot} method.
savePlot.default <- function(...) {
  
  ##seealso<< \code{\link[grDevices]{savePlot}} help page on the default \code{savePlot}
  ## method in the \pkg{grDevices} package
  
  grDevices::savePlot(...)
  
  ##value<< See \code{\link[grDevices]{savePlot}} for details.

}

#' @describeIn savePlot Save plot objects of class \code{plotSummaryAlphaPart} on the disk for
#' permanent storage.
#'
#' @export


# Save plots on the disk
savePlot.plotSummaryAlphaPart <- function(
  x,                                           ##<< plotSummaryAlphaPart, output object from
                                               ## \code{\link[AlphaPart]{plot.summaryAlphaPart}} function
  filename=paste("Rplot", type, sep="."),      ##<< character, filename to save to
  type=c("pdf", "png", "jpeg", "tiff", "bmp"), ##<< character, file/device type
  device=dev.cur(),                            ##<< device, the device to save from (not used for this method)
  pre.hook=NULL,                               ##<<
  traitsAsDir=FALSE,                           ##<<
  ...                                          ##<<
) {


  if (length(filename) > 1) stop("'filename' argument must be of length one")
  if (!("plotSummaryAlphaPart" %in% class(x))) stop("'x' must be of a 'plotSummaryAlphaPart' class")
  filenameOrig <- sub(pattern=paste(".", type, "$", sep=""), replacement="", x=filename)
  ret <- c()
  lT <- names(x)

  for (i in 1:length(x)) {
    if (traitsAsDir) {
      dir.create(path=file.path(dirname(filenameOrig), lT[i]), recursive=TRUE, showWarnings=FALSE)
      filename <- file.path(dirname(filenameOrig), lT[i], basename(filenameOrig))
    }
    fileA <- paste(filename, paste(lT[i], ".", type, sep=""), sep="_")
    ret  <- c(ret, fileA)
    cat(fileA, "\n")
    print(list(file=fileA, ...))
    do.call(what=type, args=list(file=fileA, ...))
    if (!is.null(pre.hook) && is.function(pre.hook)) pre.hook()
    print(x[[i]])
    dev.off()
  }
  
  ## --- Return ---
  
  invisible(ret)

}

