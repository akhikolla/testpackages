# undocumented.
# stores a predictor independent from the C++ library
# Ok for OKrig but not for other objects
# possible nondefault args are predictor=fitobject, CovFnParam=GCVblob$CovFnParam
savePredictor <- function(file,predictor=blackbox.getOption("fitobject"),
                          CovFnParam=blackbox.getOption("CovFnParam")) {
  if (missing(file)) {
    stop.redef("'file' argument missing, with no default")
  }
  if (inherits(predictor,"OKriglistplus")) {
    stop("savePredictor() does not yet produce reusable predictors from 'OKriglistplus' objects." )
  }
  predictor$CKrigidx <- NULL
  predictor$unique_x <- as.matrix(unique(predictor$x))
  predictor$covTheta2 <- CovFnParam[colnames(predictor$x)]^2
  predictor$smoothness <- CovFnParam["smoothness"]
  predictor$x <- NULL
  predictor$y <- NULL
  predictor$fitted.values <- NULL
  #
  unlink(file) ## this *deletes* the file if it exists
  save(  predictor = predictor , file=file)
  msgg <- paste("\"", file, "\"", sep="")
  msg <- paste("Use load(", msgg, ") to load 'predictor' in a future R session", sep="")
  cat(msg) ## with print() all quotation mark characters are expressed as "
  invisible(predictor)
}
