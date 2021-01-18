#' Deprecated LocalControl functions
#'
#' These functions are provided for compatibility with previous versions of LocalControl.
#' They may eventually be completely removed.
#' @rdname LocalControl-deprecated
#' @name LocalControl-deprecated
#' @docType package
#' @export  localControlNearestNeighbors localControlCompetingRisks plotLocalControlCIF plotLocalControlLTD
#' @aliases localControlNearestNeighbors localControlCompetingRisks plotLocalControlCIF plotLocalControlLTD
#' @usage NULL
#' @section Details:
#' \tabular{rl}{
#'   \code{localControlNearestNeighbors} \tab Now called using \code{\link{LocalControl}} with the outcomeType = "cross-sectional". \cr
#'   \code{localControlCompetingRisks}   \tab Now called using \code{\link{LocalControl}} with the outcomeType = "survival". \cr
#'   \code{plotLocalControlCIF}   \tab Now called using \code{\link{plot.LocalControlCR}}. \cr
#'   \code{plotLocalControlLTD}   \tab Now called using \code{\link{plot.LocalControlCS}}. \cr
#' }
#'
localControlNearestNeighbors <-function(data,
                                        treatmentColName,
                                        treatmentCode="",
                                        outcomeColName,
                                        clusterVars,
                                        labelColName="",
                                        numThreads=1,
                                        radiusLevels = numeric(),
                                        radStepType = "exp",
                                        radDecayRate = 0.8,
                                        radMinFract = 0.01,
                                        normalize = TRUE, verbose = FALSE)
{
  .Deprecated(old = "localControlNearestNeighbors", new = "LocalControl", package = "LocalControl")
  x = LocalControl( data = data, treatmentColName = treatmentColName, treatmentCode = treatmentCode,
                    outcomeColName = outcomeColName, clusterVars = clusterVars,
                    labelColName = labelColName, numThreads = numThreads,
                    radiusLevels = radiusLevels, radStepType = radStepType,
                    radDecayRate = radDecayRate, radMinFract = radMinFract,
                    normalize = normalize, verbose = verbose, outcomeType = "cross-sectional")
  x$summary = summary(x)
  x
}
localControlCompetingRisks <- function(data,
                                       outcomeColName,
                                       timeColName,
                                       treatmentColName,
                                       clusterVars,
                                       cenCode = 0,
                                       treatmentCode="",
                                       labelColName="",
                                       radStepType = "exp",
                                       radDecayRate = 0.8,
                                       radMinFract = 0.01,
                                       radiusLevels = numeric(),
                                       normalize = TRUE,
                                       verbose = FALSE,
                                       numThreads=1)
{
  .Deprecated(old = "localControlCompetingRisks", new = "LocalControl", package = "LocalControl")
  x = LocalControl( data = data, treatmentColName = treatmentColName, treatmentCode = treatmentCode,
                    outcomeColName = outcomeColName, clusterVars = clusterVars,
                    labelColName = labelColName, numThreads = numThreads,
                    radiusLevels = radiusLevels, radStepType = radStepType, cenCode = cenCode,
                    radDecayRate = radDecayRate, radMinFract = radMinFract, timeColName = timeColName,
                    normalize = normalize, verbose = verbose, outcomeType = "survival")
  x$summary = summary(x)
  x
}
plotLocalControlCIF = function(lccrResults, rad2plot, xlim, ylim = c(0,1),
                               col1 = "blue", col0 = "red",
                               xlab = "Time", ylab = "Cumulative incidence",
                               legendLocation = "topleft", main = "",
                               group1 = "Treatment 1", group0 = "Treatment 0", ...){
  .Deprecated(old = "plotLocalControlCIF", new = "plot.LocalControlCR", package = "LocalControl")
  plot(lccrResults, rad2plot = rad2plot, xlim = xlim, ylim = ylim,
       col1 = col1, col0 = col0, xlab = xlab, ylab = ylab,
       legendLocation = legendLocation, main = main, group1 = group1, group0 = group0, ...)
}
plotLocalControlLTD = function(lcnnResults, nnConfidence, ylim,
                               legendLocation = "bottomleft", ylab = "LTD",
                               xlab = "Fraction of maximum radius", main =""){
  .Deprecated(old = "plotLocalControlLTD", new = "plot.LocalControlCS", package = "LocalControl")
  plot(lcnnResults, nnConfidence = nnConfidence, ylim = ylim,
       legendLocation = legendLocation, ylab = ylab,  xlab = xlab, main = main)
}
NULL
