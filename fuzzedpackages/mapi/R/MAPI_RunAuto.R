#' @title Function MAPI_RunAuto
#' @export
#' 
#' @description This function is a wrapper allowing to run a complete MAPI analysis. 
#' 
#' @param samples a data.frame with names and geographical coordinates of samples. 
#'   Column names must be: 'ind', 'x', 'y'.  Optional column 'errRad' with an error radius for sample locations (eg. GPS uncertainty). 
#'   Coordinates must be projected (not latitude/longitude).
#' @param metric a data.frame or a square matrix with the pairwise metric computed for all pairs of samples. 
#'   If data.frame, column names must be: 'ind1', 'ind2', 'value'.  
#'   If matrix, sample names must be the row- and column names.
#' @param crs coordinate reference system: integer with the EPSG code, or character with proj4string. 
#'   When using dummy coordinates (eg. simulation output) you may use EPSG:3857 for example. 
#'   This allows computation but, of course, has no geographical meaning.
#' @param isMatrix Boolean. Depends on the 'metric' data:\cr
#'   TRUE if 'metric' is a square matrix with column names = row names and standing for sample names.\cr
#'   FALSE if 'metric is a three columns data.frame ('ind1', 'ind2', 'value'). \cr
#'   The default value is determined using a "matrix" class detection for 'metric' as well as identity between row and column number.
#' @param beta A value depending on spatial regularity of sampling: 0.5 for regular sampling, 0.25 for random sampling (Hengl, 2006).
#' @param ecc ellipse eccentricity value (0.975 by default).
#' @param buf optional. This parameter allows to expand or shrink the grid by a number of units in 
#'   the same reference system as the sample geographical coordinates (0 by default).
#' @param errRad global error radius for sample locations (same radius for all samples, 10 by default). 
#'   Units are in the same reference system as the sample geographical coordinates.
#'   To use different error radius values for sample locations, add a column 'errRad' in the 'sample' data (see \code{\link{mapi}}).
#' @param nbPermuts number of permutations of sample locations (0 by default).
#' @param dMin minimum distance between individuals. 0 by default.
#' @param dMax maximal distance between individuals. +Inf by default.
#' @param nbCores number of CPU cores you want to use during parallel computation. 
#'   The default value is estimated as the number of available cores minus 1, suitable for a personal computer. 
#'   On a cluster you might have to set it to a reasonable value (eg. 8) in order to keep resources for other tasks. 
#' @param N number of points used per quarter of ellipse, 8 by default. 
#'   Don't change it unless you really know what you are doing.
#' 
#' @return a spatial object of class 'sf' providing for each cell: \cr
#' - gid: Cell ID \cr
#' - x and y coordinates of cell center \cr
#' - nb_ell: number of ellipses used to compute the weighted mean \cr
#' - avg_value: weighted mean of the pairwise metric \cr
#' - sum_wgts: sum of weights of ellipses used to compute the weighted mean \cr
#' - w_stdev: weighted standard deviation of the pairwise metric \cr
#' - swQ: percentile of the sum of weights \cr
#' - geometry \cr
#' When permutations are performed: \cr
#' - permuts: list of the weighted mean values obtained from all permutations \cr
#' - proba: proportion of the permuted weighted means below the observed weighted mean \cr
#' - ltP: lower-tail p-value adjusted using the FDR procedure of Benjamini and Yekutieli \cr
#' - utP: upper-tail p-value adjusted using the FDR procedure of Benjamini and Yekutieli \cr
#'
#' @details
#' Following functions are called by \code{MAPI_RunAuto} in following order:
#' - \code{\link{MAPI_CheckData}} cleans the dataset;
#' - \code{\link{MAPI_GridAuto}} generates a grid of hexagons by calling \code{\link{MAPI_EstimateHalfwidth}} then \code{\link{MAPI_GridHexagonal}};
#' - \code{\link{MAPI_RunOnGrid}} performs the MAPI analysis.
#' 
#' NOTE: The call to \code{\link{MAPI_Tails}} is not included.
#' It should be done afterwards on the object returned by \code{MAPI_RunAuto}.
#'
#' @references
#' Benjamini, Y. & Yekutieli, D. (2001) The control of the false discovery rate in multiple testing under dependency. Annals of Statistics, 29, 1165–1188.
#' 
#' @examples
#' \dontrun{
#' data("metric")
#' data("samples")
#' # Run a MAPI analysis without permutation
#' my.results <- MAPI_RunAuto(samples, metric, crs=3857, beta=0.5, nbPermuts=0)
#' 
#' # eg. Export results to shapefile "myFirstMapiResult" in current directory
#' # to further visualize and customize the MAPI plot in SIG software.
#' library(sf)
#' st_write(my.results, dsn=".", layer="myFirstMapiResult", driver="ESRI Shapefile")
#' }
#' 

MAPI_RunAuto <- function(samples, metric, crs, isMatrix=all(class(metric)=="matrix", nrow(metric)==ncol(metric)), beta=0.25, ecc=0.975, buf=0, errRad=10, nbPermuts=0, dMin=0, dMax=Inf, nbCores=ifelse(requireNamespace("parallel", quietly=TRUE), parallel::detectCores()-1, 1), N=8) {
	data <- MAPI_CheckData(samples, metric, isMatrix=isMatrix)
	my.samples <- data[[1]]
	my.metric <- data[[2]]
	grid <- MAPI_GridAuto(samples=my.samples, crs=crs, beta=beta, buf=buf)
	resu <- MAPI_RunOnGrid(samples=my.samples, metric=my.metric, grid=grid, ecc=ecc, errRad=errRad, nbPermuts=nbPermuts, dMin=dMin, dMax=dMax, nbCores=nbCores, N=N)
	return(resu)
}
