#' @title Function MAPI_Tails
#' @export
#' 
#' @description Determine significant continuous and discontinuous areas from the result 
#'   of a MAPI analysis when run with permutations.
#' 
#' @param resu A spatial object of class 'sf' resulting from a MAPI analysis done using either
#'    \code{\link{MAPI_RunAuto}} or \code{\link{MAPI_RunOnGrid}}.
#' @param minQ Threshold under which cells with the smallest sum-of-weights percentile (range 1 .. 100) are discarded (default value = 0). 
#'   This parameter allows to discard cells for which the average value of the pairwise metric is computed 
#'   using either a small number and/or only long-distance ellipses.
#' @param alpha Significance level (default=0.05)
#' 
#' @return a spatial object of class 'sf' with the area and geometry of the polygons delineating the significant areas. 
#'    A column provides the tail for each polygon (upper or lower).
#'
#' @details When permutations are performed, in \code{\link{MAPI_RunOnGrid}} for each cell, the proportion of permuted values that are smaller or greater 
#'    than the observed value provides a lower-tailed (ltP) and upper-tailed (utP) test p-value.
#'    A false discovery rate (FDR) procedure (Benjamini and Yekutieli, 2001) is applied to account for multiple 
#'    testing (number of cells) under positive dependency conditions (spatial autocorrelation). An adjusted
#'    p-value is computed for each cell using the function \code{p.adjust} from the 'stats' package with the method 'BY'.
#'    The significance level at which FDR is controlled is set through the parameter alpha. For example, when alpha is
#'    set to 0.05, this means that 5\% of the cells detected as significant can be false positives.
#'    
#'    Significant cells belonging to the lower (or upper) tail that are spatially connected are aggregated 
#'    together to form the significant areas with the lowest (or greater) average values of the pairwise metric analyzed.
#' 
#' @examples
#' \dontrun{
#' data("metric")
#' data("samples")
#' # Run MAPI computation
#' resu <- MAPI_RunAuto(samples, metric, crs=3857, nbPermuts=1000)
#' # Discards the 10% cells with the smallest sum-of-weights 
#' #    and aggregates adjacent cells belonging to the same tail 
#' #    at a 5% significance level
#' tails <- MAPI_Tails(resu, minQ=10, alpha=0.05)
#' }
#'
#' @references 
#' Benjamini, Y. and Yekutieli, D. (2001). The control of the false discovery rate in multiple testing under dependency. Annals of Statistics 29, 1165–1188.


MAPI_Tails <- function(resu, minQ=0, alpha=0.05) {
	
	# check data availability
	if (any(colnames(resu)=="ltP") || any(colnames(resu)=="utP")) {
		message(sprintf("Significant areas aggregation, with percentile filter minQ > %d and alpha = %f ...", minQ, alpha))
	} else {
		stop("Error: no adjusted probabiblites for lower- and upper-tail columns (\"ltP\", \"utP\") found in results table provided.")
	}
	
	# Check confidence level value
	if (alpha < 0.0 || alpha > 0.5) {
		stop(sprintf("Incorrect value for parameter alpha: %f (Range: [0, 0.5])", minQ))
	}
	
	# Allow sum-of-weights percentile filtering
	if (minQ > 0) {
		my.resu <- resu[ resu$swQ > minQ , ]
	} else if (minQ <= 100) {
		my.resu <- resu
	} else {
		stop(sprintf("Incorrect value for parameter minQ: %d (Range: [0, 100])", minQ))
	}
	
	# TODO: check for consistence between alpha and nbPermuts!
	
	# get upper tail
	anyUpper <- any(c(my.resu$utP <= alpha, FALSE), na.rm=TRUE)
	if (anyUpper) { 
		tails.up.g <- sf::st_cast(sf::st_buffer(sf::st_union(my.resu$geometry[my.resu$utP <= alpha]), 0.0001), 'POLYGON')
		tails.up <- data.frame(tail=rep("upper", length(tails.up.g)))
		sf::st_geometry(tails.up) <- tails.up.g
		tails.up$area <- sf::st_area(tails.up)
	}
	# get lower tail
	anyLower <- any(c(my.resu$ltP <= alpha, FALSE), na.rm=TRUE)
	if (anyLower) { 
		tails.low.g <- sf::st_cast(sf::st_buffer(sf::st_union(my.resu$geometry[my.resu$ltP <= alpha]), 0.0001), 'POLYGON')
		tails.low <- data.frame(tail=rep("lower", length(tails.low.g)))
		sf::st_geometry(tails.low) <- tails.low.g
		tails.low$area <- sf::st_area(tails.low)
	}
	# merge and returns tails (if any)
	if (!anyUpper && !anyLower) { 
		# geometries empty (ie. no tails)
		tails <- sf::st_cast(sf::st_buffer(sf::st_union(my.resu$geometry[my.resu$utP <= alpha]), 0.0001), 'POLYGON') # returns an empty geometry
		message("... no significant area")
	} else if (anyUpper && anyLower) { 
		# both geometries exists
		tails <- rbind(tails.up, tails.low)
		message(sprintf("... %d upper-tail and %d lower-tail significant areas returned", nrow(tails.up), nrow(tails.low)))
	} else if (anyUpper) { 
		# only upper tail exists
		tails <- tails.up
		message(sprintf("... %d upper-tail significant areas returned", nrow(tails.up)))
	} else if (anyLower) { 
		# only lower tail exists
		tails <- tails.low
		message(sprintf("... %d lower-tail significant areas returned", nrow(tails.low)))
	}
	# add an id if needed
	if (anyUpper || anyLower) {
		tails$id <- 1:nrow(tails)
	}
	return(tails)
}

