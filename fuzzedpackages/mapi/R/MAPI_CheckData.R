#' @title Function MAPI_CheckData
#' @export
#' @description Check the validity of the 'samples' and 'metric' data loaded.\cr
#' Missing data are removed from 'metric', samples with missing coordinates are removed and samples that are not present in both dataset ('samples' and 'metric') are discarded.
#' 
#' @param samples a data.frame with names and geographical coordinates of samples. Column names must be: 'ind', 'x', 'y'.  
#'   Optional column 'errRad' with an error radius for sample locations (eg. GPS uncertainty). Coordinates must be projected (not latitude/longitude).
#' @param metric a data.frame or a square matrix with the pairwise metric computed for all pairs of samples. 
#'   If data.frame, column names must be: 'ind1', 'ind2', 'value'.  If matrix, sample names must be the row- and column names.
#' @param isMatrix Boolean. Depends on the 'metric' data:\cr
#'   TRUE if 'metric' is a square matrix with column names = row names and standing for sample names.\cr
#'   FALSE if 'metric is a three columns data.frame ('ind1', 'ind2', 'value'). \cr
#'   The default value is determined using a "matrix" class detection for 'metric' as well as identity between row and column number.
#' 
#' @return a list of two data.table objects corresponding to 'samples' and 'metric' after cleaning.
#' 
#' @examples
#' \dontrun{
#' data("samples")
#' data("metric")
#' # remove first sample in order to force warning
#' samples <- samples[-c(1) , ]
#' clean.list <- MAPI_CheckData(samples, metric)
#' checked.samples <- clean.list[['samples']]
#' checked.metric <- clean.list[['metric']]
#' }
#' 

MAPI_CheckData <- function(samples, metric, isMatrix=all((class(metric)=="matrix"), (nrow(metric)==ncol(metric)))) {
  
	message("Checking data...")
	
	my.samples <- data.table::as.data.table(samples)
	my.samples$ind <- as.character(my.samples$ind) # Avoids factors...
	data.table::setkey(my.samples, "ind")
	
	if (isMatrix) {
		if (class(metric) != "matrix") {stop("'isMatrix' is true, so 'metric' parameter should be a matrix class with colnames and rownames")}
		if (sum( colnames(metric) != rownames(metric) ) > 0) {
			stop("Row names and column names are different") 
		} else {
			message("`metric` matrix symmetrical ... OK")
		}
		my.metric <- data.table::as.data.table(metric)
		my.metric$tmpId <- rownames(metric)
		my.metric <- data.table::melt(my.metric, id.vars=c("tmpId"), variable.factor=FALSE, value.factor=FALSE)
	} else {
		message("`metric` provided as a vector")
		my.metric <- data.table::as.data.table(metric)
	}
	colnames(my.metric) <- c("ind1", "ind2", "value")
	my.metric <- my.metric[ !is.na(my.metric$value) , ]
	
	my.metric$ind1 <- as.character(my.metric$ind1) # Avoids factors...
	my.metric$ind2 <- as.character(my.metric$ind2) # Avoids factors...
	my.metric$value <- as.numeric(as.character(my.metric$value)) # Avoids factors...
	
	if (anyDuplicated(my.samples$ind)) {
		stop("Names of individuals are not unique.") 
	} else {
		message("Unique names ... OK")
	}
	
	my.samples <- my.samples[!(is.na(my.samples$x)|is.na(my.samples$y)) , ]
	tri <- as.vector(my.samples[ !(my.samples$ind %in% my.metric$ind1) & !(my.samples$ind %in% my.metric$ind2) ,  "ind"]$ind)
	my.samples <- my.samples[ (my.samples$ind %in% my.metric$ind1) | (my.samples$ind %in% my.metric$ind2) , ]
	if (length(tri) != 0) {
		warning(sprintf("Following samples have no metric and have been removed: %s", paste(tri, collapse=", ")))
	} else {
		message("All samples used in metric ... OK")
	}

	tri_ind1 <- my.metric[ !(my.metric$ind1 %in% my.samples$ind) , ]$ind1
	tri_ind2 <- my.metric[ !(my.metric$ind2 %in% my.samples$ind) , ]$ind2
	tri <- c(tri_ind1, tri_ind2)
	tri <- tri[ !duplicated(tri) ]
	if (length(tri) != 0) {
		my.metric <- my.metric[ (my.metric$ind1 %in% my.samples$ind)&(my.metric$ind2 %in% my.samples$ind) , ] 
		warning(sprintf("Following samples have been removed from metric: %s", paste(tri, collapse=", ")))
	}

	# in case of half-matrix, make a full one
	Mr <- my.metric[ my.metric$ind1 != my.metric$ind2 , c("ind2", "ind1", "value")]
	colnames(Mr) <- c("ind1", "ind2", "value")
	my.metric <- unique(rbind(my.metric, Mr))
	data.table::setkey(my.metric, "ind1", "ind2")
	my.metric <- my.metric[order(as.character(my.metric$ind1),as.character(my.metric$ind2)),]

	message(paste("...",  nrow(my.samples), "samples and",  nrow(my.metric), "metric left."))
	
	return( list(samples=my.samples, metric=my.metric) )
}
