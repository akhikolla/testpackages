#' @title Function MAPI_RunOnGrid
#' @export
#' @description Launch a MAPI analysis for a given grid computed with \code{\link{MAPI_GridAuto}} 
#'   or \code{\link{MAPI_GridHexagonal}} or provided by users.
#' 
#' @param samples a data.frame with names and geographical coordinates of samples. 
#'   Column names must be: 'ind', 'x', 'y'.  Optional column 'errRad' with an error radius for sample locations (eg. GPS uncertainty). 
#'   Coordinates must be projected (not latitude/longitude).
#' @param metric a data.frame or a square matrix with the pairwise metric computed for all pairs of samples. 
#'   If data.frame, column names must be: 'ind1', 'ind2', 'value'.  
#'   If matrix, sample names must be the row- and column names.
#' @param grid a spatial object of class 'sf' with the geometry of each cell. 
#'   When using your own grid, please check that the object structure is the same as returned by 
#'   \code{\link{MAPI_GridAuto}} or \code{\link{MAPI_GridHexagonal}}.
#' @param isMatrix Boolean. Depends on the 'metric' data:\cr
#'   TRUE if 'metric' is a square matrix with column names = row names and standing for sample names.\cr
#'   FALSE if 'metric is a three columns data.frame ('ind1', 'ind2', 'value'). \cr
#'   The default value is determined using a "matrix" class detection for metric as well as identity between row and column number.
#' @param ecc ellipse eccentricity value (0.975 by default).
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
#'  - gid: Cell ID \cr
#'  - x and y coordinates of cell center \cr
#'  - nb_ell: number of ellipses used to compute the weighted mean \cr
#'  - avg_value: weighted mean of the pairwise metric \cr
#'  - sum_wgts: sum of weights of ellipses used to compute the weighted mean \cr
#'  - w_stdev: weighted standard deviation of the pairwise metric \cr
#'  - swQ: percentile of the sum of weights \cr
#'  - geometry \cr
#'  When permutations are performed: \cr
#'  - permuts: list of the weighted mean values obtained from all permutations \cr
#'  - proba: proportion of the permuted weighted means below the observed weighted mean \cr
#'  - ltP: lower-tail p-value adjusted using the FDR procedure of Benjamini and Yekutieli \cr
#'  - utP: upper-tail p-value adjusted using the FDR procedure of Benjamini and Yekutieli \cr
#' 
#' @details
#' To test whether the pairwise metric values associated with the ellipses are independent of the sample locations, those are permuted 'nbPermuts' times. 
#'   At each permutation, new cell values are computed and stored to build a cumulative null distribution for each cell of the grid. 
#'   Each cell value from the observed data set is then ranked against its null distribution.
#'   For each cell, the proportion of permuted values that are smaller or greater 
#'   than the observed value provides a lower-tailed (ltP) and upper-tailed (utP) test p-value.
#'   
#' A false discovery rate (FDR) procedure (Benjamini and Yekutieli, 2001) is applied to account for multiple 
#'    testing (number of cells) under positive dependency conditions (spatial autocorrelation).  An adjusted
#'    p-value is computed for each cell using the function \code{p.adjust} from the 'stats' package with the method 'BY'.
#'    
#' @references 
#' Benjamini, Y. and Yekutieli, D. (2001). The control of the false discovery rate in multiple testing under dependency. Annals of Statistics 29, 1165–1188.
#' 
#' @examples
#' \dontrun{
#' data(metric)
#' data(samples)
#' my.grid <- MAPI_GridHexagonal(samples, crs=3857, 500) # 500m halfwidth
#'
#' # Note: 10 permutations is only for test purpose, increase to >=1000 in real life!
#' my.results <- MAPI_RunOnGrid(samples, metric, grid=my.grid, nbPermuts=10, nbCores=1)
#' 
#' # eg. Export results to shapefile "myFirstMapiResult" in current directory
#' library(sf)
#' st_write(my.results, dsn=".", layer="myFirstMapiResult", driver="ESRI Shapefile")
#' }
#' 

MAPI_RunOnGrid <- function(samples, metric, grid, isMatrix=FALSE, ecc=0.975, errRad=10, nbPermuts=0, dMin=0, dMax=Inf, nbCores=ifelse(base::requireNamespace("parallel", quietly=TRUE), parallel::detectCores()-1, 1), N=8) {
	message("MAPI COMPUTATION STARTED")
	tot <- system.time({
		data <- MAPI_CheckData(samples, metric, isMatrix=isMatrix)
		my.samples <- data[[1]]
		my.metric <- data[[2]]
		# prepare index
		data.table::setkey(my.metric, "ind1", "ind2")
		
		# fill missing errRad with parameter value
		if (any(colnames(my.samples)=="errRad")) {
			my.samples[is.na(my.samples$errRad), "errRad"] <- errRad
		} else {
			my.samples$errRad <- errRad
		}
		# get crs value from grid geometry
		crs <- sf::st_crs(grid)[[1]]
		# add geometry (point) to samples
		my.samples.geom <- sf::st_as_sf(my.samples, coords=c("x", "y"), crs=crs, remove=FALSE)
		
		## Compute ellipses
		t <- system.time({
			message("Building ellipse polygons...")
			# create locality code from hex representation of geometry and error radius
			my.samples.geom$locCode <- as.character(sf::st_as_binary(my.samples.geom$geometry, hex=TRUE))
			my.samples.geom$locCode <- apply(my.samples.geom, 1, function(r) { paste(r["locCode"], r["errRad"]) })
			
			# get distinct localities with geometry
			distinct.locations <- unique(as.data.frame(my.samples.geom[, c("locCode", "geometry", "errRad", "x", "y")]))
			
			# prepare pairs of points as cartesian product of distinct localities ...
			ellipses <- data.table::as.data.table(expand.grid(loc1=distinct.locations$locCode, loc2=distinct.locations$locCode))
			# ... and keep only the half-matrix, including diagonal as two distinct individuals may share a same locality
			ellipses <- ellipses[as.character(ellipses$loc1) <= as.character(ellipses$loc2) , ]
			
			# merge half-matrix with locations
			ellipses <- merge(ellipses, distinct.locations, by.x="loc2", by.y="locCode")
			ellipses <- merge(ellipses, distinct.locations, by.x="loc1", by.y="locCode", suffixes=c("2","1"))
			
			# computes distance between the two localities (Pythagore)
			ellipses$dist <- sqrt( (ellipses$x1 - ellipses$x2)^2 + (ellipses$y1 - ellipses$y2)^2 )
			
			# distance filtering out (if any filter set)
			if (!is.na(dMin) && dMin > 0)                       { ellipses <- ellipses[ellipses$dist >= dMin, ] }
			if (!is.infinite(dMax) && !is.na(dMax) && dMax > 0) { ellipses <- ellipses[ellipses$dist <= dMax, ] }
			
			# sort the completed half-matrix [NOTE: why?] and add rowid
			ellipses <- ellipses[order(ellipses$loc1, ellipses$loc2) , ]
			ellipses$rowid <- 1:nrow(ellipses)

			# récupération des localisations distinctes : coordonnées x,y et rayon d'erreur des 2 points formant l'ellipse in a numeric matrix
			ellipses.matrix <- as.matrix(ellipses[,c("x1", "y1", "x2", "y2", "errRad1", "errRad2")])
			
			# call to c++ function mkP4st_cpp for building the polygon (convex hull of ellipse and error circles)
			ePolys <- apply(ellipses.matrix, 1, function(r) {
				sf::st_polygon(list(mkP4st_cpp(r, N, ecc)))
			})
			
			# simplify ellipses table
			ellipses <- ellipses[,c("rowid", "loc1", "loc2", "dist")]
			# set geometry to ellipses table (same order !)
			sf::st_geometry(ellipses) <- sf::st_sfc(ePolys, crs=crs)
			
			# computes area and weight
			ellipses$area <- sf::st_area(ellipses$geometry) # area in "units"
			ellipses$weight <- 1.0 / as.numeric(ellipses$area)
			
			# free memory
			rm(ellipses.matrix)
			
			# For each ellipse we keep a simplified table with its rowid, weight and the two locality codes
			ells <- data.table::data.table(loc1=ellipses$loc1, loc2=ellipses$loc2, rowid=ellipses$rowid, weight=ellipses$weight)
			data.table::setkey(ells, "rowid")
		})
		tv <- as.vector(t) ; tv[is.na(tv)] <- 0.0
		# NOTE: %s used for stringified integers due to overflow for very large datasets (thanks to Simon Dellicour)
		message(sprintf("... %s ellipses.    [user: %0.3f, system: %0.3f, elapsed: %0.3f seconds]", as.character(nrow(ellipses)), tv[1]+tv[4], tv[2]+tv[5], tv[3]))
		
		
		## Spatial intersect between ellipses and grid
		message("Computing spatial intersection between grid cells and ellipses...")
		t <- system.time({
			inter0 <- sf::st_intersects(grid, ellipses, sparse=TRUE, prepared=TRUE)
		})
		tv <- as.vector(t) ; tv[is.na(tv)] <- 0.0
		message(sprintf("... done.  [user: %0.3f, system: %0.3f, elapsed: %0.3f seconds]", tv[1]+tv[4], tv[2]+tv[5], tv[3]))
		# We are finished with ellipses geometry, let's free some memory
		rm(ellipses) # free memory
		
		
		## Build sample pairs
		message("Building sample pairs...")
		t <- system.time({
			# build a simple table with sample code and locality code only
			my.sampleCode.locCode <- data.table::data.table(ind=my.samples.geom$ind, locCode=my.samples.geom$locCode)
			# cartesian product of samples
			sampX <- data.table::as.data.table(expand.grid(ind1=my.sampleCode.locCode$ind,ind2=my.sampleCode.locCode$ind))
			# complete matrix *without* diagonal
			sampX <- sampX[as.character(sampX$ind1) != as.character(sampX$ind2) , ]
			sampX <- sampX[order(sampX$ind1, sampX$ind2) , ]
			# merge for getting locality codes for both samples
			sampX <- merge(sampX, my.sampleCode.locCode, by.x="ind2", by.y="ind")
			sampX <- merge(sampX, my.sampleCode.locCode, by.x="ind1", by.y="ind", suffixes=c("2","1"))
			# merge with ellipses for weight
			sampX_12 <- merge(sampX, ells, by.x=c("locCode1", "locCode2"), by.y=c("loc1", "loc2"))
			sampX_21 <- merge(sampX, ells, by.x=c("locCode1", "locCode2"), by.y=c("loc2", "loc1"))
			sampX <- unique(rbind(sampX_12, sampX_21))
			# merge with metric (already symmetrized) for getting value
			sampX <- merge(sampX, my.metric, by=c("ind1", "ind2"), all.x=TRUE)
			# sort and compute a rowid
			sampX <- sampX[order(sampX$ind1, sampX$ind2) , ]
			sampX$id <- 1:nrow(sampX)
			rm(my.sampleCode.locCode, sampX_12, sampX_21) # free memory
		})
		tv <- as.vector(t) ; tv[is.na(tv)] <- 0.0
		# NOTE: %s used for stringified integers due to overflow for very large datasets (thanks to Simon Dellicour)
		message(sprintf("... %s sample pairs.  [user: %0.3f, system: %0.3f, elapsed: %0.3f seconds]", as.character(nrow(sampX)), tv[1]+tv[4], tv[2]+tv[5], tv[3]))
		
		
		## Expand results from intersection between grid and ellipses into a clean list of vectors of sample pairs ids
		message("Matching grid cells with sample pairs...")
		t <- system.time({
			# sample pairs with symmetrized localities for faster merges
			sampX2 <- rbind(sampX[ , c("id", "locCode1", "locCode2")], sampX[ , c("id", "locCode2", "locCode1")])
			data.table::setkey(sampX2, "locCode1", "locCode2")
			# iterate on intersection
			inter <- lapply(inter0, function(r) {
				# get matched ellipses row numbers
				ids <- as.vector(r)
				# get matched ellipses
				myElls <- ells[ids, ]
				# join with symmetrized sample pairs
				tmp1 <- merge(myElls, sampX2, by.x=c("loc1", "loc2"), by.y=c("locCode1", "locCode2"))
				# extract unique sample pair ids
				v <- as.vector(unique(tmp1$id))
				return(v)
			})
			# Fast count
			nbMatches <- countMatches_cpp(inter)
		})
		tv <- as.vector(t) ; tv[is.na(tv)] <- 0.0
		if (is.na(nbMatches)) { nbMatches <- 0 }
		# NOTE: %s used for stringified integers due to overflow for very large datasets (thanks to Simon Dellicour)
		message(sprintf("... %s matches found.  [user: %0.3f, system: %0.3f, elapsed: %0.3f seconds]", as.character(nbMatches), tv[1]+tv[4], tv[2]+tv[5], tv[3]))
		
		
		## Get results for unpermuted MAPI analysis
		# For each cell, we compute the number of intersecting ellipses, the sum of their weightd, the weighted mean and weighted standard deviation.
		message("Computing values for grid cells...")
		t <- system.time({
			# Call to C++ function parseInter_cpp which iterates on grid, intersetion and sample pairs for computing values
			resu <- data.table::as.data.table(parseInter_cpp(grid$gid, inter, sampX$weight, sampX$value))
			colnames(resu) <- c("gid", "nb_ell", "avg_value", "sum_wgts", "w_stdev")
			data.table::setkey(resu, "gid")
			# computes sum-of-weights percentile
			resu$swQ <- unclass(.bincode(resu$sum_wgts, stats::quantile(resu$sum_wgts, 0:100/100.0, na.rm=TRUE), include.lowest=TRUE))
			# merge to grid for getting geometry
			resu <- merge(grid, resu, by="gid", all.x=TRUE)
		})
		tv <- as.vector(t) ; tv[is.na(tv)] <- 0.0
		message(sprintf("... done.  [user: %0.3f, system: %0.3f, elapsed: %0.3f seconds]", tv[1]+tv[4], tv[2]+tv[5], tv[3]))
		

		
		## Let's start permutations block (if any)
		##########################################################################################################################################################
		if(nbPermuts > 0) {
			t <- system.time({
				message(sprintf("Starting %d permutations...", nbPermuts))
				
				# prepare permuted samples table
				my.samples_perm <- data.table::data.table(ind=my.samples$ind)
				data.table::setkey(my.samples_perm, "ind")
				
				# prepare permuted results table
				resu2 <- data.table::as.data.table(matrix(nrow=nrow(resu), ncol=nbPermuts))
				
				# prepare permutation function
				doPerm <- function(permut) {
					# shuffle individuals 
					my.samples_perm$indP <- sample(my.samples$ind)
					
					# replace unpermuted individuals by permuted ones
					sampX.P <- merge(sampX,   my.samples_perm, by.x="ind1", by.y="ind")
					sampX.P <- merge(sampX.P, my.samples_perm, by.x="ind2", by.y="ind", suffixes=c("1","2"))
					
					# merge permuted samplepairs with metric
					sampX.P <- merge(sampX.P, my.metric, all.x=TRUE, by.x=c("indP1","indP2"), by.y=c("ind1","ind2"))
					data.table::setkey(sampX.P, "id")
					
					# compute permuted result
					outPerm <- parseInterPerm_cpp(grid$gid, inter, sampX.P$weight, sampX.P$value.y)
					
					# tick progress bar if any
					if (exists("pb")) {  pb$tick()  }
					
					return(outPerm)
				}
				
				# garbage collecting in order to release memory before forks may occur
				gc(verbose=FALSE)
				
				# run permutation function (parallel or not)
				if (base::requireNamespace("pbapply", quietly=TRUE) && nbCores > 1) {
					message(sprintf("... parallelized over %d cores ...", nbCores))
					# iterate (parallelized) on permuted results table columns
					resu2 <- pbapply::pbapply(resu2, 2, doPerm, cl=nbCores)
				} else {
					message(sprintf("... unparallelized ..."))
					# try to set up a progress bar
					if (base::requireNamespace("progress", quietly=TRUE)) {
						pb <- progress::progress_bar$new(total=nbPermuts)
					}
					# iterate (un parallelized) on permuted results table columns
					resu2 <- apply(resu2, 2, doPerm)
				}
				
				message("Computing probabilities...")
				# add permuted values to results table as a list of numeric values
				resu$permuts <- NA
				for(i in 1:nrow(resu)){
					resu$permuts[i] <- list(resu2[i,])
				}
				
				# Probability computation: rank of the unpermuted value among permuted values
				resu$proba <- apply(resu,1,function(cell){
					val <- as.numeric(cell["avg_value"])
					perms <- as.numeric(cell["permuts"][[1]])
					return( sum(perms < val) / length(perms) )
				})
				
				# Benjamini-Yekutieli correction of probabilities for upper and lower tails
				resu$ltP <- stats::p.adjust(resu$proba, method="BY")
				resu$utP <- stats::p.adjust((1.0 - resu$proba), method="BY")
				
			})
			tv <- as.vector(t) ; tv[is.na(tv)] <- 0.0
			message(sprintf("... done.  [user: %0.3f, system: %0.3f, elapsed: %0.3f seconds]", tv[1]+tv[4], tv[2]+tv[5], tv[3]))
		}
		##########################################################################################################################################################
	})
	tv <- as.vector(tot) ; tv[is.na(tv)] <- 0.0
	message(sprintf("MAPI COMPUTATION ENDED.  [TOTAL: user: %0.3f, system: %0.3f, elapsed: %0.3f seconds]", tv[1]+tv[4], tv[2]+tv[5], tv[3]))
	return(resu)
}
