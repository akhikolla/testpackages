#' Conversion among data and phylogeny objects
#'
#' Function to generate a "rateData" object containing the discrete explanatory variable, continuous response variable and set of variance co-variance matrices. It loads the trait data and removes species with missing data from the data and vcv matrix. 
#' \code{as.rateData} requires either a set of matrices in rateMatrix format created using \code{as.rateMatrix} or, if no rateMatrix object is input then it requires a phylogeny in "phylo" format. If  a "phylo" object is used  \code{as.rateData} will call \code{as.rateMatrix} internally.
#' \code{as.rateMatrix} calls the "ape" function \code{vcv.phylo} multiple times and this can be slow for large phylogenies. It will often be more efficient to use \code{as.rateMatrix} first to create a "rateMatrix" object to pass to \code{as.rateData}, particularly if there are many response traits of interest to be fitted to the same phylogeny and set of reconstructed ancestral states.
#' @param y The response variable - typically a continuous trait. Specified as a column name or index
#' @param x The explanatory (discrete) variable used to define the hypothesised rate categories. Specified as a column name or index.
#' @param rateMatrix A "rateMatrix" object or NULL
#' @param phy An object of class "phylo" (see \pkg{ape}).
#' @param data A data frame containing (minimally) the x and y variables as columns with species names as rownames.
#' @param meserr.col Column name or index containing measurement error for species means.
#' @param meserr.propn Single value specifying the proportional measurement to be applied across all species.
#' @param log.y Logical, natural log transform response variable.
#' @param report_prune Logical. Prints a list of dropped species if TRUE
#' @return rateData An object of class "rateData" which is a list containing the response (y) and explanatory (x) variable along with a list of variance-covaraince matrices.
#' @references Thomas GH, Meiri S, & Phillimore AB. 2009. Body size diversification in Anolis: novel environments and island effects. Evolution 63, 2017-2030.
#' @author Gavin Thomas
#' @export
#' @examples
#' ## Read in phylogeny and data from Thomas et al. (2009)
#' data(anolis.tree)
#' data(anolis.data)
#'
#' ## Convert data to class rateData with a rateMatrix object as input
#' anolis.rateMatrix <- as.rateMatrix(phy=anolis.tree, x="geo_ecomorph", data=anolis.data)
#' 
#' anolis.rateData <- as.rateData(y="Female_SVL", x="geo_ecomorph", 
#' rateMatrix = anolis.rateMatrix, phy=NULL, data=anolis.data, log.y=TRUE)  
#'
#' ## Convert data to class rateData with a phylo object as input 
#' anolis.rateData <- as.rateData(y="Female_SVL", x="geo_ecomorph", 
#' rateMatrix = NULL, phy=anolis.tree, data=anolis.data, log.y=TRUE)
#' @export

as.rateData <-
function(y, x, rateMatrix=NULL, phy=NULL, data, meserr.col=NULL, meserr.propn=NULL, log.y=FALSE, report_prune=FALSE) {
	
		 if (is.numeric(y) && length(y) == 1) { y <- colnames(data)[y] } else { y <- y }
		 if (is.numeric(x) && length(x) == 1) { x <- colnames(data)[x] } else { x <- x }
		if (is.numeric(meserr.col) ) { meserr.col <- colnames(data)[meserr.col] } else { meserr.col <- meserr.col }

		if (!is.null(meserr.col) & !is.null(meserr.propn) ) {
				stop("Only one source of measurement error should be specified. Use messerr.col to specify a column of measurement error 
						OR use meserr.propn to specify a single numeric value as a proportion error (e.g. enter 0.05 for 5% measurement error.")
						}
	
		if (is.null(rateMatrix))	{ rateMatrix <- as.rateMatrix(phy=phy, x=x, data=data) } else { rateMatrix <- rateMatrix }
	
		if (is.null(meserr.col)) { dat <- data.frame(x=data[,x], y=data[,y], row.names = rownames(data)) }
					 
		if (!is.null(meserr.col)) { dat <- data.frame(x=data[,x], y=data[,y], meserr=data[,meserr.col], row.names = rownames(data)) }
		if (!is.null(meserr.propn)) { dat <- data.frame(x=data[,x], y=data[,y], meserr=rep(meserr.propn, dim(data)[1]), row.names = rownames(data)) }
					 
			 
	
	
		if(sum(sort(rownames(dat)) != sort(rownames(rateMatrix[[1]]))) > 0){warning("Length or names of phenotypic and of phylogenetic data do not match - non-matching taxa will be dropped")}	# check names

		
		# get names of species common to data frame and phylogeny
		sharedSpecies <- intersect(rownames(rateMatrix[[1]]), rownames(dat))
	
	
		# subset data to only those species in phylogeny and alphabetise
		dat <- dat[match(sharedSpecies ,rownames(dat)),] 
		dat<- dat[sort(rownames(dat), index.return = TRUE)$ix, ]
		
		# index missing data for y
		ccdat <- complete.cases(dat)
		idx <- which(ccdat)
		missing.dat <- rownames(dat)[which(!ccdat)]

	
		# Suppress dropped species list
		if(report_prune==TRUE) {
			if(length(missing.dat > 0)) {
				cat("Dropping species due to missing data (x, y, or meserr):", "\n")
				cat("Dropped species: ", missing.dat, "\n")
			}
			}
			
		# alphabetise and prune rate matrices to remove data with missing species
		Vmat <- vector(mode="list", length = length(rateMatrix))
		sumMat <- matrix(0, dim(rateMatrix[[1]])[1], dim(rateMatrix[[1]])[2])
		for(i in 1:length(rateMatrix)) {
			Vmatrix <- rateMatrix[[i]]
			sumMat <- sumMat + Vmatrix
			nms <- rownames(Vmatrix)
			snms <- sort(nms, index.return = TRUE)
			Vmatrix <- Vmatrix[snms$ix, snms$ix]
			Vmatrix <- Vmatrix[idx, idx]
			Vmat[[i]] <- Vmatrix
			
			}
		
	
			xdat <- as.matrix(dat[,1])
			if (log.y==TRUE) {ydat <- log(dat[idx, "y"]) } else {ydat <- dat[idx, "y"]}
			xdat <- xdat[idx, ]
			
			if (is.null(meserr.col) & is.null(meserr.propn)) {traits <- list(y = ydat, x = xdat, Vmat = Vmat, meserr=rep(0,length(ydat)))} 
			if (!is.null(meserr.col)) {	meserr <- dat[idx,"meserr"]^2
				traits <- list(y = ydat, x = xdat, Vmat = Vmat, meserr=meserr)}
			
			if (!is.null(meserr.propn)) {	meserr <- (dat[idx,"meserr"]*ydat)^2
				traits <- list(y = ydat, x = xdat, Vmat = Vmat, meserr=meserr)}
	
	
			attr(traits, "class") <- "rateData"
			return(traits)
			}




