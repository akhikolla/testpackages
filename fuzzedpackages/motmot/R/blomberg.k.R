#' Blomberg's K
#' Estimate Blomberg's K (Blomberg et al. 2003)
#' @param phy An object of class \code{phylo} (see \pkg{ape}).
#' @param y A matrix of trait values.
#' @return The estimate of the K statistic
#' @useDynLib motmot
#' @importFrom Rcpp sourceCpp
#' @seealso \code{\link{transformPhylo.ML}}, the Picante package
#' @references Blomberg SP, Garland T, & Ives AR. 2003. Testing for phylogenetic signal in comparative data: behavioral traits are more labile. Evolution 57, 717-745.

blomberg.k <- function(phy, y) {
	bm.phy <- transformPhylo.ML(phy, model="bm", y=y )
	ntip <- Ntip(phy)
	nodetimes <- nodeTimes(phy)
	vcv.int <- vcv(phy)
	vcv.phy <- sum(solve(vcv.int))
	diag.vcv <- sum(diag(vcv.int))
	mean.sq.error.0  <- bm.phy$brownianVariance[[1]]
	mat.state <- y[,1] - bm.phy$root.state
	transpose.mat <- t(mat.state)
	mean.sq.error <- transpose.mat %*% mat.state / (ntip - 1)
	mean.sq.error.div <- 1 / (ntip - 1) * ((diag.vcv - ntip / vcv.phy))

	obs <- mean.sq.error / mean.sq.error.0
	exp <- mean.sq.error.div
	k.out <- obs / exp
	return(k.out)
	}
	
	
	