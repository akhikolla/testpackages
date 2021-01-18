#' Log-likelihood estimation for traits and phylogenies
#'
#' This function calculates the log-likelihood and Brownian (co)variance for a trait(s) and a phylogeny using phylogenetically independent contrasts
#' Note that \code{as.rateMatrix} calls the CAIC function \code{vcv.array} multiple times and this can be slow for large phylogenies (though faster than using the "ape" equivalent \code{vcv.phylo}).
#' @param y A matrix of trait data. Rownames must be included and match the taxon names in the phylogeny. Can accept single traits (calculates variance) or multiple traits (calculates variance-covariance matrix).
#' @param phy An object of class "phylo" (see \pkg{ape}).
#' @param covPIC Logical - allow for covariance between multivariate traits (\code{TRUE}), or assume not covariance (\code{FALSE}). Only applicable to multivariate traits
#' @param brCov If \code{NULL} (the default), Brownian covariance is analytically estimated. If a user-supplied numerical value is suppied the likelihood is calculate given this value
#' @details The \code{phylo} object must be rooted and fully dichotomous
#' @return brownianVariance Brownian variance (or covariance for multiple traits) given the data and phylogeny
#' @return logLikelihood The log-likelihood of the model and data
#' @references Felsenstein J. 1973. Maximum-likelihood estimation of evolutionary trees from continuous characters. Am. J. Hum. Genet. 25, 471-492.
#' @references Felsenstein J. 1985. Phylogenies and the comparative method. American Naturalist 125, 1-15.
#' @references Freckleton RP & Jetz W. 2009. Space versus phylogeny: disentangling phylogenetic and spatial signals in comparative data. Proc. Roy. Soc. B 276, 21-30. 
#' @author Gavin Thomas, Rob Freckleton
#' @examples 
#' data(anolis.tree)
#' data(anolis.data)
#' ## calculate Brownian variance log-likelihood of female SVL
#' female.svl <- matrix(anolis.data[,"Female_SVL"], 
#' dimnames=list(rownames(anolis.data)))
#' input.data <- sortTraitData(phy=anolis.tree, y=female.svl, log.trait=TRUE)
#' likTraitPhylo(phy = input.data$phy, y=input.data$trait)
#' @export

likTraitPhylo<-function (y, phy, covPIC = TRUE, brCov=NULL)
{
    if (is.matrix(y) == FALSE) {
        stop("Trait data must be a matrix with taxon names as row names")
    }
    n <- length(phy$tip.label)
    k <- ncol(y)
    phy <- reorder(phy, order = "pruningwise")
    y <- as.matrix(y[phy$tip.label, ])
    contrasts <- apply(y, 2, pic.motmot, phy = phy)
    rawVariances <- c(contrasts[[1]]$contr[, 2], contrasts[[1]]$V)   
    rawContrasts <- sapply(contrasts, function(k) c(k$contr[, 1], 0))
  	t.mat <- apply(rawContrasts, 2, function(x) x / sqrt(rawVariances))
    if(is.null(brCov)) {
    		brCov <- crossprod(t.mat, t.mat) / (n-1)
    } else {
    		brCov <- as.matrix(brCov)
    }
    if(covPIC == FALSE) {
    		brCov[upper.tri(brCov)] <- 0
    		brCov[lower.tri(brCov)] <- 0
    }
    iW <- solve(brCov)
    addCon.mat <- apply(rawContrasts, 1, function(con) crossprod(con, iW %*% con)) / rawVariances
    addCons <- sum(addCon.mat)
	logLikelihood <- -0.5 * (n * k * log(2 * pi) + n * log(det(brCov)) + k * sum(log(rawVariances)) + addCons)
    return(list(brownianVariance = brCov, logLikelihood = logLikelihood))
}
