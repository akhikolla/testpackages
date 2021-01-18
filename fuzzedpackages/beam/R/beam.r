#' Bayesian inference in large Gaussian graphical models
#'
#' @param X n by p data matrix
#' @param type character. Either "marginal", "conditional" or "both". See Details.
#' @param return.only character. Either "cor", "BF", "prob". See details.
#' @param verbose logical. Whether information on progress should be be printed.
#' @param D matrix. Prior marginal correlation matrix. Must be positive definite, well-conditioned and have unit variance.
#'
#' @description 
#' This function carries out covariance and inverse-covariance estimation within the Gaussian conjugate model.
#' The scale matrix parameter of the inverse-Wishart is, by default, set to the identity, whereas the
#' degree of freedom parameter is estimated by maximization of the marginal likelihood.
#' The function also computes Bayes factors and tail probabilities (p-values) to recover the marginal and/or
#' conditional independence structure between variables.
#'
#' @details
#' The arguments \code{type} and \code{return.only} have essentially been introduced for computational and memory savings.
#' Using argument \code{type} the user may indicate whether the marginal dependencies ("marginal"), conditional dependencies
#' ("conditional") or both ("both") are to be inferred. On the other hand, the argument \code{return.only} is used to indicate
#' whether the correlations ("cor"), Bayes factors ("BF") or tail probabilities ("prob") should be returned.
#' Default is to return all three quantities both for marginal and conditional dependencies.
#'
#' @return An object of class \code{\link{beam-class}}
#'
#' @author Gwenael G.R. Leday and Ilaria Speranza
#'
#' @references
#' Leday, G.G.R. and Richardson, S. (2019). Fast Bayesian inference in large Gaussian graphical models. \emph{Biometrics}. 75(4), 1288--1298.
#'
#' @examples
#' # Load data
#' data(TCPAprad)
#' 
#' # beam
#' fit <- beam(X = TCPAprad, type="both") 
#' 
#' # Print summary
#' summary(fit)
#' 
#' # Extract matrix of marginal correlations
#' mcor(fit)[1:5, 1:5]
#' 
#' # Extract matrix of partial correlations
#' pcor(fit)[1:5, 1:5]
#' 
#' # Plot log-marginal likelihood of the Gaussian conjugate model
#' plotML(fit)
#' 
#' # Plot heatmap of marginal (upper triangle) and/or
#' # partial (lower triangle) correlation estimates
#' plotCor(fit)
#' 
#' @export

beam <- function(X, type = "conditional", return.only = c("cor", "BF", "prob"), verbose=TRUE, D=NULL){

	time0 <- proc.time()

	###########################################
	#              PREPROCESSING              #
	###########################################

	if(verbose){
	   cat("Preprocessing... ")
	}

	# Check input argument X
	assert_that(is.matrix(X))
	assert_that(not_empty(X))
	assert_that(noNA(X))
	assert_that(all(is.finite(X)))
  if(!is.null(colnames(X))){
    varlabs <- colnames(X)
  }else{
    varlabs <- character()
  }

  # Check input argument type
	assert_that(is.character(type))
	assert_that(not_empty(type))
	assert_that(length(type)==1)
	assert_that(noNA(type))
	assert_that(type%in%c("both", "marginal", "conditional"), msg="type is not recognized")
	
	# Check input argument return.only
	assert_that(is.character(return.only))
	assert_that(not_empty(return.only))
	assert_that(noNA(return.only))
	return.only <- unique(return.only)
	assert_that(all(return.only%in%c("cor", "BF", "prob")), msg="return.only is not recognized")

	# Check input argument D
	if(is.null(D)){
	  D <- matrix(0, ncol(X), ncol(X))
	}else{
  	assert_that(is.matrix(D))
  	assert_that(not_empty(D))
  	assert_that(nrow(D)==ncol(D), msg="D is not square")
  	assert_that(isSymmetric(D), msg="D is not symmetric")
  	assert_that(noNA(D))
  	assert_that(all(is.finite(D)))
	}
	
	#########################
	#         BEAM          #
	#########################
	
	# Run beam
	ind <- c("cor", "BF", "prob") %in% return.only
	res <- .beam(X=X, type=type, ronly=ind*1, D = D, verbose = verbose)

	# Add column labels
	labs <- NULL
	if(type=="marginal" || type=="both"){
	  labs <- c(labs, paste0("m_", c("cor", "logBF", "tail_prob")[ind]))
	}
	if(type=="conditional" || type=="both"){
	  labs <- c(labs, paste0("p_", c("cor", "logBF", "tail_prob")[ind]))
	}
	colnames(res$table) <- labs
	
	time1 <- proc.time() - time0

	#########################
	#        OUTPUT         #
	#########################

	# List
	out <- new("beam",
	           "table" = res$table,
	           "deltaOpt" = res$deltaOpt,
	           "alphaOpt"= res$alphaOpt,
	           "dimX"= dim(X),
	           "type"= type,
	           "varlabs" = varlabs,
	           "gridAlpha" = res$gridAlpha,
	           "valOpt" = res$valOpt,
	           "return.only" = return.only,
	           "time" = time1[3],
	           "TinvStdev" = res$TinvStdev[,1],
	           "s" = as.vector(res$s))

	return(out)

}

