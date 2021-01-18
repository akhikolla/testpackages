#' Fast inference of a conditional independence graph
#'
#' @param X n by p data matrix
#' @param thres numeric. Significance threshold to be applied on adjusted tail probabilities.
#' @param method character. Method to use for multiple comparison adjustment of tail probabilities.
#' @param verbose logical. Whether information on progress should be be printed.
#'
#' @description 
#' Fast and memory efficient reconstruction of large conditional independence networks.
#' 
#' @details
#' This function is a wrapper for \code{\link{beam}} and \code{\link{beam.select}}. It is faster and outputs a lighter object.
#' Note that the choice of shrinkage target D is the identity matrix and other choices are not supported yet.
#'
#' @return An object of class \code{\link{dgCMatrix-class}}
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
#' # lightbeam
#' res <- lightbeam(X = TCPAprad, thres=0.1) 
#' 
#' @export

lightbeam <- function(X, thres = 0.1, method = "holm", verbose=TRUE){
  
  time0 <- proc.time()
  
  ###########################################
  #              PREPROCESSING              #
  ###########################################
  
  if(verbose){
    cat("--> Preprocessing... ")
  }
  
  # Check input argument X
  assert_that(is.matrix(X))
  assert_that(not_empty(X))
  assert_that(noNA(X))
  assert_that(all(is.finite(X)))
  if(!is.null(colnames(X))){
    varlabs <- colnames(X)
  }else{
    varlabs <- paste0("V", 1:ncol(X))
  }
  
  # Check input argument thres
  assert_that(is.numeric(thres))
  assert_that(not_empty(thres))
  assert_that(noNA(thres))
  assert_that(is.finite(thres), msg="thres is not finite")
  assert_that((thres>0) & (thres<1), msg="thres must be between 0 and 1")
  
  # Check input argument method
  assert_that(is.character(method))
  assert_that(not_empty(method))
  assert_that(length(method)==1)
  assert_that(noNA(method))
  assert_that(method%in%c("holm", "bonferroni", "BH", "BY", "HC"), msg="method is not recognized")
  

  #########################
  #       LIGHTBEAM       #
  #########################
  
  # Dimension of data
  n <- nrow(X)
  p <- ncol(X)
  
  # Compute quantile (threshold on rzij2 rather than tail probability)
  th <- qbeta(thres, shape1=0.5, shape2=(n-1)/2, lower.tail=FALSE)
  
  # Run fast beam
  smat <- .lightbeam(X=X, thres=th, verbose = verbose)
  
  # Adjust tail probabilities
  nbtests <- 0.5*p*(p-1)
  smat@x <- p.adjust(smat@x, method=method, n=nbtests)
  
  # Set entries above the threshold to 0
  smat@x[smat@x>thres] <- 0
  smat <- drop0(smat)

  # Add names
  colnames(smat) <- rownames(smat) <- varlabs

  time1 <- proc.time() - time0
  
  cat("estimated computing time (in seconds): ", time1[3], "\n")
  
  #########################
  #        OUTPUT         #
  #########################
  
  return(smat)
  
}

