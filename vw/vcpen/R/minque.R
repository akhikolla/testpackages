#' MINQUE estimation of variance components
#'
#' Estimate variance components by MINQUE method, allowing multiple iterations
#'
#' @param y Numeric vector of traits. Only continuous trait currently allowed.
#' @param X Matrix of covariates (columns) for subjects (rows), matching subjects in the trait (y) vector. 
#' @param Kerns List of kernel matrices: a kernel matrix for each  variance compenent. The last kernel matrix in the list (an identity matrix) is for the residual variance component.
#' @param n.iter Number of  minque iterations
#' @param eps Default small positive value for non-positive vc estimates within iterations.
#' @return List with estimates of variance components (vc), covariate regression coefficients (beta), and residuals of model fit.
#' @examples
#' data(vcexample)
#' nvc <- 1+length(unique(doseinfo[,2]))
#' id <- 1:nrow(dose)
#' ## vcs for genetic kernel matrices
#' Kerns <- vector("list", length=nvc)
#' for(i in 1:(nvc-1)){
#'   Kerns[[i]] <- kernel_linear(dose[,grep(i, doseinfo[,2])])
#'   rownames(Kerns[[i]]) <- id
#'   colnames(Kerns[[i]]) <- id
#' }
#' ## vc for residual variance
#' Kerns[[nvc]] <- diag(nrow(dose))
#' rownames(Kerns[[nvc]]) <- id
#' colnames(Kerns[[nvc]]) <- id
#' prefit  <- minque(response, covmat, Kerns, n.iter=2)
#' prefit[1]
#' prefit[2]
#' fit <- vcpen(response, covmat, Kerns, vc_init = prefit$vc)
#' 
#' @author JP Sinnwell, DJ Schaid
## methods for minque
#' @name minque
#' @rdname minque
#' @export
minque <- function(y, X, Kerns, n.iter=1, eps=0.001){
  ## init values of vc
  vc <- rep(.5, length(Kerns))
  vc[length(Kerns)] <- 1
  
  ## eps <- 0.001

  for(i in 1:n.iter){
    fit <- minque_Rcpp(y, X, Kerns, vc)
    vc <- fit$vc
    ## for negative vc, set to small pos value for next iter
    vc <- ifelse(vc < 0, eps, vc)
  }

  fit$vc <- ifelse(fit$vc < 0, 0, fit$vc)
  return(fit)
}
