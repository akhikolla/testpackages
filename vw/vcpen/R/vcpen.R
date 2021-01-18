#' Penalized Variance Components
#'
#' Penalized Variance Component analysis 
#'
#' @param y Numeric vector of traits. Only continuous trait currently allowed.
#' @param X Matrix of covariates (columns) for subjects (rows), matching subjects in the trait (y) vector. 
#' @param Kerns List of kernel matrices: a kernel matrix for each  variance compenent. The last kernel matrix in the list (an identity matrix) is for the residual variance component.
#' @param frac1 Fraction of penalty imposed on L1 penalty, between 0 and 1 (0 for only L2; 1 for only L1 penalty).
#' @param lambda_factor Weight for each vc (values between 0 and 1) for how much it should be penalized: 0 means no penalty. Default value of NULL implies weight of 1 for all vc's.
#' @param lambda_grid Vector of lambda penalties for fitting the penalized model. Best to order values from largest to smallest so parameter estimates from a large penalty can be used as initial values for the next smaller penalty. Default value of NULL implies initial values  of seq(from=.10, to=0, by=-0.01).
#' @param maxiter Maximum number of iterations allowed during penalized fitting.
#' @param vc_init Numeric vector of initial values for variance components. Default value of NULL implies initial values determined by 2 iterations of minque estimation.
#' @param print_iter Logical:  if TRUE, print the iteration results (mainly for refined checks)
#' @param object Fitted vcpen object (used in summary method)
#' @param \dots Optional arguments for summary method
#' @param digits Signficant digits for summary method
#' 
#' @return object with S3 class vcpen
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
#' fit  <- vcpen(response, covmat, Kerns, frac1 = .6)
#' summary(fit)
#' 
#' @author JP Sinnwell, DJ Schaid
#' @name vcpen
NULL
#> NULL
#' @rdname vcpen
#' @export
vcpen <- function(y, X, Kerns, frac1=0.8, lambda_factor=NULL, lambda_grid=NULL,
                  maxiter=1000, vc_init=NULL, print_iter=FALSE){

  nvc <- length(Kerns)
  
  ## lambda_factor is the factor assigned for each vc, with values between
  ## 0 and 1; 0 means don't penalize, and 1 means give full weight
  ## for penalizing

  if( is.null(lambda_factor)){
    lambda_factor=rep(1, (nvc-1) )
  }


  if(is.null(lambda_grid)){
    lambda_grid <- seq(from=.10, to=0, by=-0.01)
  }
  
  ## vc_init is a vector of initial starting values of vcs

    if(is.null(vc_init)){
        fit.minque <- minque(y, X, Kerns, n.iter=2)
        vc_init <- fit.minque$vc
        vc_init <- ifelse(vc_init < 0.01, 0.01, vc_init)
    }
    
 
    fit <- vcpen_Rcpp(y, X, Kerns, lambda_factor, lambda_grid,
                      frac1, vc_init, maxiter,  print_iter)
  
  if(is.null(colnames(X))){
    xnames <- paste("x", 1:nrow(fit$beta_grid), sep=".")
  } else {
    xnames <- colnames(X)
  }
    
  dimnames(fit$beta_grid) <- list(xnames, fit$lambda_grid)

  df <- data.frame(t(fit$vc_grid))
  names(df) <- paste0("vc", 1:ncol(df))
  df <- cbind(lambda=fit$lambda_grid, df)
     
  fit$vc_grid <- df

  ## define eps for deciding non-zero VCs
  eps <- .001
  npar <- apply(fit$vc_grid[,-1] >  eps, 1, sum) + apply(abs(fit$beta_grid) > eps, 2, sum)
  bic_grid <- as.vector(-2*fit$logl_grid+ log(fit$n_subj)*npar)

  index <- 1:length(bic_grid)
  ## if ties, choose bic with larger lambda penalty
  is.min.bic <- bic_grid == min(bic_grid)
  index <- index[is.min.bic][1]
    
  fit$vc <-  as.vector(fit$vc_grid[index,-1])

  fit$beta <- as.vector(fit$beta_grid[, index])
  
  fit$grid_info <- data.frame(lambda = fit$lambda_grid, iter=fit$iter+1,
                          logl=fit$logl_grid, loglpen=fit$logllasso_grid,
                          bic=bic_grid, min_bic = is.min.bic)
   
  ## remove redundant info
  fit$iter_grid <- NULL
  fit$logl_grid <- NULL
  fit$logllasso_grid <- NULL
  fit$lambda_grid <- NULL
   
  class(fit) <- c("vcpen", "list")
  return(fit)
}

#' @name summary.vcpen
#' @rdname vcpen
#' @export
summary.vcpen <- function(object, ..., digits=4) {
  cat("vcpen object\n")
  cat(paste0("  N-subjects = ", object$n_subj, "\n"))
  cat(paste0("  N-VC = ", object$n_vc, "\n"))
  cat("\n  Model fits over lambda penalty grid:\n\n")
  print(object$grid_info, digits=digits, ...)
  cat("\n  VC estimates by lambda penalties:\n\n")
  print(object$vc_grid, digits=digits, ...)
  cat("\nEstimates with min BIC:\n")
  cat("beta:\n")
  print(object$beta, digits=digits)
  cat("VC estimates:\n")
  print(object$vc, digits=digits)
  invisible()    
}
