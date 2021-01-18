#' Coefficient Methods for flexreg Objects
#'
#' @param object an object of class \code{`flexreg`}, usually the result of \code{\link{flexreg}}.
#' @param ... additional arguments. Currently not used.
#' @rdname summary.flexreg
#'
#' @export
#'
coef.flexreg <- function(object, ...){
  summ <- summary(object)
  mu.model <- summ$Summary.mu[,1]
  phi.model <- summ$Summary.phi[,1]
  names(phi.model) <- rownames(summ$Summary.phi)

  l <- list("mean_model"=mu.model,"precision_model"=phi.model)

  if(!is.null(summ$Summary.add)){
    additional.par <- summ$Summary.add[,1]
    l[[3]] <- additional.par
    names(l)[3] <- "additional_par"
  }
  return(l)
}

