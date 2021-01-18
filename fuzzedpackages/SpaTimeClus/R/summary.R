#'
#' Summary function.
#' 
#' This function gives the summary of an instance of \code{\linkS4class{STCresults}}.
#' 
#' @param object instance of \code{\linkS4class{STCresults}}.
#' 
#' @name summary
#' @rdname summary-methods
#' @docType methods
#' @exportMethod summary
#' 
#' 
NULL

#' @rdname summary-methods
#' @aliases summary summary,STCresults-method
#'
setMethod(
  f="summary",
  signature = c("STCresults"),
  definition = function(object){
    cat("****************************************************************************************\n")
    cat("Data: \n")
    cat("  sample size:      ", object@data@n, "\n")
    cat("  time grid size:   ", object@data@TT, "\n")
    cat("  number of sites:  ", object@data@JJ, "\n")  

    cat("\n")
    cat("****************************************************************************************\n")
    cat("Model: \n")  
    cat("  number of components:             ", object@model@G, "\n")
    cat("  number of polynoms per component: ", object@model@K, "\n")
    cat("  degrees of polynoms:              ", object@model@Q, "\n")
    cat("  number of parameters:             ", object@model@nbparam, "\n")
    if (object@model@spatial) cat("  Spatial dependencies are modelled for each component.\n")
    if (object@model@spatial==0) cat("  Spatial dependencies are not modelled for each component.\n")
    
    cat("\n")
    cat("****************************************************************************************\n")
    cat("Criteria: \n")  
    cat("  loglikelihood  :", object@criteria@loglike, "\n")
    cat("  AIC            :", object@criteria@AIC, "\n")
    cat("  BIC            :", object@criteria@BIC, "\n")
    cat("  ICL            :", object@criteria@ICL, "\n")
    cat("  degeneracy rate:", object@criteria@degeneracy, "\n")
    cat("\n")
    cat("****************************************************************************************\n")
  }
)




#'
#' Summary function.
#' 
#' This function prints the elements  of an instance of \code{\linkS4class{STCresults}}..
#' 
#' @param x an instance of \code{\linkS4class{STCresults}}.
#' 
#' @name print
#' @rdname print-methods
#' @docType methods
#' @exportMethod print
#' 
#' 
NULL

#' @rdname print-methods
#' @aliases print print,STCresults-method
#'
setMethod(
  f="print",
  signature = c("STCresults"),
  definition = function(x){
    summary(x)    
    cat("Parameter detail\n")
    for (g in 1:x@model@G){
      cat("\n")
      cat("****************************************************************************************\n")
      cat("Parameters of component", g, "\n") 
      cat("proportion", x@param@proportions[g], "\n")
      cat("logistic coefficients \n")
      print(x@param@lambda[[g]])
      cat("\n")
      cat("polynom coefficients \n")
      deg <- max(x@model@Q)+2
      coeff <- matrix(0, x@model@K, deg)
      rownames(coeff) <- paste("polynom",1:x@model@K, sep=".")
      colnames(coeff) <- c(paste("coeff",(0:x@model@Q), sep="."), "variance")
      for (h in 1:x@model@K)  coeff[h,1:(x@model@Q+1)] <- x@param@beta[[g]][h,]
      coeff[,ncol(coeff)] <- x@param@sigma[g,]
      print(coeff)
      cat("\n")
    }
  }
)