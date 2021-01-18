#' @title Compiles \code{SimBIID_model} object
#'
#' @description Compiles an object of class \code{SimBIID_model} into an
#'              \code{XPtr} object for use in Rcpp functions, or an
#'              object of class \code{function} for calling directly from R.
#'
#' @export
#'
#' @param model An object of class \code{SimBIID_model}.
#'
#' @return An object of class \code{XPtr} that points to the compiled function, or
#'         an R \code{function} object for calling directly from R.
#'         
#' @seealso \code{\link{mparseRcpp}}
#'         
#' @examples 
#' \donttest{
#' ## set up SIR simulationmodel
#' transitions <- c(
#'     "S -> beta * S * I -> I", 
#'     "I -> gamma * I -> R"
#' )
#' compartments <- c("S", "I", "R")
#' pars <- c("beta", "gamma")
#' model <- mparseRcpp(
#'     transitions = transitions, 
#'     compartments = compartments,
#'     pars = pars
#' )
#' 
#' ## compile model to be run directly
#' model <- compileRcpp(model)
#' model
#' 
#' ## set initial states (1 initial infection 
#' ## in population of 120)
#' iniStates <- c(S = 119, I = 1, R = 0)
#' 
#' ## set parameters
#' pars <- c(beta = 0.001, gamma = 0.1)
#' 
#' ## run compiled model
#' model(pars, 0, 100, iniStates)
#' }
#' 

compileRcpp <- function(model) {
    if(missing(model)) {
        stop("'model' object missing")
    }
    if(class(model) != "SimBIID_model") {
        stop("'model' object not of class 'SimBIID_model'")
    }
    ## write to temporary file
    filename <- tempfile()
    writeLines(model$code, filename)
    
    ## compile into external pointer or function
    source(filename)
    
    ## return pointer or function
    Rcpp_object
}

