#' Set the seed from a fitted evmSim object.
#' 
#' Set the seed from a fitted \code{evmSim} object to ensure reproducibility of
#' output.
#' 
#' Sets the seed to the value used to fit the model.
#' 
#' @usage evmSimSetSeed(x)
#' @param x An object of class \code{evmSim}, as returned by \code{evm} using
#' \code{method = "simulate"}.
#' @author Harry Southworth
#' @seealso \code{\link{evm}}
#' @examples
#'  \donttest{ 
#'  data <- rnorm(1000)
#'  mod <- evm(data, qu=.7, method="simulate")
#'  evmSimSetSeed(mod)
#'  mod1 <- evm(data, qu=.7, method="simulate") # this produces the same MCMC output as mod
#' }
#' 
#' @export evmSimSetSeed
evmSimSetSeed <- function(x){
  if (oldClass(x) != "evmSim"){
      stop("This function expects an object of class \'evmSim\'.")
  }

  assign(".Random.seed", x$seed, envir=.GlobalEnv)

  invisible(x$seed)
}
