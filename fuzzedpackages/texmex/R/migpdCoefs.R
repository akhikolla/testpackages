#' Change values of parameters in a migpd object
#' 
#' Change the values of parameters in a \code{migpd} object. You might want to
#' do this after modelling marginal distributions as functions of covariates.
#' 
#' 
#' @usage migpdCoefs(object, which, coefs)
#' @param object An object of class \code{migpd}.
#' @param which Which models in the \code{migpd} object you want to change.
#' @param coefs The coefficients that you want to change to. If \code{which}
#' has length 1, \code{coefs} can be a vector of parameters.  Otherwise, it
#' should be a list of vectors, and the list should have the same length as
#' \code{which}
#' @return A \code{migpd} object. See the help for \code{\link{migpd}}.
#' @author Harry Southworth
#' @seealso \code{\link{migpd}}
#' @keywords multivariate
#' @examples
#' 
#' library(MASS)
#' liver <- liver
#' liver$ndose <- as.numeric(liver$dose)
#' d <- data.frame(alt = resid(rlm(log(ALT.M) ~ log(ALT.B) + ndose, data=liver)),
#'                 ast = resid(rlm(log(AST.M) ~ log(AST.B) + ndose, data=liver)),
#'                 alp = resid(rlm(log(ALP.M) ~ log(ALP.B) + ndose, data=liver)),
#'                 tbl = resid(rlm(log(TBL.M) ~ log(TBL.B) + ndose, data=liver)))
#' 
#' Dgpds <- migpd(d[liver$dose == "D", 1:4], mqu=.7)
#' 
#' d$ndose <- liver$ndose
#' galt <- evm(alt, data=d, qu=.7, xi = ~ ndose)
#' gast <- evm(ast, data=d, qu=.7, xi = ~ ndose)
#' galp <- evm(alp, data=d, qu=.7, xi = ~ ndose)
#' 
#' altco <- predict(galt,type="lp",newdata=data.frame(ndose=4))$obj$link[1:2]
#' astco <- predict(gast,type="lp",newdata=data.frame(ndose=4))$obj$link[1:2]
#' alpco <- predict(galp,type="lp",newdata=data.frame(ndose=4))$obj$link[1:2]
#' 
#' Dgpd <- migpdCoefs(Dgpds, which=c("alt", "ast", "alp"),
#'                    coefs=list(altco, astco, alpco))
#' 
#' summary(Dgpd)
#' summary(Dgpds)
#' 
#' 
#' @export migpdCoefs
migpdCoefs <-
  # Coerce coefficients of a migpd object.
  # Useful when coefficients are coming from
  # a model with covariates and you want to
  # learn something about the dependence between
  # margins.
function(object, which, coefs){
  if (!inherits(object, "migpd")){
    stop("object must be of class \'migpd\'")
  }

  if (length(which) != length(coefs)){
    stop("which and coefs should have the same length")
  }

  if (length(which) == 1){
    object$models[[which]]$coefficients <- coefs
  }

  for (i in 1:length(which)){
    object$models[[i]]$coefficients <- coefs[[i]]
  }
  object
}

