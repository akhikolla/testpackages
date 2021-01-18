#' @importFrom methods setOldClass
methods::setOldClass("pandemicEstimated")

traceplot_pandemicEstimated = function(object, ...){
  if (object$model_name == "poisson: static generalized logistic") p = c("a","b","c","f") else
  if (object$model_name == "negbin: static generalized logistic") p = c("a","b","c","f","phi") else
  if (object$model_name == "poisson: static seasonal generalized logistic"){
    p = c("a","b","c","f",paste0("d_",1:length(object$seasonal_effect)))
  } else
  if (object$model_name == "poisson: multi_waves(2)") p = c("a1","b1","c1","alpha1","delta1",
                                                    "a2","b2","c2","alpha2","delta2")
  rstan::traceplot(object$fit,pars=p,...)
}

#' Draw traceplot of the parameters for the pandemic model
#'
#' Uses stan's traceplot function to draw the traceplots for the relevant parameters of the estimated model.
#' @param object Output of the \code{\link{pandemic_model}} function
#' @param ... Aditional parameters passed on to the \code{\link[rstan]{traceplot}} function
#' @seealso \code{\link{pandemic_model}} and \code{\link[rstan]{traceplot}}
#' @examples
#' \dontrun{
#' dataMG = load_covid("Brazil","MG")
#' estimMG = pandemic_model(dataMG)
#' traceplot(estimMG)}
#' @importFrom rstan traceplot
#' @exportMethod traceplot
setMethod("traceplot","pandemicEstimated",traceplot_pandemicEstimated)

#' Draw estimated density of the parameters for the pandemic model
#'
#' Uses stan's stan_dens function to draw the marginal posterior for the relevant parameters of the estimated model.
#' Defined as a method for the stats::density generic function.
#' @param x Output of the \code{\link{pandemic_model}} function
#' @param ... Additional parameters passed on the \code{\link[rstan]{stan_dens}} function
#' @seealso \code{\link{pandemic_model}} and \code{\link[rstan]{stan_dens}}.
#' @examples
#' \dontrun{
#' dataMG = load_covid("Brazil","MG")
#' estimMG = pandemic_model(dataMG)
#' density(estimMG)}
#' @importFrom stats density
#' @importFrom rstan stan_dens
#' @method density pandemicEstimated
#' @export
density.pandemicEstimated = function(x, ...){
  if (x$model_name == "poisson: static generalized logistic") p = c("a","b","c","f") else
    if (x$model_name == "negbin: static generalized logistic") p = c("a","b","c","f","phi") else
      if (x$model_name == "poisson: static seasonal generalized logistic"){
        p = c("a","b","c","f",paste0("d_",1:length(x$seasonal_effect)))
      } else
        if (x$model_name == "poisson: multi_waves(2)") p = c("a1","b1","c1","alpha1","delta1",
                                                          "a2","b2","c2","alpha2","delta2")
        rstan::stan_dens(x$fit,pars=p,...)
}
