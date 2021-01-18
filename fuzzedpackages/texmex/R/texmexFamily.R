#' Create families of distributions
#' 
#' Create families of distributions for use with extreme value
#' modelling.
#' 
#' The \code{density}, \code{rng}, \code{prob} and \code{quant}
#' functions can be simple wrappers for the usual d, r, p and q
#' functions. They should take a matrix with number of columns equal
#' to the number of parameters, and a fitted model object even if the
#' model object is not used by the function.
#' 
#' Examples of "texmexFamily" objects are \code{gpd}, \code{gev}, \code{glo},
#' \code{gpdIntCensored}, \code{weibull}, \code{gumbel} and \code{egp3}.  Take a look at
#' those objects to see how the functions should be constructed.
#' 
#' The functions are used by the modelling functions to create
#' diagnostic plots, predictions, etc..
#' 
#' @aliases texmexFamily print.texmexFamily summary.texmexFamily gpd gpdIntCensored glo
#'     gev egp3 cgpd weibull gumbel print.summary.texmexFamily
#' @usage texmexFamily(name, log.lik, param, info = NULL, sandwich =
#'     NULL, start = NULL, resid = NULL, rl, delta, endpoint, density,
#'     rng, prob, quant)
#'     \method{print}{texmexFamily}(x,...)
#'     \method{summary}{texmexFamily}(object,...)
#'     \method{print}{summary.texmexFamily}(x,...)
#' @param name The name of the distribution.
#' @param log.lik The distribution's log-likelihood function.
#' @param param The names of the parameters in the model.
#' @param info Function to compute the information matrix. If not
#'     provided, the modelling functions will work with a numerical
#'     approximation.
#' @param sandwich Function to compute the filling in the Huber
#'     sandwich estimator of the covariance matrix of parameter
#'     estimates, used for dependent data. Only implemented in family
#'     \code{gpd}.
#' @param start Function to compute starting parameters for the
#'     model. If not provided, the modelling functions will try to
#'     guess.
#' @param resid Function to compute residuals for the model.
#' @param rl Function to compute return levels.
#' @param delta Function to compute adjustments for covariance for
#'     return levels.
#' @param endpoint Function to compute the upper or lower endpoint of
#'     the fitted distribution.
#' @param density Function to compute the density.
#' @param rng Function for random number generation.
#' @param prob Function to compute cumulative probabilities.
#' @param quant Function to compute quantiles.
#' @param ... Additional arguments to the print and summary methods.
#' @param x,object An object of class 'texmexFamily'.
#' @return A object of class "texmexFamily", which is essentially a
#'     list containing the input arguments. If \code{info},
#'     \code{sandwich}, \code{start}, \code{resid} are not provided,
#'     they default to \code{NULL}.
#' @note The \code{gpd}, \code{gev}, \code{weibull}, generalised logistic (\code{glo}), 
#' \code{gumbel}, \code{gpdIntCensored} and \code{egp3} families are provided. 
#' The \code{\link{evm}} function defaults to using the \code{gpd} family.
#' @author Harry Southworth
#' @seealso \code{\link{evm}}
#' @keywords models
#' @export texmexFamily
texmexFamily <-
    # Create an object of class 'texmexFamily'. It is allowable to have
    # info, start and resid as NULL, but all other information must
    # be provided by the user.
function(name, log.lik, param, info=NULL, sandwich = NULL, start=NULL, resid=NULL,
                         rl, delta, endpoint, density, rng, prob, quant){
    res <- list(name=name, log.lik=log.lik, param=param, info=info, sandwich = sandwich, start=start,
                resid=resid, rl=rl, delta=delta, endpoint=endpoint, density=density,
                rng=rng, prob=prob, quant=quant)

    oldClass(res) <- 'texmexFamily'
    res
}


#' @export
print.texmexFamily <- function(x, ...){
    cat('Family:      ', x$name, '\n')
    invisible(x)
}

#' @export
summary.texmexFamily <- function(object, ...){
    res <- list(object=object)
    if (is.null(object$info)){ 
        res$info <- 'Numerical approximation' 
    } else { 
        res$info <- 'Closed form' 
    }
    oldClass(res) <- "summary.texmexFamily"
    res
}

#' @export
print.summary.texmexFamily <- function(x, ...){
    print.texmexFamily(x$object, ...)
    cat('Parameters:  ', x$object$param, '\n')
    cat('Information: ', x$info, '\n')
    invisible(x)
}
