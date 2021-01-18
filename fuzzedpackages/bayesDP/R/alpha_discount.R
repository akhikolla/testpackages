
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @title Bayesian Discount Prior: Historical Data Weight (alpha)
#' @description \code{alpha_discount} can be used to estimate the weight
#'   applied to historical data in the context of a one- or two-arm
#'   clinical trial. \code{alpha_discount} is not used internally but is
#'   given for educational purposes.
#' @param p_hat scalar. The posterior probability of a stochastic comparison.
#'   This value can be the output of \code{posterior_probability} or a value
#'   between 0 and 1.
#' @param discount_function character. Specify the discount function to use.
#'   Currently supports \code{weibull}, \code{scaledweibull}, and
#'   \code{identity}. The discount function \code{scaledweibull} scales
#'   the output of the Weibull CDF to have a max value of 1. The \code{identity}
#'   discount function uses the posterior probability directly as the discount
#'   weight. Default value is "\code{weibull}".
#' @param alpha_max scalar. Maximum weight the discount function can apply.
#'   Default is 1.
#' @param weibull_shape scalar. Shape parameter of the Weibull discount function
#'   used to compute alpha, the weight parameter of the historical data. Default
#'   value is 3.
#' @param weibull_scale scalar. Scale parameter of the Weibull discount function
#'   used to compute alpha, the weight parameter of the historical data. Default
#'   value is 0.135.
#' @details
#'   This function is not used internally but is given for educational purposes.
#'   Given inputs \code{p_hat}, \code{discount_function}, \code{alpha_max},
#'   \code{weibull_shape}, and \code{weibull_scale} the output is the weight
#'   that would be applied to historical data in the context of a one- or
#'   two-arm clinical trial.
#'
#' @return \code{alpha_discount} returns an object of class "alpha_discount".
#'
#' An object of class \code{alpha_discount} contains the following:
#' \describe{
#'  \item{\code{alpha_hat}}{
#'    scalar. The historical data weight.
#'   }
#'  }
#'
#' @references
#' Haddad, T., Himes, A., Thompson, L., Irony, T., Nair, R. MDIC Computer
#'   Modeling and Simulation working group.(2017) Incorporation of stochastic
#'   engineering models as prior information in Bayesian medical device trials.
#'   \emph{Journal of Biopharmaceutical Statistics}, 1-15.
#'
#' @examples
#' alpha_discount(0.5)
#'
#' alpha_discount(0.5, discount_function="identity")
#'
#' @rdname alpha_discount
#' @import methods
#' @importFrom stats sd density is.empty.model median model.offset model.response pweibull quantile rbeta rgamma rnorm var vcov
#' @aliases alpha_discount,ANY-method
#' @export alpha_discount
alpha_discount <- setClass("alpha_discount")

setGeneric("alpha_discount",
           function(p_hat             = NULL,
                    discount_function = "weibull",
                    alpha_max         = 1,
                    weibull_scale     = 0.135,
                    weibull_shape     = 3){
             standardGeneric("alpha_discount")
           })

setMethod("alpha_discount",
          signature(),
          function(p_hat             = NULL,
                   discount_function = "weibull",
                   alpha_max         = 1,
                   weibull_scale     = 0.135,
                   weibull_shape     = 3){

  # Check that discount_function is input correctly
  all_functions <- c("weibull", "scaledweibull", "identity")
  function_match <- match(discount_function, all_functions)
  if(is.na(function_match)) {
    stop("discount_function input incorrectly.")
  }


  # Compute alpha discount based on distribution
  if(discount_function == "weibull"){
    alpha_hat <- pweibull(p_hat, shape=weibull_shape,
                               scale=weibull_scale)*alpha_max
  } else if(discount_function == "scaledweibull"){
    max_p <- pweibull(1, shape=weibull_shape, scale=weibull_scale)

    alpha_hat <- pweibull(p_hat, shape=weibull_shape,
                               scale=weibull_scale)*alpha_max/max_p
  } else if(discount_function == "identity"){
    alpha_hat <- p_hat
  }

  return(alpha_hat)
})
