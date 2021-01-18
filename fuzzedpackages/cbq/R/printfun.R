#' Extract CBQ Coefficients
#'
#' Create a table of coefficient results from a \code{cbq} object.
#'
#' @param object A \code{cbq} object.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A table of coefficients with their corresponding lower and upper bounds.
#' @export
#'
#' @method coef cbq
#'
coef.cbq <- coefficients.cbq <- function(object, ...) {
  coefmat <-
    cbind(matrix(object$means, nrow = object$npars), t(object$ulbs))
  row.names(coefmat) <- c(object$xnames)
  colnames(coefmat) <- c("Estimate", "LB", "UB")
  return(coefmat)
}


#' Predictions based on the fitted parameter values
#'
#' Create a vector of predictions from a \code{cbq} object.
#'
#' @param object A \code{cbq} object.
#' @param data Data used for prediction.
#' @param ci Confidence interval. The default is 0.95.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A vector of predictions.
#' @export
#'
#' @method predict cbq
#'
predict.cbq <- function(object,data, ci = 0.95, ...){
    if (missing(data)){
        dat = object$data
    } else {
        dat = data
    }
    f <- Formula::Formula(object$formula)
    x  <- stats::model.matrix(f, dat)
    d = dim(x)[2]
    estimates = as.data.frame(object$stanfit)
    betas = as.matrix(estimates[,1:d])
    
    xb = x %*% t(betas)
    prob = c((1- ci)/2, 1 - (1- ci)/2)
    pred_mean = apply(xb, 1, mean)
    pred_lower = apply(xb, 1, quantile , probs = prob[1])
    pred_upper = apply(xb, 1, quantile , probs = prob[2])
    pred_prob_mean = pald(pred_mean, mu = 0, p = object$q, sigma = 1)
    pred_prob_lower = pald(pred_lower, mu = 0, p = object$q, sigma = 1)
    pred_prob_upper = pald(pred_upper, mu = 0, p = object$q, sigma = 1)
    
    out = cbind(pred_prob_mean, pred_prob_lower, pred_prob_upper)
    colnames(out) = c("mean", "lower", "upper")
    
    return(out)
    
}

#' Print cbq object
#'
#' General print function for \code{cbq} objects, which dispatches the chosen type
#' of printing to the corresponding function.
#'
#' @param x A \code{cbq} object to be printed.
#' @param type Character string giving the type of printing, such as
#'   \code{"text"}, \code{"mcmc"}, \code{"coef"}.
#' @param ... Additional arguments to be passed to print functions.
#'
#' @export
#' @return None.
#'
#'
#'
print.cbq <- function(x, type = "text", ...) {
  printFunName <- paste0("print_", type, ".cbq")
  do.call(printFunName, args = c(list(object = x), list(...)))
}


#' Print the main results from a \code{cbq} object.
#'
#' @param object A \code{cbq} object.
#' @param digits Number of digits to display.
#'
#' @export
#'
#'
#'
print_text.cbq <- function(object, digits = 3) {
  cat("Conditional binary quantile regression \n")
  cat("\nCall:\n",
      paste(deparse(object$Call), sep = "\n", collapse = "\n"),
      "\n\n",
      sep = "")
      if (object$vi == FALSE){
          cat(
          "MCMC run for",
          object$nsim,
          "iterations, with",
          object$stanfit@sim$warmup2,
          "used. \n\n"
          )
      } else {
          cat(
          "Variational approximation with",
          object$output_samples,
          "iterations. \n\n"
          )
      }
  
  cat("Coefficients:\n")
  print(round(coef(object), digits))
  cat("\n")
}


#' Print the mcmc results from a cbq object
#'
#' This prints a number of diagnostics about the results of a \code{cbq} objects
#'
#'
#' @param object A \code{cbq} object.
#' @param ... Additional arguments to be passed to the print function.
#'
#' @export
#' @return None.
#'
#'
print_mcmc.cbq <- function(object, ...) {
  print(object$stanfit, ...)
}



#' Print cbq coefficients
#'
#' @param object A \code{cbq} object.
#' @param digits Number of digits to display.
#'
#' @export
#' @return None.
#'
#'
print_coef.cbq <- function(object, digits = 3) {
  print(round(coef(object), digits))
}
