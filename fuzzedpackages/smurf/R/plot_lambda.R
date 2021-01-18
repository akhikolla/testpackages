###############################################
#
# Plot goodness-of-fit statistics or information criteria
#
###############################################


#' @export
#' @title Plot Goodness-of-Fit Statistics or Information Criteria
#' 
#' @description Function to plot the goodness-of-fit statistics or information criteria
#'              as a function of lambda when lambda is selected in-sample, out-of-sample or using cross-validation. 
#'             
#' @param x An object for which the extraction of goodness-of-fit statistics or information criteria is meaningful. 
#'          E.g. an object of class '\code{\link[=glmsmurf-class]{glmsmurf}}', typically the result of a call to \code{\link{glmsmurf}} or \code{\link{glmsmurf.fit}}.
#' @param xlab Label for the x-axis. The default value is \code{NULL} which means that \code{substitute(log(lambda))} is used when \code{log.lambda=TRUE}
#' and \code{substitute(lambda)} when \code{log.lambda=FALSE}.
#' @param ylab Label for the y-axis. The default value is \code{NULL} which means that the y-axis label is 
#'             determined based on method that was used to select lambda.
#' @param lambda.opt Logical indicating if the optimal value of lambda should be indicated on the plot 
#'                   by a vertical dashed line. Default is \code{TRUE}.
#' @param cv1se Logical indicating if the standard errors should be indicated on the plot 
#'               when cross-validation with the one standard error rule is performed (e.g. "cv1se.dev"). Default is \code{TRUE}.            
#' @param log.lambda Logical indicating if the logarithm of lambda is plotted on the x-axis, default is \code{TRUE}.
#' @param ... Additional arguments for the \code{\link[graphics:plot.default]{plot}} function.
#' 
#' @details This plot can only be made when lambda is selected in-sample, out-of-sample or using cross-validation (possibly with the one standard error rule), 
#'          see the \code{lambda} argument of \code{\link{glmsmurf}}.
#' 
#' @seealso \code{\link{glmsmurf}}, \code{\link{glmsmurf-class}} 
#' 
#' @example /inst/examples/Rent_example2.R
plot_lambda <- function(x, ...) UseMethod("plot_lambda", x)

#' @export
#' @rdname plot_lambda
plot_lambda.glmsmurf <- function(x, xlab = NULL, ylab = NULL, lambda.opt = TRUE, cv1se = TRUE, log.lambda = TRUE, ...) {
  
  # Check if wanted objects exist
  if (all(sapply(c("lambda.vector", "lambda.method", "lambda.measures"), exists, where = x))) {
    
    # Name of method used to select lambda
    lambda.method = x$lambda.method
    
    
    # New x-axis label, only used when xlab=NULL (i.e. no user input for xlab)
    if (is.null(xlab)) {
      
      if (log.lambda) {
        
        xlab <- substitute(log(lambda))
        
      } else {
        
        xlab <- substitute(lambda)
      }
    }
    
    
    # Get validation scores and potential new y-axis label
    if (lambda.method %in% c("cv.dev", "cv1se.dev", "oos.dev")) {
      
      if (lambda.method %in% c("cv.dev", "cv1se.dev")) {
        ylab.new <- "Deviance (CV)"
        
      } else {
        ylab.new <- "Deviance (out-of-sample)"
      }
     
      lambda.measures <- x$lambda.measures$dev
      
    } else if (lambda.method %in% c("cv.mse", "cv1se.mse", "oos.mse")) {
      
      if (lambda.method %in% c("cv.mse", "cv1se.mse")) {
        ylab.new <- "Mean Squared Error (CV)"
        
      } else {
        ylab.new <- "Mean Squared Error (out-of-sample)"
      }
      
      lambda.measures <- x$lambda.measures$mse
      
    } else if (lambda.method %in% c("cv.dss", "cv1se.dss", "oos.dss")) {
      
      if (lambda.method %in% c("cv.dss", "cv1se.dss")) {
        ylab.new <- "Dawid-Sebastiani Score (CV)"
        
      } else {
        ylab.new <- "Dawid-Sebastiani Score (out-of-sample)"
      }
      
      lambda.measures <- x$lambda.measures$dss

    } else if (lambda.method == "is.aic") {
      
      ylab.new <- "AIC (in-sample)"
      
      lambda.measures <- x$lambda.measures$aic
      
    } else if (lambda.method == "is.bic") {
     
      ylab.new <- "BIC (in-sample)"
      
      lambda.measures <- x$lambda.measures$bic
      
    } else if (lambda.method == "is.gcv") {
      
      ylab.new <- "GCV score (in-sample)"
      
      lambda.measures <- x$lambda.measures$gcv
      
    } else {
      stop("Invalid method to select lambda.")
    }

    # New y-axis label, only used if ylab=NULL (i.e. no user input for ylab)
    if (is.null(ylab)) {
      ylab <- ylab.new
    }
    
    
    if (log.lambda) {

      # Take logarithm of optimal value of lambda
      lambda.opt.num <- log(x$lambda)
      # Take logarithm of vector of lambda values
      lambda.vector <- log(x$lambda.vector)
      
    } else {
      
      # Optimal value of lambda
      lambda.opt.num <- x$lambda
      # Vector of lambda values
      lambda.vector <- x$lambda.vector
    }
    
    
    # Plot validation scores vs. values of lambda or log(lambda)
    plot(lambda.vector, rowMeans(lambda.measures), xlab = xlab, ylab = ylab, ...)
    
    
    # Add vertical dashed line indicating optimal value of lambda (or its logarithm)
    if (lambda.opt) {
      abline(v = lambda.opt.num, lty = 2)
    }
    
    # One standard error rule
    if (cv1se & lambda.method %in% c("cv1se.dev", "cv1se.mse", "cv1se.dss")) {
    
      # Index of lambda with minimum average error measure
      ind0 <- .which.min.last(rowMeans(lambda.measures))
      
      if (lambda.opt) {
        # Add vertical line indicating value of lambda (or its logarithm) with minimum average error measure
        abline(v = lambda.vector[ind0], lty = 3)
      }
      
      # Standard error of error measure per value of lambda
      cv.se <- apply(lambda.measures, 1L, sd) / sqrt(ncol(lambda.measures))
      
      # Add error bars
      segments(lambda.vector, rowMeans(lambda.measures) - cv.se, 
               lambda.vector, rowMeans(lambda.measures) + cv.se)
      
      if (lambda.opt) {
        # Add horizontal line indicating minimum average error measure + corresponding SE
        abline(h = mean(lambda.measures[ind0, ]) + cv.se[ind0], lty = 3)
      }
    }
    
  } else {
    warning("Validation scores are not available since a user-specified value of lambda was used.")
  }
}


