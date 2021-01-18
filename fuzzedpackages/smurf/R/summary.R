###############################################
#
# Model summary
#
###############################################


# Print function for glmsmurf object which calls summary directly
#' @export
#' @method print glmsmurf
print.glmsmurf <- function(x, ...) {
  
  return(summary.glmsmurf(object = x, digits = 3L))
}



#' @export
#' @title Summary of a Multi-Type Regularized GLM Fitted Using the SMuRF Algorithm
#' 
#' @description Function to print a summary of a \code{glmsmurf}-object.
#' 
#' @param object An object of class '\code{\link[=glmsmurf-class]{glmsmurf}}', typically the result of a call to \code{\link{glmsmurf}} or \code{\link{glmsmurf.fit}}.
#' @param digits The number of significant digits used when printing, default is 3.
#' @param ... Additional arguments which are currently ignored.
#' 
#' @seealso \code{\link[stats]{summary.glm}}, \code{\link{glmsmurf}}, \code{\link{glmsmurf-class}}
#' 
#' @examples ## See example(glmsmurf) for examples
#' 
#' @method summary glmsmurf
summary.glmsmurf <- function(object, digits = 3L, ...) {
  
  # Handle reest
  reest <-  exists("coefficients.reest", object)
  
  # Print call, deparse turns everything into strings
  cat("\nCall:  ", paste(deparse(object$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  
  # Deviance residuals
  if (reest) {
    cat("Deviance residuals of estimated model:\n")
  } else {
    cat("Deviance residuals:\n")
  }
  # Summary of deviance residuals
  print(format(summary(residuals.glmsmurf(object = object, type = "deviance"))[-4], digits = digits), 
        quote = FALSE)
  
  if (reest) {
    cat("\n")
    cat("Deviance residuals of re-estimated model:\n")
    # Summary of deviance residuals of re-estimated model
    print(format(summary(residuals_reest.glmsmurf(object = object, type = "deviance"))[-4], digits = digits), 
          quote = FALSE)
  }
  cat("\n\n")
  
  # Coefficients
  cat("Coefficients:\n")
  
  # Convert to matrix
  mat <- as.matrix(coef(object))
  colnames(mat) <- ""
  
  # Re-estimated coefficients
  if (reest) {
    mat <- cbind(mat, coef_reest(object))
    colnames(mat) <- c("Estimated", "Re-estimated")
  } 
  
  # Format coefficients
  tmp <- format(mat, digits = digits)
  # Replace all zero elements by *
  tmp[mat == 0] <- " *"
  # Print coefficients, remove quotes
  print(tmp, quote = FALSE)
  
  # Detail meaning of '*'
  cat("---", "\n '*' indicates a zero coefficient", "\n\n")
  
  # Null deviance
  cat("\nNull deviance:", format(object$null.deviance, nsmall = digits, digits = digits), "on", 
      object$df.null, "degrees of freedom")
  ##########################
  # Estimated model
  if (reest) cat("\n-------------------", "\nEstimated model:\n")
  # Residual deviance
  cat("\nResidual deviance:", format(object$deviance, nsmall = digits, digits = digits), "on", 
      object$df.residual, "degrees of freedom")
  # AIC, BIC, GCV score
  cat("\nAIC: ", format(object$aic, nsmall = digits, digits = digits), 
      "; BIC: ", format(object$bic, nsmall = digits, digits = digits), 
      "; GCV score: ", format(object$gcv, nsmall = digits, digits = digits), sep="")
  # Penalized minus log-likelihood
  cat("\nPenalized minus log-likelihood:", object$obj.fun)
  ##########################
  # Re-estimated model (when present)
  if (reest) {
    cat("\n-------------------", "\nRe-estimated model:\n")
    # Residual deviance
    cat("\nResidual deviance:", format(object$deviance.reest, nsmall = digits, digits = digits), "on", 
        object$df.residual.reest, "degrees of freedom")
    # AIC, BIC, GCV score
    cat("\nAIC: ", format(object$aic.reest, nsmall = digits, digits = digits), 
        "; BIC: ", format(object$bic.reest, nsmall = digits, digits = digits), 
        "; GCV score: ", format(object$gcv.reest, nsmall = digits, digits = digits), sep="")
    # Objective function
    cat("\nObjective function:", object$obj.fun.reest)
    cat("\n-------------------")
  }
  
  cat("\nlambda: ", object$lambda, "; lambda1: ", object$lambda1, "; lambda2: ", object$lambda2, sep = "")
  cat("\nNumber of iterations:", object$iter)
  cat("\nFinal step size:", object$final.stepsize)
  
  conv <- object$converged
  conv_message <- character(length(conv))
  conv_message[conv == 0] <- "Succesful convergence"
  conv_message[conv == 1] <- "Maximum number of iterations reached"
  conv_message[conv == 2] <- "Two subsequent restarts"
  conv_message[conv == 3] <- "Low step size"
  
  cat("\nConvergence:", paste(conv_message, collapse = "; "))
}
