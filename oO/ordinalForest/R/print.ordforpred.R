#' @keywords internal
#' @export
print.ordforpred <-
  function(x, ...) {
    
    cat("\n")
    cat(paste("Predicted values of ", length(x$ypred), " observation(s).", sep=""), "\n")
    
    cat("\n")
    cat("Classes of ordinal target variable:", "\n")
    cat(paste(paste("\"", levels(x$ypred), "\"", sep=""), collapse=", "), "\n")
    
  }
