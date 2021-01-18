#' Generate the control variables used by the model 
#' 
#' @param tol.P : tolerance for parameters (numeric)
#' @param tol.L  : tolerance for log-likelihood  (numeric)
#' @param max.iter : maximum number of EM iterations (numeric)
#' @param SE.method : standard error estimation method  (char)
#' @param delta : increment used by forward difference  (numeric)
#' @param nknot : number of Gauss-Hermite quadrature knots  (numeric)
#' @return controlVals list

GenerateControlList <- function( control, ndim ) {

  callDepth = length(sys.calls());
  cdIndex = ifelse( callDepth > 1, callDepth - 1, 1)

  if (as.character(sys.calls()[[cdIndex]])[1] == "jmodelTM") {
    controlVals <- list(tol.P = 10 ^ (-3), tol.L = 10 ^ (-6), max.iter = 250, SE.method = 'PRES', delta = 10 ^ (- 5), 
                        nknot = 9)
  } else if ( as.character(sys.calls()[[cdIndex]])[1] == "jmodelMult") {
    controlVals <- list(tol.P = 10 ^ (-3), tol.L = 10 ^ (-6), max.iter = 250, SE.method = 'PRES', delta = 10 ^ (- 5), 
                        nknot = 11)
  } else if ( as.character(sys.calls()[[cdIndex]])[1] == "GenerateControlList") {
    warning("\n You called this function outside jmodelTM()/jmodelMult(); using jmodelTM() default values.\n") 
    controlVals <- list(tol.P = 10 ^ (-3), tol.L = 10 ^ (-6), max.iter = 250, SE.method = 'PRES', delta = 10 ^ (- 5), 
                        nknot = 9)
  } else {
    stop("GenerateControlList() is confused! Check your call stack.")
  }
  
  control <- c(control)
  namec <- names(control)
  if (length(uname <- namec[!namec %in% names(controlVals)]) > 0){
    warning("\n unknown names in 'control': ", paste(uname, collapse = ", "))
  }
  controlVals[namec] <- control
  if (controlVals$SE.method == 'PLFD' | controlVals$SE.method == 'PFDS') {
    controlVals$delta <- 10^(-3)
  }
  if (ndim == 2) {
    controlVals$nknot <- 7
  } else if (ndim >= 3) {
    controlVals$nknot <- 5
  }
  controlVals[namec] <- control

  return(controlVals)
  
}
