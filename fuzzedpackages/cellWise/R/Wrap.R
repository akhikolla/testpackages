

wrap <- function(X, locX = NULL, scaleX = NULL,
                 precScale = 1e-12, imputeNA = TRUE,
                 checkPars = list()) {
  # Transforms multivariate data X using wrapping function
  # with b = 1.5 and c = 4 and loc&scale in
  # the input arguments locX and scaleX
  #
  # args:
  #   X: data matrix
  #   locX: vector with location estimates for every column of X
  #   scaleX: vector with scale estimates for every column of X
  #   precScale: precision scale used throughout the algorithm
  #   imputeNA: whether or not to impute the NAs
  # Returns: 
  #   Xw: a matrix with wrapped data
  #   colInWrap: vector of indices of wrapped columns
  #
  
  # Check inputs
  
  if (!is.data.frame(X) & !is.matrix(X) & !is.vector(X)) {
    stop("The input data must be a vector, matrix or a data frame")
  }
  isVec = FALSE
  if (is.vector(X)) {
    X = matrix(X, ncol = 1)
    isVec = TRUE
  } else {
    X <- as.matrix(X)
  }
  
  
  if (precScale < 1e-12) {
    precScale <- 1e-12
    warning("The \"precScale\" argument was changed to 1e-12.")
  }
  
  # parameters for checkDataSet
  if (!"coreOnly" %in% names(checkPars)) {
    checkPars$coreOnly <- FALSE
  }
  if (!"silent" %in% names(checkPars)) {
    checkPars$silent <- FALSE
  }
  if (!"numDiscrete" %in% names(checkPars)) {
    checkPars$numDiscrete <- 5
  }
  if (!"precScale" %in% names(checkPars)) {
    checkPars$precScale <- 1e-12
  }
  
  if (!is.null(locX) & !is.null(scaleX)) { # custom loc/scale provided
    if (length(locX) != dim(X)[2] || length(scaleX) != dim(X)[2]) {
      stop(paste0("The arguments \"locX\" and \"scaleX\" should both be of length ",dim(X)[2]))
    }
    colInWrap <- which(scaleX > precScale)
    if (length(colInWrap) == 0) {
      stop(paste("No columns with scale > precScale were found."))
    }
  } else {
    # Start with checkDataSet to avoid transforming too discrete variables
    if (!checkPars$coreOnly) {
      # Check the data set and set aside columns and rows that do
      # not satisfy the conditions:
      out <- checkDataSet(X,
                          fracNA = 1,
                          numDiscrete = checkPars$numDiscrete,
                          precScale = checkPars$precScale,
                          silent = checkPars$silent)
      
      colInWrap <- out$colInAnalysis 
    } 
    locScale.out <- estLocScale(X[, , drop = FALSE], silent = TRUE)
    locX   <- locScale.out$loc
    scaleX <- locScale.out$scale
  }
 
  
  res <- tryCatch(.Call("_cellWise_Wrap_cpp", X[, colInWrap, drop = FALSE], 
                        locX[colInWrap], scaleX[colInWrap], precScale,
                        PACKAGE = "cellWise"), 
                  `std::range_error` = function(e) {
                    conditionMessage(e)
                  })
  if (length(colInWrap) < dim(X)[2]) {
    warning(paste(dim(X)[2] - length(colInWrap), " variable(s) had a scale <= precScale and were left out.\n\n                    The indices of the remaining columns are in $colInWrap."))
  }
  if (!imputeNA) { # if DON'T impute (by default the C++ does impute):
    if (any(is.na(X[, colInWrap]))) {
      naInds <- which(is.na(X[, colInWrap]))
      res$Xw[naInds] <- NA
    }
  }
  if (isVec) { res$Xw = as.vector(res$Xw) }
  return(list(Xw = res$Xw, colInWrap = colInWrap, loc = locX, scale = scaleX))
}