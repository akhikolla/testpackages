
DDC <- function(X, DDCpars = list()){ 
  # X is the input data, and must be a matrix or a data frame.
  #   It must always be provided.
  #
  # DDCpars contains all the options chosen by the user, but has
  # a default:
  #
  
  if (!is.data.frame(X) & !is.matrix(X)) {
    stop("The input data must be a matrix or a data frame.")
  }
  if (ncol(X) < 2) stop("The input data must have at least 2 columns.")
  
  
  if (is.null(DDCpars)) {
    DDCpars <- list()
  }
  if (!is.list(DDCpars)) {
    stop("DDCpars must be a list")
  }

  # parameters for checkDataSet
  if (!"fracNA" %in% names(DDCpars)) {
    DDCpars$fracNA <- 0.5
  }
  if (!"numDiscrete" %in% names(DDCpars)) {
    DDCpars$numDiscrete <- 3
  }
  
  # general parameters
  if (!"precScale" %in% names(DDCpars)) {
    DDCpars$precScale <- 1e-12
  }
  if (!"cleanNAfirst" %in% names(DDCpars)) {
    DDCpars$cleanNAfirst = "automatic"
  }
  if (!"tolProb" %in% names(DDCpars)) {
    DDCpars$tolProb <- 0.99
  }
  if (!"corrlim" %in% names(DDCpars)) {
    DDCpars$corrlim <- 0.5
  }
  if (!"combinRule" %in% names(DDCpars)) {
    DDCpars$combinRule <- "wmean"
  }
  if (!"returnBigXimp" %in% names(DDCpars)) {
    DDCpars$returnBigXimp <- FALSE
  }
  if (!"silent" %in% names(DDCpars)) {
    DDCpars$silent <- FALSE
  }
  if (!"nLocScale" %in% names(DDCpars)) {
    DDCpars$nLocScale = 25000
  }
  
  if (!"fastDDC" %in% names(DDCpars)) {
    if (dim(X)[2] <= 750) {
      fastDDC <- 0
      DDCpars$fastDDC <- FALSE
    } else {
      fastDDC <- 1
      DDCpars$fastDDC <- TRUE
    }
  } else {
    fastDDC <- DDCpars$fastDDC + 0
  }
  if (DDCpars$fastDDC) {
    if (dim(X)[2] < 500) {
      fastDDC <- 0
    }
  } else {
    if (dim(X)[2] > 2000) {
      cat("Consider using the option 'fastDDC == TRUE' which runs much faster on datasets with many variables.")
    }
  }
  
  if (!"standType" %in% names(DDCpars)) {
    DDCpars$standType <- "1stepM"
  }
  if (!"corrType" %in% names(DDCpars)) {
    DDCpars$corrType <- "gkwls"
  }
  
  if (!"transFun" %in% names(DDCpars)) {
    DDCpars$transFun <- "wrap"
  }
  if (!"nbngbrs" %in% names(DDCpars)) {
    DDCpars$nbngbrs <- 100
  }
  
  # Other parameters
  if (!"coreOnly" %in% names(DDCpars)) {
    DDCpars$coreOnly <- FALSE
  }
  if (!"tolProbCell" %in% names(DDCpars)) {
    DDCpars$tolProbCell <- DDCpars$tolProb
  }
  if (!"tolProbRow" %in% names(DDCpars)) {
    DDCpars$tolProbRow <- DDCpars$tolProb
  }
  if (!"tolProbReg" %in% names(DDCpars)) {
    DDCpars$tolProbReg <- DDCpars$tolProb
  }
  if (!"tolProbCorr" %in% names(DDCpars)) {
    DDCpars$tolProbCorr <- DDCpars$tolProb
  }
  if (!"includeSelf" %in% names(DDCpars)) {
    DDCpars$includeSelf <- 1
  }
  if (!"numiter" %in% names(DDCpars)) {
    DDCpars$numiter <- 1
  }
  if (!"qdim" %in% names(DDCpars)) {
    DDCpars$qdim <- 30
  }
  if (!"nCorr" %in% names(DDCpars)) {
    DDCpars$nCorr = 1000
  }
  
 # Put options in convenient order
  
  DDCpars <- DDCpars[c("fracNA", "numDiscrete", "precScale", "cleanNAfirst",
                       "tolProb", "corrlim", "combinRule", "returnBigXimp",
                       "silent", "nLocScale", "fastDDC", "standType", 
                       "corrType", "transFun", "nbngbrs",
                       "coreOnly", "tolProbCell", "tolProbRow", "tolProbReg",
                       "tolProbCorr","includeSelf", "numiter", "qdim", "nCorr")]
  
  
  # Convert option parameters to integers
  transFun   <- match(DDCpars$transFun, c("Huber", "wrap", "rank"),
                             nomatch = 2)
  standType  <- match(DDCpars$standType, 
                              c("1stepM","hubhub","wrap","mcd", "rawmcd", "wrapmedmad"),
                             nomatch = 1) - 1
  corrType   <- match(DDCpars$corrType, c("wrap", "rank", "gkwls"),
                             nomatch = 3)
  combinRule <- match(DDCpars$combinRule,
                              c("wmean","wmedian","mean","median"),
                               nomatch = 1)
  
  # Retrieve parameters from the list:
  #
  fracNA        <- DDCpars$fracNA
  numDiscrete   <- DDCpars$numDiscrete
  precScale     <- DDCpars$precScale
  returnBigXimp <- DDCpars$returnBigXimp
  silent        <- DDCpars$silent
  cleanNAfirst  <- DDCpars$cleanNAfirst

  if (!DDCpars$coreOnly) {
    # Check the data set and set aside columns and rows that do
    # not satisfy the conditions:
    out <- checkDataSet(X, fracNA, numDiscrete, precScale, silent, cleanNAfirst)
    
    checkedData <- out$remX + 0 # +0 puts it in new memory, so cpp doesn't change out$remX
  } else {# execute core only, for use with MacroPCA function
    checkedData <- X + 0
    out <- c()
  }
  
  if (fastDDC) {
    goodCols   <- which(colSums(is.na(checkedData)) / dim(checkedData)[1] <= 0.5) - 1 #-1 to convert to C++ indexing
    if (length(goodCols) < dim(checkedData)[2]) {
      cat(paste("Note that ", dim(checkedData)[2] - length(goodCols), " out of ", dim(checkedData)[2], 
                " columns in remX have over 50% of NAs. These columns will be considered standAlone variables."))
    }
  } else {
    goodCols <- (1:dim(checkedData)[2]) - 1
  }
  
  # Check nCorr  params
  if (DDCpars$nCorr == 0) {
    DDCpars$nCorr = dim(checkedData)[1]
  } else {
    DDCpars$nCorr = min(DDCpars$nCorr, dim(checkedData)[1])
  }
  
  # Carry out the actual DetectDeviatingCells algorithm on
  # the remaining dataset out1$remX :
  
  res <- tryCatch( .Call('_cellWise_DDC_cpp', checkedData, DDCpars$tolProbCell,
                         DDCpars$tolProbRow, DDCpars$tolProbReg, DDCpars$tolProbCorr,
                         DDCpars$corrlim, combinRule,
                         DDCpars$includeSelf,
                         fastDDC, DDCpars$qdim, transFun,
                         DDCpars$nbngbrs, DDCpars$numiter,
                         DDCpars$precScale, standType, corrType,
                         DDCpars$nCorr, DDCpars$nLocScale, goodCols,
                         PACKAGE = 'cellWise'),
                   "std::range_error" = function(e){
                     conditionMessage( e ) })
  
  
  loc.out <- drop(res$locX)
  scale.out <- drop(res$scaleX)
  deshrinkage.out <- drop(res$deshrinkage)
  scalestres.out <- drop(res$scalestres)
  # 
  names(loc.out) <- names(scale.out) <- colnames(checkedData)
    names(scalestres.out)  <- names(deshrinkage.out) <- colnames(checkedData)[res$colConnected]
  
  returned.result <- list(locX = loc.out, scaleX = scale.out, Z = res$Z,
                          nbngbrs = res$k, ngbrs = res$ngbrs, robcors = res$robcors,
                          robslopes = res$robslopes, deshrinkage = deshrinkage.out,
                          Xest = res$Xest,scalestres = scalestres.out,
                          stdResid = res$stdResid,
                          indcells = drop(res$indcells), Ti = res$Ti, 
                          medTi = res$medTi, madTi = res$madTi, 
                          indrows = drop(res$indrows), indall = res$indall,
                          indNAs = res$indNAs, 
                          Ximp = res$Ximp)
  if (!DDCpars$coreOnly) {
    if (returnBigXimp) {
      Ximp <- data.matrix(X, rownames.force = TRUE)
      Ximp[out$rowInAnalysis, out$colInAnalysis] <- out$Ximp
      out$Ximp <- Ximp
    }
  }
  
  DDCpars <- DDCpars[1:15] 
  return(c(list(DDCpars = DDCpars), out, returned.result))
}





