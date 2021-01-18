#' Simulation from dependence models
#'
#' Simulate Monte Carlo sample from a collection of fitted conditional
#' dependence models.
#'
#' Generates a Monte Carlo sample of the required size from a collection of
#' conditional multivariate extreme values model of Heffernan and Tawn, 2004.
#' For each marginal variable, the model that conditions on that margin is used
#' to simulate values in the part of the sample space for which that margin is
#' the largest of all marignal variables (measured on a quantile scale).
#'
#' @usage mexMonteCarlo(nSample,mexList,mult=10)
#' @param nSample Required sample size.
#' @param mexList List of fitted dependence models (returned by
#' \code{\link{mexAll}}).
#' @param mult Integer specifying what multiple of the total number of points
#' should be generated for rejection sample
#' @return A list with the following components:
#'
#' \item{nR}{For each margin, number of original Monte Carlo points replaced by
#' points generated under the corresponding conditional model.}
#' \item{MCsample}{Matrix contiaining the Monte Carlo sample, dimension
#' \code{nSample} by dimension of original dataset.} \item{whichMax}{Vector of
#' indices indicating which variable is largest (on the quantile scale)}
#' \item{whichMaxAboveThresh}{Logical vector indicating which of the variables
#' identified by \code{whichMax} are additionally above the corresponding
#' threshold for dependence estimation.}
#' @author Harry Southworth, Janet E. Heffernan
#' @references J. E. Heffernan and J. A. Tawn, A conditional approach for
#' multivariate extreme values, Journal of the Royal Statistical society B, 66,
#' 497 -- 546, 2004
#' @keywords models multivariate
#' @examples
#' \donttest{
#'   mAll <- mexAll(winter,mqu=0.7,dqu=c(0.7,0.7,0.7,0.7,0.7))
#'   mexMC <- mexMonteCarlo(5000,mAll)
#'   pairs(mexMC$MCsample)
#' }
#' @export mexMonteCarlo
mexMonteCarlo <- function(nSample, mexList, mult=10){
#set up
  d <- length(mexList)
  data <- mexList[[1]]$margins$data
  margins <- mexList[[1]]$dependence$margins
  nData <- dim(data)[1]
#Generate our Monte Carlo sample from the original dataset.
  which <- sample(1:nData, size=nSample, replace = TRUE)
  MCsampleOriginal <- data[which,]
#Transform the original dataset to the Laplace scale.
  dataLaplace <- mexTransform(mexList[[1]]$margins, margins = margins, method = "mixture")$transformed
  MCsampleLaplace <- dataLaplace[which,]
#Identify maximum components.
  whichMax <- apply(MCsampleLaplace,1,which.max)
# Now identify which of the maximal components lie above their associated conditional dependence model thresholds.
  dth <- sapply(mexList,function(l)l$dependence$dth)
  dqu <- sapply(mexList,function(l)l$dependence$dqu)
  whichMaxAboveThresh <- sapply(1:nSample,function(i)MCsampleLaplace[i,whichMax[i]] >= dth[whichMax[i]])
# generate large samples from each of the Conditional models,
  mexKeep <- lapply(1:d,function(i){
    mc <- predict.mex(mexList[[i]],pqu=dqu[i],nsim=nSample*d*mult)
    mc$data$simulated[mc$data$CondLargest,order(c(i,c(1:d)[-i]))]
  })
# Replace original sample by samples from conditional models.
  nR <- rep(0,d)
  names(nR) <- names(data)
  for(i in 1:d){
      replace <- whichMax == i & whichMaxAboveThresh
      nReplace <- sum(replace)
      if(nReplace > 0){
         nR[i] <- nReplace
         MCsampleOriginal[replace,] <- as.matrix(mexKeep[[i]])[1:nReplace,]
      }
  }

  res <- list(nR=nR, MCsample=MCsampleOriginal, whichMax=whichMax, whichMaxAboveThresh=whichMaxAboveThresh)
  oldClass(res) <- "mexMC"
  res
}

#' @rdname mex
#' @export
mexAll <- function(data,mqu,dqu) {
    d <- dim(data)[2]
    res <- lapply(1:d,function(i){mex(data,which=i,mqu=mqu,dqu=dqu[i])})
    names(res) <- names(data)
    oldClass(res) <- "mexList"
    res
}

#' @rdname mex
#' @export
print.mexList <- function(x, ...){
    print(x[[1]])
    for(i in 2:length(x)){
        cat("\n______\n")
        print(x[[i]][[2]])
    }
    invisible(x)
}
