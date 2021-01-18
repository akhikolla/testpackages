#' @title Find density modes with HMM-VB
#'
#' @description This function finds the density modes with HMM-VB. First, for each data point it finds an optimal state sequence
#' using Viterbi algorithm. Next, it uses Modal Baum-Welch algorithm (MBW) to find the modes
#' of distinct Viterbi state sequences. Data points associated the same modes form clusters. 
#' @param data A numeric vector, matrix, or data frame of observations.
#' Categorical values are not allowed. If a matrix or data frame, rows 
#' correspond to observations and columns correspond to variables.
#' @param model An object of class 'HMMVB' that contains trained HMM-VB obtained 
#' by the call to function \code{hmmvbTrain}. 
#' @param nthread An integer specifying the number of threads used in clustering.
#' @param bicObj An object of class 'HMMVBBIC' which stores results of model selection. 
#' If provided, argument \code{model} is ignored.
#' @return An object of class 'HMMVBclust'.
#' @seealso \code{\link{HMMVB-class}}, \code{\link{HMMVBclust-class}}, \code{\link{hmmvbTrain}}
#' @examples
#' # find modes using trained HMM-VB
#' Vb <- vb(1, dim=4, numst=2)
#' set.seed(12345)
#' hmmvb <- hmmvbTrain(iris[,1:4], VbStructure=Vb)
#' modes <- hmmvbFindModes(iris[,1:4], model=hmmvb)
#' show(modes)
#' 
#' \donttest{
#' # find modes using HMMVBBIC object obtained in model selection
#' Vb <- vb(1, dim=4, numst=1)
#' set.seed(12345)
#' modelBIC <- hmmvbBIC(iris[,1:4], VbStructure=Vb)
#' modes <- hmmvbClust(iris[,1:4], bicObj=modelBIC)
#' show(modes)}



hmmvbFindModes <- function(data, model=NULL, nthread=1, bicObj=NULL){
  data = as.matrix(data)
  
  nthread = as.integer(nthread)
  
  if ((length(nthread) != 1) || nthread < 1)
    stop('nthread should be a positive integer >= 1!\n')
  
  if (!is.null(bicObj) && (!isS4(bicObj) || !is(bicObj, "HMMVBBIC")))
    stop('If provided, bicObj should be an instance of class HMMVBBIC!\n')
  
  if (!is.null(bicObj)){
    
    res <- rcpp_findModes(t(data), getOptHMMVB(bicObj), nthread)
    
  }  else{
    
    if (is.null(model) || !isS4(model) || !is(model, "HMMVB"))
      stop('If argument bicObj is not provided, model argument should be given and it must be an instance of class HMMVB!\n')
    
    res <- rcpp_findModes(t(data), model, nthread)
  }
  
  res@data <- data
  res@clsid = res@clsid + 1
  
  return(res)
  
}
