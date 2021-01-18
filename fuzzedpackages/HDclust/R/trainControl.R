#' @title Parameters for HMM-VB training.
#' @name trainControl
#' @description This function creates a list with parameters for estimating an HMM-VB,
#' which is used as an argument for \code{hmmvbTrain}.
#' @param ninit0 The number of initializations for default scheme 0, under
#' which the k-means clustering for entire dataset is used to initialize the model.
#' @param ninit1 The number of initializations for default scheme 1, under
#' which the k-means clustering for a subset of data is used to initialize the model.
#' @param ninit2 The number of initializations for default scheme 2, under
#' which a random subset of data is used as cluster centroids to initialize the model.
#' @param epsilon Stopping criteria for Baum-Welch algorithm. Should be a small number in 
#' range (0,1).
#' @param diagCov A logical value indicating whether or not variable block covariance matrices
#' will be diagonal.
#' @return The named list with parameters.
#' @seealso \code{\link{hmmvbTrain}}
#' @examples 
#' # setting up multiple initialization schemes
#' Vb <- vb(1, dim=4, numst=2)
#' set.seed(12345)
#' hmmvb <- hmmvbTrain(iris[,1:4], VbStructure=Vb, 
#'           trControl=trainControl(ninit0=2, ninit1=2, ninit2=2))
#' show(hmmvb)
#' 
#' # forcing diagonal covariance matrices
#' Vb <- vb(1, dim=4, numst=2)
#' set.seed(12345)
#' hmmvb <- hmmvbTrain(iris[,1:4], VbStructure=Vb, 
#'           trControl=trainControl(diagCov=TRUE))
#' show(hmmvb)
trainControl <- function(ninit0=1, ninit1=0, ninit2=0, 
                              epsilon=1.0e-4, diagCov=FALSE){
  ninit0 = as.integer(ninit0)
  ninit1 = as.integer(ninit1)
  ninit2 = as.integer(ninit2)
  
  if ((length(ninit0)!=1) || (ninit0 < 0))
    stop('Invalid ninit0 argument provided! ninit0 should be a scalar non-negative integer\n')
  
  if ((length(ninit1)!=1) || (ninit1 < 0))
    stop('Invalid ninit0 argument provided! ninit0 should be a scalar non-negative integer\n')
  
  if ((length(ninit2)!=1) || (ninit2 < 0))
    stop('Invalid ninit0 argument provided! ninit0 should be a scalar non-negative integer\n')
  
  if ((!is.numeric(epsilon)) || (length(epsilon)!=1) || (epsilon <= 0) || (epsilon >= 1))
    stop('Invalid epsilon argument provided! epsilon should be a scalar positive double 
         in range (0,1)\n')
  
  if (ninit0 == 0 && ninit1 == 0 && ninit2 == 0){
    message("All initialization parameters ninit0, ninit1 and ninit2 are 0! Thus ninit0 is set to 1\n")
    ninit0 = 1
  }
  
  if (!is.logical(diagCov) || length(diagCov)!=1)
    stop('Invalid diagCov argument provided! diagCov should be a scalar logical\n')
  
  return(list(ninit0=ninit0, ninit1=ninit1, ninit2=ninit2,
              epsilon=epsilon, diagCov=diagCov))
}
