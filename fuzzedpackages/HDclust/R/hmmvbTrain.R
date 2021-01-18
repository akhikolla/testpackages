#' @title Train HMM-VB
#'
#' @description This function estimates parameters for HMM-VB
#' using the Baum-Welch algorithm. If the variable block structure
#' is not provided, the function will first find the structure by
#' a greedy search algorithm that minimizes BIC.
#' @param data A numeric vector, matrix, or data frame of observations.
#' Categorical values are not allowed. If a matrix or data frame, rows 
#' correspond to observations and columns correspond to variables.
#' @param VbStructure An object of class 'VB'. If supplied, variable block 
#' structure stored in VbStructure is used to train HMM-VB. If not 
#' provided, a search algorithm will be perfomed to find a variable block
#' structure with minimal BIC.
#' @param searchControl A list of control parameters for variable block structure
#' search. This parameter is ignored if variable block structure VbStructure is provided.
#' The defaults are set by the call \code{vbSearchControl()}.
#' @param trControl A list of control parameters for HMM-VB training algorithm.
#' The defaults are set by the call \code{hmmvbTrainControl()}.
#' @param nthread An integer specifying the number of threads used in searching and 
#' training routines.
#' @return An object of class 'HMMVB' providing estimation for HMM-VB.
#' The details of output components are as follows:
#' \item{VbStructure}{An object of class 'VB' with variable block structure for 
#' HMM-VB}
#' \item{HmmChain}{A list of objects of class 'HMM' with trained Hidden Markov Models
#' for each variable block.}
#' \item{diagCov}{A logical value indicating whether or not covariance matrices 
#' for mixture models are diagonal.}
#' \item{BIC}{BIC value for provided variable block structure or optimal BIC value
#' for found variable block structure.}
#' @seealso \code{\link{VB}}, \code{\link{vb}}, \code{\link{vbSearchControl}}, 
#' \code{\link{trainControl}}
#' @examples 
#' # Train HMM-VB with known variable block structure
#' data("sim3")
#' Vb <- vb(2, dim=40, bdim=c(10,30), numst=c(3,5), varorder=list(c(1:10),c(11:40)))
#' set.seed(12345)
#' hmmvb <- hmmvbTrain(sim3[,1:40], VbStructure=Vb)
#' show(hmmvb)
#' 
#' \donttest{
#' # Train HMM-VB with unknown variable block structure using default parameters
#' data("sim2")
#' set.seed(12345)
#' hmmvb <- hmmvbTrain(sim2[,1:5])
#' show(hmmvb)
#' 
#' # Train HMM-VB with unknown variable block structure using with ten permutations
#' # and several threads
#' data("sim2")
#' set.seed(12345)
#' hmmvb <- hmmvbTrain(sim2[1:5], searchControl=vbSearchControl(nperm=10), nthread=3)
#' show(hmmvb)}
#' @export 
hmmvbTrain <- function(data, VbStructure=NULL,
                       searchControl=vbSearchControl(), trControl=trainControl(),
                       nthread=1){
  data <- as.matrix(data)
  
  if (!is.null(VbStructure) && (!isS4(VbStructure) || !is(VbStructure, "VB")))
    stop('if provided VbStructure should be an instance of class VB!\n')
  
  nthread = as.integer(nthread)
  
  if ((length(nthread) != 1) || nthread < 1)
    stop('nthread should be a positive integer >= 1!\n')
  
  if (!is.null(VbStructure))
    if (dim(data)[2] != VbStructure@dim)
      stop('Dimensionality of data is different from dimensionality of variable block structure!')
  
  HmmVb <- rcpp_trainHmmVb(t(data), VbStructure, searchControl, trControl, nthread, vb, hmm, mkhmmvb)
  
  
  return(HmmVb)
  #newVbStructure, HmmChain,
  
}
