## documentation in ./man/VB-class.Rd
setClass("VB", slots=c(nb="integer",
                       dim="integer",
                       bdim="integer",
                       numst="integer",
                       varorder="list"))


#' @title Make an instance of class "VB"
#' @description This function creates a variable block structure.
#' @param nb The number of variable blocks.
#' @param dim Dimensionality of the data.
#' @param bdim An integer vector specifying dimensionality of each variable block.
#' This argument can be omitted if the variable block structure has a single
#' block (case of GMM).
#' @param numst An integer vector specifying the number of mixture models in each 
#' variable block.
#' @param varorder A list of integer vectors specifying the variable order in 
#' each variable block. This argument can be omitted if variable structure has a single
#' variable block (GMM).
#' @return An object of class "VB".
#' @seealso \code{\link{VB}}
#' @examples
#' # variable block structure for GMM with 3 dimensions and 2 mixture states
#' Vb <- vb(1, dim=3, numst=2)
#' 
#' # variable block structure with 2 variable blocks
#' Vb <- vb(2, dim=10, bdim=c(4,6), numst=c(3,11), varorder=list(c(1:4),c(5:10)))
#' @export
vb <- function(nb, dim, bdim=NULL,
               numst, varorder=NULL){
  
  nb = as.integer(nb)
  dim = as.integer(dim)
  numst = as.integer(numst)
  
  if ((nb <= 0) || (length(nb)!=1))
    stop("nb should be a scalar positive integer!\n")
  
  if ((dim <= 0) || (length(dim)!=1))
    stop("dim should be a scalar positive integer!\n")
  
  if (any(numst <= 0) || (length(numst) != nb))
    stop("numst should be should be a vector of length nb with positive integers!\n")
  
  # if variable block structure is GMM, simplify definition
  if (nb==1){
    bdim <- c(dim)
    varorder=list(1:dim)
  } else {
    if (is.null(bdim))
      stop("bdim should be provided if nb > 1!\n")
    
    bdim = as.integer(bdim)
    
    if (any(bdim <= 0) || (length(bdim) != nb))
      stop("bdim should be should be a vector of length nb with positive integers!\n")
  
    if (sum(bdim) != dim)
      stop("Total dimensionality in bdim doesn't match dim!\n")
  
    if (is.null(varorder))
      stop("varorder should be provided if nb > 1!\n")
    
    if (!is.list(varorder) || (length(varorder) != nb))
      stop("varorder should be a list with number of elements equal to nb!\n")
  
    # check if variable order is valid
    for (i in 1:nb){
      varorder[[i]] = as.integer(varorder[[i]])
      if (any(varorder[[i]] <= 0) || (length(varorder[[i]]) != bdim[i]))
        stop("Variable block ",i," is invalid! Each variable block should contain positive
             integers and have dimension equal to the value in bdim vector\n")
    }
  
    if (length(unique(unlist(varorder))) != length(unlist(varorder)))
      stop("Variable order list contains repeats!\n")
  
    if (any(unlist(varorder) > dim))
      stop("Variable order list has values larger than total dimensionality!\n")
  }

  return (new("VB", nb=nb, dim=dim, bdim=bdim, numst=numst, varorder=varorder))
}

## documentation in ./man/HMM-class.Rd
setClass("HMM", slots=c(dim="numeric", numst="numeric",
                                prenumst="numeric", a00="numeric",
                                a="matrix", mean="matrix", sigma="list",
                                sigmaInv="list", sigmaDetLog="numeric"))

#' @title Make an instance of "HMM" class.
#' @name hmm
#' @description This function creates a Hidden Markov Model on a variable block.
#' @param dim Dimensionality of the data in HMM.
#' @param numst An integer vector specifying the number of HMM states.
#' @param prenumst An integer vector specifying the number of states of previous 
#' variable block HMM.
#' @param a00 Probabilities of HMM states.
#' @param a Transition probability matrix from states in the previous variable block
#' to the states in the current one.  
#' @param mean A numerical matrix with state means. \emph{k}th row corresponds to the
#' \emph{k}th state.
#' @param sigma A list containing state covariance matrices.
#' @param sigmaInv A list containing state inverse covariance matrices.
#' @param sigmaDetLog A vector with \eqn{log(|sigma|)} for each state.
#' @return An object of class 'HMM'.
#' @keywords internal
hmm <- function(dim, numst, prenumst,
                a00, a, mean, sigma, sigmaInv, sigmaDetLog){
  if (!is.numeric(dim) || (dim <= 0) || !isTRUE(all.equal(dim, as.integer(dim))) || (length(dim)!=1))
    stop("dim should be a scalar positive integer!\n")
  
  if (!is.numeric(numst) || (numst <= 0) || !isTRUE(all.equal(numst, as.integer(numst))) || (length(numst)!=1))
    stop("numst should be a scalar positive integer!\n")
  
  if (!is.numeric(prenumst) || (prenumst <= 0) || !isTRUE(all.equal(prenumst, as.integer(prenumst))) || (length(prenumst)!=1))
    stop("prenumst should be a scalar positive integer!\n")
  
  if (!is.numeric(sigmaDetLog) || (length(sigmaDetLog)!=numst))
    stop("sigmaDetLog should be a vector of doubles with length numst!\n")
  
  if (!is.matrix(mean))
    stop("mean should be a matrix!\n")
  
  if (dim(mean)[1] != numst)
    stop("Number of rows in mean doesn't match with numst!\n")
  
  if (dim(mean)[2] != dim)
    stop("Number of cols in mean doesn't match with dim!\n")
  
  if (!is.matrix(a) || any(a < 0) || any(a > 1))
    stop("a should be a matrix with each element in range (0,1) !\n")
  
  for (i in 1:dim(a)[1]){
    if (!all.equal(sum(a[i,]), 1.0))
      message("Transition probabilities in matrix a are not properly normalized!\n")
  }
  
  if (dim(a)[1] != prenumst || dim(a)[2] != numst)
    stop("Transition probability matrix a has wrong dimensionality!\n")
  
  if (!is.numeric(a00) || any(a00 < 0) || any(a00 > 1))
    stop("Prior probability in a00 should be doubles in range (0,1) !\n")
  
  if (!all.equal(sum(a00), 1.0))
    message("Prior probabilities in a00 are not properly normalized!\n")
  
  if (length(a00) != numst)
    stop("Length of prior probabilities a00 doesn't match with numst!\n")
  
  if (!is.list(sigma) || (length(sigma) != numst))
    stop("sigma should be a list with number of elements in equal to numst!\n")
  
  for (i in 1:length(sigma)){
    if (!is.matrix(sigma[[i]]))
      stop("sigma", i, "should be a numeric matrix!\n")
    
    if (dim(sigma[[i]])[1] != dim || dim(sigma[[i]])[2] != dim)
      stop("Dimensions of sigma for state ", i, " don't match with dim!\n")
  }
  
  if (!is.list(sigmaInv) || (length(sigmaInv) != numst))
    stop("sigmaInv should be a list with number of elements in equal to numst!\n")
  
  for (i in 1:length(sigmaInv)){
    if (!is.matrix(sigmaInv[[i]]))
      stop("sigmaInv", i, "should be a numeric matrix!\n")
    
    if (dim(sigmaInv[[i]])[1] != dim || dim(sigmaInv[[i]])[2] != dim)
      stop("Dimensions of sigmaInv for state ", i, " don't match with dim!\n")
  }
  
  return (new("HMM", dim=dim, numst=numst,
               prenumst=prenumst, a00=a00,
               a=a, mean=mean, sigma=sigma, sigmaInv=sigmaInv, sigmaDetLog=sigmaDetLog))
}


## documentation in ./man/HMMVB-class.Rd
setClass("HMMVB", slots=c(VbStructure="VB", HmmChain="list", BIC="numeric", diagCov="logical", Loglikehd="numeric"))

#' @title Make an instance of "HMMVB" class.
#' @name mkhmmvb
#' @description This function creates an instance of "HMMVB" class. The function is called inside 
#' C code during HMMVB model training. It is not meant to be called by users.
#' @param VbStructure An object of class 'VB' that contains variable block structure.
#' @param HmmChain A list of objects of class 'HMM' with trained Hidden Markov Models
#' for each variable block.
#' @param diagCov A logical value indicating whether or not covariance matrices 
#' for mixture models are diagonal.
#' @param BIC BIC value for provided variable block structure or optimal BIC value
#' for found variable block structure.
#' @param Loglikehd Loglikelihood value for each data point
#' @return An instance of class "HMMVB"
#' @keywords internal
mkhmmvb <- function(VbStructure, HmmChain, BIC, diagCov, Loglikehd){
  
  if (length(HmmChain) != VbStructure@nb)
    stop("Number of HMM models doesn't match nb in variable block structure!\n")
  
  for (i in 1:length(HmmChain)){
    if (HmmChain[[i]]@dim != VbStructure@bdim[i])
      stop("Dimensionality of HMM model",i,"doesn't match bdim in variable block structure!\n")
    
    if (HmmChain[[i]]@numst != VbStructure@numst[i])
      stop("Number of states of HMM model",i,"doesn't match numst in variable block structure!\n")
  }
  
  if (!is.numeric(BIC) ||  (length(BIC)!=1))
    stop("BIC should be a scalar double")
  
  if (!is.logical(diagCov) ||  (length(diagCov)!=1))
    stop("diagCov should be a scalar logical")
  
  return (new("HMMVB", VbStructure=VbStructure, HmmChain=HmmChain, BIC=BIC, diagCov=diagCov, Loglikehd=Loglikehd))
}

## documentation in ./man/HMMVBBIC-class.Rd
setClass("HMMVBBIC", slots=c(numst="numeric",
                               BIC="numeric",
                               optHMMVB="HMMVB"))

## documentation in ./man/HMMVBclust-class.Rd
setClass("HMMVBclust", slots=c(data="matrix",
                                              clustParam="list",
                                              clsid="numeric",
                                              size="numeric",
                                              Loglikehd="numeric"))
