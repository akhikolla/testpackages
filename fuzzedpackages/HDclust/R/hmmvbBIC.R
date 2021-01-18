#' @title BIC for HMM-VB
#'
#' @description This function finds an optimal number of mixture components (states) for 
#' HMM-VB using the Bayesian Information Criterion (BIC). The variable block 
#' structure is provided as input and then BIC is estimated for HMM-VB with 
#' different configurations of states for the variable blocks.
#' @param data A numeric vector, matrix, or data frame of observations.
#' Categorical values are not allowed. If a matrix or data frame, rows 
#' correspond to observations and columns correspond to variables.
#' @param VbStructure An object of class 'VB'. Variable block 
#' structure stored in VbStructure is used to train HMM-VB model. \code{numst} parameter
#' of the variable block structure is ignored.
#' @param configList A list of integer vectors specifying number of states in each variable
#' block for which BIC is to be calculated.
#' @param numst An integer vector specifying the numbers of mixture components (states) in 
#' each variable block for which BIC is to be calculated. Number of states is the same for 
#' all variable blocks. The argument is ignored if \code{configList} argument is provided.
#' @param trControl A list of control parameters for HMM-VB training algorithm.
#' The defaults are set by the call \code{hmmvbTrainControl()}.
#' @param nthread An integer specifying the number of threads used in searching and 
#' training routines.
#' @return A named list with estimated BIC values and the number of states or state configurations
#' for which BIC was calculated.
#' @seealso \code{\link{VB}}, \code{\link{vb}}, \code{\link{trainControl}}
#' @examples 
#' \donttest{
#' # Default search for the optimal number of states for HMM-VB model 
#' data("sim3")
#' Vb <- vb(2, dim=40, bdim=c(10,30), numst=c(1,1), varorder=list(c(1:10),c(11:40)))
#' set.seed(12345)
#' hmmvbBIC(sim3[1:40], VbStructure)
#' 
#' # Search for the optimal number of states for HMM-VB model using 
#' # provided values for the number of states 
#' data("sim3")
#' Vb <- vb(2, dim=40, bdim=c(10,30), numst=c(1,1), varorder=list(c(1:10),c(11:40)))
#' set.seed(12345)
#' hmmvbBIC(sim3[1:40], VbStructure=Vb, numst=c(2L, 4L, 6L))
#'
#' # Search for the optimal number of states for HMM-VB model using 
#' # provided configurations of the number of states 
#' data("sim3")
#' Vb <- vb(2, dim=40, bdim=c(10,30), numst=c(1,1), varorder=list(c(1:10),c(11:40)))
#' set.seed(12345)
#' configs = list(c(1,2), c(3,5), c(6,7))
#' hmmvbBIC(sim3[1:40], VbStructure=Vb, configList=configs)}
#' @export 
hmmvbBIC <- function(data, VbStructure, configList=NULL, numst = 1:10,
                       trControl=trainControl(),
                       nthread=1){
  data <- as.matrix(data)
  
  if ((!isS4(VbStructure) || !is(VbStructure, "VB")))
    stop('VbStructure should be an instance of class VB!\n')
  
  if (dim(data)[2] != VbStructure@dim)
    stop('Dimensionality of data is different from dimensionality of variable block structure!')
  
  nthread = as.integer(nthread)
  
  if ((length(nthread) != 1) || nthread < 1)
    stop('nthread should be a positive integer >= 1!\n')
  
  if (!is.null(configList)){
    if (!is.list(configList))
      stop('if provided, configList should be a list!\n')
    
    for (i in 1:length(configList)){
      if (!is.numeric(configList[[i]]) || !isTRUE(all.equal(configList[[i]], as.integer(configList[[i]]))) 
          || (any(configList[[i]] <= 0)) || length(configList[[i]])!=VbStructure@nb)
        stop("configList should contain integer vectors with positive numbers of states. All vectors
             should have number of elements equal to the number of variable blocks in variable blocks 
             structure provided!\n")
    }
    
    BIC <- rep(0, length(configList))
    
    pb <- txtProgressBar(min=0, max=length(configList), style=3)
    ipbar <- 1
    
    minBIC = 1e100
    optHMMVB = NULL
    
    for (i in 1:length(configList)){
      setTxtProgressBar(pb, ipbar)
      
      NewVb = vb(VbStructure@nb, VbStructure@dim, VbStructure@bdim,
                 configList[[i]], VbStructure@varorder)
      
      result <- tryCatch({
        
        if (i == length(configList))
          HmmVb <- rcpp_trainHmmVb(t(data), NewVb, vbSearchControl(), trControl, nthread, vb, hmm, mkhmmvb, TRUE)
        else
          HmmVb <- rcpp_trainHmmVb(t(data), NewVb, vbSearchControl(), trControl, nthread, vb, hmm, mkhmmvb, FALSE)
        
        
        if (HmmVb@BIC < minBIC){
          minBIC = HmmVb@BIC
          optHMMVB = HmmVb
        }
        
        HmmVb@BIC
      }, error = function(err) {
        cat(paste("Error for numst = ",numst[i],":",err))
        
        return(NA)
        
      } 
      ) # END tryCatch
      
      BIC[i] = result
      
      ipbar <- ipbar + 1
    }
    
    return(new("HMMVBBIC", BIC=BIC, numst=numeric(0), optHMMVB=optHMMVB))
  }
  
  else{
    if (!is.numeric(numst) || !isTRUE(all.equal(numst, as.integer(numst))) || (any(numst <= 0)))
      stop("numst should be a vector with positive integers!\n")
    
    BIC <- rep(0, length(numst))
    pb <- txtProgressBar(min=0, max=length(numst), style=3)
    ipbar <- 1
    
    minBIC = 1e100
    optHMMVB = NULL
    
    for (i in 1:length(numst)){
      setTxtProgressBar(pb, ipbar)
      newnumst = rep(numst[i], VbStructure@nb)
      
      NewVb = vb(VbStructure@nb, VbStructure@dim, VbStructure@bdim,
                 newnumst, VbStructure@varorder)
      
      result <- tryCatch({
        
        if (i == length(numst))
          HmmVb <- rcpp_trainHmmVb(t(data), NewVb, vbSearchControl(), trControl, nthread, vb, hmm, mkhmmvb, TRUE)
        else
          HmmVb <- rcpp_trainHmmVb(t(data), NewVb, vbSearchControl(), trControl, nthread, vb, hmm, mkhmmvb, FALSE)
        
        if (HmmVb@BIC < minBIC){
          minBIC = HmmVb@BIC
          optHMMVB = HmmVb
        }
        
        HmmVb@BIC
      }, error = function(err) {
        cat(paste("Error for numst = ",numst[i],":",err))
        
        return(NA)
        
      } 
      ) # END tryCatch
      
      BIC[i] = result
      
     
      ipbar <- ipbar + 1
    }
    
    return(new("HMMVBBIC", BIC=BIC, numst=numst, optHMMVB=optHMMVB))
  }
}
