#'BoltzMM: A package for probability computation, data generation, and model estimation of fully-visible Boltzmann machines.
#'
#'The BoltzMM package allows for computation of probability mass functions of fully-visible Boltzmann machines via \code{pfvbm} and \code{allpfvbm}.
#'Random data can be generated using \code{rfvbm}. Maximum pseudolikelihood estimation of parameters via the MM algorithm can be conducted using \code{fitfvbm}.
#'Computation of partial derivatives and Hessians can be performed via \code{fvbmpartiald} and \code{fvbmHessian}.
#'Covariance estimation and normal standard errors can be computed using \code{fvbmcov} and \code{fvbmstderr}.
#'
#'@author Andrew T. Jones and Hien D. Nguyen
#'@references H.D. Nguyen and I.A. Wood (2016), Asymptotic normality of the maximum pseudolikelihood estimator for fully-visible Boltzmann machines, IEEE Transactions on Neural Networks and Learning Systems, vol. 27, pp. 897-902.
#'
#' H.D. Nguyen and I.A. Wood (2016), A block successive lower-bound maximization algorithm for the maximum pseudolikelihood estimation of fully visible Boltzmann machines, Neural Computation, vol 28, pp. 485-492.
#'
#'@docType package
#'@name BoltzMM
NULL

#'@importFrom stats pnorm
NULL

#'Standard errors for the parameter elements of a fitted fully-visible Boltzmann machine.
#'@description Computes the normal approximation standard errors from the sandwich estimator of the covariance matrix for a maximum pseudolikelihood estimated fully-visible Boltzmann machine.
#'@param data An N by n matrix, where each of the N rows contains a length n string of spin variables  (i.e. each element is -1 or 1).
#'@param covarmat A covariance matrix generated from \code{fvbmcov}.
#'@return A list containing 2 objects: a vector containing the standard errors corresponding to the bias parameters \code{bvec_se}, and a matrix containing the standard errors corresponding to the interaction parameters \code{Mmat_se}.
#'@references H.D. Nguyen and I.A. Wood (2016), Asymptotic normality of the maximum pseudolikelihood estimator for fully-visible Boltzmann machines, IEEE Transactions on Neural Networks and Learning Systems, vol. 27, pp. 897-902.
#'@author Andrew T. Jones and Hien D. Nguyen
#'@examples
#'# Generate num=1000 random strings of n=3 binary spin variables under bvec and Mmat.
#'num <- 1000
#'bvec <- c(0,0.5,0.25)
#'Mmat <- matrix(0.1,3,3) - diag(0.1,3,3)
#'data <- rfvbm(num,bvec,Mmat)
#'# Fit a fully visible Boltzmann machine to data, starting from parameters bvec and Mmat.
#'model <- fitfvbm(data,bvec,Mmat)
#'# Compute the sandwich covariance matrix using the data and the model.
#'covarmat <- fvbmcov(data,model,fvbmHess)
#'# Compute the standard errors of the parameter elements according to a normal approximation.
#'fvbmstderr(data,covarmat)
#'@export
fvbmstderr <- function(data,covarmat) {
  N <- dim(data)[1]
  D <- dim(data)[2]

  stderr <- sqrt(diag(covarmat))/sqrt(N)

  bvec <- stderr[c(1:D)]
  Mmat <- matrix(0,D,D)
  Mmat[lower.tri(Mmat)] <- stderr[-c(1:D)]
  Mmat <- Mmat + t(Mmat)
  return(list(bvec_se = bvec, Mmat_se = Mmat))
}


#'Hessian of the log-pseudolikelihood function for a fitted fully-visible Boltzmann machine.
#'@description Computes the Hessian with respect to all unique parameter elements of the bias vector and interaction matrix of a fully-visible Boltzmann machine, for some random length n string of spin variables (i.e. each element is -1 or 1) and some fitted parameter values.
#'@param data An N by n matrix, where each of the N rows contains a length n string of spin variables  (i.e. each element is -1 or 1).
#'@param model List generated from \code{fitfvbm}.
#'@return The n+choose(n,2) by n+choose(n,2) Hessian matrix, summed over the N rows of \code{data} and evaluated at the fitted parameter values provided in \code{model}. Each row (column) is a unique element of the bias vector and interaction matrix. The rows are arranged in lexicographical order with the bias elements first, followed by the interaction elements. For example, if n=3, the order would be bias[1], bias[2] bias[3], interaction[1,2], interaction[1,3], and interaction[2,3].
#'@references H.D. Nguyen and I.A. Wood (2016), Asymptotic normality of the maximum pseudolikelihood estimator for fully-visible Boltzmann machines, IEEE Transactions on Neural Networks and Learning Systems, vol. 27, pp. 897-902.
#'@author Andrew T. Jones and Hien D. Nguyen
#'@examples # Generate num=1000 random strings of n=3 binary spin variables under bvec and Mmat.
#'num <- 1000
#'bvec <- c(0,0.5,0.25)
#'Mmat <- matrix(0.1,3,3) - diag(0.1,3,3)
#'data <- rfvbm(num,bvec,Mmat)
#'# Fit a fully visible Boltzmann machine to data, starting from parameters bvec and Mmat.
#'model <- fitfvbm(data,bvec,Mmat)
#'# Compute the Hessian matrix summed over all num rows of data.
#'fvbmHess(data,model)
#'@export
fvbmHess <- function(data, model) {
  bvec <- model[[2]]
  Mmat <- model[[3]]

  N <- dim(data)[1]
  D <- length(bvec)

  HessComps <- list()
  HessComps[[1]] <- matrix(0,D+1,D+1)
  for (jj in 1:D) {
    HessComps[[jj]] <- matrix(0,D+1,D+1)
    for (ii in 1:N) {
      x_bar <- as.matrix(c(1,data[ii,]),D+1,1)
      HessComps[[jj]] <- HessComps[[jj]] - x_bar%*%t(x_bar)/
        cosh(sum(Mmat[jj,]*data[ii,])+bvec[jj])^2
    }
  }

  Index <- matrix(0,D,D)
  Index[lower.tri(Index)] <- 1:(D*(D-1)/2)
  Index <- Index + t(Index)

  BigHess <- matrix(0,D+D*(D-1)/2,D+D*(D-1)/2)
  for (jj in 1:D) {
    WHICH <- which(Index[lower.tri(Index)]%in%Index[jj,])
    #Index[lower.tri(Index)] is 1:(D*(D-1)/2)
    #Index[jj,] is row j of Index


    NonZero <- HessComps[[jj]][-c(jj+1),]
    NonZero <- NonZero[,-c(jj+1)]
    BigHess[c(jj,D+WHICH),c(jj,D+WHICH)] <- BigHess[c(jj,D+WHICH),c(jj,D+WHICH)] +
      NonZero
  }

  return(BigHess)
}


#'Hypothesis testing for a fully-visible Boltzmann machine.
#'@description Tests the hypothesis that the true bias and interaction parameter values are those in \code{nullmodel}, given \code{data} and \code{model}.
#'@param data An N by n matrix, where each of the N rows contains a length n string of spin variables  (i.e. each element is -1 or 1).
#'@param model List generated from \code{fitfvbm}.
#'@param nullmodel A list containing two elements: a vector of length n \code{bvec}, and an n by n matrix \code{Mmat}. A list generated by \code{fitfvbm} is also sufficient.
#'@return A list containing 4 objects: a vector containing the z-scores corresponding to the bias parameters \code{bvec_z},a vector containing the p-values corresponding to the bias parameters \code{bvec_p},a matrix containing the z-scores corresponding to the interaction parameters \code{Mmat_z},  and a matrix containing the standard errors corresponding to the interaction parameters \code{Mmat_p}.
#'@references H.D. Nguyen and I.A. Wood (2016), Asymptotic normality of the maximum pseudolikelihood estimator for fully-visible Boltzmann machines, IEEE Transactions on Neural Networks and Learning Systems, vol. 27, pp. 897-902.
#'@author Andrew T. Jones and Hien D. Nguyen
#'@examples
#'# Generate num=1000 random strings of n=3 binary spin variables under bvec and Mmat.
#'num <- 1000; bvec <- c(0,0.5,0.25); Mmat <- matrix(0.1,3,3) - diag(0.1,3,3);
#'data <- rfvbm(num,bvec,Mmat)
#'# Fit a fully visible Boltzmann machine to data, starting from parameters bvec and Mmat.
#'model <- fitfvbm(data,bvec,Mmat)
#'
#'#Propose a null hypothesis model
#'nullmodel <- list(bvec = c(0,0,0), Mmat = matrix(0,3,3))
#'
#'# Compute z-scores
#'fvbmtests(data,model,nullmodel)
#'@export
fvbmtests <- function(data,model,nullmodel) {
  #Compute the z-scores
  #get Hessian
  Hess <- fvbmHess(data,model)
  #get Covaraince matrix
  covmat <- fvbmcov(data,model,fvbmHess)
  #get standard errors
  stderr <- fvbmstderr(data,covmat)
  #z-scores for bias parameters
  zb <- (model$bvec-nullmodel$bvec)/stderr$bvec_se
  #z-scores for interaction parameters
  zM <- (model$Mmat-nullmodel$Mmat)/stderr$Mmat_se
  diag(zM)<-NA
  #p-values from z-scores
  pvalb <- 2*pnorm(-abs(zb))
  pvalM <- 2*pnorm(-abs(zM))
  diag(pvalM)<-NA
  #return list
  return(list(bvec_z = zb, bvec_p = pvalb, Mmat_z = zM, Mmat_p = pvalM))
}

#'Marginal probability function for a fully-visible Boltzmann machine.
#'@description Computes the marginal probabilities (for values = +1 in each coordinate) under under some specified bias vector and interaction matrix, specified by \code{bvec} and \code{Mmat}, respectively.
#'@param bvec Vector of length n containing real valued bias parameters.
#'@param Mmat Symmetric n by n matrix, with zeros along the diagonal, containing the interaction parameters.
#'@return Vector of length n containing the marginal probabilities of +1 in each coordinate.
#'@references H.D. Nguyen and I.A. Wood (2016), Asymptotic normality of the maximum pseudolikelihood estimator for fully-visible Boltzmann machines, IEEE Transactions on Neural Networks and Learning Systems, vol. 27, pp. 897-902.
#'@author Andrew T. Jones and Hien D. Nguyen
#'@examples
#'#Compute the marginal probabilities under bvec and Mmat.
#'# Set the parameter values
#'bvec <- c(0,0.5,0.25)
#'Mmat <- matrix(0.1,3,3) - diag(0.1,3,3)
# Compute the marginal probabilities
#'marginpfvbm(bvec,Mmat)
#'@export
marginpfvbm <- function(bvec, Mmat) {
  # Get dimension of vector
  n <- length(bvec)

  ## Get all strings of length n and probabilities
  strings <- expand.grid(rep(list(0:1),n))
  allprob <-allpfvbm(bvec,Mmat)

  #sum over allprob
  margins <- array(NA,n)
  for (i in seq_len(n)){
    margins[i] <- sum(allprob[which(strings[,i]==1)])
  }
  #return marginal probability vector
  return(margins)
}

#'@title Senate voting data from the 45th Australian Parliament.
#'
#'@description A dataset he data from the first sitting of the Senate of the 45th
#' Australian Parliament, until the final sitting of the year 2016. The first division during
#'this period was conducted on the 31st of August 2016, and the last division was performed
#'on the 1st of December 2016. In total, 147 divisions were performed during this period.
#'
#'Each row represents a division(vote), each column is a party or independent.
#'Data is either "Yes" or "No" depending on the vote. Absences and abstentions are left as NA.
#'See \url{https://hal.archives-ouvertes.fr/hal-01927188v1} for details of data preparation.
#'
#'@source \url{www.aph.gov.au/Parliamentary_Business/Statistics/Senate_StatsNet/General/divisions}
#'@docType data
#'@keywords datasets
#'@name senate
#'@usage data(senate)
#'@format A data frame with 147 rows (votes) and 9 variables (parties).
#'@author Jessica J. Bagnall
#'@examples
#'dim(senate)
NULL


