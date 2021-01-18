###################################################################################
##' Mixture model of Hidden Markov Models.
##'
##' @description
##' This function performs maximum likelihood inference of mixture of hidden Markov models
##'
##' @param y list. Observations (one element of the list per subject, missing values should be encoded by NA).
##' @param K numeric. Number of mixture components.
##' @param M numeric. Number of states.
##' @param smartinit boolean. If TRUE, MLE of the mixture neglecting the time dependency is used for the initialization of the algorithm.
##' @param nbinit numeric.Number of initializations.
##' @param tol numeric. It indicates the maximal gap between two successive iterations of EM algorithm which stops the algorithm
##' @param nbKeep numeric. It indicates the number of chains used for the final EM algorithm
##' @param iterSmall numeric. It indicates  the number of iterations for each SmallEM algorithm
##' @param nbcores numeric.  It defines the numerber of cores used by the algorithme
##'
##' @return Returns an instance of \linkS4class{mhmmresults}.
##' @references  Du Roy de Chaumaray, M. and Marbac, M. and Navarro, F. (2019). Mixture of hidden Markov models for accelerometer data. arXiv preprint arXiv:1906.01547
##' @export mhmm
##'
##' @examples
##' data(accelero)
##' # To make the estimation <5
##' res <- mhmm(accelero, K = 2, M = 4, nbcores = 1, nbinit = 5, iterSmall = 2)
##' plot(res, 1)
##' 
##'  \donttest{
##' data(accelero)
##' # It is better to increase the number of random initializations
##' res <- mhmm(accelero, K = 2, M = 4, nbcores = 1)
##' plot(res, 1)
##' }
##'
mhmm <- function(y,
                 K,
                 M,
                 smartinit=TRUE,
                 nbinit = 100,
                 tol = 10**(-4),
                 nbKeep = min(20, nbinit),
                 iterSmall = 10,
                 nbcores = 1) {
  checkinputs(y, K, M, smartinit, nbinit, tol, nbKeep, iterSmall, nbcores)
  dataobs <- mhmmData(y)
  lambda <- NULL
  if (smartinit) lambda <- smartinitparam(y, M, nbcores)
  all.starts <- replicate(nbcores, replicate(ceiling(nbinit/nbcores), init.mhmm(y, K, M, lambda), simplify = FALSE), simplify = FALSE)
  out <- mclapply(all.starts,
                  EMmhmmCPP,
                  yiR = dataobs,
                  tolR = tol,
                  nbKeepR = ceiling(nbKeep/nbcores),
                  itersmallR = iterSmall,
                  mc.cores = nbcores)
  all.loglike <- sapply(out, function(u) u$loglike)
  if (any(is.nan(all.loglike))){
    all.loglike[which(is.nan(all.loglike))] <- -Inf
  }
  out <- out[[which.max(all.loglike)]]
  tuneOutput(out, dataobs, y)
}
 
probamhmm <- function(y, param){
  dataobs <- mhmmData(y)
  tuneOutput(EMmhmmCPP(yiR = dataobs, paramR = list(param), tolR = Inf, nbKeepR = 1, itersmallR = 0), dataobs, y)
}

tuneOutput <- function(out, dataobs, y){
  out$param@A <- lapply(out$param@A, function(u) matrix(u, out$param@M, out$param@M))
  out$param@lambda$eps <- as.numeric(out$param@lambda$eps)
  out$param@lambda$a <- as.numeric(out$param@lambda$a)
  out$param@lambda$b <- as.numeric(out$param@lambda$b)
  out$param@delta <- as.numeric(out$param@delta)
  K <- length(out$param@delta)
  nbobs <- length(na.omit(unlist(y)))
  out <- new("mhmmresults",
             param = out$param,
             data = dataobs,
             partitions = list(clusters = apply(out$probaClusters, 1, which.max),
                               states = lapply(out$probaStates, function(u) lapply(u, function(v) apply(matrix(v, ncol = out$param@M), 1, which.max)))),
             probabilities = list(clusters = out$probaClusters,
                                  states = lapply(out$probaStates, function(u) lapply(u, function(v) matrix(v, ncol = out$param@M)))),
             loglike = out$loglike,
             bic = out$loglike - 0.5 * (3 * out$param@M + K * K + (K-1) + K) * log(nbobs),
             meanvalueperstates = (1-out$param@lambda$eps) * out$param@lambda$a / out$param@lambda$b
  )
  for (i in 1:length(y)){
    tmp <- y[[i]]
    tmp[which(!is.na(tmp))] <- unlist(out@partitions$states[[i]])
    out@partitions$states[[i]] <- tmp
    tmp <- out@probabilities$states[[i]]
    
    out@probabilities$states[[i]] <- matrix(NA, length(y[[i]]), length(out@param@lambda$eps))
    if (length(out@param@lambda$eps)>1){
      for (h in 1:length(out@param@lambda$eps))
        out@probabilities$states[[i]][which(!is.na(y[[i]])),h] <- unlist(sapply(tmp, function(i) i[,h]))
    }else{
      out@probabilities$states[[i]][which(!is.na(y[[i]])),] <- 1
    }
    
  }
  out@meantimesperstates <-   t(sapply(out@probabilities$states, function(u) colMeans(u, na.rm = T)))
  out@partitions$MAPstates <- out@partitions$states
  for (i in 1:length(out@data@yi)) out@partitions$MAPstates[[i]][which(!is.na(out@partitions$MAPstates[[i]]))] <- unlist(lapply(out@data@yi[[i]], Viterbi, z=out@partitions$clusters[i], p=out@param))
  out@data@yi <- y
  or <- order(out@meanvalueperstates)
  out@meanvalueperstates <- out@meanvalueperstates[or]
  out@meantimesperstates <- out@meantimesperstates[,or]
  out@param@lambda$eps <- out@param@lambda$eps[or]
  out@param@lambda$a <- out@param@lambda$a[or]
  out@param@lambda$b <- out@param@lambda$b[or]
  for (k in 1:length(out@param@A)){
    out@param@A[[k]] <- out@param@A[[k]][,or]
    out@param@A[[k]] <- out@param@A[[k]][or,]    
    ref <- out@param@A[[k]]
    for (it in 1:200) ref <- ref %*% out@param@A[[k]]
    out@param@pi[k,] <- ref[1,] 
  }

  out@probabilities$states <- lapply(out@probabilities$states, function(u) u[,or])  
  tmp <- cbind(or, 1:length(or))
  tmp2 <- tmp[order(tmp[,1]),]
  out@partitions$states <- lapply(out@partitions$states, function(u) tmp2[u,2])
  out@partitions$MAPstates <- lapply(out@partitions$MAPstates, function(u) tmp2[u,2])
  out
}