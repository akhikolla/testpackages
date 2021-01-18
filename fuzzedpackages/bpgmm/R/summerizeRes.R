#' summerizePgmmRJMCMC
#'
#' @param pgmmResList result list from pgmmRJMCMC
#' @param trueCluster true cluster allocation
#' @importFrom mclust adjustedRandIndex
#' @export
#'
#' @examples
#' library("fabMix")
#' library("mclust")
#' library("pgmm")
#' library("mvtnorm")
#' library("mcmcse")
#' library("MASS")
#' library("gtools")
#' n <- 50
#' p <- 10
#' q <- 4
#' K <- 10
#' syntheticDataset <- simData(
#'   sameLambda = TRUE, sameSigma = TRUE, K.true = K, n = n, q = q, p = p,
#'   sINV_values = 1 / ((1:p))
#' )
#' nsim <- 5
#' burn <- 0
#' X <- t(syntheticDataset$data)
#' qnew <- 4
#' Mstep <- 1
#' Vstep <- 1
#' constraint <- c(0, 0, 0)
#' mInit <- 20
#' mVec <- c(1, 20)
#' \donttest{
#' res <- pgmmRJMCMC(X, mInit, mVec, qnew,
#'   niter = nsim, burn = burn, constraint = constraint,
#'   Mstep = Mstep, Vstep = Vstep
#' )
#' }
#'
#' \donttest{
#' summerizePgmmRJMCMC(res, syntheticDataset$class)
#' }
#' @export
summerizePgmmRJMCMC <- function(pgmmResList, trueCluster = NULL) {
  Zalloc <- sumerizeZ(res$ZmatList)

  nCluster <- table(sapply(res$ZmatList, function(x) {
    length(unique(x))
  }))

  nConstraint <- res$constraintList
  nConstraint <- listToStrVec(nConstraint)
  nConstraint <- table(nConstraint, dnn = "")



  sumRes <- list(Zalloc = Zalloc, nCluster = nCluster, nConstraint = nConstraint)

  if (!is.null(trueCluster)) {
    ari <- adjustedRandIndex(trueCluster, Zalloc)
    sumRes$ari <- ari
  }

  sumRes
}
