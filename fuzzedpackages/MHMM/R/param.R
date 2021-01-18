########################################################################################################################
## Classe S4 mhmmparam
########################################################################################################################
###################################################################################
##' Constructor of \code{\linkS4class{mhmmparam}} class
##'
##'  
##' \describe{
##'   \item{K}{numeric. Number of classes.}
##'   \item{M}{numeric. Number of latent states (activity levels).}
##'   \item{A}{list. Matrices of the transition probailities for each class.}
##'   \item{delta}{numeric. Proportions of the classes.}
##'   \item{pi}{matrix. Probabilities of the latent states per class (stationary distribution of the Markov chains).}
##'   \item{lambda}{list. Parameters of the emission laws (eps: proportions of the zero-inflated component, a: shapes of the gamma distributions, b: rates of the gamma distributions)}
##' }
##'
#' @examples
#'   getSlots("mhmmparam")
#'
#' @name mhmmparam-class
#' @rdname mhmmparam-class
#' @exportClass mhmmparam
##' @export
setClass(
  Class = "mhmmparam",
  representation = representation(
    K="numeric",
    M="numeric",
    A="list",
    delta="numeric",
    pi="matrix",
    lambda="list"
  ),
  prototype = prototype(
    K=numeric(),
    M=numeric(),
    A=list(),
    delta=numeric(),
    pi=matrix(),
    lambda=list()
  )
)


init.mhmm <- function(y, K, M, lambda=NULL){
  delta <- runif(K)
  A <- replicate(K,  {tmp <- matrix(runif(M * M), M, M); for (l in 1:M) tmp[l, ] <- tmp[l,]/sum(tmp[l,]);tmp;}, simplify = FALSE)
  eps <- runif(M)
  if (!is.null(lambda)){
    lambda$eps <- eps
  }else{
    tmpech <- NULL
    while (length(tmpech)<(2*M)){
      tmpech <- na.omit(y[[sample(1:length(y), 1)]])
      tmpech <- tmpech[which(tmpech!=0)]
    }
    val <- sample(tmpech, M)
    b <- sqrt( val / var(tmpech, na.rm = TRUE) ) + .1
    lambda <- list(eps = eps, a = val * b + 0.1, b = b)
  }
  delta <- runif(K)
  new("mhmmparam", K = K, M = M, A = A, delta = delta/sum(delta), pi =  matrix(1/M, K, M),lambda = lambda)
}