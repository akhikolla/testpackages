##' @title Generate toy data for running GADAG
##' @description This function generates toy data that can be used to run GADAG: the adjacency matrix of a DAG with p nodes and the design matrix with n observations of the distribution of the p nodes.
##' @param n Number of samples in the design matrix.
##' @param p Number of nodes of the DAG.
##' @param edgemin Minimal value for the non-null edges of the DAG (between 0 and 1).
##' @param type Form of the DAG. It can be chosen between 7 alternatives: \code{"star"}, \code{"bistar"}, \code{"full"}, \code{"path"}, \code{"quadristar"}, \code{"sixstar"} (see details below).
##' @param seed Fix the seed.
##' @return A list containing the design nxp matrix X (with samples in rows and variables in columns) and the adjacency matrix G associated to the DAG with p nodes.
##' @author \packageAuthor{GADAG}.
##' @details One of the following seven alternatives can be chosen for the DAG form:
##' \itemize{
##' \item{\code{"star"}}{ star-shaped DAG (all active edges start from node 1),}
##' \item{\code{"bistar"}}{ half of the edges start from node 1 and the other half from node 2,}
##' \item{\code{"full"}}{ full DAG (all the edges are active),}
##' \item{\code{"path"}}{ path-shaped DAG (all the nodes are connected by a single path),}
##' \item{\code{"quadristar"}}{ node 1 is connected to nodes 2 to 4, each being connected to 1/3 of the rest of the nodes,}
##' \item{\code{"sixstar"}}{ same as \code{"quadristar"}, with 6 nodes.}
##' }
##' @seealso \code{\link{GADAG}}, \code{\link{GADAG_Run}}.
##' @examples
##'  #############################################################
##'  # Generating toy data
##'  #############################################################
##'  toy_data <- generateToyData(n=100, p=10)
##'
##'  # toy_data is a list of two matrices corresponding to a "star"
##'  # DAG (node 1 activates all other nodes):
##'  # - toy_data$X is a 100x10 design matrix
##'  # - toy_data$G is the 10x10 adjacency matrix (ground trough)
##'
##'  \dontrun{
##'  # generate another type of data: a DAG with 100 nodes in a path form
##'  toy_data <- generateToyData(n=100, p=100,type="path")
##'  }
##'
##'  \dontrun{
##'  # set the minimal edge value to 1
##'  toy_data <- generateToyData(n=100, p=10, edgemin=1) # all edges are set to 1
##'  }

generateToyData <- function(n, p, edgemin=0, type="star", seed=42) {
  #-----------------------------------------------------
  # Generates DAGs and associated data
  # The form of the DAG can be chosen between 7 alternatives (see below)
  # The intensity of the active edges is random (uniformly distributed)
  #
  # INPUTS
  # n: number of observations
  # p: number of nodes
  # edgemin: minimal value for the non-null edges. Scalar in [0,1]
  # type: "star": star-shaped DAG (all active edges start from node 1)
  #       "bistar": half of the edges start from node 1 and the other half from node 2
  #       "full": full DAG (all the edges are active)
  #       "path": path-shaped DAG (all the nodes are connected by a single path)
  #       "quadristar": node 1 is connected to nodes 2 to 4, each being connected to 1/3 of the rest of the nodes
  #       "sixstar": same as "quadristar", with 6 nodes
  # seed: to fix the seed

  # OUTPUTS
  # X: matrix of observations (n*p)
  # G: matrix of the DAG (p*p)
  #-----------------------------------------------------
  if (n < 0 || p < 0){
    stop("Please choose p and n as non-negative numbers.")
  }

  if (n==1){
    stop("You need more than one sample to run GADAG.")
  }

  if (edgemin >1 || edgemin < 0){
    stop("edgemin has to taken between O and 1.")
  }

  set.seed(seed)

  G <- matrix(0,p,p)

  if (type=="star") {
    G[1,] <- c(0, edgemin+runif(p-1)*(1-edgemin))
  } else if (type=="bistar") {
    p1 <- round((p-1)/2)
    p2 <- p - p1 - 1
    G[1,] <- c(0, edgemin+runif(p1)*(1-edgemin), rep(0, p2))
    G[2,] <- c(rep(0, p1+1), runif(p2))
  } else if (type=="full") {
    G <- matrix(edgemin+runif(p*p)*(1-edgemin), p, p)
    G[lower.tri(G, diag=TRUE)] <- 0
  } else if (type=="path") {
    for (i in 2:p) {
      G[i-1,i] <- edgemin+runif(1)*(1-edgemin)
    }
  } else if (type=="quadristar") {
    p1 <- round((p-4)/3)
    G[1,2:4] <- edgemin+runif(3)*(1-edgemin)
    G[2,] <- edgemin+runif(p)*(1-edgemin)
    G[3,] <- edgemin+runif(p)*(1-edgemin)
    G[4,] <- edgemin+runif(p)*(1-edgemin)
    G[2,c(1:4, (p1+5):p)] <- 0
    G[3,c(1:(p1+4), (2*p1+5):p)] <- 0
    G[4, 1:(2*p1+4)] <- 0
  } else if (type=="sixstar") {
    p1 <- round((p-6)/5)
    G[1,2:6] <- edgemin+runif(5)*(1-edgemin)
    G[2,] <- edgemin+runif(p)*(1-edgemin)
    G[3,] <- edgemin+runif(p)*(1-edgemin)
    G[4,] <- edgemin+runif(p)*(1-edgemin)
    G[5,] <- edgemin+runif(p)*(1-edgemin)
    G[6,] <- edgemin+runif(p)*(1-edgemin)
    G[2,c(1:6, (p1+7):p)] <- 0
    G[3,c(1:(p1+6), (2*p1+7):p)] <- 0
    G[4,c(1:(2*p1+6), (3*p1+7):p)] <- 0
    G[5,c(1:(3*p1+6), (4*p1+7):p)] <- 0
    G[6, 1:(4*p1+6)] <- 0
  }

  eps <- matrix(rnorm(n*p,mean=0,sd=1),n,p)
  X <- eps %*% ginv(diag(1,p,p)-G)
  X = X - matrix(1,n,1) %*% colMeans(X)
  X = X / sqrt( (1/n) * matrix(1,n,1) %*% colSums(X^2))
  return(list(X=X,G=G))
}
