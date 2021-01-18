#################################################################
##' @title Compute the fitness of a potential solution
##' @description Internal function of the genetic algorithm that evaluates the fitness (penalized log-likelihood) of a potential solution, given as a pair of a permutation (P) and a triangular matrix (T).
##' @usage fitness(P,X,T,lambda)
##' @param P A permutation from [1,p] in a matrix form.
##' @param X Design matrix, with samples (n) in rows and variables (p) in columns.
##' @param T A pxp lower-triangular matrix.
##' @param lambda Parameter of penalization (>0).
##' @return A numeric value corresponding to the fitness of the potential solution.
##' @author \packageAuthor{GADAG}
##' @seealso \code{\link{GADAG}}, \code{\link{GADAG_Run}}, \code{\link{evaluation}}.
##' @rawNamespace export(fitness)
##' @examples
##'  #############################################################
##'  # Loading toy data
##'  #############################################################
##'  data(toy_data)
##'  # toy_data is a list of two matrices corresponding to a "star"
##'  # DAG (node 1 activates all other nodes):
##'  # - toy_data$X is a 100x10 design matrix
##'  # - toy_data$G is the 10x10 adjacency matrix (ground trough)
##'
##'  ############################################################
##'  # Creating a candidate solution
##'  ############################################################
##'  # define parameters
##'  p <- ncol(toy_data$G)
##'
##'  # permutation matrix
##'  Perm <- sample(p) # permutation in a vector form
##'  P <- matrix(0,p,p)
##'  P[p*0:(p-1) + Perm] <- 1 # Perm is tranformed into a matrix form
##'
##'  # lower-triangular matrix
##'  T <- matrix(rnorm(p),p,p)
##'  T[upper.tri(T,diag=TRUE)] <- 0
##'
##'  ########################################################
##'  # Computing the fitness of the potential solution
##'  ########################################################
##'  Fitness <- fitness(P=P, X=toy_data$X, T=T, lambda=0.1)
##'  print(Fitness) # here is the fitness of the candidate solution (P,T)

fitness <- function(P,X,T,lambda){
  # INPUTS
  # P: permutation matrix (p*p)
  # X: observation matrix (n*p)
  # T: triangular matrix (p*p)
  # lambda: penalization term (scalar)

  # OUTPUTS
  # f: fitness value

  n <- dim(X)[1]
  f <- (1/n)*norm(X - X %*% P %*% T %*% t(P),'f')^2 + lambda*(sum(abs(T)))
  return(f)
}
