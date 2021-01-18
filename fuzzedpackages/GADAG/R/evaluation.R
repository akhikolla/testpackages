#################################################################
#### evaluating a solution ####
##' @title Evaluate the fitness of a population
##' @description Internal function of the genetic algorithm that evaluates the fitness (penalized log-likelihood) of a set (population) of permutations. It internally computes the best triangular matrix associated to each permutation of the population.
##' @param Pop Population of permutations from [1,p]: matrix with \code{pop.size} rows and p columns, each row corresponding to one permutation of the population.
##' @param X Design matrix, with samples (n) in rows and variables (p) in columns.
##' @param XtX (optional) Cross-product of X; computed if not provided.
##' @param lambda Parameter of penalization (>0).
##' @param grad.control A list containing the parameters for controlling the inner optimization, i.e. the gradient descent
##' \itemize{
##' \item{\code{tol.obj.inner}}{ tolerance (>0),}
##' \item{\code{max.ite.inner}}{ maximum number of iterations (>0).}
##' }
##' @param ncores Number of cores (>1, depending on your computer).
##' @return A list with the following elements:
##' \itemize{
##' \item{Tpop}{ Matrix with pxp columns, each column corresponding to the best triangular matrix (in a vector form) associated to each permutation of the population.}
##' \item{f}{ Fitness of the population.}
##' }
##' @rawNamespace export(evaluation)
##' @seealso \code{\link{GADAG}}, \code{\link{GADAG_Run}}, \code{\link{fitness}}.
##' @return A list with the following elements:
##' \itemize{
##' \item{\code{Tpop}}{ Matrix with p rows and pop.size columns, each column corresponding to the best triangular matrix (in a vector form) associated to each permutation of the population.}
##' \item{\code{f}}{ Fitness of the population.}
##' }
##' @author \packageAuthor{GADAG}
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
##'  ########################################################
##'  # Creating a population of permutations
##'  ########################################################
##'  # first, define the parameters
##'  p <- ncol(toy_data$G) # number of nodes
##'  pop.size <- 10 # population size
##'
##'  # then create your population of permutations
##'  Pop <- matrix(data = 0, nrow = pop.size, ncol = p)
##'  for(i in 1:pop.size){
##'      Pop[i,] = sample(p)
##'  }
##'
##'  ########################################################
##'  # Evaluating the fitness of the population
##'  ########################################################
##'  # evaluation of the whole population
##'  Evaluation <- evaluation(Pop=Pop,X=toy_data$X,lambda=0.1)
##'  print(Evaluation$f) # here is the fitness of the whole population
##'
##'  # evaluation of one of the permutation
##'  Perm <- Pop[1,]
##'  Evaluation <- evaluation(Pop=Perm,toy_data$X,lambda=0.1)
##'
##'  # optimal matrix T associated to Perm:
##'  T <- matrix(Evaluation$Tpop,p,p)

evaluation <- function(Pop, X, XtX=NULL, lambda, grad.control = list(tol.obj=1e-6, max.ite=50), ncores=1){
  # Pop: population of permutations (pop.size*p)
  # X: observation matrix (n*p)
  # XtX: t(X)%*%X matrix, precomputed for speed
  # lambda: penalization term
  # tol.obj: tolerance for the gradient descent (on the norm of T)
  # max.ite: maximum number of iterations for the gradient descent
  #
  # OUTPUTS
  # List of two with
  # Tpop: optimal T values for the population (pop.size*p^2 matrix, one row for each individual)
  # f: fitness of the population
  if (is.null(grad.control$tol.obj)){
    tol.obj <- 1e-6
  } else {
    tol.obj <- grad.control$tol.obj
  }
  if (is.null(grad.control$max.ite)){
    max.ite <- 50
  } else {
    max.ite <- grad.control$max.ite
  }

  if (is.null(XtX)) {
    XtX <- crossprod(X)
  }

  n <- dim(X)[1]
  p <- dim(X)[2]
  L <- (2/n) * norm(XtX,'f') # Lispchitz constant
  if (max.ite < 0){
    stop("max.ite should be non-negative.")
  }
  if (is.vector(Pop)==TRUE){
    Pop <- t(as.matrix(Pop,ncol=length(Pop),nrow=1))
  }
  pop.size <- dim(Pop)[1]
  if (ncol(Pop)!=ncol(X)){
    stop("The number of variables of Pop does not correspond to the number of variables of X.")
  }
  my.gradientdescent <- function(i){
    gradientdescent(P=chrom(Pop[i,]), n=n, XtX=XtX, L=L, lambda=lambda, maxite=max.ite, tolobj=tol.obj)
  }

  Tpop <- matrix(unlist(mclapply(X=1:pop.size, FUN=my.gradientdescent, mc.cores = ncores, mc.preschedule=FALSE)), nrow=pop.size, byrow=TRUE)

  my.fitness <- function(i){
    fitness(P=chrom(Pop[i,]), X, matrix(Tpop[i,], p, p), lambda=lambda)
  }

  f <- unlist(mclapply(X=1:pop.size, FUN=my.fitness))

  return(list(Tpop = Tpop, f = f))
}
