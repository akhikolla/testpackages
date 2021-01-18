\name{GADAG-package}
\alias{GADAG}
\docType{package}
\title{
  \packageTitle{GADAG}
}
\seealso{\code{\link{GADAG_Run}}, \code{\link{GADAG_Analyze}}}

\description{
  \packageDescription{GADAG}
}
\details{

  GADAG aims at recovering the structure of an unknow DAG G, whose edges represent the interactions that exist between p nodes, using n noisy observations
  of these nodes (design matrix X).
  GADAG is more precisely based on a l1-penalized (to make the estimated graph sparse enough) maximum log-likelihood estimation procedure, with the constraint that the estimated graph is a DAG.
  This DAG learning problem is particularly critical in the high-dimensional setting, the exploration of
  the whole of set of DAGs being a NP-hard problem.
  GADAG proposes an original formulation for the estimated DAG, splitting the
  initial problem into two sub-problems: node ordering and graph topology search.
  The node ordering, modelled as a permutation of [1,p] or the associated pxp matrix P, represents the importance of the p nodes of the graph,
  from the node with the smallest number of children to the node with the largest number of children.
  The topological structure of the graph, which is given as a lower triangular matrix T,
  then sets the graph edges weights (including 0, equivalent to no edges).
  GADAG works as follows:  it efficiently looks for the best permution in an outer loop with a genetic algorithm,
  while a nested loop is used to find the optimal T associated to each given P. The latter internal optimization
  problem is solved by a steepest gradient descent approach.

  The DESCRIPTION file:
    \packageDESCRIPTION{GADAG}
}
\author{
  \packageAuthor{GADAG}

  Maintainer: \packageMaintainer{GADAG}
}
\keyword{Directed Acyclic Graphs - Convex program - Optimization}
\examples{
  #############################################################
  # Loading toy data
  #############################################################
  data(toy_data)
  # toy_data is a list of two matrices corresponding to a "star"
  # DAG (node 1 activates all other nodes):
  # - toy_data$X is a 100x10 design matrix
  # - toy_data$G is the 10x10 adjacency matrix (ground trough)

  #############################################################
  # Running GADAG
  #############################################################
  # Simple run, with only the penalty term specified
  GADAG_results <- GADAG_Run(X=toy_data$X, lambda=0.1)
  print(GADAG_results$G.best) # optimal adjacency matrix graph

  # Expensive run with many evaluations if we refine the
  # termination conditions
  \dontrun{
  n.gen <- 1e10 # we allow a very large number of iterations
  tol.Shannon <- 1e-10 # the entropy of Shannon of the population
                       # has to be very small
  pop.size <- 5*ncol(toy_data$G) # this is usually a good
                                 # population size
  max.eval <- n.gen * pop.size # maximal number of nested
                               # evaluation
  GADAG_results <- GADAG_Run(X=toy_data$X, lambda=0.1,
       GADAG.control=list(n.gen=n.gen, tol.Shannon=tol.Shannon,
                          pop.size = pop.size, max.eval=max.eval))
  print(GADAG_results$G.best) # optimal adjacency matrix graph         }

  # Expensive run if we also increase the population size
  \dontrun{
  pop.size <- 10*ncol(toy_data$G)
  GADAG_results <- GADAG_Run(X=toy_data$X, lambda=0.1,
      GADAG.control=list(pop.size=pop.size))
  }

  # You can have more information about the evolution of the
  # algorithm by turning return.level on
  \dontrun{
  return.level <- 1
  GADAG_results <- GADAG_Run(X=toy_data$X, lambda=0.1, return.level = return.level)
  print(GADAG_results$f.best.evol) # this shows the evolution of the fitness
                                   # across the iterations
  }
}
