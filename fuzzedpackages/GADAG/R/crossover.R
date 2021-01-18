###################################################################
###### cross-over operator ######

##' @keywords internal
# @title Crossover permutations
# @description Internal function of the genetic algorithm (crossover operator) to combine elements of a population of permutations.
# @param Pop Population of permutations from [1,p] (output from create.population() function for example).
# @param p.xo Crossover probability.
# @return A list with the following elements:
# \itemize{
# \item{Children}{ Population of permutations created after crossovering permutations of Pop.}
# \item{I.cross}{ Vector of integers from [1,p] corresponding to the elements of Pop that have a non-zero probability for crossovering.}
# }
# @seealso \code{\link{GADAG}}, \code{\link{GADAG_Run}}.
# @author \packageAuthor{GADAG}
# @examples
#  ########################################################
#  # Creating a population of permutations
#  ########################################################
#  Population <- create.population(p=10, pop.size=20)
#
#  ########################################################
#  # Crossovering the population
#  ########################################################
#  Population_Crossover <- crossover(Pop=Population)

crossover <- function(Pop,p.xo=.25){
  # INPUTS
  # Pop: population of permutations (pop.size*p)
  # p.xo: probability of crossover
  #
  # OUTPUTS
  # Children: population of permutations (pop.size.p)
  # I.cross: T/F indicating which permutations have a non-zero probability for crossover

  if (is.vector(Pop)==TRUE){
    stop("Pop should have at least two elements to be run.")
  }
  pop.size <- dim(Pop)[1]
  p        <- dim(Pop)[2]

  # Select subpopulation for cross-over
  do.xo <- runif(pop.size) < p.xo
  I.cross <- c(1:pop.size)[do.xo]
  if (length(I.cross)>1){
    if (length(I.cross)%%2 !=0){
      I.cross <- I.cross[-1]
    }
    Pop <- Pop[I.cross,]
    pop.size.crossover <- length(I.cross)

    # Number of swaps
    n.swaps <- sample.int(p, size=pop.size.crossover/2, replace=TRUE)

    # All the swaps
    mysample.int <- function(j){return(sample.int(p, size=max(n.swaps), replace=FALSE))}
    all.swaps    <- t(apply(matrix(1:(pop.size.crossover/2)), 1, mysample.int))

    if (max(n.swaps)==1) all.swaps <- t(all.swaps)

    Children <- matrix(0,ncol=p,nrow=pop.size.crossover)
    for (i in seq(1,pop.size.crossover-1,2)) {
      swaps.nodes <- all.swaps[(i+1)/2, 1:n.swaps[(i+1)/2]]
      Children[i:(i+1),swaps.nodes] <- Pop[i:(i+1),swaps.nodes]
      Children[i,-swaps.nodes] <- Pop[i+1,which(Pop[i+1,] %in% (1:p)[-Pop[i,swaps.nodes]])]
      Children[i+1,-swaps.nodes] <- Pop[i,which(Pop[i,] %in% (1:p)[-Pop[i+1,swaps.nodes]])]
    }
  } else {
    Children <- c("No children")
  }
  return(list(Children=Children,I.cross=I.cross))
}
