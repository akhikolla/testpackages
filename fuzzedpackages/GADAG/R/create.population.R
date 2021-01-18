###################################################################
#### creating a population of permutation vectors ####

##' @keywords internal
# @title Create a population of permutations
# @description Internal function of the genetic algorithm to create a population of permutations.
# @usage create.population(p, pop.size)
# @param p Length of each permutation (>0), corresponding to the number of nodes of the DAG to recover for the DAG learning problem.
# @param pop.size Length (number of permutations) of the population (>0).
# @return A pop.size x p matrix corresponding to the population of permutations.
# @author \packageAuthor{GADAG}
# @seealso \code{\link{GADAG}}, \code{\link{GADAG_Run}}.
# @examples
#  ########################################################
#  # Creating a population of permutations
#  ########################################################
#  Population <- create.population(p=10, pop.size=20)
#

create.population = function(p, pop.size){
  # INPUTS
  # p: number of variables
  # pop.size: initial population size
  #
  # OUTPUTS
  # Pop: population of permutations (pop.size*p)
  if (p < 0 || pop.size < 0){
    stop('p and pop.size should be non-negative.')
  }
  Pop <- matrix(data = 0, nrow = pop.size, ncol = p)
  for(i in 1:pop.size){
    Pop[i,] = sample(p)
  }
  return(Pop)
}
