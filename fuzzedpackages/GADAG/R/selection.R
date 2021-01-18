#################################################################
#### Selection operator ####

##' @keywords internal
# @title Select the permutations for crossovering
# @description Internal function of the genetic algorithm to select the permutations that will probably crossover within a population of permutations.
# @usage selection(Pop, f.pop)
# @param Pop Population of permutations from [1,p] (output from create.population() function for example).
# @param f.pop Fitness of the population of permutations ($f element from evaluation() list output).
# @return A vector of integers between 1 and p corresponding to the elements of the population that will probably crossover.
# @author \packageAuthor{GADAG}
# @seealso \code{\link{GADAG}}, \code{\link{GADAG_Run}}.
# @examples
#  ########################################################
#  # Loading toy data
#  ########################################################
#  data(toy_data)
#
#  ########################################################
#  # Creating a population of permutations
#  ########################################################
#  Population <- create.population(p=ncol(toy_data$X), pop.size=20)
#
#  ########################################################
#  # Evaluating the fitness of the population
#  ########################################################
#  Evaluation <- evaluation(Pop=Population, X=toy_data$X,
#                           XtX=crossprod(toy_data$X),lambda=0.1)
#
#  ########################################################
#  # Selecting the elements for crossovering
#  ########################################################
#  Selection <- selection(Pop=Population, f.pop=Evaluation$f)

selection = function(Pop, f.pop){
  # INPUTS
  # Pop: population of permutations (pop.size*p)
  # f.pop: fitness of the population
  # type: type of selection; only proportional selection is available for now
  #
  # OUTPUTS
  # Index vector of the selected elements

  pop.size <- dim(Pop)[1]

  fmin <- min(f.pop)
  fmax <- max(f.pop)
  if (fmax > fmin) {
    weights <- abs((f.pop-fmax)/(fmax-fmin))*0.98 + 0.02 # weights which are associated to the population
  } else {
    weights <- rep(1, pop.size)
  }
  V  <- c(0, cumsum(weights)) / sum(weights)
  return(findInterval(runif(pop.size), V))
}
