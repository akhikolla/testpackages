#################################################################
##### Mutation operator #####

##' @keywords internal
# @title Mutate permutations
# @description Internal function (mutation operator) to mutate elements of a population of permutations.
# @param Pop Population of permutations from [1,p] (output from create.population() function for example).
# @param p.mut Mutation probability.
# @return Population of permutations created after mutating permutations of Pop.
# @author \packageAuthor{GADAG}
# @seealso \code{\link{GADAG}}, \code{\link{GADAG_Run}}.
# @examples
#  ########################################################
#  # Creating a population of permutations
#  ########################################################
#  Population <- create.population(p=10, pop.size=20)
#
#  ########################################################
#  # Mutating the population
#  ########################################################
#  Population_Mutation <- mutation(Pop=Population)

mutation <- function(Pop, p.mut=0.05){
  # INPUTS
  # Pop: population of permutations (pop.size*p)
  # p.mut: probability of mutation
  #
  # OUTPUTS
  # newPop: population of permutations (pop.size*p)

  if (is.vector(Pop)==TRUE){
    Pop <- t(as.matrix(Pop,ncol=length(Pop),nrow=1))
  }

  pop.size <- dim(Pop)[1]
  p <- dim(Pop)[2]

  newPop <- Pop
  I.mut  <- which(runif(pop.size) < p.mut)

  if (length(I.mut)>0) {
    Popmut <- Pop[I.mut,]
    rank.mut1 <- sample.int(n=(p-1), size=length(I.mut), replace=TRUE)
    rank.mut2 <- sample.int(n=(p-1), size=length(I.mut), replace=TRUE)
    rank.mut <- cbind(rank.mut1,rank.mut2)
    if (length(I.mut)>1){
      for (j in (1:length(I.mut))){
        Popmut[j,rank.mut[j,]] <- Pop[I.mut[j],rev(rank.mut[j,])]
        newPop[I.mut,] <- Popmut
      }
    } else {
      Popmut[rank.mut] <- Pop[I.mut,rev(rank.mut)]
      newPop[I.mut,] <- Popmut
    }
  }
  return(newPop)
}
