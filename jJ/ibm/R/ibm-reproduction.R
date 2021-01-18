# Reproduction process ----------------------------------------------------


#' @title Reproduction Process
#' @description This functions performs the 'reproduction' process over an object,
#' increasing the number of individuals. It is a generic, S3 methods can be specified
#' for a particular specification of the population. 
#' @param object The population object, containing the information about individuals.
#' @param rates The reproduction rate or rates.
#' @param \dots Additional arguments for different methods.
#' @details The rate can be a single value or a value for each individual calculated
#' externally. No recycling is allowed.
reproduction = function(object, rates, ...) {
  UseMethod("reproduction")
}

#' @export
#' @method reproduction default
reproduction.default = function(object, rates, random=TRUE, newborns=1, ...) {
  n = object
  if(any(length(n)!=1, !is.numeric(n))) stop("You must provide a integer population size.")
  pop = seq_len(n)
  rates = .checkRates(rates, n)
  if(isTRUE(random)) {
    pop = sample(pop, replace=TRUE)
    rates = rates[pop]
  }
  new = rbinom(n = n, size = newborns, prob=rates)
  new = rep(pop, times=new)
  return(new)
}

#' @export
#' @method reproduction matrix
reproduction.matrix = function(object, rates, random=TRUE, newborns=1, n=NULL,
                               index.return=FALSE, ...) {
  if(any(is.null(n), n>length(object))) n = nrow(object)
  new = reproduction(object=n, rates=rates, random=random, newborns=newborns)
  if(isTRUE(index.return)) return(new)
  if(length(new)!=0) object[seq(from=n+1, to=n+length(new)), ] = object[new, ]
  return(object)
}

# #' @export
# #' @method reproduction populationMatrix
# reproduction.populationMatrix = function(object, rates, random=TRUE, newborns=1) {
#   n = object
#   pop = seq_len(n)
#   rates = .checkRates(rates, n)
#   if(isTRUE(random)) {
#     pop = sample(pop, replace=TRUE)
#     rates = rates[pop]
#   }
#   new = rbinom(n = n, size = newborns, prob=rates)
#   new = rep(pop, times=new)
#   return(new)
# }



