# Mortality process -------------------------------------------------------

#' @title Mortality Process
#' @description This functions performs the 'mortality' process over an object,
#' decreasing the number of individuals. It is a generic, S3 methods can be specified
#' for a particular specification of the population. 
#' @param object The population object, containing the information about individuals.
#' @param rates The mortality rate or rates.
#' @param \dots Additional arguments for different methods.
#' @details The rate can be a single value or a value for each individual calculated
#' externally. No recycling is allowed.
mortality = function(object, rates, ...) {
  UseMethod("mortality")  
}


#' @export
#' @method mortality default
mortality.default = function(object, rates, survivors=FALSE, ...) {
  n = object
  if(any(length(n)!=1, !is.numeric(n))) stop("You must provide a integer population size.")
  pop = seq_len(n)
  rates = .checkRates(rates, n)
  dead = which(rbinom(n = n, size = 1, prob=rates)==1)
  if(isTRUE(survivors)) return(setdiff(pop, dead)) else return(dead) 
}

#' @export
#' @method mortality matrix
mortality.matrix = function(object, rates, n=NULL, ...) {
  if(any(is.null(n), n>length(object))) n = nrow(object)
  dead = mortality(object=n, rates=rates, survivors=FALSE)
  if(length(dead)!=0) object[dead, ] = NA
  return(object)
}

# #' @export
# #' @method mortality populationMatrix
# mortality.populationMatrix = function(object, rates) {
#   n = object
#   if(is.null(n)) n = tip
#   pop = seq_len(n)
#   rates = .checkRates(rates, n)
#   deaths = rbinom(n = n, size = 1, prob=rates)
#   ndeath = sum(deaths)
#   if(ndeath==0) return(object)
#   if(ndeath==n) {
#     object[, ] = NA
#     attr(object, which="n") = 0
#     return(object)
#   }
#   die = which(deaths==1)
#   object[die, ] = NA
#   attr(object, which="dead") = die
#   return(object)
# }
