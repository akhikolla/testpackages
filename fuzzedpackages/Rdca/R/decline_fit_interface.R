#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @useDynLib Rdca, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom Rdpack reprompt
#' @import minpack.lm
## usethis namespace: end
NULL


# ********** Optimization Models **********


exponential_fit <- function(decline_fit_lst, time_lst) {

   time <- as.vector(time_lst$t)
   input_unit <- as.character(decline_fit_lst$input_unit)
   output_unit <- as.character(decline_fit_lst$output_unit)
   fluid <- as.character(decline_fit_lst$fluid)
   model <- as.character(decline_fit_lst$model)
   fit_data <- as.character(decline_fit_lst$fit_data)
   prod <- as.vector(decline_fit_lst$prod_data)
   par <- as.vector(decline_fit_lst$initial_parameters)
   lower <- as.vector(decline_fit_lst$lower)
   upper <- as.vector(decline_fit_lst$upper)
   control <- decline_fit_lst$control
   if (length(time) != length(prod)) {
      stop("'prod_data' and 'time/date' must have equal length")
   }
   error_fun <- function(par, model, input_unit, output_unit, fluid, fit_data, prod, time) {

      decline_lst <- list(input_unit = input_unit, output_unit = output_unit, fluid = fluid, qi = par[1], Di = par[2], b = 0)
      dca <- exponential(decline_lst, time)
      if (fit_data == "rate") {
         SSE <- (dca[,2] - prod) ^ 2
      } else {
         SSE <- (dca[,3] - prod) ^ 2
      }
      return(SSE)
   }
   nls.out <- nls.lm(par = par, fn = error_fun, model = model, input_unit = input_unit,
                     output_unit = output_unit, fluid = fluid, fit_data = fit_data,
                     prod = prod, time = time, lower = lower, upper = upper, control = control)
   par_estimated <- nls.out$par
   lst <- list(input_unit = input_unit, output_unit = output_unit, fluid = fluid, qi = par_estimated[1],
               Di = par_estimated[2], b = 0)
   attr(lst, "nls.out") <- nls.out
   class(lst) <- c("exponential", "decline")
   return(lst)
}


harmonic_fit <- function(decline_fit_lst, time_lst) {

   time <- as.vector(time_lst$t)
   input_unit <- as.character(decline_fit_lst$input_unit)
   output_unit <- as.character(decline_fit_lst$output_unit)
   fluid <- as.character(decline_fit_lst$fluid)
   model <- as.character(decline_fit_lst$model)
   fit_data <- as.character(decline_fit_lst$fit_data)
   prod <- as.vector(decline_fit_lst$prod_data)
   par <- as.vector(decline_fit_lst$initial_parameters)
   lower <- as.vector(decline_fit_lst$lower)
   upper <- as.vector(decline_fit_lst$upper)
   control <- decline_fit_lst$control
   if (length(time) != length(prod)) {
      stop("'prod_data' and 'time/date' must have equal length")
   }

   error_fun <- function(par, model, input_unit, output_unit, fluid, fit_data, prod, time) {

      decline_lst <- list(input_unit = input_unit, output_unit = output_unit, fluid = fluid, qi = par[1], Di = par[2], b = 1)
      dca <- harmonic(decline_lst, time)
      if (fit_data == "rate") {
         SSE <- (dca[,2] - prod) ^ 2
      } else {
         SSE <- (dca[,3] - prod) ^ 2
      }
      return(SSE)
   }
   nls.out <- nls.lm(par = par, fn = error_fun, model = model, input_unit = input_unit,
                     output_unit = output_unit, fluid = fluid, fit_data = fit_data,
                     prod = prod, time = time, lower = lower, upper = upper, control = control)

   par_estimated <- nls.out$par
   lst <- list(input_unit = input_unit, output_unit = output_unit, fluid = fluid, qi = par_estimated[1],
               Di = par_estimated[2], b = 1)
   attr(lst, "nls.out") <- nls.out
   class(lst) <- c("harmonic", "decline")
   return(lst)
}



hyperbolic_fit <- function(decline_fit_lst, time_lst) {

   time <- as.vector(time_lst$t)
   input_unit <- as.character(decline_fit_lst$input_unit)
   output_unit <- as.character(decline_fit_lst$output_unit)
   fluid <- as.character(decline_fit_lst$fluid)
   model <- as.character(decline_fit_lst$model)
   fit_data <- as.character(decline_fit_lst$fit_data)
   prod <- as.vector(decline_fit_lst$prod_data)
   par <- as.vector(decline_fit_lst$initial_parameters)
   lower <- as.vector(decline_fit_lst$lower)
   upper <- as.vector(decline_fit_lst$upper)
   control <- decline_fit_lst$control
   if (length(time) != length(prod)) {
      stop("'prod_data' and 'time/date' must have equal length")
   }
   error_fun <- function(par, model, input_unit, output_unit, fluid, fit_data, prod, time) {

      decline_lst <- list(input_unit = input_unit, output_unit = output_unit, fluid = fluid, qi = par[1], Di = par[2], b = par[3])
      dca <- hyperbolic(decline_lst, time)
      if (fit_data == "rate") {
         SSE <- (dca[,2] - prod) ^ 2
      } else {
         SSE <- (dca[,3] - prod) ^ 2
      }
      return(SSE)
   }
   nls.out <- nls.lm(par = par, fn = error_fun, model = model, input_unit = input_unit,
                     output_unit = output_unit, fluid = fluid, fit_data = fit_data,
                     prod = prod, time = time, lower = lower, upper = upper, control = control)
   par_estimated <- nls.out$par
   lst <- list(input_unit = input_unit, output_unit = output_unit, fluid = fluid, qi = par_estimated[1],
               Di = par_estimated[2], b = par_estimated[3])
   attr(lst, "nls.out") <- nls.out
   class(lst) <- c("hyperbolic", "decline")
   return(lst)
}



modified_hyperbolic_fit <- function(decline_fit_lst, time_lst) {

   time <- as.vector(time_lst$t)
   input_unit <- as.character(decline_fit_lst$input_unit)
   output_unit <- as.character(decline_fit_lst$output_unit)
   fluid <- as.character(decline_fit_lst$fluid)
   model <- as.character(decline_fit_lst$model)
   fit_data <- as.character(decline_fit_lst$fit_data)
   prod <- as.vector(decline_fit_lst$prod_data)
   par <- as.vector(decline_fit_lst$initial_parameters)
   lower <- as.vector(decline_fit_lst$lower)
   upper <- as.vector(decline_fit_lst$upper)
   control <- decline_fit_lst$control
   if (length(time) != length(prod)) {
      stop("'prod_data' and 'time/date' must have equal length")
   }
   error_fun <- function(par, model, input_unit, output_unit, fluid, fit_data, prod, time) {

      decline_lst <- list(input_unit = input_unit, output_unit = output_unit, fluid = fluid,
                          qi = par[1], Di = par[2], b = par[3], Dt = par[4])
      dca <- modified_hyperbolic(decline_lst, time)
      if (fit_data == "rate") {
         SSE <- (dca[,2] - prod) ^ 2
      } else {
         SSE <- (dca[,3] - prod) ^ 2
      }
      return(SSE)
   }
   nls.out <- nls.lm(par = par, fn = error_fun, model = model, input_unit = input_unit,
                     output_unit = output_unit, fluid = fluid, fit_data = fit_data,
                     prod = prod, time = time, lower = lower, upper = upper, control = control)
   par_estimated <- nls.out$par
   lst <- list(input_unit = input_unit, output_unit = output_unit, fluid = fluid, qi = par_estimated[1],
               Di = par_estimated[2], b = par_estimated[3], Dt = par_estimated[4])
   attr(lst, "nls.out") <- nls.out
   class(lst) <- c("modified_hyperbolic", "decline")
   return(lst)
}
