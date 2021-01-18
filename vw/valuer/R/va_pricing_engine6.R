#Copyright 2017 Ivan Zoccolan

#This file is part of valuer.

#Valuer is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#Valuer is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#A copy of the GNU General Public License is available at
#https://www.R-project.org/Licenses/ and included in the R distribution
#(in directory share/licenses).


#' Variable Annuity  pricing engine with GBM and generic mortality
#' @description
#' Class providing a variable annuity pricing engine with the underlying
#' reference risk neutral fund modeled as a Geometric Brownian Motion and the
#' intensity of mortality process specified by a generic SDE
#' (stochastic mortality).
#' The value of the VA contract is estimated by means of the Monte Carlo
#' method if the policyholder cannot surrender (the so called "static"
#' approach), and by means of Least Squares Monte Carlo in case the
#' policyholder can surrender the contract (the "mixed" approach).\cr
#' See \bold{References} -\code{[BMOP2011]} for a description of the mixed
#' and static approaches and the algorithm implemented by this class,
#' \code{[LS2001]} for Least Squares Monte Carlo.
#' @docType class
#' @importFrom orthopolynom laguerre.polynomials
#' @importFrom RcppEigen fastLmPure
#' @export
#' @return Object of \code{\link{R6Class}}
#' @format \code{\link{R6Class}} object.
#' @section Methods:
#' \describe{
#'  \item{\code{new}}{Constructor method with arguments:
#'   \describe{
#'    \item{\code{product}}{\code{\link{va_product}} object}
#'    \item{\code{interest}}{\code{\link{constant_parameters}} object
#'    with the interest rate}
#'    \item{\code{spot}}{\code{numeric} scalar with the initial fund price}
#'    \item{\code{volatility}}{\code{\link{constant_parameters}} object with
#'    the volatility}
#'    \item{\code{dividends}}{\code{\link{constant_parameters}} object with
#'    the dividend rate}
#'    \item{\code{mortality_parms}}{A list of parameters
#'    specifying the demographic processes.
#'    See \code{\link{mortality_BMOP2011}} for an example.}
#'   }
#'  }
#'  \item{\code{death_time}}{Returns the time of death index. If the
#'   death doesn't occur during the product time-line it returns the
#'   last index of the product time-line}
#'  \item{\code{simulate_financial_paths}}{Simulates \code{npaths} paths
#'   of the underlying fund of the VA contract and the discount factors
#'   (interest rate) and saves them into private fields for later use.}
#'  \item{\code{simulate_mortality_paths}}{Simulates \code{npaths} paths
#'   of the intensity of mortality and saves them into private fields
#'   for later use.}
#'  \item{\code{get_fund}}{Gets the \code{i}-th path of the underlying fund
#'   where \code{i} goes from 1 to \code{npaths}}
#'  \item{\code{do_static}}{Estimates the VA contract value by means of
#'   the static approach (Monte Carlo), see \bold{References}. It takes as
#'   arguments:
#'  \describe{
#'   \item{\code{the_gatherer}}{\code{gatherer} object to hold
#'    the point estimates}
#'   \item{\code{npaths}}{positive integer with the number of paths to
#'    simulate}
#'   \item{\code{simulate}}{boolean to specify if the paths should be
#'    simulated from scratch, default is TRUE.}
#'   }
#'  }
#'  \item{\code{do_mixed}}{Estimates the VA contract by means of
#'  the mixed approach (Least Squares Monte Carlo), see \bold{References}.
#'  It takes as arguments:
#'  \describe{
#'   \item{\code{the_gatherer}}{\code{gatherer} object to hold
#'    the point estimates}
#'   \item{\code{npaths}}{positive integer with the number of paths to
#'    simulate}
#'   \item{\code{degree}}{positive integer with the maximum degree of
#'    the weighted Laguerre polynomials used in the least squares by LSMC}
#'   \item{\code{freq}}{string which contains the frequency of the surrender
#'    decision. The default is \code{"3m"} which corresponds to deciding
#'    every three months if surrendering the contract or not.}
#'   \item{\code{simulate}}{boolean to specify if the paths should be
#'    simulated from scratch, default is TRUE.}
#'   }
#'  }
#'  \item{\code{get_discount}}{Arguments are \code{i,j}.
#'  Gets the \code{j}-th discount factor corresponding to the \code{i}-th
#'  simulated path of the discount factors. This method must be implemented
#'  by sub-classes.}
#'  \item{\code{fair_fee}}{Calculates the fair fee for a contract using the
#'  bisection method. Arguments are:
#'   \describe{
#'    \item{\code{fee_gatherer}}{\code{\link{data_gatherer}} object to hold
#'    the point estimates}
#'    \item{\code{npaths}}{\code{numeric} scalar with the number of MC
#'    simulations to run}
#'    \item{\code{lower}}{\code{numeric} scalar with the lower fee corresponding
#'    to positive end of the bisection interval}
#'    \item{\code{upper}}{\code{numeric} scalar with the upper fee corresponding
#'    to the negative end of the bisection interval}
#'    \item{\code{mixed}}{\code{boolean} specifying if the mixed method has
#'    to be used. The default is \code{FALSE}}
#'    \item{\code{tol}}{\code{numeric} scalar with the tolerance of the
#'    bisection algorithm. Default is \code{1e-4}}
#'    \item{\code{nmax}}{positive \code{integer} with the maximum number of
#'    iterations of the bisection algorithm}
#'   \item{\code{simulate}}{boolean specifying if financial and mortality
#'    paths should be simulated.}
#'   }
#'  }
#' }
#' @references
#' \enumerate{
#'  \item{[BMOP2011]}{ \cite{Bacinello A.R., Millossovich P., Olivieri A.
  #'  ,Pitacco  E., "Variable annuities: a unifying valuation approach."
  #'  In: Insurance: Mathematics and Economics 49 (2011), pp. 285-297.}}
#'  \item{[LS2001]}{ \cite{Longstaff F.A. e Schwartz E.S. Valuing
  #'  american options by simulation: a simple least-squares approach.
  #'  In: Review of Financial studies 14 (2001), pp. 113-147}}
#' }
#'@examples
#'\dontrun{
#'#Sets up the payoff as a roll-up of premiums with roll-up rate 1%
#'
#'rate <- constant_parameters$new(0.01)
#'
#'premium <- 100
#'rollup <- payoff_rollup$new(premium, rate)
#'
#'#Ten years time-line
#'begin <- timeDate::timeDate("2016-01-01")
#'end <- timeDate::timeDate("2025-12-31")
#'
#'#Age of the policyholder.
#'age <- 60
#'# A constant fee of 4% per year (365 days)
#'fee <- constant_parameters$new(0.04)
#'
#'#Barrier for a state-dependent fee. The fee will be applied only if
#'#the value of the account is below the barrier
#'barrier <- Inf
#'#Withdrawal penalty applied in case the insured surrenders the contract
#'#It is a constant penalty in this case
#'penalty <- penalty_class$new(type = 1, 0.01)
#'#Sets up the contract with GMAB guarantee
#'contract <- GMAB$new(rollup, t0 = begin, t = end, age = age, fee = fee,
#'barrier = barrier, penalty = penalty)
#'
#'#Interest rate
#'r <- constant_parameters$new(0.03)
#'#Initial value of the underlying fund
#'spot <- 100
#'#Volatility
#'vol <- constant_parameters$new(0.2)
#'#Dividend rate
#'div <- constant_parameters$new(0.0)
#'#Gatherer for the MC point estimates
#'the_gatherer <- mc_gatherer$new()
#'#Number of paths to simulate
#'no_of_paths <- 10
#'
#'#Sets up the pricing engine specifying the va_contract, the interest rate
#'#the parameters of the Weibull intensity of mortality, the initial fund
#'#value, the volatility and dividends rate
#'engine <- va_bs_engine2$new(contract, r, spot,
#'volatility=vol, dividends=div, mortality_BMOP2011)
#'
#'#Estimates the contract value by means of the static approach.
#'
#'engine$do_static(the_gatherer, no_of_paths)
#'the_gatherer$get_results()
#'
#'#Estimates the contract value by means of the mixed approach.
#'#To compare with the static approach we won't simulate the underlying
#'#fund paths again.
#'
#'the_gatherer_2 <- mc_gatherer$new()
#'
#'engine$do_mixed(the_gatherer_2, no_of_paths, degree = 3,
#'freq = "3m", simulate = FALSE)
#'the_gatherer_2$get_results()
#'}


va_bs_engine2 <- R6::R6Class("va_bs_engine", inherit = va_bs_engine,
 public = list(
  initialize = function(product, interest, spot, volatility, dividends, mortality_parms, ...){
   if(!missing(product))
    if (inherits(product, "va_product")){
     private$the_product <- product
     private$times <- product$get_times()
     private$no_time_intervals <- length(private$times) - 1
     private$drifts <- numeric(private$no_time_intervals)
     private$standard_deviations <- numeric(private$no_time_intervals)
     private$variates <- numeric(private$no_time_intervals)
    } else stop(error_msg_1_("product", "va_product"))
   else stop(error_msg_1_("product", "va_product"))

   if(!missing(interest))
    if(inherits(interest, "parameters")){
     private$r <- interest
    } else stop(error_msg_1_("interest", "parameters"))
   else stop(error_msg_1_("interest", "parameters"))

   if(!missing(spot))
    if (is_positive_scalar(spot)){
     private$spot <- spot
    } else stop(error_msg_3_("spot"))
   else stop(error_msg_3_("spot"))

   if(!missing(volatility) & !missing(dividends)){
    if(inherits(volatility, "parameters") & inherits(dividends, "parameters")){
     for(j in seq(private$no_time_intervals)){
      this_variance <- volatility$integral_square(private$times[j],
                                                  private$times[j+1])
      drift <- interest$integral(private$times[j], private$times[j+1])
      drift <- drift - dividends$integral(private$times[j], private$times[j+1])
      private$drifts[j] <- drift - 0.5 * this_variance
      private$standard_deviations[j] <- sqrt(this_variance)
     }
    } else stop(error_msg_2("parameters"))
   } else stop(error_msg_2("parameters"))

   private$mortality_parms <- mortality_parms
   private$mortality_model <- do.call(yuima::setModel,
                                       mortality_parms[[2]])
   #
   no_time_intervals <- length(private$the_product$get_times()) - 1
   private$samp <- yuima::setSampling(
    Terminal = tail(private$the_product$times_in_yrs(), 1),
    n = no_time_intervals)
  },
  simulate_mortality_paths = function(npaths){

    ind <- private$mortality_parms[[3]]
    #Sets the arguments to call yuima::simulate
    parms <- list(
      object = private$mortality_model,
      xinit = private$mortality_parms[[1]]$xinit,
      sampling = private$samp,
      true.parameter = private$mortality_parms[[1]]
    )
    #Sets up storage for the intensity of mortality
    #and int of mort integrals paths
    len <-  length(private$the_product$get_times())
    private$mu <- matrix(NA, npaths, len)
    private$mu_integrals <- matrix(NA, npaths, len)
    dt <- private$samp@delta
    #Sets up the intensity of mortality and
    #calculates the integrals of the intensity of mortality
    #Will use foreach to run the simulations in parallel
    #Imported in NAMESPACE
    #yuima::simulate
    #yuima::get.zoo.data
    if(requireNamespace("foreach", quietly = TRUE)){
      loop <- foreach::foreach(iterators::icount(npaths))
      data_paths <- foreach::`%dopar%`(loop, {
        zoo_paths <- do.call(simulate, parms)
        get.zoo.data(zoo_paths)
      })
      for(i in seq(npaths)){
        private$mu[i, ] <- data_paths[[i]][[ind]]
        mu_ <- head(private$mu[i, ], -1)
        private$mu_integrals[i, ] <- cumsum(c(0, mu_ * dt))
      }
    } else for(i in seq(npaths)){
      zoo_paths <- do.call(simulate, parms)
      private$mu[i, ] <- get.zoo.data(zoo_paths)[[ind]]
      mu_ <- head(private$mu[i, ], -1)
      private$mu_integrals[i, ] <- cumsum(c(0, mu_ * dt))

    }
  },
  death_time = function(i){
    ind <- which(private$mu_integrals[i, ] > rexp(1))
    if (length(ind) != 0)
      res <- min(ind)
    else res <- length(private$the_product$get_times()) + 1
    res
  }
 ),
 private = list(
  mortality_model = "yuima.model-class",
  #Stores the demographic parameters needed to
  #set the model above by yuima::setModel and
  #run yuima::simulate
  mortality_parms ="list",
  #Stores the times of a simulate path
  samp = "yuima.sampling-class",
  #Matrix to hold the simulated paths of the
  #stochastic intensity of mortality
  mu = "matrix",
  #Matrix to hold the integrals of the simulated paths
  #of stochastic intensity of mortality
  mu_integrals = "matrix",
  #Method to get Laguerre polynomials of state variables.
  #Arguments are:
  #paths - numeric vector of indexes of the paths to consider
  #time - numeric scalar with the time index
  #degree - positive scalar with the max degree of the Laguerre polynomials
  bases = function(paths, time, degree){
    #orthopolynom::laguerre.polynomials it's imported in NAMESPACE
    #for performance reasons.
    res <- laguerre.polynomials(degree, normalized = TRUE)
    x <- private$fund[paths, time]
    mu <- private$mu[paths, time]
    #Normalizes to avoid underflows in calculating
    #the exponential below.
    x <- x / private$the_product$get_premium()

    x <- sapply(seq_along(res), function(i){
      exp(-0.5 * x) * (as.function(res[[i]])(x))
    })
    mu <- sapply(seq_along(res), function(i){
      exp(-0.5 * mu) * (as.function(res[[i]])(mu))
    })
    cbind(x, mu)
  }
 )
)



