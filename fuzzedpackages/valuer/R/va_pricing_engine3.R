#Copyright 2016 Ivan Zoccolan

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


#' Variable Annuity  pricing engine with GBM and Makeham
#' @description
#' Class providing a variable annuity pricing engine with the underlying
#' reference risk neutral fund modeled as a Geometric Brownian Motion and the
#' intensity of mortality  modeled by the Makeham intensity of mortality.
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
#'    \item{\code{A}}{\code{numeric} scalar argument of the intensity
#'    of mortality function \code{\link{makeham}}}
#'    \item{\code{B}}{\code{numeric} scalar argument of the intensity
#'    of mortality function \code{\link{makeham}}}
#'    \item{\code{spot}}{\code{numeric} scalar with the initial fund price}
#'    \item{\code{volatility}}{\code{\link{constant_parameters}} object with
#'    the volatility}
#'    \item{\code{dividends}}{\code{\link{constant_parameters}} object with
#'    the dividend rate}
#'    \item{\code{c}}{\code{numeric} scalar argument of the intensity
#'    of mortality function \code{\link{makeham}}}
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
#'    \item{\code{simulate}}{boolean specifying if financial and mortality
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
#'no_of_paths <- 1e2
#'
#'#Sets up the pricing engine specifying the va_contract, the interest rate
#'#the parameters of the Makeham intensity of mortality, the initial fund
#'#value, the volatility and dividends rate
#'engine <- va_mkh_engine$new(contract, r, A = 0.0001, B = 0.00035, spot,
#'volatility = vol, dividends = div, c = 1.075)
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


va_mkh_engine <- R6::R6Class("va_mkh_engine", inherit = va_bs_engine,
 public = list(
  initialize = function(product, interest, A, B, spot, volatility, dividends, c){

    ifelse(!missing(A),  c1 <- 0.0001, c1 <- A)
    ifelse(!missing(B), c2 <- 0.00035, c2 <- B)
    if(missing(c)) c <- 1.075

    super$initialize(product, interest, c1, c2, spot, volatility, dividends)

    if(!missing(c))
     if(is_positive_scalar(c))
       private$mu_3 <- c
     else stop(error_msg_5("c"))
    else stop(error_msg_5("c"))

  },
  simulate_mortality_paths = function(npaths){
   age <- private$the_product$get_age()
   A <- private$mu_1
   B <- private$mu_2
   c <- private$mu_3
   t_yrs  <- head(private$the_product$times_in_yrs(), -1)
   #Makeham intensity of mortality
   integrand <- function(x) {A + B * (c ^ (age + x))}

   mu_integrals <- sapply(t_yrs, function(u) {
     stats::integrate(integrand, 0, u)$value
   })
   mu_integrals
  }
 ),
 private = list(
  #Third parameter of the Makeham function
  mu_3 = 1.075
 )
)
