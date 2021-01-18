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

##########################DESIGN  COMMENTS ####################################
#Design notes of the variable annuities pricing engine.
#
#The overall design of the VA engine it's inspired by the ideas described in
#M. S. Joshi. 2008. C++ Design Patterns and Derivatives Pricing (2nd ed.).
#Cambridge University Press, New York, NY, USA.
#My thanks to the Author for his coincise and elegant description of a
#design for a pricing engine of a path dependent derivative.
#
#The pricing engine for a variable annuity needs to extend the cited design in
#two directions. First of all, it has to take into account not only the
#financial risk implicit in the underlying fund but also the demographic risk.
#Secondly it needs a method  for pricing in case the policyholder is
#allowed to surrender the contract. This is similar to price a path dependent
#American derivative via Least Squares Monte Carlo (LSMC).
#Without entering into much details this is called the "mixed" approach by
#Bacinello,  Millossovich, Olivieri and Pitacco in their article:
#"Variable annuities: a unifying valuation approach."
#Insurance: Mathematics and Economics 49 #(2011), pp. 285-297.
#while the "static" approach is the one which  doesn't allow the policyholder
#to surrender.
#The engine returns both the static and mixed approach calculations
#by means of its methods do_static and do_mixed.
#The method do_mixed implements the algorithm described by Bacinello et All.
#The method do_static implements a standard Monte Carlo (MC) routine.
#My thanks to Bacinello et All for their brilliant framework to price
#VA contracts.
#All mistakes with this implementation of their ideas is attributable
#to me only.
#To compare the values of the contracts obtained from both approaches
#do_mixed and do_static can be applied on the same simulated paths.
#Both do_mixed and do_static methods are passed a statistics gatherer object
#to collect the point estimates. The implementation of a statistic gatherer
#quickly suggested the use of classes with reference semantics and that's
#mainly why I begun to use R6.
#
#The basic design principle is that the engine knows how to simulate the
#financial and demographic processes but does not know anything about the
#va contract itself.
#Info about the contract such as the type of guarantee, expiry, dates
#of due payments, age of the insured are saved into a va_product object.
#This object is passed to the pricing engine at initialization.
#The pricing engine pulls the info about the contract through the va_product
#public methods. See va_products.R
#This enables the pricing of any VA contract we want since we just need to
#write the appropriate va_product subclass and pass the corresponding
#object to the engine.

#The pricing algorithms are defined in the base class va_engine.
#Vise-versa the simulation of the source of risks is demanded to subclasses of
#this base class, such as  va_bs_engine and va_sde_engine.
#va_bs_engine simulates the underlying fund with a GBM and the demographic risk
#by means of the deterministic Weibull intensity of mortality.
#The class va_sde_engine leverages on the beautiful yuima package to specify
#both the financial and demographic processes by means of
#a system of stochastic differential equations.
#So it is pretty general and flexible in that respect.
#However bear in mind the cost of this flexibility is far longer execution
#time and assumes some fair knowledge of yuima.
#Please see the help pages for financials_BMOP2011, mortality_BMOP2011,
#financials_BBM2010 and mortality_BBM2010.
#To improve the execution time of the simulation of paths in the va_sde_engine,
#parallel execution of the yuima simulate function is implemented via
#foreach. However paraller execution proved to be voracious for memory
#in Windows (SOCK clusters). I recommend to use it in Linux as fork clusters
#are more efficient in terms of memory utilization.
#To make things easier for the user I've provided  example parameters
#to set up the engine.
#The simulated paths of the underlying fund, interest rate, discount factors,
#intensity of mortality and VA cash flows are saved into private matrices.
#As a consequence the  memory foot print of the engine object increases
#with the number of simulations and may get conspicuous.
#
#The public method death_time simulates the death times.
#This method gets from a public method simulate_mortality_path() the
#simulated paths  of the integrals of the intensity of mortality
#over the product time-line and uses them  to simulate the death time.
#simulate_mortality_paths() is a standard interface of the base VA engine
#which will be implemented in sub classes depending we want  stochastic
#mortality as in va_sde_engine or deterministic mortality as in va_bs_engine
#With the deterministic choice the intensity of mortality
#and intensity of mortality integral paths don't change so they will be
#calculated during initialization and saved into a private field that can be
#used later on ( e.g: by  death_time).
#So a sub class for deterministic mortality has a private field to hold
#the mortality path and the initialize method is overridden to call
#simulate_financial_paths() and save the path once for all.
#The implementation of death_time is overridden to get the calculated
#mortality integrals from the private field.
#Sub classes implementing stochastic mortality calls simulate_financial_paths()
#to get a number of simulated path and save them into the private matrix.
#Method death_time corresponding to a particular i-th path will have to fetch
#the i-th path from the private matrix.
#The pricer passes the death time to the cash_flows method of the VA product
#so that it will return the cash_flows up to death time. See the design
#comments in va_product.R
#
######################END DESIGN  COMMENTS#####################################

#' Generic Variable Annuity  pricing engine
#' @description
#' Class providing an interface for a generic VA pricing engine.\cr
#' This class shouldn't be instantiated but used as base class for
#' variable annuity pricing engines.
#' The value of the VA contract is estimated by means of the Monte Carlo
#' method if the policyholder cannot surrender (the so called "static"
#' approach), and by means of Least Squares Monte Carlo in case the
#' policyholder can surrender the contract (the "mixed" approach).\cr
#' See \bold{References} -\code{[BMOP2011]} for a description of the mixed
#' and static approaches and the algorithm implemented by this class,
#' \code{[LS2001]} for Least Squares Monte Carlo.
#' @docType class
#' @importFrom R6 R6Class
#' @importFrom orthopolynom laguerre.polynomials
#' @importFrom RcppEigen fastLmPure
#' @return Object of \code{\link{R6Class}}
#' @format \code{\link{R6Class}} object.
#' @section Methods:
#' \describe{
#'  \item{\code{new}}{Constructor method}
#'  \item{\code{death_time}}{Returns the time of death index. If the
#'  death doesn't occur during the product time-line it returns the
#'  last index of the product time-line plus one.}
#'  \item{\code{simulate_financial_paths}}{Simulates \code{npaths} paths
#'   of the underlying fund of the VA contract and the discount factors
#'   (interest rate) and saves them into private fields for later use.}
#'  \item{\code{simulate_mortality_paths}}{Simulates \code{npaths} paths
#'   of the intensity of mortality and saves them into private fields
#'   for later use.}
#'  \item{\code{get_fund}}{Gets the \code{i}-th path of the underlying fund
#'   where \code{i} goes from 1 to \code{npaths}}
#'  \item{\code{do_static}}{Estimates the VA contract value by means of
#'  the static approach (Monte Carlo), see \bold{References}. It takes as
#'  arguments:
#'   \describe{
#'     \item{\code{the_gatherer}}{\code{gatherer} object to hold
#'     the point estimates}
#'     \item{\code{npaths}}{positive integer with the number of paths to
#'     simulate}
#'     \item{\code{simulate}}{boolean to specify if the paths should be
#'     simulated from scratch, default is TRUE.}
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


va_engine <- R6::R6Class("va_engine",
 public = list(
  initialize = function(product, ...){
   if(!missing(product))
    if (inherits(product, "va_product")) {
     private$the_product <- product
    } else stop(error_msg_1_("product", "va_product"))
   else stop(error_msg_1_("product", "va_product"))
  },
  simulate_mortality_paths = function(npaths){

  },
  death_time = function(i){

  },
  simulate_financial_paths = function(npaths){

  },
  get_fund = function(i) {

  },
  do_static = function(the_gatherer, npaths, simulate = TRUE){
   #Estimates the VA by means of the static approach
   #implemented with Monte Carlo
   times_len = length(private$the_product$get_times())
   ind <- seq(npaths)

   if(simulate){
    #Simulates financials
    self$simulate_financial_paths(npaths)
    #Simulates mortality paths
    self$simulate_mortality_paths(npaths)
  }

  #Simulates times of death
  private$tau <- sapply(ind, function(i){
    self$death_time(i)})

   #Initial cash flow matrix
   cash <- matrix(0, nrow = npaths, ncol = times_len)
   #Adjusted time of death to avoid 'subscript out of bounds'
   #errors with the cash matrix
   adj_tau <- sapply(ind, function(i) ifelse(private$tau[i] > times_len,
                                             times_len, private$tau[i]))
   #With the static method we set the penalty as 1
   #since we  want the surrender values to be all zeros.
   old_penalty <- private$the_product$get_penalty_object()
   penalty_1 <- penalty_class$new(type =1, const = 1)
   private$the_product$set_penalty_object(penalty_1)
   for (i in ind){
     #Discount factors from time of death
     #This is used by GMIB (Ib, Ic) and GMWB (Wa) riders for
     #the death benefit at death_time.
    discounts <- self$get_discount(i) / self$get_discount(i, adj_tau[i])

    cash[i, ] <- private$the_product$cash_flows(self$get_fund(i),
                                                private$tau[i], discounts)
   }
   res <- sapply(ind, function(i){
     sum(cash[i, 1:adj_tau[i]] *
           self$get_discount(i, 1:adj_tau[i]))
   })
   private$the_product$set_penalty_object(penalty = old_penalty)
   the_gatherer$dump_result(res)
  },
  do_mixed = function(the_gatherer, npaths, degree = 3, freq = "3m", simulate = TRUE, adjust = TRUE){
   #Estimates the VA by means of the mixed approach
   #implemented with Least Squares Monte Carlo
   times_len <- length(private$the_product$get_times())
   ind <- seq(npaths)
   if(simulate){
   #Simulates financials
   self$simulate_financial_paths(npaths)
   #Simulates mortality paths
   self$simulate_mortality_paths(npaths)
   }

   #Simulates times of death
   private$tau <- sapply(ind, function(i){
     self$death_time(i)})

   surrender_times <- private$the_product$surrender_times(freq)
   survival_times <- private$the_product$survival_benefit_times()
   tt <- c(surrender_times, survival_times)

   #Initial cash flow matrix
   cash <- matrix(0, nrow = npaths, ncol = times_len)

   #Adjusted time of death to avoid 'subscript out of bounds'
   #errors with the cash matrix
   adj_tau <- sapply(ind, function(i) ifelse(private$tau[i] > times_len,
                                             times_len, private$tau[i]))
   #Initial surrender times vector
   sur_ts <- adj_tau
   ############################

   for (i in ind){
     #Discount factors from time of death
     #This is used by GMIB (Ib, Ic) and GMWB (Wa) riders for
     #the death benefit at death_time.
     discounts <- self$get_discount(i) / self$get_discount(i, adj_tau[i])
     cash[i, ] <- private$the_product$cash_flows(self$get_fund(i),
                                                 private$tau[i], discounts)

     #Sets all zero in the cash flow matrix but for times corresponding
     #to surrender decision, living benefit times and death time.
     cash[i, -sort(unique(c(tt, private$tau[i])))] <- 0
   }

   #In case the penalty is  decreasing with time
   #and the financial paths are continuous (e.g: no jump diffusion),
   #recent studies proved that the rational policyholder would decide
   #to surrender only when the account value is below the barrier.
   #Therefore in these hypotheses for the LSMC regression we want to
   #take into consideration only paths where the insured is alive and
   #the  account is below the barrier.
   #The bound parameter will then be the barrier by the penalty so that
   #it can be compared with the cash matrix.

   if(adjust)
     bound <- private$the_product$get_barrier()
   else bound <- Inf

   bound <- rep(bound, length(surrender_times))

   if(inherits(private$the_product, "GMWB")) {
     surv_times <- survival_times[survival_times %in% surrender_times]
     bound[surv_times] <- bound[surv_times] - cash[1, 1]
   }

   penalty <- private$the_product$get_penalty()

   if (length(penalty) != 1) {
      bound <- (1 - penalty[surrender_times]) * bound
    } else bound <- (1 - penalty[1]) * bound

   for(t in rev(surrender_times)){

    t_idx <- which(t == surrender_times)
    #See comment above about bound
    h_t <- which(private$tau > t & cash[, t] < bound[t_idx])
    not_h_t <- which(private$tau > t & cash[, t] >= bound[t_idx])
    if(length(h_t) != 0){
    #Continuation value at time t
     c_t <- sapply(h_t, function(i){
       sum(cash[i, (t+1):sur_ts[i]] *
        self$get_discount(i, (t+1):sur_ts[i]) / self$get_discount(i, t))
       })
     #### Regression ####
     #Regressors
     x_t <- private$bases(h_t, t, degree)
     #Estimated continuation values
     #RcppEigen::fastLmPure it's imported in NAMESPACE
     chat_t <- fastLmPure(x_t, c_t)$fitted.values
     #### Comparison between surrender values and estimated
     #### continuation values ####
     for (i in seq_along(h_t)){
      surv_ben <- private$the_product$survival_benefit(self$get_fund(h_t[i]),
                                                       private$tau[h_t[i]], t)

      surrender <-  cash[h_t[i], t] -  surv_ben

      if (surrender > chat_t[i])
       sur_ts[h_t[i]] <- t
      else cash[h_t[i], t] <- surv_ben
     }
    }
    #We need to set the cashflow at time t back to just the survival benefit
    #(so zero clear the surrender value) given we're not taking the decision
    #to surrender on all paths but those indexed by h_t
     for (i in not_h_t){
      surv_ben <- private$the_product$survival_benefit(self$get_fund(i),
                                                       private$tau[i], t)
      cash[i, t] <- surv_ben
     }
   }
   res <- sapply(seq_along(sur_ts), function(i) {
     sum(cash[i, 1:sur_ts[i]] *
       self$get_discount(i, 1:sur_ts[i]))
     })
   the_gatherer$dump_result(res)
  },
  get_discount = function(i, j) {

    },
  fair_fee = function(fee_gatherer, npaths, lower, upper, mixed = FALSE, tol = 1e-4, nmax = 100, simulate = TRUE){
    #Argument checks
    if(!missing(fee_gatherer))
     if(inherits(fee_gatherer, "data_gatherer")){}
     else stop(error_msg_1_("fee_gatherer", "data_gatherer"))
    else stop(error_msg_1_("fee_gatherer", "data_gatherer"))

    if(!is_positive_integer(npaths)) stop(error_msg_4("npaths"))
    if(!is_not_negative_scalar(lower)) stop(error_msg_7("lower"))
    if(!is_positive_scalar(upper)) stop(error_msg_3_("upper"))
    if(upper <= lower) stop("Argument upper should be greater than lower")
    if(!is.logical(mixed)) stop("Argument mixed should be a logical")
    if(!is_between(tol, 0, 1)) stop(error_msg_8("tol"))
    if(!is_positive_integer(nmax)) stop(error_msg_4("nmax"))
    #######

    P <- private$the_product$get_premium()

    #Finds the contract value corresponding to the lower fee
    private$the_product$set_fee(constant_parameters$new(lower))

    the_gatherer <- mc_gatherer$new()

    if(mixed) self$do_mixed(the_gatherer, npaths, simulate)
    else self$do_static(the_gatherer, npaths, simulate)

    old_value <- the_gatherer$get_results()[[1]]  - P

    if(old_value < 0) stop("Please decrease lower\n")

    n <- 1
    while(n <= nmax){
     cat("#")
     current <- (lower + upper) / 2
     private$the_product$set_fee(constant_parameters$new(current))

     fee_gatherer$dump_result(current)

     if(mixed) self$do_mixed(the_gatherer, npaths, simulate = FALSE)
     else self$do_static(the_gatherer, npaths, simulate = FALSE)

     current_value <- the_gatherer$get_results()[[1]] - P
     se <- the_gatherer$get_results()[[2]]
     #Debug
     cat("Old value: ", old_value, "\n")
     cat("Current value", current_value, "\n")

     if(abs(current_value) <= se | (upper - lower) / 2 <= tol){
      cat("\nFair fee is", current)
      break
     } else {
      n <- n + 1
      if (sign(current_value) == sign(old_value)){
       lower <- current
       old_value <- current_value
      } else upper <- current
     }
    }
   if(n > nmax){
     warning("The bisection method did not converge.\nCurrent fee value: ", current)
   }
  }
 ),
 private = list(
  #A va_product object with the VA contract.
  the_product = "va_product",
  #A numeric vector with the simulated times of death
  tau = "numeric",
  #Method to get Laguerre polynomials of state variables.
  #Arguments are:
  #paths - numeric vector of indexes of the paths to consider
  #time - numeric scalar with the time index
  #degree - positive scalar with the max degree of the Laguerre polynomials
  #To be implemented in sub-classes
  bases = function(paths, time, degree){
    NULL
  }
 )
)



#' Variable Annuity  pricing engine with GBM
#' @description
#' Class providing a variable annuity pricing engine with the underlying
#' reference risk neutral fund modeled as a Geometric Brownian Motion and the
#' intensity of mortality  modeled by the Weibull intensity of mortality.
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
#'    \item{\code{c1}}{\code{numeric} scalar argument of the intensity
#'    of mortality function \code{\link{mu}}}
#'    \item{\code{c2}}{\code{numeric} scalar argument of the intensity
#'    of mortality function \code{\link{mu}}}
#'    \item{\code{spot}}{\code{numeric} scalar with the initial fund price}
#'    \item{\code{volatility}}{\code{\link{constant_parameters}} object with
#'    the volatility}
#'    \item{\code{dividends}}{\code{\link{constant_parameters}} object with
#'    the dividend rate}
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
#'#the parameters of the Weibull intensity of mortality, the initial fund
#'#value, the volatility and dividends rate
#'engine <- va_bs_engine$new(contract, r, c1=90.43, c2=10.36, spot,
#'volatility=vol, dividends=div)
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




va_bs_engine <- R6::R6Class("va_bs_engine", inherit = va_engine,
 public = list(
  initialize = function(product, interest, c1, c2, spot, volatility, dividends, ...){
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

   if(!missing(c1))
    if (is_positive_scalar(c1))
      private$mu_1 <- c1
    else stop(error_msg_5("c1"))
   else private$mu_1 <- 90.43

   if(!missing(c2))
    if(is_positive_scalar(c2))
      private$mu_2 <- c2
    else stop(error_msg_5("c2"))
   else private$mu_2 <- 10.36

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
    }else stop(error_msg_2("parameters"))
   } else stop(error_msg_2("parameters"))

   private$mu_integrals <- self$simulate_mortality_paths()
  },
  simulate_financial_paths = function(npaths){
    cf_times <- private$times
    private$fund <- t(vapply(seq(npaths), function(index) {
      private$simulate_financial_path()
    }, FUN.VALUE = vector("numeric", length(cf_times))))
    t0 <- cf_times[1]
    log_discounts <- vector(mode="numeric", length = length(cf_times))
    for (i in seq_along(cf_times))
      log_discounts[i] <- -private$r$integral(t0, cf_times[i])
    private$discounts <- exp(log_discounts)
  },
  get_fund = function(i) private$fund[i, ],
  get_discount = function(i, j) private$discounts[j],
  simulate_mortality_paths = function(npaths){
    age <- private$the_product$get_age()
    c1 <- private$mu_1
    c2 <- private$mu_2
    t_yrs  <- head(private$the_product$times_in_yrs(), -1)
    #Deterministic intensity of mortality (Weibull)
    integrand <- function(t) {
      (c1 ^ (-c2)) * c2 * ((age + t) ^ (c2 - 1))
    }
    mu_integrals <- sapply(t_yrs, function(u) {
      stats::integrate(integrand, 0, u)$value
    })
    mu_integrals
  },
  death_time = function(i){
    ind <- which(private$mu_integrals > rexp(1))
    if (length(ind) != 0)
      res <- min(ind)
    else res <- length(private$times) + 1
    res
  }
 ),
 private = list(
  #Product times
  times = "timeDate",
  #parameters object with the risk-neutral interest rate
  r = "parameters",
  #numeric vector with the drifts of the underlying fund
  drifts = "numeric",
  #numeric vector with the standard_deviations of the
  #underlying fund
  standard_deviations = "numeric",
  #A numeric vector holding the simulated standard normal values
  variates = "numeric",
  #numeric scalar with the initial value of the underlying fund
  spot = "numeric",
  no_time_intervals = "numeric",
  #numeric vector with the integrals along the product
  #time-line of the intensity of mortality
  mu_integrals = "numeric",
  #matrix with  simulated paths of the underlying fund
  fund = "matrix",
  #numeric vector with the discount factors
  discounts = "numeric",
  #Intensity of mortality parameters
  mu_1 = 90.43,
  mu_2 = 10.36,
  #Simulates one path of the underlying fund
  simulate_financial_path = function(){
    private$variates <- rnorm(private$no_time_intervals)
    current_log_spot <- private$drifts +
      private$standard_deviations * private$variates
    current_log_spot <- cumsum(current_log_spot)
    c(private$spot, private$spot*exp(current_log_spot))
  },
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
    #Normalizes to avoid underflows in calculating
    #the exponential below.
    x <- x / private$the_product$get_premium()
    sapply(seq_along(res), function(i){
      exp(-0.5 * x) * (as.function(res[[i]])(x))
    })
  }
 )
)


