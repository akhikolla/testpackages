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


######################### DESIGN COMMENTS ######################################
#Implementation of VA products
#A base class for riders is needed to make a standard interface as each product
#has to store the same info such as fee, barrier (for state-dependent fees),
#penalties for withdrawals, type of guarantee payoff (rollup / ratchet).
#So a base class for VA products called va_product is defined. This class
#provides the interface only and should not be instantiated.
#A new va_product subclass is needed for each VA contract rider
#For example we have a specialized GMAB va_product class for VA with
#GMAB rider, etc.
#A roll-up payoff or ratchet payoff object will be passed into the initialize
#of the product.
#The cash_flows method returns all possible  cash_flows, so withdrawal in case
#the insured surrenders the contract or if there's a GMWB rider
#as well as other living / death benefits depending on the riders.
#The engine will use the cash_flows depending if we're doing static or mixed
#The cash flows will be saved in private field of the engine.
#To calculate the cash flows we need first to calculate the account of the
#insured. The formula which calculates the account cannot be vectored since
#each value depends on the previous one.
#So for performance reasons it is implemented in C++ and interfaced with Rcpp
#This is the calc_account function.
#A simple state-dependent fee structure with a single barrier is implemented
#in the calc_account function.
#Public methods will return the times of survival benefit payments
#and the possible surrender times.
#A public method returns the survival benefit at any given time.
#The surrender penalty can be either constant or a decreasing function of time.
##########################DESIGN COMMENTS END###################################


#Defines a base class for a product which will be inherited by
#specialized classes. This is  exported but should not be instantiated.

#' Generic Variable Annuity  product class
#' @description  Class providing an interface for a generic VA product object.
#' This class shouldn't be instantiated but used as base class for
#' implementing  products with contract riders such as GMAB, GMIB, etc.
#' It supports a simple state-dependent fee structure with a single barrier.\cr
#' See \bold{References} for a description of variable annuities life
#' insurance products, their guarantees and fee structures.
#' @docType class
#' @return Object of \code{\link{R6Class}}
#' @format \code{\link{R6Class}} object.
#' @section Methods:
#'  \describe{
#'   \item{\code{new}}{Constructor method with arguments:
#'    \describe{
#'     \item{\code{payoff}}{\code{payoff} object of the GMAB guarantee}
#'     \item{\code{t0}}{\code{\link{timeDate}} object with
#'     the issue date of the contract}
#'     \item{\code{t}}{\code{timeDate} object with the end date of the
#'     accumulation period}
#'     \item{\code{t1}}{\code{timeDate} object with the end date of the
#'      life benefit payment}
#'     \item{\code{age}}{\code{numeric} positive scalar with the age
#'      of the policyholder}
#'     \item{\code{fee}}{\code{\link{constant_parameters}} object with
#'      the fee}
#'     \item{\code{barrier}}{\code{numeric} positive scalar with the
#'      state-dependent fee barrier}
#'     \item{\code{penalty}}{\code{\link{penalty_class}} object with the
#'     penalty}
#'    }
#'   }
#'   \item{\code{get_times}}{get method for the product time-line.
#'    Returns a \code{\link{timeDate}} object}
#'   \item{\code{get_age}}{get method for the age of the insured}
#'   \item{\code{set_age}}{set method for the age of the insured}
#'   \item{\code{get_barrier}}{get method for the state-dependent fee barrier.
#'    Returns a positive scalar with the barrier}
#'   \item{\code{set_barrier}}{set method for the state-dependent fee barrier.
#'     Argument must be a positive scalar.}
#'   \item{\code{set_penalty_object}}{the argument \code{penalty} is a
#'   \code{\link{penalty_class}} object which is stored in a private field.}
#'   \item{\code{get_penalty_object}}{gets the \code{\link{penalty_class}} object.}
#'   \item{\code{set_penalty}}{set method for the penalty applied in case of
#'    surrender. The argument must be a scalar between 0 and 1.}
#'   \item{\code{get_penalty}}{get method for the surrender penalties. It can be
#'   a scalar between 0 and 1 in case the penalty is constant or a numeric vector
#'   in case the penalty varies with time.}
#'   \item{\code{set_fee}}{set method for the contract fee. The argument is
#'      a \code{\link{constant_parameters}} object with the fee.}
#'   \item{\code{set_payoff}}{set method for the \code{\link{payoff_guarantee}}
#'    object.}
#'   \item{\code{survival_benefit_times}}{returns a \code{numeric} vector with
#'    the survival benefit time indexes.}
#'   \item{\code{surrender_times}}{returns a \code{numeric} vector with the
#'    surrender time indexes. Takes as argument a string with the frequency
#'    of the decision if surrendering the contract,  e.g. "3m"
#'    corresponds to a surrender decision taken every 3 months.}
##'  \item{\code{times_in_yrs}}{returns the product time-line in
#'    fraction of year}
#'   \item{\code{cash_flows}}{returns a \code{numeric} vector with the
#'    cash flows of the product. It takes as argument \code{spot_values} a
#'    \code{numeric} vector which holds the values of the underlying fund this
#'    method will calculate the cash flows from}
#'   \item{\code{survival_benefit}}{Returns a numeric scalar corresponding to
#'    the survival benefit.
#'    The arguments are \code{spot_values} vector which holds the values of
#'    the underlying fund and \code{t} the time index of the survival benefit.
#'    The function will return 0 if there's no survival benefit at the
#'    specified time}
#'   \item{\code{get_premium}}{Returns the premium as non negative scalar}
#' }
#' @references
#' \enumerate{
#' \item{[BMOP2011]}{ \cite{Bacinello A.R., Millossovich P., Olivieri A.,
#'  Pitacco  E., "Variable annuities: a unifying valuation approach."
#' In: Insurance: Mathematics and Economics 49 (2011), pp. 285-297.
#' }}
#' \item{[BHM2014]}{ \cite{Bernard C., Hardy M. and Mackay A. "State-dependent
#' fees for variable annuity guarantees." In: Astin Bulletin 44 (2014),
#' pp. 559-585.}}
#' }

va_product <- R6::R6Class("va_product",
 public = list(
  initialize = function(payoff, t0, t, t1, age, fee, barrier, penalty, ...){
   if (!missing(payoff))
    if (inherits(payoff, "payoff_guarantee")) private$the_payoff <- payoff
    else stop(error_msg_1("payoff_guarantee"))
   else stop("Please provide a guarantee payoff object\n")

   #Initializes the issue (start) date t0
   if (!missing(t0))
    if (is_date(t0))
      private$t0 <- t0
    else stop(error_msg_1_("t0", "timeDate"))
   else stop(error_msg_1_("t0", "timeDate"))

   #Initializes the end of accumulation date t
   #Guarantees it's a date greater than t0 if not missing
   if (!missing(t))
    if (is_date(t))
     if (isTRUE(t0 < t)){}
     else stop(error_msg_11("t", "t0"))
    else  stop(error_msg_1_("t", "timeDate"))


   #Initializes the end of guaranteed benefit date t1
   #Guarantees it's a date greater than t0 if not missing
   if (!missing(t1))
    if(is_date(t1))
      if (isTRUE(t0 < t1)){}
      else stop(error_msg_11("t1", "t0"))
    else  stop(error_msg_1_("t1", "timeDate"))

   if(missing(t) & missing(t1))
    stop(error_msg_1_("t", "timeDate"))

   if (!missing(t) & !missing(t1))
    if(isTRUE(t <= t1)){
      private$t <- t
      private$t1 <- t1
    } else stop(error_msg_11("t1", "t"))

   if (missing(t) & !missing(t1)) {
    private$t1 <- t1
    private$t <- t1
   }

   if (!missing(t) & missing(t1)){
    private$t <- t
    private$t1 <- t
   }



   private$times <- timeDate::timeSequence(private$t0, private$t1)
   # Normalizes the product time line into year fractions
   private$times_yrs <- yr_fractions(private$times)


   if(!missing(age))
    if (is_positive_integer(age))
     private$the_age <- age
    else stop(error_msg_4("age"))
   else private$the_age <- 60
   if (!missing(fee))
    if(inherits(fee, "constant_parameters")) private$the_fee <- fee
    else stop(error_msg_1_("fee", "constant_parameters"))
   else private$the_fee <- constant_parameters$new(0.02, 365)

   if (!missing(barrier))
    if (is_positive_scalar(barrier)) private$the_barrier <- barrier
    else stop(error_msg_3_("barrier"))
   else private$the_barrier <- Inf

   if (!missing(penalty))
    if (inherits(penalty, "penalty_class"))
      private$the_penalty <- penalty
    else stop(error_msg_1_("penalty", "penalty_class"))
   else private$the_penalty <- penalty_class$new(type = 1, const = 1)

   private$surv_times <- length(private$times)

   if (identical(private$the_penalty$get_type(), 1))
     private$penalty <- private$the_penalty$get()
   else for (i in seq_along(private$times_yrs))
     private$penalty[i] <- private$the_penalty$get(private$times_yrs[i])

  },
  set_payoff = function(payoff){
    if (!missing(payoff))
      if (inherits(payoff, "payoff_guarantee")) private$the_payoff <- payoff
      else stop(error_msg_1("payoff_guarantee"))
    else stop("Please provide a guarantee payoff object\n")
  },
  get_age = function() private$the_age,
  set_age = function(age){
   if(!missing(age))
    if (is_positive_integer(age))
     private$the_age <- age
    else stop(error_msg_4("age"))
   else private$the_age <- 60
  },
  get_barrier = function() private$the_barrier,
  set_barrier = function(barrier) {
   if (!missing(barrier))
    if (is_positive_scalar(barrier)) private$the_barrier <- barrier
    else stop(error_msg_3_("barrier"))
   else private$the_barrier <- Inf
  },
  set_penalty_object = function(penalty){
   if (!missing(penalty))
    if (inherits(penalty, "penalty_class"))
     private$the_penalty <- penalty
    else stop(error_msg_1_("penalty", "penalty_class"))
   else private$the_penalty <- penalty_class$new(type = 1, const = 1)

   if (identical(private$the_penalty$get_type(), 1))
     private$penalty <- private$the_penalty$get()
   else for (i in seq_along(private$times_yrs))
     private$penalty[i] <- private$the_penalty$get(private$times_yrs[i])
  },
  get_penalty_object = function(time) private$the_penalty,
  set_penalty = function(penalty) {
    private$the_penalty$set(penalty)
    if (identical(private$the_penalty$get_type(), 1))
      private$penalty <- private$the_penalty$get()
    else for (i in seq_along(private$times_yrs))
      private$penalty[i] <- private$the_penalty$get(private$times_yrs[i])
  },
  get_penalty = function( ) private$penalty,
  set_fee = function(fee){
    if (!missing(fee))
      if(inherits(fee, "constant_parameters")) private$the_fee <- fee
      else stop(error_msg_1_("fee", "constant_parameters"))
      else private$the_fee <- constant_parameters$new(0.02, 365)
  },
  get_premium = function() private$the_payoff$get_premium(),
  get_times = function() private$times,
  times_in_yrs = function() private$times_yrs,
  survival_benefit_times = function(){},
  surrender_times = function(){},
  cash_flows = function(spot_values, ...) spot_values
 ),
 private = list(
  #Issue date of the contract
  t0 = "timeDate",
  #End of the accumulation period date
  t = "timeDate",
  #End of guaranteed benefit date
  t1 = "timeDate",
  #timeDate object to store the product time-line
  times = "timeDate",
  #payoff_guarantee object which stores the type of payoff
  the_payoff = "payoff_guarantee",
  #A posite scalar with the age of the insured
  the_age = "numeric",
  #A positive scalar with the annual VA contract fee.
  the_fee = "constant_parameters",
  #A positive scalar with the barrier for state-dependent fees.
  the_barrier = "numeric",
  #A scalar or numeric vector  with the withdrawal penalty.
  penalty = numeric(0),
  #A penalty object
  the_penalty = "penalty_class",
  #A numeric vector with the product time-line
  #in fraction of years
  times_yrs = "numeric",
  #Survival benefit time index
  surv_times = "numeric"
 )
)


#' Variable Annuity with GMAB guarantee
#' @description
#' Class for VA with Guaranteed Minimum Accumulation Benefit (GMAB).
#' It supports a simple state-dependent fee structure with a single barrier.\cr
#' See \bold{References} for a description of variable annuities life
#' insurance products, their guarantees and fee structures.
#' @docType class
#' @export
#' @return Object of \code{\link{R6Class}}
#' @format \code{\link{R6Class}} object.
#' @section Methods:
#'  \describe{
#'   \item{\code{new}}{Constructor method with arguments:
#'    \describe{
#'     \item{\code{payoff}}{\code{payoff} object of the GMAB guarantee}
#'     \item{\code{t0}}{\code{\link{timeDate}} object with
#'     the issue date of the contract}
#'     \item{\code{t}}{\code{timeDate} object with the end date of the
#'     accumulation period}
#'     \item{\code{t1}}{\code{timeDate} object with the end date of the
#'      life benefit payment}
#'     \item{\code{age}}{\code{numeric} positive scalar with the age
#'      of the policyholder}
#'     \item{\code{fee}}{\code{\link{constant_parameters}} object with
#'      the fee}
#'     \item{\code{barrier}}{\code{numeric} positive scalar with the
#'      state-dependent fee barrier}
#'     \item{\code{penalty}}{\code{\link{penalty_class}} object with the
#'      penalty}
#'    }
#'   }
#'   \item{\code{get_times}}{get method for the product time-line.
#'    Returns a \code{\link{timeDate}} object}
#'   \item{\code{get_age}}{get method for the age of the insured}
#'   \item{\code{set_age}}{set method for the age of the insured}
#'   \item{\code{get_barrier}}{get method for the state-dependent fee barrier.
#'    Returns a positive scalar with the barrier}
#'   \item{\code{set_barrier}}{set method for the state-dependent fee barrier.
#'    Argument must be a positive scalar.}
#'   \item{\code{set_penalty_object}}{the argument \code{penalty} is a
#'   \code{\link{penalty_class}} object which is stored in a private field.}
#'   \item{\code{get_penalty_object}}{gets the \code{\link{penalty_class}} object.}
#'   \item{\code{set_penalty}}{set method for the penalty applied in case of
#'    surrender. The argument must be a scalar between 0 and 1.}
#'   \item{\code{get_penalty}}{get method for the surrender penalties. It can be
#'   a scalar between 0 and 1 in case the penalty is constant or a numeric vector
#'   in case the penalty varies with time.}
#'   \item{\code{set_fee}}{set method for the contract fee. The argument is
#'      a \code{\link{constant_parameters}} object with the fee.}
#'   \item{\code{set_payoff}}{set method for the \code{\link{payoff_guarantee}}
#'    object.}
#'   \item{\code{survival_benefit_times}}{returns a \code{numeric} vector with
#'    the survival benefit time indexes.}
#'   \item{\code{surrender_times}}{returns a \code{numeric} vector with the
#'    surrender time indexes. Takes as argument a string with the frequency
#'    of the decision if surrendering the contract,  e.g. "3m"
#'    corresponds to a surrender decision taken every 3 months.}
#'   \item{\code{times_in_yrs}}{returns the product time-line in
#'    fraction of year}
#'   \item{\code{cash_flows}}{returns a \code{numeric} vector with the
#'    cash flows of the product. It takes as argument \code{spot_values} a
#'    \code{numeric} vector which holds the values of the underlying fund and
#'     \code{death_time} a time index with the time of death}
#'   \item{\code{survival_benefit}}{Returns a numeric scalar corresponding to
#'    the survival benefit.
#'    The arguments are \code{spot_values} vector which holds the values of
#'    the underlying fund and \code{t} the time index of the survival benefit.}
#'   \item{\code{get_premium}}{Returns the premium as non negative scalar}
#' }
#' @references
#' \enumerate{
#' \item{[BMOP2011]}{ \cite{Bacinello A.R., Millossovich P., Olivieri A.,
  #'  Pitacco  E., "Variable annuities: a unifying valuation approach."
  #' In: Insurance: Mathematics and Economics 49 (2011), pp. 285-297.
  #' }}
#' \item{[BHM2014]}{ \cite{Bernard C., Hardy M. and Mackay A. "State-dependent
  #' fees for variable annuity guarantees." In: Astin Bulletin 44 (2014),
  #' pp. 559-585.}}
#' }
#'@examples
#'#Sets up the payoff as a roll-up of premiums with roll-up rate 1%
#'
#'rate <- constant_parameters$new(0.01)
#'
#'premium <- 100
#'rollup <- payoff_rollup$new(premium, rate)
#'
#'#Five years time-line
#'begin <- timeDate::timeDate("2016-01-01")
#'end <- timeDate::timeDate("2020-12-31")
#'
#'age <- 60
#'# A constant fee of 2% per year (365 days)
#'fee <- constant_parameters$new(0.02)
#'
#'#Barrier for a state-dependent fee. The fee will be applied only if
#'#the value of the account is below the barrier
#'barrier <- 200
#'
#'#Withdrawal penalty applied in case the insured surrenders the contract
#'#It is a constant penalty in this case
#'penalty <- penalty_class$new(type = 1, 0.01)
#'
#'
#'#Sets up a VA contract with GMAB guarantee. The guaranteed miminum
#'#is the roll-up of premiums with rate 1%
#'contract <- GMAB$new(rollup, t0 = begin, t = end, age = age,  fee = fee,
#'barrier = barrier, penalty = penalty)



GMAB <- R6::R6Class("GMAB", inherit = va_product,
 public = list(
  survival_benefit_times = function() private$surv_times,
  surrender_times = function(freq){
    #Check on freq units
    units <- c("m", "w", "d")
    freq_unit = gsub("[ 0-9]", "", freq, perl = TRUE)
    if (!(freq_unit %in% units)) stop(error_msg_10())
    #
    surr_dates <- timeDate::periods(private$times, freq, freq)$to
    surr_idx <- vector(mode = "numeric", length = length(surr_dates))
    for (i in seq_along(surr_dates))
      surr_idx[i] <- which(surr_dates[i] == private$times)
    head(surr_idx, -1)
  },
  cash_flows = function(spot_values, death_time, ...){
   fee <- private$the_fee$get()
   barrier <- private$the_barrier
   penalty <- private$penalty
   len <- length(spot_values)

   if (death_time <= length(private$times)){
    ben <- rep(0, death_time)
    out <- calc_account(spot_values[1:death_time], ben, fee, barrier, penalty)
    if(death_time < length(private$times)){
      out <- rep(out, length.out=len)
      out[(death_time+1):len] <- 0
    }
   } else {
    time_int <- c(private$t0, private$t)
    ben <- rep(0, len)
    out <- calc_account(spot_values, ben, fee, barrier, penalty)
    #GMAB living benefit
    last <- length(out)
    out[last] <- private$the_payoff$get_payoff(out[last], time_int, out)
   }
   out
  },
  survival_benefit = function(spot_values, death_time, t){
    last <- length(private$times)
    penalty <- private$penalty
    if (t == last & t != death_time){
      fee <- private$the_fee$get()
      barrier <- private$the_barrier
      t0 <- head(private$times, 1)
      t1 <- tail(private$times, 1)
      ben <- rep(0, last)
      out <- calc_account(spot_values, ben, fee, barrier, penalty)
      out <- private$the_payoff$get_payoff(out[last], c(t0, t1), out)
    } else out <- 0
    out
  }
 )
)


#' Variable Annuity with GMAB and GMDB guarantees
#' @description
#' Class for a VA with Guaranteed Minimum Accumulation Benefit (GMAB)
#' and Guaranteed Minimum Accumulation Benefit (GMDB).
#' It supports a simple state-dependent fee structure with a single barrier.\cr
#' See \bold{References} for a description of variable annuities life
#' insurance products, their guarantees and fee structures.
#' @docType class
#' @export
#' @return Object of \code{\link{R6Class}}
#' @format \code{\link{R6Class}} object.
#' @section Methods:
#'  \describe{
#'   \item{\code{new}}{Constructor method with arguments:
#'    \describe{
#'     \item{\code{payoff}}{\code{payoff} object of the GMAB guarantee}
#'     \item{\code{t0}}{\code{\link{timeDate}} object with
#'     the issue date of the contract}
#'     \item{\code{t}}{\code{timeDate} object with the end date of the
#'     accumulation period}
#'     \item{\code{t1}}{\code{timeDate} object with the end date of the
#'      life benefit payment}
#'     \item{\code{age}}{\code{numeric} positive scalar with the age
#'     of the policyholder}
#'     \item{\code{fee}}{\code{\link{constant_parameters}} object with
#'     the fee}
#'     \item{\code{barrier}}{\code{numeric} positive scalar with the
#'     state-dependent fee barrier}
#'     \item{\code{penalty}}{\code{\link{penalty_class}} object with the
#'      penalty}
#'     \item{\code{death_payoff}}{\code{payoff} object with the payoff
#'     of the GMDB guarantee}
#'    }
#'   }
#'   \item{\code{get_times}}{get method for the product time-line.
#'    Returns a \code{\link{timeDate}} object}
#'   \item{\code{get_age}}{get method for the age of the insured}
#'   \item{\code{set_age}}{set method for the age of the insured}
#'   \item{\code{get_barrier}}{get method for the state-dependent fee barrier.
#'    Returns a positive scalar with the barrier}
#'   \item{\code{set_barrier}}{set method for the state-dependent fee barrier.
#'     Argument must be a positive scalar.}
#'   \item{\code{set_penalty_object}}{the argument \code{penalty} is a
#'   \code{\link{penalty_class}} object which is stored in a private field.}
#'   \item{\code{get_penalty_object}}{gets the \code{\link{penalty_class}} object.}
#'   \item{\code{set_penalty}}{set method for the penalty applied in case of
#'    surrender. The argument must be a scalar between 0 and 1.}
#'   \item{\code{get_penalty}}{get method for the surrender penalties. It can be
#'   a scalar between 0 and 1 in case the penalty is constant or a numeric vector
#'   in case the penalty varies with time.}
#'   \item{\code{set_fee}}{set method for the contract fee. The argument is
#'      a \code{\link{constant_parameters}} object with the fee.}
#'   \item{\code{set_payoff}}{set method for the \code{\link{payoff_guarantee}}
#'   object of the GMAB rider}
#'   \item{\code{set_death_payoff}}{set method for the
#'   \code{\link{payoff_guarantee}}  object of the GMDB rider}
#'   \item{\code{survival_benefit_times}}{returns a \code{numeric} vector with
#'    the survival benefit time indexes.}
#'   \item{\code{surrender_times}}{returns a \code{numeric} vector with the
#'    surrender time indexes. Takes as argument a string with the frequency
#'    of the decision if surrendering the contract,  e.g. "3m"
#'    corresponds to a surrender decision taken every 3 months.}
#'   \item{\code{times_in_yrs}}{returns the product time-line in
#'    fraction of year}
#'   \item{\code{cash_flows}}{returns a \code{numeric} vector with the
#'    cash flows of the product. It takes as argument \code{spot_values} a
#'    \code{numeric} vector which holds the values of the underlying fund and
#'     \code{death_time} a time index with the time of death}
#'   \item{\code{survival_benefit}}{Returns a numeric scalar corresponding to
#'    the survival benefit.
#'    The arguments are \code{spot_values} vector which holds the values of
#'    the underlying fund and \code{t} the time index of the survival benefit.}
#'   \item{\code{get_premium}}{Returns the premium as non negative scalar}
#' }
#' @references
#' \enumerate{
#' \item{[BMOP2011]}{ \cite{Bacinello A.R., Millossovich P., Olivieri A.,
#' Pitacco  E., "Variable annuities: a unifying valuation approach."
#' In: Insurance: Mathematics and Economics 49 (2011), pp. 285-297.
#' }}
#' \item{[BHM2014]}{ \cite{Bernard C., Hardy M. and Mackay A. "State-dependent
#' fees for variable annuity guarantees." In: Astin Bulletin 44 (2014),
#' pp. 559-585.}}
#' }
#'@examples
#'#Sets up the payoff as a roll-up of premiums with roll-up rate 1%
#'
#'rate <- constant_parameters$new(0.01)
#'
#'premium <- 100
#'rollup <- payoff_rollup$new(premium, rate)
#'
#'#Five years time-line
#'begin <- timeDate::timeDate("2016-01-01")
#'end <- timeDate::timeDate("2020-12-31")
#'#Age of the insured
#'age <- 60
#'# A constant fee of 2% per year (365 days)
#'fee <- constant_parameters$new(0.02, 365)
#'
#'#Barrier for a state-dependent fee. The fee will be applied only if
#'#the value of the account is below the barrier
#'barrier <- 200
#'
#'#Withdrawal penalty applied in case the insured surrenders the contract
#'#It is a constant penalty in this case
#'penalty <- penalty_class$new(type = 1, 0.01)
#'
#'#Sets up the GMAB + GMDB with the same payoff for survival and death
#'#benefits
#'contract <- GMAB_GMDB$new(rollup, t0 = begin, t = end, age = age, fee =fee,
#'barrier = barrier, penalty = penalty, death_payoff = rollup)

GMAB_GMDB <- R6::R6Class("GMAB_GMDB", inherit = GMAB,
 public = list(
  initialize = function(payoff, t0, t, t1, age, fee, barrier, penalty,
                        death_payoff){
   super$initialize(payoff, t0, t, t1, age, fee, barrier, penalty)
   if (!missing(death_payoff))
     if (inherits(death_payoff, "payoff_guarantee"))
     private$the_death_payoff <- death_payoff
     else stop(error_msg_1_("death_payoff", "payoff_guarantee"))
   else stop("Please provide a guarantee payoff object\n")
   },
  set_death_payoff = function(death_payoff){
    if (!missing(death_payoff))
      if (inherits(death_payoff, "payoff_guarantee"))
        private$the_death_payoff <- death_payoff
      else stop(error_msg_1_("death_payoff", "payoff_guarantee"))
    else stop("Please provide a guarantee payoff object\n")
  },
  cash_flows = function(spot_values, death_time, ...){
   fee <- private$the_fee$get()
   barrier <- private$the_barrier
   penalty <- private$penalty
   len <- length(spot_values)
   t0 <- private$t0
   if (death_time <= length(private$times)){
    ben <- rep(0, death_time)
    out <- calc_account(spot_values[1:death_time], ben, fee, barrier, penalty)
    #GMDB death benefit
    last <- length(out)
    t <- private$times[death_time]
    out[last] <- private$the_death_payoff$get_payoff(out[last], c(t0, t), out)
    if(death_time < length(private$times)){
      out <- rep(out, length.out=len)
      out[(death_time+1):len] <- 0
    }
   } else {
    t <- private$t
    ben <- rep(0, len)
    out <- calc_account(spot_values, ben, fee, barrier, penalty)
    #GMAB living benefit
    last <- length(out)
    out[last] <- private$the_payoff$get_payoff(out[last], c(t0, t), out)
   }
   out
  }
 ),
 private = list(
  #A payoff_guarantee object which stores the type
  #of payoff (e.g: roll-up, ratchet, etc) for the death benefit.
  the_death_payoff = "payoff_guarantee"
 )
)


#' Variable Annuity with GMDB guarantee
#' @description
#' Class for VA with Guaranteed Minimum Death Benefit (GMDB).
#' It supports a simple state-dependent fee structure with a single barrier.\cr
#' See \bold{References} for a description of variable annuities life
#' insurance products, their guarantees and fee structures.
#' @docType class
#' @export
#' @return Object of \code{\link{R6Class}}
#' @format \code{\link{R6Class}} object.
#' @section Methods:
#'  \describe{
#'   \item{\code{new}}{Constructor method with arguments:
#'   \describe{
#'     \item{\code{payoff}}{\code{payoff} object of the GMDB guarantee}
#'     \item{\code{t0}}{\code{\link{timeDate}} object with
#'     the issue date of the contract}
#'     \item{\code{t}}{\code{timeDate} object with the end date of the
#'     accumulation period}
#'     \item{\code{t1}}{\code{timeDate} object with the end date of the
#'      life benefit payment}
#'     \item{\code{age}}{\code{numeric} positive scalar with the age
#'     of the policyholder}
#'     \item{\code{fee}}{\code{\link{constant_parameters}} object with
#'     the fee}
#'     \item{\code{barrier}}{\code{numeric} positive scalar with the
#'     state-dependent fee barrier}
#'     \item{\code{penalty}}{\code{\link{penalty_class}} object with the
#'      penalty}
#'    }
#'   }
#'   \item{\code{get_times}}{get method for the product time-line.
#'    Returns a \code{\link{timeDate}} object}
#'   \item{\code{get_age}}{get method for the age of the insured}
#'   \item{\code{set_age}}{set method for the age of the insured}
#'   \item{\code{get_barrier}}{get method for the state-dependent fee barrier.
#'    Returns a positive scalar with the barrier}
#'   \item{\code{set_barrier}}{set method for the state-dependent fee barrier.
#'    Argument must be a positive scalar.}
#'   \item{\code{set_penalty_object}}{the argument \code{penalty} is a
#'   \code{\link{penalty_class}} object which is stored in a private field.}
#'   \item{\code{get_penalty_object}}{gets the \code{\link{penalty_class}} object.}
#'   \item{\code{set_penalty}}{set method for the penalty applied in case of
#'    surrender. The argument must be a scalar between 0 and 1.}
#'   \item{\code{get_penalty}}{get method for the surrender penalties. It can be
#'   a scalar between 0 and 1 in case the penalty is constant or a numeric vector
#'   in case the penalty varies with time.}
#'   \item{\code{set_fee}}{set method for the contract fee. The argument is
#'      a \code{\link{constant_parameters}} object with the fee.}
#'   \item{\code{survival_benefit_times}}{returns a \code{numeric} vector with
#'    the survival benefit time indexes.}
#'   \item{\code{surrender_times}}{returns a \code{numeric} vector with the
#'    surrender time indexes. Takes as argument a string with the frequency
#'    of the decision if surrendering the contract,  e.g. "3m"
#'    corresponds to a surrender decision taken every 3 months.}
#'   \item{\code{times_in_yrs}}{returns the product time-line in
#'    fraction of year}
#'   \item{\code{cash_flows}}{returns a \code{numeric} vector with the
#'    cash flows of the product. It takes as argument \code{spot_values} a
#'    \code{numeric} vector which holds the values of the underlying fund and
#'     \code{death_time} a time index with the time of death}
#'   \item{\code{survival_benefit}}{Returns a numeric scalar corresponding to
#'    the survival benefit.
#'    The arguments are \code{spot_values} vector which holds the values of
#'    the underlying fund and \code{t} the time index of the survival benefit.}
#'   \item{\code{get_premium}}{Returns the premium as non negative scalar}
#' }
#' @references
#' \enumerate{
#' \item{[BMOP2011]}{ \cite{Bacinello A.R., Millossovich P., Olivieri A.,
  #'  Pitacco  E., "Variable annuities: a unifying valuation approach."
  #' In: Insurance: Mathematics and Economics 49 (2011), pp. 285-297.
  #' }}
#' \item{[BHM2014]}{ \cite{Bernard C., Hardy M. and Mackay A. "State-dependent
  #' fees for variable annuity guarantees." In: Astin Bulletin 44 (2014),
  #' pp. 559-585.}}
#' }
#'@examples
#'#Sets up the payoff as a roll-up of premiums with roll-up rate 2%
#'
#'rate <- constant_parameters$new(0.02)
#'
#'premium <- 100
#'rollup <- payoff_rollup$new(premium, rate)
#'
#'begin <- timeDate::timeDate("2016-01-01")
#'end <- timeDate::timeDate("2020-12-31")
#'
#'age <- 60
#'# A constant fee of 0.02% per year (365 days)
#'fee <- constant_parameters$new(0.02)
#'
#'#Barrier for a state-dependent fee. The fee will be applied only if
#'#the value of the account is below the barrier
#'barrier <- Inf
#'
#'#Withdrawal penalty applied in case the insured surrenders the contract
#'#It is a constant penalty in this case
#'penalty <- penalty_class$new(type = 1, 0.01)
#'
#'#Sets up a VA contract with GMDB guarantee. The guaranteed miminum
#'#is the roll-up of premiums with rate 2%
#'
#'contract <- GMDB$new(rollup, t0 = begin, t = end, age = age,  fee = fee,
#'barrier = barrier, penalty = penalty)



GMDB <- R6::R6Class("GMDB", inherit = GMAB,
 public = list(
  cash_flows = function(spot_values, death_time, ...){
   fee <- private$the_fee$get()
   barrier <- private$the_barrier
   penalty <- private$penalty
   len <- length(spot_values)
   t0 <- private$t0
   if (death_time <= length(private$times)){
    ben <- rep(0, death_time)
    out <- calc_account(spot_values[1:death_time], ben, fee, barrier, penalty)
    #GMDB death benefit
    last <- length(out)
    t <- private$times[death_time]
    out[last] <- private$the_payoff$get_payoff(out[last], c(t0, t), out)
    if(death_time < length(private$times)){
      out <- rep(out, length.out=len)
      out[(death_time+1):len] <- 0
    }
   } else {
    ben <- rep(0, len)
    out <- calc_account(spot_values, ben, fee, barrier, penalty)
   }
   out
  },
  survival_benefit = function(spot_values, death_time, t){
    last <- private$surv_times
    penalty <- private$penalty
    if (t == last & t != death_time){
      fee <- private$the_fee$get()
      barrier <- private$the_barrier
      ben <- rep(0, last)
      out <- calc_account(spot_values, ben, fee, barrier, penalty)
      out <- out[last]
    } else out <- 0
    out
  }
 )
)

