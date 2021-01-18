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


#' Variable Annuity with GMIB guarantee
#' @description
#' Class for VA with Guaranteed Minimum Income Benefit (GMIB).
#' A GMIB rider provides a lifetime annuity from a specified future time.
#' Types of GMIB supported are a whole-life annuity (Ia),
#' an annuity-certain (Ib) or annuity-certain followed by a deferred
#' whole-life annuity (Ic).
#' It supports a simple state-dependent fee structure with a single barrier.\cr
#' See \bold{References} for a description of variable annuities life
#' insurance products, their guarantees and fee structures.
#' @details
#' The annuity payment is assumed to be annual and it's calculated as
#' the annuitization rate by the roll-up or ratchet payoff at the end
#' of the accumulation period \code{t}.
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
#'    \item{\code{penalty}}{\code{\link{penalty_class}} object with the
#'      penalty}
#'     \item{\code{eta}}{\code{numeric} scalar with the market
#'      annuitisation rate}
#'     \item{\code{type}}{\code{string} with the income benefit type:
#'      it can be 'Ia' for a whole-life annuity, 'Ib' for an
#'      annuity-certain with maturity t1, 'Ic' for an annuity certain
#'      with maturity t1 followed with a deferred life-annuity if
#'      the insured is alive after t1.}
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
#'    \code{numeric} vector which holds the values of the underlying
#'     fund, \code{death_time} a time index with the time of death and
#'    \code{discounts} a \code{numeric} vector with the discount factors
#'     at time of death. These latest are used to calculate the death
#'     benefit for type Ib and Ic.}
#'   \item{\code{survival_benefit}}{Returns a numeric scalar corresponding to
#'    the survival benefit.
#'    The arguments are \code{spot_values} vector which holds the values of
#'    the underlying fund and \code{time} the time index of the survival benefit.
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
#'@examples
#'#Sets up the payoff as a roll-up of premiums with roll-up rate 1%
#'
#'rate <- constant_parameters$new(0.01)
#'
#'premium <- 100
#'rollup <- payoff_rollup$new(premium, rate)
#'
#'
#'t0 <- timeDate::timeDate("2016-01-01")
#'
#'#Five year accumulation period
#'t <- timeDate::timeDate("2020-12-31")
#'
#'#Five year annuity certain period
#'t1 <- timeDate::timeDate("2025-12-31")
#'
#'age <- 60
#'
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
#'#Sets up a VA contract with GMIB guarantee, whole-life (Ia).
#'contract <- GMIB$new(rollup, t0 = t0, t = t, age = age,  fee = fee,
#'barrier = barrier, penalty = penalty, eta = 0.04)
#'
#'#Sets up a VA contract with GMIB gurantee annuity-certain with
#'#maturity t1
#'contract <- GMIB$new(rollup, t0 = t0, t = t, t1 = t1,  age = age,
#'fee = fee, barrier = barrier, penalty = penalty, eta = 0.04, type = "Ib")


GMIB <- R6::R6Class("GMIB", inherit = va_product,
 public = list(
  initialize = function(payoff, t0, t, t1, age, fee, barrier, penalty,
                        eta = 0.05, type){

   super$initialize(payoff, t0, t, t1, age, fee, barrier, penalty)

   if(!missing(eta))
    if(is_between(eta, 0, 1))
     private$eta <- eta
    else stop(error_msg_8("eta"))

   if(!missing(type))
    if (is_identical_to_any(type, c("Ia", "Ib", "Ic")))
     private$type <- type
    else stop("type must be either \"Ia\" or \"Ib\" or \"Ic\"")
   else private$type <- "Ia"

   if(is_identical_to_any(private$type, c("Ia", "Ic"))){
    max_age <- 120
    private$times <- timeDate::timeSequence(private$t0,
                      length.out = 365 * (max_age - private$the_age))
    # Normalizes the product time line into year fractions
    private$times_yrs <- yr_fractions(private$times)
   }

   if(is_identical_to_any(private$type, "Ib") & missing(t1))
     stop(error_msg_1_("t1", "timeDate"))

   ts <- private$times
   s_id <- which(private$t == ts)
   endpoints <- timeDate::.endpoints(ts, on = "y")
   out <- unique(sort(c(s_id, endpoints)))
   private$surv_times <- out[out >= s_id]
  },
  survival_benefit_times = function() private$surv_times,
  surrender_times = function(freq = "3m"){
    #Check on freq units
    units <- c("m", "w", "d")
    freq_unit = gsub("[ 0-9]", "", freq, perl = TRUE)
    if (!(freq_unit %in% units)) stop(error_msg_10())
    #surrender is not admitted during the payout phase
    #of a GMIB rider
    tseq <- timeDate::timeSequence(from = private$t0,
                         to = private$t)
    surr_dates <- timeDate::periods(tseq, freq, freq)$to
    surr_idx <- vector(mode = "numeric", length = length(surr_dates))
    for (i in seq_along(surr_dates))
      surr_idx[i] <- which(surr_dates[i] == private$times)
    head(surr_idx, -1)
  },
  cash_flows = function(spot_values, death_time, discounts){
    #discounts is a vector with the discount factors
    #from death_time till t1. It must be generated and passed
    #by the engine.
    fee <- private$the_fee$get()
    barrier <- private$the_barrier
    penalty <- private$penalty
    len <- length(spot_values)
    t_id <- which(private$t == private$times)
    if (death_time <= t_id){
      ben <- rep(0, death_time)
      out <- calc_account(spot_values[1:death_time], ben, fee, barrier, penalty)
      out <- rep(out, length.out=len)
      out[(death_time + 1) : len] <- 0
    } else {
      time_int <- c(private$t0, private$t)
      ben <- rep(0, t_id)
      out <- calc_account(spot_values[1:t_id], ben, fee, barrier, penalty)
      last <- length(out)
      income <- private$the_payoff$get_payoff(out[last], time_int, out)
      out <- rep(out, length.out = len)
      out[t_id:len] <- 0
      out[self$survival_benefit_times()] <- private$eta * income
      switch(private$type,
       "Ia" = {#Whole life annuity
         if(death_time < len ) out[(death_time +1) : len] <- 0
        },
       "Ib" = {#Annuity certain
         t1_id <- which(private$t1 == private$times)
         if(death_time < len){
           out[death_time] <- sum(discounts[death_time : len] *
                                    out[death_time : len])
           out[(death_time + 1) : len] <- 0
         }
        },
       "Ic" = {#Annuity certain + whole life annuity
         t1_id <- which(private$t1 == private$times)
         if(death_time < t1_id){
          out[death_time] <- sum(discounts[death_time : t1_id] *
                                   out[death_time : t1_id])
          out[(death_time + 1) : len] <- 0
         } else out[(death_time + 1) : len] <- 0
       }
      )
    }
    out
  },
  survival_benefit = function(spot_values, death_time, time){
    if(time %in% private$surv_times & time != death_time){
     fee <- private$the_fee$get()
     barrier <- private$the_barrier
     penalty <- private$penalty
     t_id <- which(private$t == private$times)
     time_int <- c(private$t0, private$t)
     ben <- rep(0, t_id)
     out <- calc_account(spot_values[1:t_id], ben, fee, barrier, penalty)
     last <- length(out)
     income <- private$the_payoff$get_payoff(out[last], time_int, out)
     out <- private$eta * income
    }else out <- 0
  out
  }
 ),
 private = list(
  #Market annuitisation rate
  eta = 0.05,
  #Type of the Income Benefit
  type = "Ia"
 )
)
