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


#' Variable Annuity with GMWB guarantee
#' @description
#' Class for a VA product with Guaranteed Minimum Withdrawal Benefit (GMWB).
#' A GMWB rider allows for periodic withdrawals from the policy account.
#' Types of GMWB  supported are withdrawals up to a fixed date independent
#' of survival (Wa), withdrawals up to fixed date only if the insured is
#' alive (Wb) or whole life withdrawals (Wc).
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
#'     \item{\code{payoff}}{\code{payoff_GMWB} object with the amount of
#'     the  periodic withdrawal}
#'     \item{\code{t0}}{\code{\link{timeDate}} object with
#'     the issue date of the contract}
#'     \item{\code{t1}}{\code{timeDate} object with the end date of the contract}
#'     \item{\code{age}}{\code{numeric} positive scalar with the age
#'      of the policyholder}
#'     \item{\code{fee}}{\code{\link{constant_parameters}} object with
#'      the fee}
#'     \item{\code{barrier}}{\code{numeric} positive scalar with the
#'      state-dependent fee barrier}
#'     \item{\code{penalty}}{\code{\link{penalty_class}} object with the
#'      penalty}
#'     \item{\code{type}}{\code{string} with the GMWB contract type:
#'      it can be \code{'Wa'} for withdrawals up to \code{t1} independent
#'       of survival,\code{'Wb'} for withdrawals up to \code{t1} only
#'       if the insured is alive, \code{'Wc'} for whole life withdrawals.}
#'     \item{\code{freq}}{\code{string} with the frequency of withdrawals
#'     expressed in months (e.g. \code{'12m'} stands for yearly withdrawals).}
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
#'    cash flows of the product. It takes as argument: \code{spot_values} a
#'    \code{numeric} vector which holds the values of the underlying
#'     fund, \code{death_time} a time index with the time of death and
#'    \code{discounts} a \code{numeric} vector with the discount factors
#'     at time of death. These latest are used to calculate the death
#'     benefit for the GMWB of type Wa}
#'   \item{\code{survival_benefit}}{Returns a numeric scalar corresponding to
#'    the survival benefit.
#'    The arguments are: \code{spot_values} vector which holds the values of
#'    the underlying fund, \code{death_time} time index of the time of death
#'    and \code{time} the time index of the survival benefit.
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
#'#Sets up the periodic payment.
#'
#'premium <- 100
#'beta <- 0.1
#'GMWB_payment <- payoff_GMWB$new(premium, beta)
#'
#'#Issue date of the contract
#'t0 <- timeDate::timeDate("2016-01-01")
#'
#'#Ten years expiration of the guarantee
#'
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
#'#Sets up a VA contract with GMWB guarantee type Wa with yearly
#'#withdrawals for 10 years.
#'
#'contract <- GMWB$new(GMWB_payment, t0 = t0, t1 = t1, age = age,  fee = fee,
#'barrier = barrier, penalty = penalty, type = "Wa", freq = "12m")
#'


GMWB <- R6::R6Class("GMWB", inherit = va_product,
 public = list(
  initialize = function(payoff, t0, t1, age, fee, barrier, penalty, type, freq){

   if (!missing(payoff))
    if (inherits(payoff, "payoff_GMWB")){
     super$initialize(payoff, t0 = t0, t1 = t1, age = age, fee = fee,
                      barrier = barrier, penalty = penalty)
    } else stop(error_msg_1_("payoff", "payoff_GMWB"))
   else stop(error_msg_1_("payoff", "payoff_GMWB"))

   if(!missing(type))
    if (is_identical_to_any(type, c("Wa", "Wb", "Wc")))
      private$type <- type
     else stop("type must be either \"Wa\" or \"Wb\" or \"Wc\"")
   else private$type <- "Wa"

   if(is_identical_to_any(private$type,  "Wc")){
     max_age <- 120
     private$times <- timeDate::timeSequence(private$t0,
                       length.out = 365 * (max_age - private$the_age))
     # Normalizes the product time line into year fractions
     private$times_yrs <- yr_fractions(private$times)
     private$t <- tail(private$times, 1)
   }
   #checks units and sets the frequency private field
   if (!missing(freq)){
    units <- c("m")
    freq_unit = gsub("[ 0-9]", "", freq, perl = TRUE)
    if (!(freq_unit %in% units)) stop("Allow unit for argument freq is 'm'")
    else private$freq <- freq
   } else private$freq <- "12m"

   #Initializes the withdrawals  dates

   surv_dates <- timeDate::periods(private$times, freq, freq)$to
   surv_dates[length(surv_dates)] <- private$times[length(private$times)]
   surv_idx <- vector(mode = "numeric", length = length(surv_dates))
   for (i in seq_along(surv_dates))
     surv_idx[i] <- which(surv_dates[i] == private$times)
   private$surv_times <- c(1, surv_idx)

  },
  get_freq = function() private$freq,
  set_freq = function(freq) {
    #checks units and sets the frequency private field
    units <- c("m")
    freq_unit = gsub("[ 0-9]", "", freq, perl = TRUE)
    if (!(freq_unit %in% units)) stop("Allow unit for argument freq is 'm'")
    else private$freq <- freq

    #Initializes the withdrawals  dates

    surv_dates <- timeDate::periods(private$times, freq, freq)$to
    surv_dates[length(surv_dates)] <- private$times[length(private$times)]
    surv_idx <- vector(mode = "numeric", length = length(surv_dates))
    for (i in seq_along(surv_dates))
      surv_idx[i] <- which(surv_dates[i] == private$times)
    private$surv_times <- c(1, surv_idx)

  },
  survival_benefit_times = function( ) private$surv_times,
  surrender_times = function(freq = "3m"){
   #Check on freq units
   units <- c("m", "w", "d")
   freq_unit = gsub("[ 0-9]", "", freq, perl = TRUE)
   if (!(freq_unit %in% units)) stop(error_msg_10())
   surr_dates <- timeDate::periods(private$times, freq, freq)$to
   surr_dates[length(surr_dates)] <- private$times[length(private$times)]
   surr_idx <- vector(mode = "numeric", length = length(surr_dates))
   for (i in seq_along(surr_dates))
    surr_idx[i] <- which(surr_dates[i] == private$times)
   head(surr_idx, -1)
  },
  cash_flows = function(spot_values, death_time, discounts){
   fee <- private$the_fee$get()
   barrier <- private$the_barrier
   penalty <- private$penalty
   len <- length(spot_values)
   times <-  private$surv_times
   #Withdrawals at the beginning of each period
   ts <- head(times[times != death_time], -1)
   ben <- rep(0, len)
   #Withdrawals at the beginning of each period
   ben[ts] <-  private$the_payoff$get_payoff()
   if (death_time < len){
    out <- calc_account(spot_values[1 : death_time], ben, fee, barrier, penalty)
    out <- rep(out, length.out=len)
    if(private$type == "Wa") {
     #For type Wa the at least the present value of guaranteed payments is
     #payed in case of death
     remaining_payments <- sum(discounts[death_time : len] * ben[death_time : len])
     out[death_time] <- max(out[death_time], remaining_payments)
    }
    #To avoid counting the surrender value at t = 0 during
    #the discounting in the do_mixed method
    out[1] <- 0
    #Adds the withdrawal payment to the cash flow
    out[ts] <- out[ts] + private$the_payoff$get_payoff()
    out[(death_time + 1) : len] <- 0
   } else {
    out <- calc_account(spot_values, ben, fee, barrier, penalty)
    #To avoid counting the surrender value at t = 0 during
    #the discounting in the do_mixed method
    out[1] <- 0
    #Adds the withdrawal payment to the cash flow
    out[ts] <- out[ts] + private$the_payoff$get_payoff()
   }
   out
  },
  survival_benefit = function(spot_values, death_time, time){
   times <-   private$surv_times

   if(time %in% times & time != death_time)
     out <- private$the_payoff$get_payoff()
   else out <- 0

   if(time == tail(times, 1) & time != death_time){
     fee <- private$the_fee$get()
     barrier <- private$the_barrier
     penalty <- private$penalty
     ben <- rep(0, length(spot_values))
     ben[head(times, -1)] <-  private$the_payoff$get_payoff()
     out <- calc_account(spot_values, ben, fee, barrier, penalty)
     out <- tail(out, 1)
   }
   out
  }
 ),
 private = list(
  #Type of GMWB rider
  type = "Wa",
  #Frequency of the withdrawals
  freq = "12m"
 )
)
