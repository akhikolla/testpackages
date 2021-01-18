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



.onUnload <- function (libpath) {
  library.dynam.unload("valuer", libpath)
}


#Utility functions

is_numeric_scalar <- function(x, ...) tryCatch(inherits(x, "numeric") & isTRUE(length(x) == 1),
                                          error = function(e) FALSE)

is_positive_scalar <- function(x, ...) {
  tryCatch(inherits(x, "numeric") & isTRUE(length(x) == 1) &
    isTRUE(x > 0), error = function(e) FALSE)
}

is_not_negative_scalar <- function(x, ...){
  tryCatch(inherits(x, "numeric") & isTRUE(length(x) == 1) &
    isTRUE(x >= 0), error = function(e) FALSE)
}

is_integer <- function(x, ...){

  if (tryCatch((inherits(x, "numeric")), error = function(e) FALSE))
    return(isTRUE(x == round(x)))
  else FALSE

}

is_positive_integer <- function(x, ...){

  if (tryCatch((inherits(x, "numeric")), error = function(e) FALSE))
    return(isTRUE(x == round(x)) & isTRUE(x > 0))
  else FALSE

}

are_dates <- function(dates, ...) tryCatch(inherits(dates, "timeDate"), error = function(e) FALSE)

is_date <- function(date, ...) tryCatch((inherits(date, "timeDate") & (length(date) == 1)),
                                   error = function(e) FALSE )

same_length <- function(x,y, ...) tryCatch(isTRUE(length(x) == length(y)), error = function(e) FALSE)

is_payoff <- function(x, ...) tryCatch(inherits(x, "payoff"), error = function(e) FALSE)

is_between <- function(x,lower,upper, ...) tryCatch(isTRUE(x >= lower) & isTRUE(x <= upper), error = function(e) FALSE)

istrue <- function(x, ...) tryCatch(isTRUE(x), error = function(e) FALSE)

is_identical_to_any <-  function(x, y, ...) tryCatch(any(sapply(y, function(z) identical(x, z))),
                                                     error = function(e) FALSE)
#Error messages,

error_msg_1 <- function(object) paste("Argument must be a", object, "object")
error_msg_2 <- function(object) paste("Arguments must be", object, "objects")
error_msg_3 <- "Argument must be a positive scalar"
error_msg_3_ <- function(arg) paste("Argument", arg, "must be a positive scalar")
error_msg_1_ <- function(arg, object) paste("Argument ", arg, "must be a ", object, "object")
error_msg_4 <- function(arg) paste("Argument", arg, "must be a positive integer")
error_msg_5 <- function(arg) paste("Argument", arg, "must be a numeric scalar")
error_msg_6 <- "The delivery time should be one of the possible product dates"
error_msg_7 <- function(arg) paste("Argument", arg, "must be a non negative scalar")
error_msg_8 <- function(arg) paste("Argument", arg, "must be between 0 and 1")
error_msg_9 <- function(arg) paste("Argument", arg, "must be a numeric vector")

error_msg_10 <- function() {
                  msg1 <- "Allowed units are"
                  msg2 <- "\"m\" for 4 weeks, \"w\" for weeks,"
                  msg3 <- "\"d\" for days"
                  paste(msg1, msg2, msg3)
}

error_msg_11 <- function(arg1, arg2) paste(arg1, "must be after", arg2)


#' Normalizes a timeDate sequence into year fractions
#'@param times A \code{\link{timeDate}} sequence

yr_fractions <- function(times){

 t_periods <- timeDate::periods(times, period = "12m", by="12m")
 t_diffs <- as.numeric(timeDate::difftimeDate(t_periods$to, t_periods$from)) + 1
 ind <- seq_along(t_diffs)
 out <- sapply(ind, function(i) i - 1 + seq(t_diffs[i]) / t_diffs[i])
 c(0, unlist(out))

}


#Takes square root if positive otherwise returns zero.
#To be used with mean reverting squared root processes (CIR SDE)

#'Square root utility function
#'@description
#'Takes square root if positive otherwise returns zero.
#'To be used with mean reverting squared root processes (CIR SDE)
#'@export
#'@param x numeric scalar

sq <- function(x) ifelse(x > 0, sqrt(x), 0)

#Deterministic intensity of mortality ( Weibull )

#' Weibull intensity of mortality
#' @param t time as numeric scalar
#' @param x age as numeric scalar
#' @param c1 numeric scalar
#' @param c2 numeric scalar
#' @export

mu  <-  function(t, x, c1, c2) {(c1^(-c2))*c2*((x + t)^(c2 -1))}

#Deterministic intensity of mortality ( Makeham )

#' Makeham's intensity of mortality
#' @param t time as numeric scalar
#' @param x age as numeric scalar
#' @param A numeric scalar
#' @param B numeric scalar
#' @param c numeric scalar
#' @export

makeham <- function(t, x, A, B, c) A + B * (c ^ (x + t))


#Cumulative standard deviation
#http://stackoverflow.com/questions/2765374/
#efficient-calculation-of-matrix-cumulative-standard-deviation-in-r

cumsd <- function(x) {
  ind_na <- !is.na(x)
  nn <- cumsum(ind_na)
  x[!ind_na] <- 0
  sq(cumsum(x^2) / (nn-1) - (cumsum(x))^2/(nn-1)/nn)
}


