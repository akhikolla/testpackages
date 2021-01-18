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




#Defines a base parameter class

parameters <- R6::R6Class("parameters",
 public = list(
  integral = function(t1, t2) {},
  integral_square = function(t1, t2) {}
 )
)



#' Constant parameter class
#' @description
#' Class providing a constant parameter object with methods to calculate
#' the integral of the parameter and  the squared parameter over a time span.
#' @docType class
#' @examples
#' r <- constant_parameters$new(0.01)
#' #Over the full year (365 days) the integral should evaluate to 0.01
#' r$integral(timeDate::timeDate("2016-07-09"), timeDate::timeDate("2017-07-09"))
#' #Over the full year the integral square should evaluate to 0.001
#' r$integral_square(timeDate::timeDate("2016-07-09"), timeDate::timeDate("2017-07-09"))
#' @export
#' @return Object of \code{\link{R6Class}}
#' @format \code{\link{R6Class}} object.
#' @section Methods:
#' \describe{
#'   \item{\code{integral}}{Calculates the integral given the initial
#'    and final times. The arguments are two \code{\link{timeDate}}
#'    object with the initial  and final times.
#'    It returns a \code{numeric} scalar with the integral}
#'   \item{\code{integral_square} (\code{public})}{Calculates the integral
#'    of the squared constant parameter given the initial and final times.
#'    The arguments are two \code{\link{timeDate}} object with the initial
#'    and final times. It returns a \code{numeric} scalar with the integral}
#'   \item{\code{get} (\code{public})}{get the constant}
#' }

constant_parameters <- R6::R6Class("constant_parameters",
 inherit = parameters,
 public = list(
  initialize = function(const, t_span){
   if (!missing(t_span)){
    if (is_positive_integer(t_span)) private$time_span <- t_span
    else stop(error_msg_4("t_span"))
   } else private$time_span <- 365

   if(!missing(const)){
    if (is_numeric_scalar(const)){
     private$constant <- const / private$time_span
     private$constant_square <- const ^ 2 / private$time_span
    } else stop(error_msg_5("const"))
   } else {
    private$constant <- 0.0
    private$constant_square <- 0.0
   }
  },
  integral = function(t1, t2){
   as.numeric((t2 - t1) * private$constant)
  },
  integral_square = function(t1, t2){
   as.numeric((t2 - t1) * private$constant_square)
  },
  get = function() private$constant
 ),
 private = list(
  #A positive integer which stores the time span (in days).
  time_span = "numeric",
  #A numeric scalar which stores the constant parameter.
  constant = "numeric",
  #A numeric scalar which stores the square of the constant parameter
  constant_square = "numeric"
 )
)



#' Surrender penalty class
#' @description
#' Class providing a surrender charge. It supports a constant surrender
#' charge (type 1) and two surrender charges decreasing with time,
#' (type 2 and type 3).
#' @docType class
#' @examples
#' #Sets a constant penalty
#' penalty <- penalty_class$new(type = 1, const = 0.03)
#' penalty$get()
#' penalty$set(0.04)
#' penalty$get()
#' #Sets a time decreasing penalty of type 2
#' penalty <- penalty_class$new(type = 2, const = 0.08, T = 10)
#' penalty$get(time = 0)
#' penalty$get(time = 2)
#' penalty$set(0.05)
#' penalty$get(time = 0)
#' #Sets a time decreasing penalty of type 3
#' penalty <- penalty_class$new(type = 3, const = 0.08, T = 10)
#' penalty$get(time = 0)
#' penalty$get(time = 2)
#' penalty$set(0.05)
#' penalty$get(time = 0)
#' @export
#' @return Object of \code{\link{R6Class}}
#' @format \code{\link{R6Class}} object.
#' @section Methods:
#' \describe{
#'   \item{\code{new}}{Initialization methods with arguments:
#'   \describe{
#'    \item{\code{type}}{type of the surrender charge. It can be
#'    1 (constant) or 2 or 3 (decreasing with time).}
#'    \item{\code{const}}{positive integer between 0 and 1
#'    with the maximum surrender charge.}
#'    \item{\code{T}}{Positive integer with expiry of the VA product.}
#'   }
#'  }
#'   \item{\code{get}}{get the surrender penalty. Argument is \code{time}
#'   a scalar in [0, T].}
#'   \item{\code{set}}{set the maximum surrender penalty.}
#'   \item{\code{get_type}}{get the type of the surrender penalty}
#' }


penalty_class <- R6::R6Class("penalty_class",
 public = list(
  initialize = function(type, const = 1, T = 5, ...){

   if(is_identical_to_any(type, c(1, 2, 3)))
     private$type <- type
   else stop("Type must be either 1, 2 or 3")

   switch(type,
    #Constant penalty
    '1' = {
     if (!missing(const))
      if (is_between(const, 0, 1)) private$the_penalty <- const
      else stop(error_msg_8("const"))
     else private$the_penalty <- 1

     #
     private$T <- T

    },
    #Types of penalty decreasing with time
    '2' = {
     if (!missing(const))
      if (is_between(const, 0, 1)) private$const <- const
      else stop(error_msg_8("const"))
     else private$const <- 0.05

     if(!missing(T))
      if(is_positive_integer(T)){
        private$T <- T
        private$the_penalty <- function(t) private$const * (1 - t / T)^3
      } else stop(error_msg_4("T"))
     else stop(error_msg_4("T"))
    },
    '3' = {
     if (!missing(const))
      if (is_between(const, 0, 1)) private$const <- const
      else stop(error_msg_8("const"))
     else private$const <- 0.08

     if(!missing(T))
      if(is_positive_integer(T)){
        private$T <- T
        private$the_penalty <- function(t) 1 - exp(-(private$const / T) * (T - min(T, t)))
      } else stop(error_msg_4("T"))
     else stop(error_msg_4("T"))
    }
   )
  },
  get = function(time){
    if(missing(time)) time <- 0
    if (time > private$T | time < 0) 0
    else  if (private$type == 1)
      private$the_penalty
    else private$the_penalty(time)
  },
  set = function(const){
    if (!missing(const))
     if (is_between(const, 0, 1)){
      switch(private$type,
       '1' = {
        private$const <- const
        private$the_penalty <- const
       },
       '2' = {
        private$const <- const
       },
       '3' = {
        private$const <- const / 10
       }
      )
     } else stop(error_msg_8("const"))
    else stop(error_msg_8("const"))
  },
  get_type = function() private$type
 ),
 private = list(
  #Type of surrender penalty
  type = NULL,
  #Surrender penalty
  the_penalty= NULL,
  #Const for time decreasing functions
  const = NULL,
  #Time up to which the surrender penalty is applicable
  T = numeric(0)
 )
)
