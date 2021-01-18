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


#Defines the statistics gatherer class
#The statistics gatherer object will be changed its state by functions
#and other object methods and must retain the change at each step.
#In R we will need to use classes with reference semantics such as
#Reference Classes or R6


#This class defines the interface only.
#Implementation will be done in subclasses inheriting this one.


gatherer <-  R6::R6Class("gatherer",
 public = list(
   dump_result = function(result) {},
   get_results = function() {}
 ),
 private = list(
   values = NULL
 )
)


#'Monte Carlo gatherer
#' @description Class which defines a gatherer for the Monte Carlo simulated
#' values. It has methods to return the Monte Carlo estimate and Monte Carlo
#' Standard Error of the estimate as well as a convergence table.
#' @docType class
#' @importFrom ggplot2 ggplot
#' @export
#' @return Object of \code{\link{R6Class}}
#' @format \code{\link{R6Class}} object.
#' @section Methods:
#' \describe{
#'   \item{\code{new}}{Constructor method}
#'   \item{\code{dump_result}}{Saves the argument \code{result} which is
#'   a numeric scalar}
#'   \item{\code{get_results}}{Returns the  Monte Carlo estimate and the
#'   (estimated) Monte Carlo Standard Error of the estimate}
#'   \item{\code{convergence_table}}{Returns the convergence table}
#'   \item{\code{plot}}{Plots a Monte Carlo convergence graph at 95\% level}
#' }
#'
mc_gatherer  <- R6::R6Class("mc_gatherer", inherit = gatherer,
 public = list(
  initialize = function(){ private$values <- 0.0 },
  dump_result = function(result){
   if (inherits(result, "numeric"))
     private$values <-  result
   else stop(error_msg_9("result"))
  },
  get_results = function(){
   data.frame(mean = mean(private$values),
              se  =  sd(private$values) / sqrt(length(private$values))
   )},
  convergence_table = function(){
    cum_mean <- cumsum(private$values) / seq_along(private$values)
    cum_se <- cumsd(private$values) / sqrt(seq_along(private$values))
    private$conv_table <- data.frame(mean = cum_mean, se = cum_se)
    private$conv_table
  },
  path_done = function() length(private$values),
  plot = function(){
    table <- self$convergence_table()
    df <- data.frame(No_of_simulations = seq(self$path_done()),
                     lower = table$mean - 1.96 * table$se, mean = table$mean,
                     upper = table$mean + 1.96 * table$se)
    df <- tail(df, -1)
    p <-  ggplot2::ggplot(df, ggplot2::aes(x = No_of_simulations, y = mean,
                                           colour = "Estimates"))
    p <- p + ggplot2::labs(x = "Number of simulated paths", y = "")
    p <- p + ggplot2::geom_line(ggplot2::aes(y = lower, colour = "Bounds"))
    p <- p + ggplot2::geom_line(ggplot2::aes(y = upper, colour = "Bounds"))
    p <- p + ggplot2::geom_line() + ggplot2::labs(colour ="")
    p <- p + ggplot2::ggtitle("Monte Carlo convergence graph at 95% level")
    p <- p + ggplot2::theme_bw()
    p
  }
 ),
 private = list(
  conv_table = NULL
 )
)


#'Simple data gatherer
#' @description Class which defines a simple data gatherer to hold estimates
#' calculated in a loop.
#' @docType class
#' @export
#' @return Object of \code{\link{R6Class}}
#' @format \code{\link{R6Class}} object.
#' @section Methods:
#' \describe{
#'   \item{\code{new}}{Constructor method}
#'   \item{\code{dump_result}}{Saves the argument \code{result} which is
#'   a numeric scalar}
#'   \item{\code{get_results}}{Returns a numeric vector with the point
#'   estimates.}
#' }
#'

data_gatherer <- R6::R6Class("data_gatherer", inherit = gatherer,
  public = list(
    dump_result = function(result){
      if (inherits(result, "numeric"))
       private$values <- c(private$values, result)
      else stop(error_msg_9("result"))
    },
    get_results = function() private$values

  )
)
