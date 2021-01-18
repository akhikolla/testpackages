
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



#' Variable Annuity pricing engine with general fund processes and Weibull mortality
#' @description
#' Class providing a variable annuity pricing engine where the underlying
#' reference fund is specified by an arbitrary system of stochastic
#' differential equations. In contrast, the interest rates is constant and
#' the intensity of mortality is deterministic and given by the Weibull
#' function.
#' The fund paths are simulated by means of the
#' \href{https://CRAN.R-project.org/package=yuima}{yuima} package. \cr
#' The value of the VA contract is estimated by means of the Monte Carlo
#' method if the policyholder cannot surrender (the so called "static"
#' approach), and by means of Least Squares Monte Carlo in case the
#' policyholder can surrender the contract (the "mixed" approach).\cr
#' See \bold{References} -\code{[BMOP2011]} for a description of the mixed
#' and static approaches and the algorithm implemented by this class,
#' \code{[LS2001]} for Least Squares Monte Carlo and \code{[YUIMA2014]}
#' for \code{yuima}.
#' @docType class
#' @importFrom orthopolynom laguerre.polynomials
#' @importFrom RcppEigen fastLmPure
#' @importFrom yuima simulate get.zoo.data
#' @export
#' @return Object of \code{\link{R6Class}}
#' @format \code{\link{R6Class}} object.
#' @section Methods:
#' \describe{
#'  \item{\code{new}}{Constructor method with arguments:
#'   \describe{
#'    \item{\code{product}}{A \code{\link{va_product}}
#'    object with the VA product.}
#'    \item{\code{financial_parms}}{A list of parameters
#'    specifying the financial processes.
#'    See \code{\link{financials_BZ2016bis}} for an example.}
#'    \item{\code{interest}}{\code{\link{constant_parameters}} object with the
#'    constant interest rate}
#'    \item{\code{c1}}{\code{numeric} scalar argument of the intensity
#'    of mortality function \code{\link{mu}}}
#'    \item{\code{c2}}{\code{numeric} scalar argument of the intensity
#'    of mortality function \code{\link{mu}}}
#'   }
#'  }
#'  \item{\code{death_time}}{Returns the time of death index. If the
#'  death doesn't occur during the product time-line it returns the
#'  last index of the product time-line plus one.}
#'  \item{\code{simulate_financial_paths}}{Simulates \code{npaths} paths
#'  of the underlying fund of the VA contract and the discount factors
#'  (interest rate) and saves them into private fields for later use.}
#'  \item{\code{simulate_mortality_paths}}{Simulates \code{npaths} paths
#'  of the intensity of mortality and saves them into private fields
#'  for later use.}
#'  \item{\code{get_fund}}{Gets the \code{i}-th path of the underlying fund
#'  where \code{i} goes from 1 to \code{npaths}.}
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
#'   \describe{
#'    \item{\code{the_gatherer}}{\code{gatherer} object to hold
#'     the point estimates}
#'     \item{\code{npaths}}{positive integer with the number of paths to
#'     simulate}
#'     \item{\code{degree}}{positive integer with the maximum degree of
#'     the weighted Laguerre polynomials used in the least squares by LSMC}
#'     \item{\code{freq}}{string which contains the frequency of the surrender
#'     decision. The default is \code{"3m"} which corresponds to deciding every
#'     three months if surrendering the contract or not.}
#'     \item{\code{simulate}}{boolean to specify if the paths should be
#'     simulated from scratch, default is TRUE.}
#'   }
#'  }
#'  \item{\code{get_discount}}{Arguments are \code{i,j}.
#'  Gets the \code{j}-th discount factor corresponding to the \code{i}-th
#'  simulated path of the discount factors.}
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
#'  \item{[YUIMA2014]}{ \cite{Alexandre Brouste, Masaaki Fukasawa, Hideitsu
  #'   Hino, Stefano M. Iacus, Kengo Kamatani, Yuta Koike, Hiroki Masuda,
  #'   Ryosuke Nomura, Teppei Ogihara, Yasutaka Shimuzu, Masayuki Uchida,
  #'   Nakahiro Yoshida (2014). The YUIMA Project: A Computational
  #'   Framework for Simulation and Inference of Stochastic Differential
  #'   Equations. Journal of Statistical Software, 57(4), 1-51.
  #'   URL http://www.jstatsoft.org/v57/i04/.}}
#'  }
#'@examples
#'#Sets up the payoff as a roll-up of premiums with roll-up rate 2%
#'
#'rate <- constant_parameters$new(0.02)
#'
#'premium <- 100
#'rollup <- payoff_rollup$new(premium, rate)
#'
#'#constant interest rate
#'r <- constant_parameters$new(0.03)
#'
#'#Five years time-line
#'begin <- timeDate::timeDate("2016-01-01")
#'end <- timeDate::timeDate("2020-12-31")
#'
#'#Age of the policyholder.
#'age <- 50
#'# A constant fee of 2% per year (365 days)
#'fee <- constant_parameters$new(0.02)
#'
#'#Barrier for a state-dependent fee. The fee will be applied only if
#'#the value of the account is below the barrier
#'barrier <- 200
#'#Withdrawal penalty applied in case the insured surrenders the contract
#'#It is a constant penalty in this case
#'penalty <- penalty_class$new(type = 1, 0.02)
#'#Sets up the contract with GMAB guarantee
#'contract <- GMAB$new(rollup, t0 = begin, t = end, age = age, fee = fee,
#'barrier = barrier, penalty = penalty)
#'
#'#Sets up a gatherer of the MC point estimates
#'the_gatherer  <- mc_gatherer$new()
#'no_of_paths <- 10
#'
#'#Sets up the pricing engine
#'engine <- va_sde_engine3$new(contract, financials_BZ2016bis, interest = r)
#'
#'#Estimates the contract value by means of the static approach
#'
#'engine$do_static(the_gatherer, no_of_paths)
#'the_gatherer$get_results()
#'
#'
#'#Estimates the contract value by means of the mixed approach
#'#To compare with the static approach we don't simulate the underlying
#'#fund paths again.
#'
#'the_gatherer_2 <- mc_gatherer$new()
#'
#'engine$do_mixed(the_gatherer_2, no_of_paths, degree = 3, freq = "3m",
#'simulate = FALSE)
#'the_gatherer_2$get_results()

va_sde_engine3 <- R6::R6Class("va_sde_engine2", inherit = va_engine,
 public = list(
  initialize = function(product, financial_parms, interest, c1, c2){
   super$initialize(product)
   private$times <- product$get_times()
   no_time_intervals <- length(private$times) - 1
   private$financial_parms <- financial_parms
   private$financial_model <- do.call(yuima::setModel,financial_parms[[2]])

   if(!missing(interest))
     if(inherits(interest, "parameters")){
       private$r <- interest
     } else stop(error_msg_1_("interest", "parameters"))
   else stop(error_msg_1_("interest", "parameters"))
   #Sets up and stores the discount factor vector
   cf_times <- private$times
   t0 <- cf_times[1]
   log_discounts <- vector(mode="numeric", length = length(cf_times))
   for (i in seq_along(cf_times))
     log_discounts[i] <- -private$r$integral(t0, cf_times[i])
   private$discounts <- exp(log_discounts)

   if(!missing(c1))
    if (is_positive_scalar(c1))
     private$mu_1 <- c1
    else stop(error_msg_5("c1"))
   else private$mu_1 <- 88.14778
   if(!missing(c2))
    if(is_positive_scalar(c2))
     private$mu_2 <- c2
    else stop(error_msg_5("c2"))
   else private$mu_2 <- 10.002

   private$samp <- yuima::setSampling(
   Terminal = tail(private$the_product$times_in_yrs(), 1),
                                    n = no_time_intervals)

   private$mu_integrals <- self$simulate_mortality_paths()




  },
  simulate_financial_paths = function(npaths){
   ind <- private$financial_parms[[3]]
   #Builds parameter list for yuima::simulate
   parms <- list(object = private$financial_model,
                 xinit = private$financial_parms[[1]]$xinit,
                 sampling = private$samp,
                 true.parameter =  private$financial_parms[[1]])

   len <-  length(private$the_product$get_times())
   #Sets storage for fund
   private$fund <- matrix(NA, npaths, len)
   #Simulates the underlying fund spot prices
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
    for (i in seq(npaths)){
     private$fund[i, ] <- exp(as.numeric(data_paths[[i]][[ind[1]]]))
    }
   } else for (i in seq(npaths)){
      #Imported in NAMESPACE
      #yuima::simulate
      #yuima::get.zoo.data
      zoo_paths <- do.call(simulate, parms)
      data_paths <- get.zoo.data(zoo_paths)
      private$fund[i, ] <- exp(as.numeric(data_paths[[ind[1]]]))
    }
   },
   get_fund = function(i) private$fund[i, ],
   get_discount = function(i,j) private$discounts[j],
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
   #Stores the yuima financial model
   financial_model = "yuima.model-class",
   #Stores the financial parameters needed to
   #set the model above by yuima::setModel and
   #run yuima::simulate
   financial_parms = "list",
   #Intensity of mortality parameters
   mu_1 = 88.14778,
   mu_2 = 10.002,
   #Product times
   times = "timeDate",
   #time-line of the intensity of mortality
   mu_integrals = "numeric",
   #Stores the times of a simulate path
   samp = "yuima.sampling-class",
   #Matrix to hold the simulated paths of the
   #underlying fund
   fund = "matrix",
   #Matrix to hold the simulated paths of the
   #stochastic interest rate
   r = "vector",
   #stochastic discount factors
   discounts = "vector",
   #Method to get Laguerre polynomials of state variables.
   #Arguments are:
   #paths - numeric vector of indexes of the paths to consider
   #time - numeric scalar with the time index
   #degree - positive scalar with the max degree of
   #the Laguerre polynomials
   bases = function(paths, time, degree){
    #orthopolynom::laguerre.polynomials it's imported in NAMESPACE
    res <- laguerre.polynomials(degree, normalized = TRUE)
    x <- private$fund[paths, time]
    #Normalizes to avoid underflows in calculating
    #the exponential below.
    x <- x / private$the_product$get_premium()

    x <- sapply(seq_along(res), function(i){
            exp(-0.5 * x) * (as.function(res[[i]])(x))
            })

   }
  )
)
