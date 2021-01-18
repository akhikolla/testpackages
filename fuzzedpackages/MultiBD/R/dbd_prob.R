##########################################################
### Transition probabilities of a death/birth-death process
##########################################################

#' Transition probabilities of a death/birth-death process
#'
#' Computes the transition pobabilities of a death/birth-death process
#' using the continued fraction representation of its Laplace transform
#' @param t time
#' @param a0 total number of type 1 particles at \code{t = 0}
#' @param b0 total number of type 2 particles at \code{t = 0}
#' @param mu1 death rate of type 1 particles (a two variables function)
#' @param lambda2 birth rate of type 2 particles (a two variables function)
#' @param mu2 death rate function of type 2 particles (a two variables function)
#' @param gamma transition rate from type 2 particles to type 1 particles (a two variables function)
#' @param a lower bound for the total number of type 1 particles (default \code{a = 0})
#' @param B upper bound for the total number of type 2 particles
#' @param nblocks number of blocks
#' @param tol tolerance
#' @param computeMode computation mode
#' @param nThreads number of threads
#' @param maxdepth maximum number of iterations for Lentz algorithm
#' @return a matrix of the transition probabilities
#' @references Ho LST et al. 2016. "Birth(death)/birth-death processes and their computable transition probabilities with statistical applications". In review.
#'
#' @examples
#' \dontrun{
#' data(Eyam)
#'
#' loglik_sir <- function(param, data) {
#'   alpha <- exp(param[1]) # Rates must be non-negative
#'   beta  <- exp(param[2])
#'
#'   # Set-up SIR model
#'   drates1 <- function(a, b) { 0 }
#'   brates2 <- function(a, b) { 0 }
#'   drates2 <- function(a, b) { alpha * b     }
#'   trans12 <- function(a, b) { beta  * a * b }
#'
#'   sum(sapply(1:(nrow(data) - 1), # Sum across all time steps k
#'              function(k) {
#'                log(
#'                  dbd_prob(  # Compute the transition probability matrix
#'                    t  = data$time[k + 1] - data$time[k], # Time increment
#'                    a0 = data$S[k], b0 = data$I[k],       # From: S(t_k), I(t_k)
#'                    drates1, brates2, drates2, trans12,
#'                    a = data$S[k + 1], B = data$S[k] + data$I[k] - data$S[k + 1],
#'                    computeMode = 4, nblocks = 80         # Compute using 4 threads
#'                  )[1, data$I[k + 1] + 1]                 # To: S(t_(k+1)), I(t_(k+1))
#'                )
#'              }))
#'   }
#'
#'   loglik_sir(log(c(3.204, 0.019)), Eyam) # Evaluate at mode
#' }
#'
#' # Birth-death-shift model for transposable elements
#'
#' lam = 0.0188; mu = 0.0147; v = 0.00268; # birth, death, shift rates
#'
#' drates1 <- function(a, b) { mu * a }
#' brates2 <- function(a, b) { lam * (a + b) }
#' drates2 <- function(a, b) { mu * b }
#' trans12 <- function(a, b) { v * a }
#'
#' # Get transition probabilities
#' p <- dbd_prob(t = 1, a0 = 10, b0 = 0,
#'               drates1, brates2, drates2, trans12,
#'               a = 0, B = 50)
#'
#' @export
dbd_prob <- function(t, a0, b0, mu1, lambda2, mu2, gamma, a=0, B,
                    nblocks=256, tol=1e-12, computeMode=0, nThreads=4,
                    maxdepth=400) {

  ###################
  ### Input checking
  ###################

  ## a>=0, a<=a0, B >=a0+b0-a
  if(a < 0) stop("a cannot be negative.")
  if (a > a0) stop("a cannot be bigger than a0.")
  if (B < a0+b0-a) stop("B is too small.")

  ###########################################################
  ### store rate values in matrices
  ### rates are transformed to fit birth/birth-death process
  ###########################################################

  l1 <- function(u,v){
    if (v > B) return(0)
    return(mu1(a0-u,B-v))
  }

  l2 <- function(u,v){
    if (v > B) return(0)
    return(mu2(a0-u,B-v))
  }

  m2 <- function(u,v){
    if (v > B) return(0)
    return(lambda2(a0-u,B-v))
  }

  g <- function(u,v){
    if (v > B) return(0)
    return(gamma(a0-u,B-v))
  }

  ###########################
  ### Call bbd_prob function
  ###########################

  res = matrix(0, nrow=a0-a+1, ncol=B+1)
  res[(a0-a+1):1,(B+1):1] = bbd_prob(t, 0, B-b0, l1, l2, m2, g, A=a0-a, B,
                                     nblocks, tol, computeMode, nThreads, maxdepth)

  colnames(res) = 0:B
  rownames(res) = a:a0

  return(res)
}

# function() {
#
#   sum(
#     sapply(
#       1:(nrow(data) - 1), # Sum across all time steps k
#       function(k) {
#         log(
#           dbd_prob( # Compute the transition probability matrix
#             t = data$time[k + 1] - data$time[k], # Time increment
#             a0 = data$S[k], b0 = data$I[k], # From: S(t_k), I(t_k)
#             drates1, brates2, drates2, trans12,
#                  a = data$S[k + 1],
#                  B = data$S[k] + data$I[k] - data$S[k + 1],
#                  computeMode = 4, nblocks = 80 # Compute using 4 threads
#             )[1, data$I[k + 1] + 1] # To: S(t_(k+1)), I(t_(k+1))
#           )}))
# }
