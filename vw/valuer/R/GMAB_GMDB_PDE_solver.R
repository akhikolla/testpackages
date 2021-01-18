#Copyright 2016-2017 Ivan Zoccolan

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



#'PDE Pricing of Variable Annuity
#'
#'\code{va_pde_pricer} returns the price of a VA with GMAB and GMDB
#'guarantees. The underlying fund is a GBM and the intensity of mortality
#'is deterministic. The fee has a  state-dependent structure.
#'@param F0 numeric scalar with the initial value of the underlying fund
#'@param r numeric scalar with the constant interest rate
#'@param sigma numeric scalar with the constant volatility
#'@param x numeric integer with the age of the insured
#'@param delta numeric scalar with the roll-up rate of the GMAB and GMDB
#'@param fee numeric scalar with the state-dependent base fee
#'@param beta numeric scalar with the state-dependent barrier
#'It should be greater than \code{F0}.
#'If set to \code{Inf} the fee structure becomes constant
#'@param T numeric integer with the maturity of the contract
#'@param dt numeric scalar with the discretization step of the
#'time dimension
#'@param dF numeric scalar with the discretization step for
#'the fund dimension
#'@param lambda function with the intensity of mortality.
#'Default is \code{\link{makeham}} with parameters
#'\code{x = 50, A = 0.0001, B = 0.00035, c = 1.075}
#'@param K function with the  surrender penalty.
#'@param Fmax numeric scalar with the maximum fund value
#'@return numeric scalar with the VA price
#'
#'@details
#'This function resolves the PDE in [MK2017] by means of the  finite
#'difference implicit method. It requires the package limSolve to be
#'installed.
#'@references
#'\enumerate{
#'  \item{[MK2017]}{ \cite{A. MacKay, M. Augustyniak, C. Bernard, and M.R.
#'    Hardy.
#'  Risk management of policyholder behavior in equity-linked life insurance.
#'  The Journal of Risk and Insurance, 84(2):661-690, 2017. DOI: 10.1111/jori.12094}}
#'}
#'@examples
#'lambda <- function(t) mu(t, x = 50, c1 = 90.43 , c2 = 10.36)
#'K <- function(t, T) {0.05 * ( 1 - t / T)^3}
#'
#'va <- va_pde_pricer(lambda = lambda, K = K, Fmax = 200)
#'
#'va
#'@importFrom stats integrate
#'@export

va_pde_pricer <- function(F0 = 100,  r = 0.03,  sigma = 0.165, x = 50,
                                 delta = 0.00, fee = 0.01, beta = 150, T = 5,
                                 dt = 0.1, dF = 0.1,
                                 lambda = function(t) makeham(t, x = 50, A = 0.0001,
                                                          B = 0.00035,
                                                          c = 1.075),
                                 K = function(t, T) {0.05 * ( 1 - t / T)^3},
                                 Fmax = 500){

  if(!requireNamespace("limSolve", quietly = TRUE)){
   stop("This function requires the package limSolve.
        Please run install.packages(\"limSolve\")")
  }

  #Life probability according to the intensity of
  #mortality function lambda

  lprob <- function(t, u) {

    exp(-stats::integrate(lambda, t, u)$value)

  }


  N <- round(T / dt) #Time dimension discretized in N intervals
  M <- round(Fmax / dF) #Possible fund values discretized in M intervals


  #Functions to calculate the coefficients of the linear equations
  #to determine V[i+1, j] based on V[i,j-1], V[i,j], V[i,j+1]

  A <- function(j) {
    0.5 * j * dt *  (r - ifelse(j*dF < beta, fee, 0)) - 0.5 * sigma^2 * j^2 * dt
  }

  B <- function(i, j){
    1 + j^2 * sigma^2 * dt + (r + lambda(i*dt))*dt

  }

  C <- function(j){
    -0.5 * j^2 * sigma^2 * dt - 0.5 * j * dt * (r - ifelse(j*dF < beta, fee, 0))

  }

  D <- function(i, j) {
    -lambda(i*dt) * dt * max(F0 * exp(delta * i*dt), j*dF)
  }

  #Continuation values

  V <- matrix(0, M+1, N+1)

  #Boundary conditions

  #Right bound
  G <- F0 * exp(delta * T)

  V[1:(M+1), N+1] <- sapply(rev(0:M)*dF, function(x) max(G, x))

  #Upper bound
  V[1, ] <- Fmax

  #Lower bound

  #V0_func calculates the lower bounds by using
  #the closed formula 4

  V0_func <- function(i){

    V0_integrand <- Vectorize(function(u) {

      esp <- (r - delta) * u  - r * i * dt

      exp(-esp) * lprob(i*dt, u) * lambda(u)

    })

    esp <- ((r - delta) *  T) - r * i * dt

    res <- exp(-esp) * lprob(i*dt, T)

    res <- res + integrate(V0_integrand, i*dt, T)$value

    res <- F0 * res

    res

  }

  for (i in 0:N)  V[M+1, i+1] <- V0_func(i)

  ###########################################

  #Defines the matrix AA of the linear system
  #AAx = bb where AA is tri-diagonal.

  lower_diag <- rep(0, M-2)
  upper_diag <- lower_diag
  diag <- rep(0, M-1)

  #Sets first the lower diagonal and the upper diagonal of the
  #AA matrix since they don't depend on the time index i

  for (j in 1:(M - 1))  {

    if (j > 1) lower_diag[j - 1] <- A(j)

    if (j < M - 1) upper_diag[j] <- C(j)

  }

  for (i in (N-1):0){

    #Calculates diagonal in the AA matrix
    #and add the d(i, j) term to the bb vector

    bb <- rep(0, M-1)

    for (j in 1:(M - 1))  {
      #Sets the diagonal of the AA matrix
      diag[j] <- B(i, j)

      bb[j] <- -D(i, j)
    }

    bb <- bb + V[M:2, i + 2]
    bb[1] <- bb[1] - A(1) * V[M+1, i + 1]
    bb[M-1] <- bb[M-1] - C(M-1)*V[1, i+1]

    #Solves the linear system using  Solve.tridiag from limSolve
    #which calls the dgtsv() routine from LAPACK
    #This is way faster than the standard solve function as it requires
    #just O(M) operations.

    V[M:2, i+1] <- limSolve::Solve.tridiag(lower_diag, diag, upper_diag, bb)

    #optimal stopping condition
    #The continuation value is equal to the withdrawal value
    eps = dF / 10
    for (j in M:2){
      if ((V[j, i+1] - (1 - K(i*dt, T))*(M - j + 1)*dF) < eps){
        V[j, i+1] <- (1 - K(i*dt, T))*(M -j + 1)*dF
      }
    }

  }

  ans <-  V[M+1 - F0/dF, 1]

  return(invisible(ans))

}
