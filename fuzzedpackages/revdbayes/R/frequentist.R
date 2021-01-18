# =========================== gp_mle ===========================

gp_mle <- function(gp_data) {
  # Maximum likelihood estimation for the generalized Pareto distribution
  #
  # Performs maximum likelihood estimation for the generalized Pareto
  # distribution.  Uses the function \code{gpdmle} associated with
  # Grimshaw (1993), which returns MLEs of sigma and k = - \xi.
  #
  # Grimshaw, S. D. (1993) Computing Maximum Likelihood Estimates
  #   for the Generalized Pareto Distribution.  Technometrics, 35(2), 185-191.
  #   and Computing (1991) 1, 129-133. https://doi.org/10.1007/BF01889987.
  #
  # Args:
  #   gp_data : A numeric vector containing positive values, assumed to be a
  #             random sample from a generalized Pareto distribution.
  #
  # Returns:
  #   A list with components
  #     mle  : A numeric vector.  MLEs of GP parameters sigma and xi.
  #     nllh : A numeric scalar.  The negated log-likelihood at the MLE.
  #
  # Call Grimshaw (1993) function, note: k is -xi, a is sigma
  pjn <- grimshaw_gp_mle(gp_data)
  temp <- list()
  temp$mle <- c(pjn$a, -pjn$k)  # mle for (sigma,xi)
  sc <- rep(temp$mle[1], length(gp_data))
  xi <- temp$mle[2]
  temp$nllh <- sum(log(sc)) + sum(log(1 + xi * gp_data / sc) * (1 / xi + 1))
  return(temp)
}

# =========================== gp_pwm ===========================

#' Probability-weighted moments estimation of generalised Pareto parameters
#'
#' Uses the methodology of Hosking and Wallis (1987) to estimate the parameters
#' of the generalised Pareto (GP) distribution.
#'
#' @param gp_data A numeric vector of raw data, assumed to be a random sample
#'   from a probability distribution.
#' @param u A numeric scalar.  A threshold.  The GP distribution is fitted to
#'   the excesses of \code{u}.
#' @return A list with components
#'   \itemize{
#'     \item{\code{est}:} A numeric vector.  PWM estimates of GP parameters
#'       \eqn{\sigma} (scale) and \eqn{\xi} (shape).
#'     \item{\code{se}:} A numeric vector.  Estimated standard errors of
#'       \eqn{\sigma} and \eqn{\xi}.
#'     \item{\code{cov}:} A numeric matrix.  Estimate covariance matrix of the
#'       the PWM estimators of \eqn{\sigma} and \eqn{\xi}.
#'   }
#' @references Hosking, J. R. M. and Wallis, J. R. (1987) Parameter and
#'  Quantile Estimation for the Generalized Pareto Distribution.
#'  Technometrics, 29(3), 339-349. \url{https://doi.org/10.2307/1269343}.
#' @seealso \code{\link{gp}} for details of the parameterisation of the GP
#'   distribution.
#' @examples
#' u <- quantile(gom, probs = 0.65)
#' gp_pwm(gom, u)
#' @export
gp_pwm <- function(gp_data, u = 0) {
  # Probability weighted moments estimation for the generalized Pareto
  # distribution.
  #
  # Args:
  #   gp_data : A numeric vector of raw data, assumed to be a random sample
  #             from a probability distribution.
  #   u       : A numeric scalar.  A threshold.  The GP distribution is
  #             fitted to the excesses of u.
  # Returns:
  #   A list with components
  #     est  : A numeric vector.  PWM estimates of GP parameters sigma and xi.
  #     se   : A numeric vector.  Estimated standard errors of sigma and xi.
  #    cov   : A numeric matrix.  Estimate covariance matrix of the the PWM
  #            estimators of sigma and xi.
  #
  n <- length(gp_data)
  exceedances <- gp_data[gp_data > u]
  excess <- exceedances - u
  nu <- length(excess)
  xbar <- mean(excess)
  a0 <- xbar
  gamma <- -0.35
  delta <- 0
  pvec <- ((1:nu) + gamma) / (nu + delta)
  a1 <- mean(sort(excess) * (1 - pvec))
  xi <- 2 - a0 / (a0 - 2 * a1)
  sigma <- (2 * a0 * a1) / (a0 - 2 * a1)
  pwm <- c(sigma, xi)
  names(pwm) = c("sigma","xi")
  denom <- nu * (1 - 2 * xi) * (3 - 2 * xi)
  if (xi > 0.5) {
    denom <- NA
    warning("Asymptotic Standard Errors not available for PWM when xi>0.5.")
  }
  one <- (7 - 18 * xi + 11 * xi^2 - 2 * xi^3) * sigma ^ 2
  two <- (1 - xi) * (1 - xi + 2 * xi^2) * (2 - xi) ^ 2
  cov <-  - sigma * (2 - xi) * (2 - 6 * xi + 7 * xi ^ 2 - 2 * xi ^ 3)
  pwm_varcov <- matrix(c(one, cov, cov, two), 2)/denom
  colnames(pwm_varcov) <- c("sigma","xi")
  rownames(pwm_varcov) <- c("sigma","xi")
  se <- sqrt(diag(pwm_varcov))
  return(list(est = pwm, se = se, cov = pwm_varcov))
}

# =========================== gp_lrs ===========================

#' Linear Combinations of Ratios of Spacings estimation of generalised Pareto
#' parameters
#'
#' Uses the Linear Combinations of Ratios of Spacings (LRS) methodology of
#' (Reiss and Thomas, 2007, page 134) to estimate the parameters of the
#' generalised Pareto (GP) distribution, based on a sample of positive values.
#'
#' @param x A numeric vector containing only \strong{positive} values, assumed
#'   to be a random sample from a generalized Pareto distribution.
#' @return A numeric vector of length 2.  The estimates of the scale parameter
#'   \eqn{\sigma} and the shape parameter \eqn{\xi}.
#' @seealso \code{\link{gp}} for details of the parameterisation of the GP
#'   distribution.
#' @references Reiss, R.-D., Thomas, M. (2007) Statistical Analysis of
#'   Extreme Values with Applications to Insurance, Finance, Hydrology and
#'   Other Fields.Birkhauser.
#'   \url{https://doi.org/10.1007/978-3-7643-7399-3}.
#' @examples
#' u <- quantile(gom, probs = 0.65)
#' gp_lrs((gom - u)[gom > u])
#' @export
gp_lrs <- function(x) {
  # LRS estimation for the generalized Pareto distribution.
  #
  # Args:
  #   x : A numeric vector containing positive values, assumed to be a
  #       random sample from a generalized Pareto distribution.
  #
  # Returns:
  #   A numeric vector.  Estimates of parameters sigma and xi.
  #
  n <- length(x)                             # sample size
  x <- sort(x)                               # put data in ascending order
  q0 <- 1 / (n + 1)
  q2 <- n / (n + 1)                          # for sample minimum and maximum
  a <- sqrt((1 - q2) / (1 - q0))
  q1 <- 1 - a *(1 - q0)                      # `middle' quantile
  n0 <- 1
  n1 <- round((n + 1) * q1)
  n2 <- n                                    # corresponding order statistics
  ns <- c(n0,n1,n2)
  qs <- c(q0,q1,q2)
  xs <- x[ns]
  r_hat <- (xs[3] - xs[2]) / (xs[2] - xs[1])
  xi_hat <- -log(r_hat) / log(a)
  sigma_hat <- xi_hat * xs[3] / ((1 - q2) ^ (-xi_hat) - 1)
  return(c(sigma_hat, xi_hat))
}

# =========================== gp_obs_info ===========================

gp_obs_info <- function(gp_pars, y) {
  # Observed information for the generalized Pareto distribution
  #
  # Calculates the observed information matrix for a random sample \code{y}
  # from the generalized Pareto distribution, i.e. the negated Hessian matrix
  # of the generalized Pareto log-likelihood, evaluated at \code{gp_pars}.
  #
  # Args:
  #   gp_pars : A numeric vector. Parameters sigma and xi of the
  #   generalized Pareto distribution.
  #   y       : A numeric vector. A sample of positive data values.
  #
  # Returns:
  #   A 2 by 2 numeric matrix.  The observed information matrix.
  #
  s <- gp_pars[1]
  x <- gp_pars[2]
  i <- matrix(NA, 2, 2)
  i[1,1] <- -sum((1 - (1 + x) * y * (2 * s + x * y) / (s + x * y) ^ 2) / s ^ 2)
  i[1,2] <- i[2,1] <- -sum(y * (1 - y / s) / (1 + x * y / s) ^ 2 / s ^ 2)
  i[2,2] <- sum(2 * log(1 + x * y / s) / x ^ 3 - 2 * y / (s + x * y) / x ^ 2 -
                  (1 + 1 / x) * y ^ 2 / (s + x * y) ^ 2)
  return(i)
}

# =========================== grimshaw_gp_mle ===========================

#' Maximum likelihood estimation of generalised Pareto parameters
#'
#' Uses the methodology of Grimshaw (1993) to find the MLEs of the parameters
#' of the generalised Pareto distribution, based on a sample of positive
#' values.  The function is essentially the same as that made available with
#' Grimshaw (1993), with only minor modifications.
#' @param x A numeric vector containing only \strong{positive} values, assumed
#'   to be a random sample from a generalized Pareto distribution.
#' @return A numeric vector of length 2.  The estimates of the \strong{negated}
#'   shape parameter \eqn{k (= -\xi)} and the scale parameter
#'   \eqn{a (= \sigma)}.
#' @references Grimshaw, S. D. (1993) Computing Maximum Likelihood Estimates
#'   for the Generalized Pareto Distribution.  Technometrics, 35(2), 185-191.
#'   and Computing (1991) 1, 129-133.
#'   \url{https://doi.org/10.1080/00401706.1993.10485040}.
#' @seealso \code{\link{gp}} for details of the parameterisation of the GP
#'   distribution, in terms of \eqn{\sigma} and \eqn{\xi}.
#' @examples
#' u <- quantile(gom, probs = 0.65)
#' grimshaw_gp_mle((gom - u)[gom > u])
#' @export
grimshaw_gp_mle <- function(x) {
  #  Argument for function:
  #
  #  x     the sample values from the GPD
  #
  #
  #
  #  Returned from the function:
  #
  #  k       the mle of k
  #  a       the mle of a
  #
  n<-length(x)
  xq<-sort(x)
  xbar<-mean(x)
  sumx2<-sum(x^2)/n
  x1<-xq[1]
  xn<-xq[n]
  #
  #  Find the local maxima/minima of the likelihood.
  #
  #
  #  Initialize epsilon as the accuracy criterion
  #
  epsilon<-10^(-6)/xbar
  #
  #  The local maxima/minima must be found numerically by
  #  finding the zero(s) of h().
  #
  #  Algorithm for finding the zero(s) of h().
  #
  #
  #  Any roots that exist must be within the interval (lobnd,hibnd).
  #
  lobnd<-2*(x1-xbar)/x1^2
  if(lobnd>=0){
    lobnd<- -epsilon
  }
  hibnd<-(1/xn)-epsilon
  if(hibnd<=0){
    hibnd<-epsilon
  }
  #
  #  If h''(0)>0, look for one negative and one positive zero of h().
  #  If h''(0)<0, look for two negative and two positive zeros of h().
  #
  secderiv<-sumx2-2*xbar^2  #{ Evaluate h''(0). }
  if(secderiv>0){
    #
    #
    #  Look for one negative and one positive zero of h().
    #
    #
    thzeros<-cbind(c(0,0),c(0,0))
    nzeros<-2
    #
    #  Begin with the initial value at lobnd.
    #
    hlo<-(1+sum(log(1-lobnd*x))/n)*(sum(1/(1-lobnd*x))/n)-1
    if(hlo<0){
      thlo<-lobnd       #{  Orient the search so h(thlo)<0  }
      thhi<- -epsilon
    }
    else{
      thlo<- -epsilon
      thhi<-lobnd
    }
    thzero<-lobnd    #{  Initial value for modified Newton-Raphson is lobnd. }
    dxold<-abs(thhi-thlo)
    dx<-dxold
    temp1<-sum(log(1-thzero*x))/n
    temp2<-sum(1/(1-thzero*x))/n
    temp3<-sum(1/(1-thzero*x)^2)/n
    h<-(1+temp1)*(temp2)-1
    hprime<-(temp3-temp2^2-temp1*(temp2-temp3))/thzero

    #
    #  Newton-Raphson Algorithm to find the zero of the function h()
    #  for a given initial starting point.
    j<-1
    maxiter<-100  #{Maximum number of mod. Newton-Raphson iterations}
    while(j<=maxiter){
      #
      #  Determine whether it is better to use Bisection (if N-R is
      #  out of range or not decreasing fast enough) or Newton-Raphson.
      #
      c1<-(((thzero-thhi)*hprime-h)*((thzero-thlo)*hprime-h)>=0)
      c2<-(abs(2*h)>abs(dxold*hprime))
      if(c1+c2>=1){
        dxold<-dx
        dx<-(thhi-thlo)/2
        thzero<-thlo+dx
        if(thlo==thzero){       #{Change in root is negligible}
          j<-1000
        }
      }
      else{
        dxold<-dx
        dx<-h/hprime
        temp<-thzero
        thzero<-thzero-dx
        if(temp==thzero){       #{Change in root is negligible}
          j<-1001
        }
      }
      #
      #  Determine if convergence criterion is met.
      #
      if(abs(dx)<epsilon*abs(thlo+thhi)/2){
        j<-999
      }
      temp1<-sum(log(1-thzero*x))/n
      temp2<-sum(1/(1-thzero*x))/n
      temp3<-sum(1/(1-thzero*x)^2)/n
      h<-(1+temp1)*(temp2)-1
      hprime<-(temp3-temp2^2-temp1*(temp2-temp3))/thzero
      if(h<0){            #{Maintain the bracket on the root}
        thlo<-thzero
      }
      else{
        thhi<-thzero
      }
      j<-j+1
    }

    if(j>maxiter+1){
      thzeros[1,]<-cbind(thzero,j)
    }

    #
    #  Begin with the initial value at hibnd.
    #
    hlo<-(1+sum(log(1-epsilon*x))/n)*(sum(1/(1-epsilon*x))/n)-1
    if(hlo<0){
      thlo<-epsilon       #{  Orient the search so h(thlo)<0  }
      thhi<-hibnd
    }
    else{
      thlo<-hibnd
      thhi<-epsilon
    }
    thzero<-hibnd    #{  Initial value for modified Newton-Raphson is hibnd. }
    # 14/8/2018, PJN: Try a different initial value, away from hibnd
    thzero <- (hibnd + epsilon) / 2
    #
    dxold<-abs(thhi-thlo)
    dx<-dxold
    temp1<-sum(log(1-thzero*x))/n
    temp2<-sum(1/(1-thzero*x))/n
    temp3<-sum(1/(1-thzero*x)^2)/n
    h<-(1+temp1)*(temp2)-1
    hprime<-(temp3-temp2^2-temp1*(temp2-temp3))/thzero

    #
    #  Newton-Raphson Algorithm to find the zero of the function h()
    #  for a given initial starting point.
    j<-1
    maxiter<-100  #{Maximum number of mod. Newton-Raphson iterations}
    while(j<=maxiter){
      #
      #  Determine whether it is better to use Bisection (if N-R is
      #  out of range or not decreasing fast enough) or Newton-Raphson.
      #
      c1<-(((thzero-thhi)*hprime-h)*((thzero-thlo)*hprime-h)>=0)
      c2<-(abs(2*h)>abs(dxold*hprime))
      if(c1+c2>=1){
        dxold<-dx
        dx<-(thhi-thlo)/2
        thzero<-thlo+dx
        if(thlo==thzero){       #{Change in root is negligible}
          j<-1000
        }
      }
      else{
        dxold<-dx
        dx<-h/hprime
        temp<-thzero
        thzero<-thzero-dx
        if(temp==thzero){       #{Change in root is negligible}
          j<-1001
        }
      }
      #
      #  Determine if convergence criterion is met.
      #
      if(abs(dx)<epsilon*abs(thlo+thhi)/2){
        j<-999
      }
      temp1<-sum(log(1-thzero*x))/n
      temp2<-sum(1/(1-thzero*x))/n
      temp3<-sum(1/(1-thzero*x)^2)/n
      h<-(1+temp1)*(temp2)-1
      hprime<-(temp3-temp2^2-temp1*(temp2-temp3))/thzero
      if(h<0){            #{Maintain the bracket on the root}
        thlo<-thzero
      }
      else{
        thhi<-thzero
      }
      j<-j+1
    }


    if(j>maxiter+1){
      thzeros[2,]=cbind(thzero,j)
    }
  }

  else{
    #
    #
    #  Look for two negative and two positive zeros of h().
    #
    #
    thzeros<-matrix(rep(0,8),ncol=2)
    nzeros<-4
    #
    #  Begin with the initial value at lobnd.
    #
    hlo<-(1+sum(log(1-lobnd*x))/n)*(sum(1/(1-lobnd*x))/n)-1
    if(hlo<0){
      thlo<-lobnd       #{  Orient the search so h(thlo)<0  }
      thhi<- -epsilon
    }
    else{
      thlo<- -epsilon
      thhi<-lobnd
    }
    thzero<-lobnd    #{  Initial value for modified Newton-Raphson is lobnd. }
    dxold<-abs(thhi-thlo)
    dx<-dxold
    temp1<-sum(log(1-thzero*x))/n
    temp2<-sum(1/(1-thzero*x))/n
    temp3<-sum(1/(1-thzero*x)^2)/n
    h<-(1+temp1)*(temp2)-1
    hprime<-(temp3-temp2^2-temp1*(temp2-temp3))/thzero

    #
    #  Newton-Raphson Algorithm to find the zero of the function h()
    #  for a given initial starting point.
    j<-1
    maxiter<-100  #{Maximum number of mod. Newton-Raphson iterations}
    while(j<=maxiter){
      #
      #  Determine whether it is better to use Bisection (if N-R is
      #  out of range or not decreasing fast enough) or Newton-Raphson.
      #
      c1<-(((thzero-thhi)*hprime-h)*((thzero-thlo)*hprime-h)>=0)
      c2<-(abs(2*h)>abs(dxold*hprime))
      if(c1+c2>=1){
        dxold<-dx
        dx<-(thhi-thlo)/2
        thzero<-thlo+dx
        if(thlo==thzero){       #{Change in root is negligible}
          j<-1000
        }
      }
      else{
        dxold<-dx
        dx<-h/hprime
        temp<-thzero
        thzero<-thzero-dx
        if(temp==thzero){       #{Change in root is negligible}
          j<-1001
        }
      }
      #
      #  Determine if convergence criterion is met.
      #
      if(abs(dx)<epsilon*abs(thlo+thhi)/2){
        j<-999
      }
      temp1<-sum(log(1-thzero*x))/n
      temp2<-sum(1/(1-thzero*x))/n
      temp3<-sum(1/(1-thzero*x)^2)/n
      h<-(1+temp1)*(temp2)-1
      hprime<-(temp3-temp2^2-temp1*(temp2-temp3))/thzero
      if(h<0){            #{Maintain the bracket on the root}
        thlo<-thzero
      }
      else{
        thhi<-thzero
      }
      j<-j+1
    }

    if(j>maxiter+1){
      thzeros[1,]<-cbind(thzero,j)
    }
    #
    #  Look at the derivative to determine where the second root lies.
    #   If h'(0)>0, second root lies between thzero and -epsilon.
    #   If h'(0)<0, second root lies between lobnd and thzero.
    #
    temp1<-sum(log(1-thzero*x))/n
    temp2<-sum(1/(1-thzero*x))/n
    temp3<-sum(1/(1-thzero*x)^2)/n
    hprime<-(temp3-temp2^2-temp1*(temp2-temp3))/thzero
    if(hprime>0){
      #
      #  h'(0)>0, so the second zero lies between thzero and -epsilon.
      #
      #
      #  Establish Initial Values.
      #
      thlo<-thzero
      thhi<- -epsilon
      thzero<-thhi
      dx<-thlo-thhi

      j<-1
      maxiter<-100  #{Maximum number of bisection iterations}
      while(j<=maxiter){
        dx<-.5*dx
        thmid<-thzero+dx
        hmid<-(1+sum(log(1-thmid*x))/n)*(sum(1/(1-thmid*x))/n)-1
        if(hmid<0){
          thzero<-thmid
        }
        if(hmid==0){   #{Zero of h() has been found}
          j<-999
        }
        if(abs(dx)<epsilon*abs(thlo+thhi)/2){
          j<-999
        }
        j<-j+1
      }

      if(j>maxiter+1){
        thzeros[2,]<-cbind(thzero,j)
      }
    }
    else{
      #
      #  h'(0)<0, so the second zero lies between lobnd and thzero.
      #
      #
      #  Establish Initial Values.
      #
      thlo<-lobnd
      thhi<-thzero
      thzero<-thlo
      dx<-thhi-thlo

      j<-1
      maxiter<-100  #{Maximum number of bisection iterations}
      while(j<=maxiter){
        dx<-.5*dx
        thmid<-thzero+dx
        hmid<-(1+sum(log(1-thmid*x))/n)*(sum(1/(1-thmid*x))/n)-1
        if(hmid<0){
          thzero<-thmid
        }
        if(hmid==0){   #{Zero of h() has been found}
          j<-999
        }
        if(abs(dx)<epsilon*abs(thlo+thhi)/2){
          j<-999
        }
        j<-j+1
      }

      if(j>maxiter+1){
        thzeros[2,]<-cbind(thzero,j)
      }
    }
    #
    #  Begin with the initial value at hibnd.
    #
    hlo<-(1+sum(log(1-epsilon*x))/n)*(sum(1/(1-epsilon*x))/n)-1
    if(hlo<0){
      thlo<-epsilon       #{  Orient the search so h(thlo)<0  }
      thhi<-hibnd
    }
    else{
      thlo<-hibnd
      thhi<-epsilon
    }
    thzero<-hibnd    #{  Initial value for modified Newton-Raphson is hibnd. }
    # 14/8/2018, PJN: Try a different initial value, away from hibnd
    thzero <- (hibnd + epsilon) / 2
    #
    dxold<-abs(thhi-thlo)
    dx<-dxold
    temp1<-sum(log(1-thzero*x))/n
    temp2<-sum(1/(1-thzero*x))/n
    temp3<-sum(1/(1-thzero*x)^2)/n
    h<-(1+temp1)*(temp2)-1
    hprime<-(temp3-temp2^2-temp1*(temp2-temp3))/thzero

    #
    #  Newton-Raphson Algorithm to find the zero of the function h()
    #  for a given initial starting point.
    j<-1
    maxiter<-100  #{Maximum number of mod. Newton-Raphson iterations}
    while(j<=maxiter){
      #
      #  Determine whether it is better to use Bisection (if N-R is
      #  out of range or not decreasing fast enough) or Newton-Raphson.
      #
      c1<-(((thzero-thhi)*hprime-h)*((thzero-thlo)*hprime-h)>=0)
      c2<-(abs(2*h)>abs(dxold*hprime))
      if(c1+c2>=1){
        dxold<-dx
        dx<-(thhi-thlo)/2
        thzero<-thlo+dx
        if(thlo==thzero){       #{Change in root is negligible}
          j<-1000
        }
      }
      else{
        dxold<-dx
        dx<-h/hprime
        temp<-thzero
        thzero<-thzero-dx
        if(temp==thzero){       #{Change in root is negligible}
          j<-1001
        }
      }
      #
      #  Determine if convergence criterion is met.
      #
      if(abs(dx)<epsilon*abs(thlo+thhi)/2){
        j<-999
      }
      temp1<-sum(log(1-thzero*x))/n
      temp2<-sum(1/(1-thzero*x))/n
      temp3<-sum(1/(1-thzero*x)^2)/n
      h<-(1+temp1)*(temp2)-1
      hprime<-(temp3-temp2^2-temp1*(temp2-temp3))/thzero
      if(h<0){            #{Maintain the bracket on the root}
        thlo<-thzero
      }
      else{
        thhi<-thzero
      }
      j<-j+1
    }

    if(j>maxiter+1){
      thzeros[3,]<-cbind(thzero,j)
    }
    #
    #  Look at the derivative to determine where the second root lies.
    #   If h'(0)>0, second root lies between thzero and hibnd.
    #   If h'(0)<0, second root lies between epsilon and thzero.
    #
    temp1<-sum(log(1-thzero*x))/n
    temp2<-sum(1/(1-thzero*x))/n
    temp3<-sum(1/(1-thzero*x)^2)/n
    hprime<-(temp3-temp2^2-temp1*(temp2-temp3))/thzero
    if(hprime>0){
      #
      #  h'(0)>0, so the second zero lies between thzero and hibnd.
      #
      #
      #  Establish Initial Values.
      #
      thlo<-thzero
      thhi<-hibnd
      thzero<-thhi
      dx<-thlo-thhi

      j<-1
      maxiter<-100  #{Maximum number of bisection iterations}
      while(j<=maxiter){
        dx<-.5*dx
        thmid<-thzero+dx
        hmid<-(1+sum(log(1-thmid*x))/n)*(sum(1/(1-thmid*x))/n)-1
        if(hmid<0){
          thzero<-thmid
        }
        if(hmid==0){   #{Zero of h() has been found}
          j<-999
        }
        if(abs(dx)<epsilon*abs(thlo+thhi)/2){
          j<-999
        }
        j<-j+1
      }

      if(j>maxiter+1){
        thzeros[4,]<-cbind(thzero,j)
      }
    }
    else{
      #
      #  h'(0)<0, so the second zero lies between epsilon and thzero.
      #
      #
      #  Establish Initial Values.
      #
      thlo<-epsilon
      thhi<-thzero
      thzero<-thlo
      dx<-thhi-thlo

      j<-1
      maxiter<-100  #{Maximum number of bisection iterations}
      while(j<=maxiter){
        dx<-.5*dx
        thmid<-thzero+dx
        hmid<-(1+sum(log(1-thmid*x))/n)*(sum(1/(1-thmid*x))/n)-1
        if(hmid<0){
          thzero<-thmid
        }
        if(hmid==0){   #{Zero of h() has been found}
          j<-999
        }
        if(abs(dx)<epsilon*abs(thlo+thhi)/2){
          j<-999
        }
        j<-j+1
      }

      if(j>maxiter+1){
        thzeros[4,]<-cbind(thzero,j)
      }
    }
  }
  #
  #  Of the candidate zero(s) of h(), determine whether they correspond
  #  to a local maximum or minimum of the log-likelihood.
  #
  #  Eliminate any non-convergent roots}
  thetas<-thzeros[thzeros[,2]>maxiter+1,]
  nzeros<-nrow(thetas)
  proll<-rep(0,nzeros)
  mles<-matrix(rep(0,4*nzeros),ncol=4)
  i<-1
  while(i<=nzeros){
    temp1<-sum(log(1-thetas[i,1]*x))
    mles[i,1]<- -temp1/n
    mles[i,2]<-mles[i,1]/thetas[i,1]
    mles[i,3]<- -n*log(mles[i,2])+(1/mles[i,1]-1)*temp1
    mles[i,4]<-999
    i<-i+1
  }
  ind <- seq_along(mles[, 4])
  ind<-ind[mles[,4]==999]
  if(sum(ind)==0){   #{ Check to see if there are any local maxima. }
    nomle<-0          #{ If not, there is no mle. }
  }
  else{
    nomle<-1
  }
  if(nomle!=0){
    mles<-mles[ind,]
    # 6/7/2018: Start of codeprovided by Leo  Belzile, to fix a bug
    # Discard values of k that are greater than 1
    outside <- which(mles[,1]> 1+1e-10)
    if (length(outside) > 0){ #at most 2 such values
        mles[outside,3] <- -Inf #replace by hard bound
    }
    # End of code provided by Leo
    nmles<-nrow(mles)
    #
    #  Add the boundary value where k=1 to the candidates for the
    #  maximum of the log-likelihood.
    #
    mles<-rbind(mles,c(1,xn,-n*log(xn),999))
    nmles<-nmles+1
    #
    #  Choose of the candidate mles whichever has the largest
    #  log-likelihood.
    #
    maxlogl<-max(mles[,3])
    ind<-order(mles[,3])
    ind<-ind[nmles]
    k<-mles[ind,1]
    #  label(k,'GPD mle of k')
    a<-mles[ind,2]
    #  label(a,'GPD mle of a')
  }
  else{
    #
    #  No Maximum Likelihood Estimators were found.
    #
    k<-NA
    a<-NA
  }
  return(list(k=k,a=a))
}

###############################################################################

gev_mle <- function(init, ...){
    x <- stats::optim(init, gev_loglik, ..., control = list(fnscale = -1),
               hessian = FALSE)
    temp <- list()
    temp$mle <- x$par
    temp$nllh <- -x$value
    temp$hessian <- try(stats::optimHess(temp$mle, gev_loglik, ...),
                        silent = TRUE)
    if (inherits(temp$hessian, "try-error")) {
      temp$cov <- temp$se <- NULL
      return(temp)
    }
    temp$cov <- try(-solve(temp$hessian), silent = TRUE)
    if (inherits(temp$cov, "try-error")) {
      temp$cov <- temp$se <- NULL
    } else {
      temp$se <- sqrt(diag(temp$cov))
    }
    return(temp)
}

gev_fish <- function(theta) {
  #
  # Calculates the Fisher information for a single observation from a
  # GEV(mu, sigma, xi) distribution, evaluated at xi = 0
  #
  # Args:
  #   theta      : vector (mu, sigma, xi)
  #
  # Returns: a 3 by 3 matrix
  #
  sigma <- theta[2]
  gam <- 0.5772157
  zet <- 1.2020569
  i_mm <- 1 / sigma^2
  i_ss <- (pi ^ 2 / 6 - (1 - gam) ^ 2 ) / sigma^2
  i_xx <- pi ^ 2 / 6 - pi ^ 2 * gam / 2 + gam ^ 2 - gam ^ 3 - 2 * zet +
          2 * gam * zet + pi ^ 2 * gam ^ 2 / 4 + gam ^ 4 / 4 + 3 * pi ^ 4 / 80
  i_ms <- (gam - 1) / sigma^2
  i_mx <- (pi ^ 2 / 6 + gam ^ 2 - 2 * gam) / (2 * sigma)
  i_sx <- (4 * gam + 4 * zet + pi ^ 2 * gam + 2 * gam ^ 3 - pi ^ 2
          - 6 * gam ^ 2) / (4 * sigma)
  col1 <- c(i_mm, i_ms, i_mx)
  col2 <- c(i_ms, i_ss, i_sx)
  col3 <- c(i_mx, i_sx, i_xx)
  return(matrix(c(col1, col2, col3), 3, 3))
}

#---------------------------------#
#    PWM		      	    #
#---------------------------------#
# Amended function from fExtremes #
#---------------------------------#

gev_pwm <- function(x){
  y <- function(x, w0, w1, w2) {
      (3 ^ x - 1) / (2 ^ x - 1) - (3 * w2 - w0) / (2 * w1 - w0)
  }
  n <- length(x)
  nmom <- 3
  moments <- rep(0, nmom)
  pwm_ests <- function(x) {
    x <- rev(sort(x))
    moments[1] <- mean(x)
    for (i in 1:n) {
        weight <- 1 / n
        for (j in 2:nmom) {
            weight <- weight * (n - i - j + 2) / (n - j + 1)
            moments[j] <- moments[j] + weight * x[i]
        }
    }
    w0 <- moments[1]
    w1 <- moments[2]
    w2 <- moments[3]
    xi <- stats::uniroot(f = y, interval = c(-5, +5), w0 = w0, w1 = w1,
                         w2 = w2)$root
    beta <- (2 * w1 - w0) * xi / gamma(1 - xi) / (2 ^ xi - 1)
    mu <- w0 + beta * (1 - gamma(1 - xi)) / xi
    c(mu, beta, xi)
  }
  pwm <- pwm_ests(x)
  xi <- pwm[3]
#
# Table 1 from Hosking, Wallis and Wood (1985)
# .. but with k replaced by xi, where xi=-k
# xi w1l w12 wl3 w22 w23 w33
  r1 <- c(0.4, 1.6637, 1.3355, 1.1405, 1.8461, 1.1628, 2.9092)
  r2 <- c(0.3, 1.4153, .8912, .5640, 1.2574, .4442, 1.4090)
  r3 <- c(0.2, 1.3322, .6727, .3926, 1.0013, .2697, .9139)
  r4 <- c(0.1, 1.2915, .5104, .3245, .8440, .2240, .6815)
  r5 <- c(0.0, 1.2686, .3704, .2992, .7390, .2247, .5633)
  r6 <- c(-0.1, 1.2551, .2411, .2966, .6708, .2447, .5103)
  r7 <- c(-0.2, 1.2474, .1177, .3081, .6330, .2728, .5021)
  r8 <- c(-0.3, 1.2438, -.0023, .3297, .6223, .3033, .5294)
  r9 <- c(-0.4, 1.2433, -.1205, .3592, .6368, .3329, .5880)
  r_mat <- rbind(r1, r2, r3, r4, r5, r6, r7, r8, r9)
#
  pwm.varcov <- se <- NA
  if (xi > 0.5) {
     warning("Asymptotic Standard Errors not available for PWM when xi>0.5.")
  }
  if (xi >= -0.4 & xi <= 0.4){
    w_row_1 <- which.min(abs(xi - r_mat[, 1]))
    w_1 <- xi - r_mat[w_row_1, 1]
    w_row_2 <- w_row_1 - sign(w_1)
    w_2 <- xi - r_mat[w_row_2, 1]
    wij <- 10 * ((0.1 - abs(w_1)) * r_mat[w_row_1, -1] +
                   (0.1 - abs(w_2)) * r_mat[w_row_2, -1])
    pwm_varcov <- matrix(wij[c(1, 2, 3, 2, 4, 5, 3, 5, 6)], 3, 3) / n
    se <- sqrt(diag(pwm_varcov))
  }
  if (xi < -0.4){
    n_sim <- 100
    my_pars <- matrix(NA, ncol = 3, nrow = n_sim)
    for (i in 1:n_sim){
      x <- rgev(n = n, loc = pwm[1], scale = pwm[2], shape = pwm[3])
      my_pars[i, ] <- pwm_ests(x)
    }
    pwm_varcov <- stats::var(my_pars)
    se <- sqrt(diag(pwm_varcov))
  }
  return(list(est = pwm, se = se, cov = pwm_varcov))
}

###############################################################################

gev_lrs <- function(x){
  n <- length(x)                          # sample size
  x <- sort(x)                            # put data in ascending order
  q0 <- 1 / (n + 1)
  q2 <- n / (n+1)                         # for sample minimum and maximum
  a <- sqrt(log(q2) / log(q0))
  q1 <- q0 ^ a                            # `middle' quantile
  n0 <- 1
  n1 <- round((n + 1) * q1)
  n2 <- n                                 # corresponding order statistics
  ns <- c(n0, n1, n2)
  qs <- c(q0, q1, q2)
  xs <- x[ns]
  r_hat <- (xs[3] - xs[2]) / (xs[2] - xs[1])
  xi_hat <- -log(r_hat) / log(a)
  fq0 <- ((-log(q0)) ^ (-xi_hat) - 1) / xi_hat
  fq2 <- ((-log(q2)) ^ (-xi_hat) - 1) / xi_hat
  sigma_hat <- (xs[3] - xs[1]) / (fq2 - fq0)
  mu_hat <- xs[1] - sigma_hat * fq0
  return(c(mu_hat, sigma_hat, xi_hat))
}

###############################################################################

os_mle <- function(init, gumbel = FALSE, ...){
  x <- stats::optim(init, os_loglik, ..., control = list(fnscale = -1),
                    hessian = FALSE, gumbel = gumbel)
  temp <- list()
  temp$mle <- x$par
  if (gumbel) {
    temp$mle <- c(temp$mle, 0)
  }
  temp$nllh <- -x$value
  temp$hessian <- try(stats::optimHess(temp$mle, os_loglik, gumbel = gumbel,
                                       ...), silent = TRUE)
  if (inherits(temp$hessian, "try-error")) {
    temp$cov <- temp$se <- NULL
    return(temp)
  }
  temp$cov <- try(-solve(temp$hessian), silent = TRUE)
  if (inherits(temp$cov, "try-error")) {
    temp$cov <- temp$se <- NULL
  } else {
    temp$se <- sqrt(diag(temp$cov))
  }
  return(temp)
}
