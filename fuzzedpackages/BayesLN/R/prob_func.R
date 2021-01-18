#' @importFrom graphics abline plot
#' @importFrom stats dgamma lm qnorm quantile var integrate uniroot rnorm sd
#'
NULL


#' @title SMNG and logSMNG Distributions
#'
#' @description Density function, distribution function, quantile function and random generator for the SMNG distribution and the logSMNG.
#' It requires the specification of a five prameters vector: \code{mu}, \code{delta}, \code{gamma}, \code{lambda} and
#' \code{beta}.
#' @param x,q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Sample size.
#' @param mu Location parameter, default set to 0.
#' @param delta Concentration parameter, must be positive.
#' @param gamma Tail parameter, must be positive.
#' @param lambda Shape parameter.
#' @param beta Skewness parameter, default set to 0 (symmetric case).
#' @param rel_tol Level of relative tolerance required for the \code{integrate} procedure or for the infinite sum convergence check.
#' Default set to \code{1e-5}.
#' @param inf_sum Logical: if FALSE (default) the integral representation of the SMNG density is used,
#' otherwise the infinite sum is employed.
#'
#' @details
#' The SMNG distribution is a normal scale-mean mixture distribution with a GIG as mixing distribution. The
#' density can be expressed as an infinite sum of Bessel K functions and it is characterized by 5 parameters.
#'
#' Moreover, if X is SMNG distributed, then \eqn{Z=exp(X)} is distributed as a log-SMNG distribution.
#'
#' @return
#' \code{dSMNG} and \code{dlSMNG} provide the values of the density function at a quantile \code{x} for, respectively
#'  a SMNG distribution and a log-SMNG.
#'
#'\code{pSMNG} and \code{plSMNG} provide the cumulative distribution function at a quantile \code{q}.
#'
#' \code{qSMNG} and \code{qlSMNG} provide the quantile corresponding to a probability level \code{p}.
#'
#' \code{rSMNG} and \code{rlSMNG} generate \code{n} independent samples from the desired distribution.
#'
#'
#' @examples
#'\donttest{
#' ### Plots of density and cumulative functions of the SMNG distribution
#' x<-seq(-10,10,length.out = 500)
#'
#' plot(x,dSMNG(x = x,mu = 0,delta = 1,gamma = 1,lambda = 1,beta= 2),
#'     type="l",ylab="f(x)")
#' lines(x,dSMNG(x = x,mu = 0,delta = 1,gamma = 1,lambda = 1,beta= -2),col=2)
#' title("SMNG density function")
#'
#' plot(x,pSMNG(q = x,mu = 0,delta = 1,gamma = 1,lambda = 1,beta= 2),
#'     type="l",ylab="F(x)")
#' lines(x,pSMNG(q = x,mu = 0,delta = 1,gamma = 1,lambda = 1,beta= -2),col=2)
#' title("SMNG cumulative function")
#'
#'
#' ### Plots of density and cumulative functions of the logSMNG distribution
#' x<-seq(0,20,length.out = 500)
#'
#' plot(x,dlSMNG(x = x,mu = 0,delta = 1,gamma = 1,lambda = 2,beta = 1),
#'     type="l",ylab="f(x)",ylim = c(0,1.5))
#' lines(x,dlSMNG(x = x,mu = 0,delta = 1,gamma = 1,lambda = 2,beta = -1),col=2)
#' title("logSMNG density function")
#'
#' plot(x,plSMNG(q = x,mu = 0,delta = 1,gamma = 1,lambda = 2,beta = 1),
#'     type="l",ylab="F(x)",ylim = c(0,1))
#' lines(x,plSMNG(q = x,mu = 0,delta = 1,gamma = 1,lambda = 2,beta = -1),col=2)
#' title("logSMNG cumulative function")
#'}
#'
#'@name SMNGdistribution

NULL


#' Integrand SMNG density function
#'
#' Function that implements the integral part that appears in the SMNG density function.
#'
#' @keywords internal

integral <- function(t, a, b, c, l) {
  t ^ (-2 * l) * exp(-0.5 * (a * t ^ 2 + b / t ^ 2 - c * t))
}


#' Addend SMNG density function
#'
#' Function that implements the j-th addend of the SMNG density function expressed as an infinite sum.
#'
#' @keywords internal



add <- function(j, delta_bar, r, lambda, delta, gamma, beta, mu){
  expBes <- exp(delta * gamma - delta_bar * gamma)
  expBes * besselK(x = delta_bar * gamma,
    nu = (-lambda + 0.5) + (j - 1) / 2, expon.scaled = T) * ((sqrt(gamma) * r) / sqrt(delta_bar)) ^ (j - 1) *
    (1 / factorial(j - 1)) / (besselK(x = delta * gamma, nu = lambda, expon.scaled = T) * sqrt(2 * pi)) *
    (gamma / delta_bar) ^ (-lambda + 0.5)
}


#' Summation SMNG density function
#'
#' Function that adds terms to the infinite sum of the SMNG density function until the relative
#' error reaches the requested threshold fixed by \code{rel_tol}.
#'
#' @keywords internal


inf_sum<-function(x,lambda, delta, gamma, beta, mu, rel_tol){
  delta_bar <- sqrt((x - mu) ^ 2 + delta ^ 2)
  r <- beta * (x - mu)
  par_sum <- add(j = 1, delta_bar = delta_bar,  r = r, lambda = lambda, delta = delta, gamma = gamma,
      beta = beta, mu = mu)
  j <- 2
  rel_err <- par_sum
  while (abs(rel_err) > rel_tol) {
    n_add <-  add(j = j, delta_bar = delta_bar, r = r, lambda = lambda, delta = delta, gamma = gamma,
        beta = beta, mu = mu)
    par_sum <- par_sum + n_add
    rel_err <- n_add / par_sum
    j <- j + 1
  }
  par_sum
}


#' @name dSMNG
#' @rdname SMNGdistribution
#' @export

dSMNG <- function(x, mu = 0, delta, gamma, lambda, beta = 0, inf_sum = FALSE, rel_tol = 1e-5) {
    if (delta <= 0)
      stop("delta must be positive")
    if (gamma <= 0)
      stop("gamma must be positive")
    if (inf_sum == FALSE) {
      k <- ((gamma / delta) ^ lambda * exp(-beta ^ 2 / 2)) / (besselK(x = delta * gamma, nu = lambda) *
                                                            sqrt(2 * pi))
      x_cen <- (x - mu)
      int <- sapply(X = x_cen,
                    function(X) integrate(integral, lower = 0, upper = Inf, a = (X ^ 2 + delta ^ 2),
                                          b = gamma ^ 2, c = 2 * beta * X, l = lambda, rel.tol = rel_tol)$value)
      dens <- k * int
    } else{
      k <- ((gamma / delta) ^ lambda * exp(-beta ^ 2 / 2))
      summ <- sapply(X = x,
                     function(X) inf_sum(x = X, lambda = lambda, delta = delta, gamma = gamma,  beta = beta,
                                         mu = mu, rel_tol = rel_tol))
      dens <- k * summ
    }
    return(dens)
}

#' @name pSMNG
#' @rdname SMNGdistribution
#' @export

pSMNG<-function(q, mu, delta, gamma, lambda, beta, rel_tol = 1e-5){
  p <- sapply(X = q, function(X) integrate(dSMNG, lower = -Inf, upper = X, lambda = lambda,
                                           delta = delta, gamma = gamma, beta = beta, mu = mu,
                                           inf_sum = F, rel_tol = rel_tol, rel.tol = rel_tol)$value)
  p
}


#' Function for finding SMNG quantiles
#'
#' Non-linear function whose zero can be fuond by the \code{uniroot} procedure in order to obtain the SMNG quantile
#' corresponding to probability \code{p}.
#'
#' @keywords internal




zerfun <- function(x, mu, delta, gamma, lambda, beta, p, rel_tol) {
  f <- pSMNG(x, mu = mu, delta = delta, gamma = gamma, lambda = lambda, beta = beta, rel_tol = rel_tol) - p
  return(f)
}

#' @name qSMNG
#' @rdname SMNGdistribution
#' @export


qSMNG <- function(p, mu, delta, gamma, lambda, beta, rel_tol= 1e-5){
  interval <- c(mu - 4 * delta, mu + 4 * delta)
  q <- uniroot(f = zerfun, interval = interval, mu = mu, delta = delta, gamma = gamma, l = lambda,
               beta = beta, p = p, rel_tol = rel_tol, extendInt = "yes", tol = rel_tol)$root
  return(q)
}

#' @name rSMNG
#' @rdname SMNGdistribution
#' @export

rSMNG <- function(n, mu, delta, gamma, lambda, beta) {
  if (delta <= 0)
    stop("delta must be positive")
  if (gamma <= 0)
    stop("gamma must be positive")
  sigma2 <- ghyp::rgig(n = n, lambda = lambda, chi = delta ^ 2, psi = gamma ^ 2)
  s <- rnorm(n = n, mean = mu + beta * sqrt(sigma2), sd = sqrt(sigma2))
  return(s)
}




#' @name dlSMNG
#' @rdname SMNGdistribution
#' @export

dlSMNG <- function(x, mu = 0, delta, gamma, lambda, beta, inf_sum = FALSE, rel_tol = 1e-5) {
    dl <- numeric(length(x))
    for (i in 1:length(x)) {
      dl[i] <- ifelse(x[i] <= 0, no = (1 / x[i]) * dSMNG(log(x[i]), mu = mu, delta = delta,
                       gamma = gamma, lambda = lambda, beta = beta, inf_sum = inf_sum, rel_tol = rel_tol),
                       yes = 0)
    }
    return(dl)
  }



#'@name plSMNG
#' @rdname SMNGdistribution
#' @export

plSMNG<-function(q, mu, delta, gamma, lambda, beta, rel_tol = 1e-5){
  p <- sapply(X = q, function(X) integrate(dlSMNG, lower = 0,upper = X,
                                           lambda = lambda, delta = delta,
                                           gamma = gamma, beta = beta, mu = mu,inf_sum = F,
                                           rel_tol = rel_tol, rel.tol = rel_tol)$value)
  p
}

#' Function for finding logSMNG quantiles
#'
#' Non-linear function whose zero can be fuond by the \code{uniroot} procedure in order to obtain the logSMNG quantile
#' corresponding to probability \code{p}.
#'
#' @keywords internal

zerfun_log <- function(x, mu, delta, gamma, lambda, beta, p, rel_tol){
  f <- plSMNG(x, mu = mu, delta = delta, gamma = gamma, lambda = lambda, beta = beta, rel_tol = rel_tol) - p
  return(f)
}

#'@name qlSMNG
#' @rdname SMNGdistribution
#' @export


qlSMNG<-function(p, mu, delta, gamma, lambda, beta, rel_tol = 1e-5){
  interval <- c(0, 40)
  q<-uniroot(f = zerfun_log, interval = interval, mu = mu, delta = delta, gamma = gamma, l = lambda, beta = beta,
             p = p, rel_tol = rel_tol, extendInt = "yes", tol = rel_tol)$root
  return(q)
}


#'@name rlSMNG
#' @rdname SMNGdistribution
#' @export

rlSMNG <- function(n, mu, delta, gamma, lambda, beta) {
  s <- exp(rSMNG(n = n, mu = mu, delta = delta, gamma = gamma, lambda = lambda, beta = beta))
  return(s)
}


#' @title SMNG Moments and Moment Generating Function
#'
#' @description Functions that implement the mean, the generic moments (both raw and centered) and the moment generating function of the SMNG distribution.
#'
#' @param j Order of the moment.
#' @param r Coefficient of the MGF. Can be viewed also as the order of the logSMNG moments.
#' @param mu Location parameter, default set to 0.
#' @param delta Concentration parameter, must be positive.
#' @param gamma Tail parameter, must be positive.
#' @param lambda Shape parameter.
#' @param beta Skewness parameter, default set to 0 (symmetric case).
#' @param rel_tol Level of relative tolerance required for the \code{integrate} procedure or for the infinite sum.
#' Default set to \code{1e-5}.
#' @param inf_sum Logical: if FALSE (default), the integral representation of the SMNG density is used,
#' otherwise the infinite sum is employed.
#' @param type String that indicate the kind of moment to comupute. Could be \code{"central"} (default) or \code{"raw"}.
#'
#'@details
#'If the mean (i.e. the first order raw moment) of the SMNG distribution is required, then the function \code{meanSMNG} could be use.
#'
#'On the other hand, to obtain the generic \emph{j}-th moment both \code{"raw"} or \code{"centered"} around the mean, the function \code{momentSMNG} could be used.
#'
#'Finally, to evaluate the Moment Generating Function (MGF) of the SMNG distribution in the point \code{r}, the function \code{SMNG_MGF} is provided.
#'It is defined only for points that are lower then the parameter \code{gamma}, and for integer values of \code{r} it could also considered as the
#'\emph{r}-th raw moment of the logSMNG distribution. The last function is implemented both in the integral form, which uses the routine \code{\link{integrate}},
#'or in the infinite sum structure.
#'
#'@examples
#'
#' ### Comparisons sample quantities vs true values
#' sample <- rSMNG(n=1000000,mu = 0,delta = 2,gamma = 2,lambda = 1,beta = 2)
#' mean(sample)
#' meanSMNG(mu = 0,delta = 2,gamma = 2,lambda = 1,beta = 2)
#'
#' var(sample)
#' SMNGmoment(j = 2,mu = 0,delta = 2,gamma = 2,lambda = 1,beta = 2,type = "central")
#' SMNGmoment(j = 2,mu = 0,delta = 2,gamma = 2,lambda = 1,beta = 2,type = "raw")-
#'                         meanSMNG(mu = 0,delta = 2,gamma = 2,lambda = 1,beta = 2)^2
#'
#' mean(exp(sample))
#' SMNG_MGF(r = 1,mu = 0,delta = 2,gamma = 2,lambda = 1,beta = 2)
#'
#'
#' @name SMNGmoments
NULL


#' Addend SMNG moment generating function
#'
#' Function that implements the j-th addend of the SMNG moment generating function expressed as an infinite sum.
#'
#'
#' @keywords internal



add_MGF <- function(j, r = r, lambda, delta, gamma, beta) {
  (((beta * r) ^ j) / factorial(j)) * ((delta / sqrt(gamma ^ 2 - r ^ 2)) ^ (j / 2)) *
    besselK(x = delta * sqrt(gamma ^ 2 - r ^ 2), nu = lambda + j / 2) / (besselK(x = delta * gamma, nu = lambda))
}

#' Summation SMNG moment generating function
#'
#' Function that adds terms to the infinite sum of the SMNG moment generating function until the relative
#' error reaches the requested threshold fixed by \code{rel_tol}.
#'
#' @keywords internal

inf_sum_MGF<-function(r,lambda, delta, gamma, beta,rel_tol){
  par_sum <- add_MGF(j = 0,r = r,lambda =lambda, delta = delta, gamma = gamma, beta = beta)
  j <- 1
  rel_err <- par_sum
  while(abs(rel_err) > rel_tol){
    n_add <- add_MGF(j = j, r = r, lambda = lambda, delta = delta, gamma = gamma, beta = beta)
    par_sum <- par_sum + n_add
    rel_err <- n_add / par_sum
    j <- j + 1
  }
  return(par_sum)
}

#' Integrand SMNG moment generating function
#'
#' Function that implements the integral part that appears in the SMNG moment generating function.
#'
#' @keywords internal


integral_MGF <- function(t, a, b, c, l){
  t ^ ( - 2 * l - 1) * exp(- 0.5 * (a * t ^ 2 + b / t ^ 2 - c / t))
}

#'@name SMNG_MGF
#' @rdname SMNGmoments
#' @export


SMNG_MGF <- function(r, mu = 0, delta, gamma, lambda, beta = 0, inf_sum = FALSE, rel_tol = 1e-5) {
  if (delta <= 0)
    stop("delta must be positive")
  if (gamma <= 0)
    stop("gamma must be positive")
  if(gamma <= r)
    stop("The gamma parameter must be greater than the order r")
  if(inf_sum == FALSE){
    k <- exp(mu * r) * (((gamma / delta) ^ lambda) / (besselK(x = delta * gamma,nu = lambda)))
    b <- gamma ^ 2 - r ^ 2
    a <- delta ^ 2
    c <- 2 * beta * r
    int <- integrate(integral_MGF, lower = 0,upper = Inf, a = a, b = b, c = c, l = lambda, rel.tol = rel_tol)$value
    MGF <- k * int
  } else {
    k <- exp(mu * r) * ((gamma / sqrt(gamma ^ 2 - r ^ 2)) ^ lambda)
    summ <- inf_sum_MGF(r = r, lambda = lambda, delta = delta, gamma = gamma, beta = beta, rel_tol = rel_tol)
    MGF <- k * summ
  }
  return(MGF)
}


#' Ratio of Bessel K functions
#'
#' Ratio of Bessel K functions with equal argument. The difference between the two orders is \code{nu_diff}.
#'
#' @keywords internal

RatioBesselK <- function (x, nu, nu_diff)
{
  besselK(x, nu + nu_diff, expon.scaled = TRUE) / besselK(x, nu, expon.scaled = TRUE)
}


#' @name meanSMNG
#' @rdname SMNGmoments
#' @export

meanSMNG <- function(mu, delta, gamma, lambda, beta){
  if (delta <= 0)
    stop("delta must be positive")
  if (gamma <= 0)
    stop("gamma must be positive")
  mean <- mu + beta * sqrt(delta / gamma) * RatioBesselK(x = delta * gamma, nu = lambda, nu_diff = 0.5)
  return(mean)
}





#' SMNG moments centered in \code{mu}
#'
#' Function that implements the formula of the SMNG moments ceneterd with respect to the location parameter \code{mu}.
#' It is used to compute all the moments through the function \code{momentRecursion}
#'
#'
#' @keywords internal

SMNGZmoment <- function(j, mu, d, g, l, beta) {
  if(j == 0) {
    mom <- 1
  } else if(j %% 2 == 0){
    mom <- (2 * (d / g)) ^ (j / 2) * RatioBesselK(x = d * g, nu = l, nu_diff = j / 2) * (gamma(j / 2 + 0.5) / sqrt(pi)) *
      Re(fAsianOptions::kummerM(x = - beta ^ 2 / 2, a = - j / 2, b = 0.5))
  }else{
    mom <- beta * (2 * (d / g)) ^ (j / 2) * RatioBesselK(x = d * g, nu = l, nu_diff = j / 2) * ((sqrt(2) * gamma( j / 2 + 1)) / sqrt(pi))*
      Re(fAsianOptions::kummerM(x = - beta ^ 2 / 2, a = 0.5 - j / 2, b = 1.5))
  }
  return(mom)
}

#' Recursion used for SMNG moments
#'
#' Recursive forumla that allowas to obtain both the raw and the centered moments starting form the
#' moments centered with respect to \code{mu}.
#'
#'
#' @keywords internal

momentRecursion <- function(center, j, mu, d, g, l, beta) {
  t <- seq(0, j, by = 1)
  add <- numeric(length(t))
  for(i in 1 : length(t)){
    add[i] <- choose(n = j, k = t[i]) * (mu - center) ^ (j - t[i])*
      SMNGZmoment(j = t[i], mu = mu, d = d, g = g, l = l, beta = beta)
  }
  cen <- sum(add)
  return(cen)
}


#' @name SMNGmoment
#' @rdname SMNGmoments
#' @export

SMNGmoment <- function(j, mu, delta, gamma, lambda, beta, type = "central") {
  if (delta <= 0)
    stop("delta must be positive")
  if (gamma <= 0)
    stop("gamma must be positive")
  if (type == "central") {
    mean <- meanSMNG(mu = mu, delta = delta, gamma = gamma, lambda = lambda, beta = beta)
    mom <- momentRecursion(center = mean, j = j, mu = mu, d = delta, g = gamma, l = lambda, beta = beta)
  } else if (type == "raw") {
    mom <- momentRecursion(center = 0, j = j, mu = mu, d = delta, g = gamma, l = lambda, beta = beta)
  } else{
    stop("type of the moment must be 'central' or 'raw'")
  }
  return(mom)
}



#'GH Moment Generating Function
#'
#' Function that implements the moment generating function of the Generalized Hyperbolyc (GH) distribution.
#'
#' @param r Coefficient of the MGF. Can be viewd also as the order of the log-GH moments.
#' @param mu Location parameter, default set to 0.
#' @param delta Concentration parameter, must be positive.
#' @param alpha Tail parameter, must be positive and greater than the modulus of \code{beta}.
#' @param lambda Shape parameter.
#' @param beta Skewness parameter, default set to 0 (symmetric case).
#'
#' @details This function allows to evaluate the moment generating function of the GH distribution in the point \code{r}.
#'It is defined only for points that are lower than the value of \eqn{\gamma}, that is defined as:
#'\eqn{\gamma^2=\alpha^2-\beta^2.}
#'For integer values of \code{r}, it could also be considered as the
#'\emph{r}-th raw moment of the log-GH distribution.
#'
#'
#'
#'
#'
#' @export


GH_MGF <- function(r, mu = 0, delta, alpha, lambda, beta = 0) {
  if (delta <= 0)
    stop("delta must be positive")
  if (alpha <= 0)
    stop("alpha must be positive")
  gamma <- sqrt(alpha ^ 2 - beta ^ 2)
  if(gamma <= 0)
    stop("alpha must be greater than beta")
  if(gamma <= r)
    stop("The value of gamma must be greater than the order r")
  gamma_bar <- alpha ^ 2 - (beta + r) ^ 2
  expBes <- exp(delta * gamma - delta * sqrt(gamma_bar))
   MGF <- exp(mu * r) * (gamma ^ 2 / gamma_bar) ^ (lambda / 2) *
     expBes * besselK(x = delta * sqrt(gamma_bar), nu = lambda, expon.scaled = T) /
     besselK(x = delta * gamma, nu = lambda, expon.scaled = T)
  return(MGF)
}


