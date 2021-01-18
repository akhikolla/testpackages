# Copyright 2015-2019 Venelin Mitov
#
# This file is part of POUMM.
#
# POUMM is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# POUMM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with POUMM  If not, see <http://www.gnu.org/licenses/>.

#' @rdname OU
#' @title Distribution of an Ornstein-Uhlenbeck Process at Time \eqn{t}, Given 
#'   Initial State at Time \eqn{0}
#'   
#' @description An Ornstein-Uhlenbeck (OU) process represents a continuous time 
#'   Markov chain parameterized by an initial state \eqn{x_0}, selection 
#'   strength \eqn{\alpha>0}, long-term mean \eqn{\theta}, and time-unit 
#'   variance \eqn{\sigma^2}. Given \eqn{x_0}, at time \eqn{t}, the state of the
#'   process is characterized by a normal distribution with mean \eqn{x_0 
#'   exp(-\alpha t) + \theta (1 - exp(-\alpha t))} and variance \eqn{\sigma^2 
#'   (1-exp(-2 \alpha t)) / (2 \alpha)}. In the limit \eqn{\alpha -> 0}, the OU 
#'   process converges to a Brownian motion process with initial state \eqn{x_0}
#'   and time-unit variance \eqn{\sigma^2} (at time \eqn{t}, this process is 
#'   characterized by a normal distribution with mean \eqn{x_0} and variance 
#'   \eqn{t \sigma^2}.
#'   
#' @param n Integer, the number of values to sample.
#' @param z Numeric value or vector of size n.
#' @param z0 Numeric value or vector of size n, initial value(s) to condition 
#'   on.
#' @param t Numeric value or vector of size n, denoting the time-step.
#' @param alpha,theta,sigma Numeric values or n-vectors, parameters of the OU 
#'   process; alpha and sigma must be non-negative. A zero alpha is interpreted 
#'   as the Brownian motion process in the limit alpha -> 0.
#' @param log Logical indicating whether the returned density should is on the logarithmic scale.
#'  
#' @details Similar to dnorm and rnorm, the functions described in this
#'   help-page support single values as well as vectors for the parameters z,
#'   z0, t, alpha, theta and sigma.
#' @name OU
NULL

#' @describeIn OU probability density
#' @return dOU returns the conditional probability density(ies) of the elements 
#'   in z, given the initial state(s) z0, time-step(s) t and OU-parameters by
#'   alpha, theta and sigma.
#' @export
dOU <- function(z, z0, t, alpha, theta, sigma, log = TRUE) {
  ett <- exp(-alpha * t)
  sd <- sigma * sqrt(t)
  a <- alpha > 0
  sd[a] <- sigma[a] * sqrt((1 - ett[a]^2) / (2 * alpha[a]))
  mean <- z0 * ett + theta * (1 - ett)
  nan <- is.infinite(t) & alpha == 0
  mean[nan] <- z0[nan]
  
  dnorm(z, mean = mean, sd = sd, log = log)
}

#' @describeIn OU random generator
#'   
#' @return rOU returns a numeric vector of length n, a random sample from the
#'   conditional distribution(s) of one or n OU process(es) given initial
#'   value(s) and time-step(s).
#' @examples 
#' z0 <- 8
#' t <- 10
#' n <- 100000
#' sample <- rOU(n, z0, t, 2, 3, 1)
#' dens <- dOU(sample, z0, t, 2, 3, 1)
#' var(sample)  # around 1/4
#' varOU(t, 2, 1) 
#' 
#' @export
rOU <- function(n, z0, t, alpha, theta, sigma) {
  ett <- exp(-alpha * t)
  sd <- sigma * sqrt(t)
  a <- alpha > 0
  
  sd[a] <- sigma[a] * sqrt((1 - ett[a]^2) / (2 * alpha[a]))
  mean <- z0 * ett + theta * (1 - ett)
  nan <- is.infinite(t) & alpha == 0
  mean[nan] <- z0[nan]
  
  rnorm(n, mean = mean, sd = sd)
}

#' @describeIn OU mean value
#'   
#' @return meanOU returns the expected value of the OU-process at time t.
#' @export
meanOU <- function(z0, t, alpha, theta) {
  ett <- exp(-alpha * t)
  mean <- z0 * ett + theta * (1 - ett)
  nan <- is.infinite(t) & alpha==0
  mean[nan] <- z0[nan]
  mean
}


#' @describeIn OU variance
#'   
#' @return varOU returns the expected variance of the OU-process at time t.
#' @export
varOU <- function(t, alpha, sigma) {
  ett <- exp(-alpha * t)
  var <- sigma^2 * t
  a <- alpha > 0
  var[a] <- sigma[a]^2 * (1 - ett[a]^2) / (2 * alpha[a])
  var
}

#' @describeIn OU standard deviation
#'   
#' @return sdOU returns the standard deviation of the OU-process at time t.
#' @export
sdOU <- function(t, alpha, sigma) {
  ett <- exp(-alpha * t)
  sd <- sigma * sqrt(t)
  a <- alpha > 0
  sd[a] <- sigma[a] * sqrt((1 - ett[a]^2) / (2 * alpha[a]))
  sd
}


#' Generation of a random trajectory of an OU process starting from a given
#' initial state (only for test purpose)
#' @inheritParams rTrajectoryOU
#' 
#' @details Used for test purpose only. This is an internal function and is
#'   appropriate for small time-steps only.
rTrajectoryOUDef <- function(z0, t, alpha, theta, sigma, steps = 1) {
  if(length(t) != steps) {
    t <- rep(t, length.out = steps)
  }
  wien <- sigma * rnorm(steps, 0, sqrt(t))
  z <- numeric(steps + 1)
  z[1] <- z0
  for(i in 1:steps)
    z[i+1] <- z[i] + alpha * (theta - z[i]) * t[i] + wien[i]
  z[-1]
}

#' Generation of a random trajectory of an OU process starting from a given
#' initial state
#' 
#' @description Generates a trajectory xt given initial state z0 according to an
#'   Ornstein-Uhlenbeck process.
#'   
#' @param z0 Numeric value, initial state.
#' @param t Numeric value or vector of size steps, denoting the time-step(s).
#' @param alpha,theta,sigma Numeric values, parameters of the OU process.

#' @param steps Integer, number of steps.
#' @return A numeric vector of length steps containing the generated values
#'   at times 0+t, 0+2t, ..., 0+steps*t.
#'   
#' @examples 
#' z0 <- 0
#' nSteps <- 100
#' t <- 0.01
#' trajectory <- rTrajectoryOU(z0, t, 2, 2, 1, steps = nSteps)
#' plot(trajectory, type = 'l')
#' 
#' @export
rTrajectoryOU <- function(z0, t, alpha, theta, sigma, steps = 1) {
  if(length(t) != steps) {
    t <- rep(t, length.out = steps)
  }
  z <- numeric(steps + 1)
  z[1] <- z0
  ett <- exp(-alpha * t)
  sd <- sigma * sqrt(t)
  a <- alpha > 0
  sd[a] <- sigma[a] * sqrt((1 - ett[a]^2) / (2 * alpha[a]))
  deltas <- rnorm(steps, mean = theta * (1 - ett), sd = sd)
  
  for(i in 1:steps) {
    z[i+1] <- z[i] * ett[i] + deltas[i]
  }
  
  z[-1]
}


#' @rdname PhylogeneticH2
#' @title Phylogenetic Heritability
#'   
#' @description The phylogenetic heritability, \eqn{H^2}, is defined as the 
#'   ratio of the genetic variance over the total phenotypic variance expected 
#'   at a given evolutionary time t (measured from the root of the tree). Thus,
#'   the phylogenetic heritability connects the parameters alpha, sigma and
#'   sigmae of the POUMM model through a set of equations. The functions
#'   described here provide an R-implementation of these equations.
#'   
#' @param t Numeric value denoting evolutionary time (i.e. distance from the 
#'   root of a phylogenetic tree).
#' @param alpha,sigma Numeric values or n-vectors, parameters of the OU process;
#'   alpha and sigma must be non-negative. A zero alpha is interpreted as the 
#'   Brownian motion process in the limit alpha -> 0.
#' @param sigmae Numeric, environmental phenotypic deviation at the tips.
#' @param H2 Phylogenetic heritability at time t.
#' @param z Numerical vector of observed phenotypes.
#' @param tree A phylo object.
#' @param tFrom,tTo Numerical minimal and maximal root-tip distance to limit the
#'   calculation.
#'   
#' @return All functions return numerical values or NA, in case of invalid 
#'   parameters
#' @examples 
#' # At POUMM stationary state (equilibrium, t=Inf)
#' H2 <- H2(alpha = 0.75, sigma = 1, sigmae = 1, t = Inf)     # 0.4
#' alpha <- alpha(H2 = H2, sigma = 1, sigmae = 1, t = Inf)    # 0.75
#' sigma <- sigmaOU(H2 = H2, alpha = 0.75, sigmae = 1, t = Inf) # 1
#' sigmae <- sigmae(H2 = H2, alpha = 0.75, sigma = 1, t = Inf) # 1
#' 
#' # At finite time t = 0.2
#' H2 <- H2(alpha = 0.75, sigma = 1, sigmae = 1, t = 0.2)     # 0.1473309
#' alpha <- alpha(H2 = H2, sigma = 1, sigmae = 1, t = 0.2)    # 0.75
#' sigma <- sigmaOU(H2 = H2, alpha = 0.75, sigmae = 1, t = 0.2) # 1
#' sigmae <- sigmae(H2  =  H2, alpha = 0.75, sigma = 1, t = 0.2) # 1
#' 
#'    
#' @name PhylogeneticH2
#' @seealso OU
NULL

#' @describeIn PhylogeneticH2 Calculate alpha given time t, H2, sigma and sigmae
#' 
#' @importFrom lamW lambertW0
#'
#' @export
alpha <- function(H2, sigma, sigmae, t = Inf) {
  # if(!all(H2>=0 & H2 <=1 & sigma >= 0 & sigmae >=0 & t>=0)) {
  #   warning(paste0("Function `alpha()` was called on invalid parameters. ", 
  #                  "Check that H2>=0 & H2 <=1 & sigma >= 0 & sigmae >=0 & t>=0.",
  #                  "Parameter values were: ", 
  #                  toString(c(H2=H2, sigma=sigma, sigmae=sigmae, t=t))))
  #   NA
  # }
  if(is.infinite(t)) {
    ifelse(H2 == 1,
      # Assume correctly defined PMM, i.e. Brownian motion with sigma > 0 and 
      # environmental deviation with sigmae >= 0
      0,
      # H2 is the phylogenetic heritability at equilibrium
      sigma^2 * (1 - H2) / (2 * H2 * sigmae^2)    
    )
  } else {
    # finite t
    y <- sigma^2 / sigmae^2 * (1 / H2 - 1)
    ifelse(is.na(y) | is.infinite(y),
      as.double(NA),
      (t * y + lambertW0((-exp(-t * y)) * t * y)) / (2 * t)
    )
  }
}

#' @describeIn PhylogeneticH2 Calculate sigma given time t, H2 at time t, alpha
#'   and sigmae
#'   
#' @note This function is called sigmaOU and not simply sigma to avoid a conflict 
#' with a function sigma in the base R-package.
#'   
#' @export
sigmaOU <- function(H2, alpha, sigmae, t=Inf) {
  res <- if(is.infinite(t)) {
    ifelse(alpha == 0,
      # BM
      sqrt(sigmae^2 * H2 / (t * (1 - H2))),
      # alpha>0, OU in equilibrium
      sqrt(2 * alpha * H2 * sigmae^2 / (1 - H2))
    )
  } else {
    # t is finite
    ifelse(alpha == 0,
      # BM
      sqrt(sigmae^2 * H2 / (t * (1 - H2))),
      # alpha>0, OU not in equilibrium
      sqrt(2 * alpha * H2 * sigmae^2 / ((1 - exp(-2 * alpha * t)) * (1 - H2)))
    ) 
  }
  names(res) <- NULL
  res
}

#' @describeIn PhylogeneticH2 Calculate sigmae given alpha, sigma, and H2 at
#'   time t
#' @details The function sigmae uses the formula H2 = varOU(t, alpha, sigma) /
#'   (varOU(t, alpha, sigma) + sigmae^2)
#'   
#' @export
sigmae <- function(H2, alpha, sigma, t = Inf) {
  sigmaG2 <- varOU(t, alpha, sigma)
  sqrt(sigmaG2 / H2 - sigmaG2)
}

#' @describeIn PhylogeneticH2 "Empirical" phylogenetic heritability estimated
#'   from the empirical variance of the observed phenotypes and sigmae
#'   
#' @importFrom stats var
#' @export
H2e <- function(z, sigmae, tree=NULL, tFrom=0, tTo=Inf) {
  if(!is.null(tree)) {
    tipTimes <- nodeTimes(tree)[1:length(tree$tip.label)]
    z <- z[which(tipTimes >= tFrom & tipTimes <= tTo)]
  }
  res <- 1 - (sigmae^2 / var(z))
  names(res) <- NULL
  res
}

#' Phylogenetic heritability estimated at time t
#' @param alpha,sigma numeric, parameters of the OU process acting on the
#'   genetic contributions
#' @param sigmae numeric, environmental standard deviation
#' @param t time from the beginning of the process at which heritability should
#'   be calculated, i.e. epidemiologic time
#' @param tm average time for within host evolution from getting infected until
#'   getting measured or passing on the infection to another host
#' @export
H2 <- function(alpha, sigma, sigmae, t = Inf, tm = 0) {
  lenoutput <- max(sapply(list(alpha, sigma, sigmae, t, tm), length))
  if(lenoutput>1) {
    alpha <- rep(alpha, lenoutput / length(alpha))
    sigma <- rep(sigma, lenoutput / length(sigma))
    sigmae <- rep(sigmae, lenoutput / length(sigmae))
    t <- rep(t, lenoutput / length(t))
    tm <- rep(tm, lenoutput / length(tm))
  }
  sigmag2 <- ifelse(alpha > 0, 
                    0.5*(1 - exp(-2 * alpha * t)) / alpha * sigma^2, 
                    t * sigma^2)  
  sigmagm2 <- ifelse(alpha > 0, 
                     0.5*(1 - exp(-2 * alpha * tm)) / alpha * sigma^2, 
                     tm * sigma^2)
  
  values <- ifelse(t > tm, 
                   (sigmag2 - sigmagm2) / (sigmag2 + sigmae^2), 
                   rep(0, length(t)))
  
  values <- ifelse(alpha == 0 & sigma > 0 & is.infinite(t) & !is.infinite(sigmae),
                   1, values)
  values
}

#' Expected covariance of two tips at given root-tip time and phylogenetic distance
#' 
#' @param alpha,sigma,sigmae POUMM parameters
#' @param t A non-negative number or vector time from the beginning of the POUMM
#'   process (root-tip distance). If a vector, the evaluation is done on each 
#'   couple (row) from cbind(t, tau).
#' @param tau A non-negative number or vector indicating the phylogenetic 
#'   distance between two tips, each of them located at time t from the root. 
#'   If a vector, the evaluation is done on each couple (row) from cbind(t, tau).
#' @param tanc A non-negative number or vector indication the root-mrca distance
#' for a couple of tips. Defaults to t-tau/2 corresponding to an ultrametric tree.
#' @param corr Logical indicating whether correlation should be returned instead
#' of covariance.
#' @param as.matrix Logical indicating if a variance-covariance matrix should be
#' returned.
#' 
#' @details The function assumes that the two tips are at equal distance t from 
#' the root. This implies that the root-tip distance of their mrca is t - tau/2.
#' 
#' @return If as.matrix == FALSE, a number. Otherwise a two by two symmetric
#'   matrix. If t or tau is a vector of length bigger than 1, then a vector of 
#'   numbers or a list of matrices.
#' 
#' @export
covPOUMM <- function(alpha, sigma, sigmae, t, tau, tanc = t - tau/2, 
                     corr = FALSE, as.matrix = FALSE) {
  if(length(t) == 1 & length(tau) == 1 & length(tanc) == 1) {
    covMat <- covVTipsGivenTreePOUMM(
      tree = NULL, alpha = alpha, sigma = sigma, sigmae = sigmae,
      tanc = rbind(c(t, tanc), c(tanc, t)), 
      tauij = rbind(c(0, tau), c(tau, 0)), corr = corr)
    if(as.matrix) {
      covMat
    } else {
      covMat[1,2]
    }
  } else {
    ttau <- cbind(t, tau, tanc)
    apply(ttau, 1, function(.) {
      t <- .[1]
      tau <- .[2]
      tanc <- .[3]
      covMat <- covVTipsGivenTreePOUMM(
        tree = NULL, alpha = alpha, sigma = sigma, sigmae = sigmae,
        tanc = rbind(c(t, tanc), c(tanc, t)), 
        tauij = rbind(c(0, tau), c(tau, 0)), corr = corr)
      if(as.matrix) {
        covMat
      } else {
        covMat[1,2]
      }
    })
  }
}

#' Variance covariance matrix of the values at the tips of a tree under an OU 
#' process
#' 
#' @param tree A phylo object.
#' @param alpha,sigma Non-negative numeric values, parameters of the OU process.
#' @param sigmae Non-negative numeric value, environmental standard deviation at
#'   the tips.
#' @param corr Logical indicating if a correlation matrix shall be returned.
#' @param tanc Numerical matrix with the time-distance from the root of the tree
#'   to the mrca of each tip-couple. If NULL it will be calculated.
#' @param tauij Numerical matrix with the time (patristic) distance between each
#'   pair of tips. If NULL, it will be calculated.
#' @return a variance covariance or a correlation matrix of the tips in tree.
#' @references (Hansen 1997) Stabilizing selection and the comparative analysis 
#'   of adaptation.
#'   
#' @export
covVTipsGivenTreePOUMM <- function(
  tree, alpha = 0, sigma = 1, sigmae = 0, tanc = NULL, tauij = NULL, corr = FALSE) {
  
  # distances from the root to the mrca's of each tip-couple.
  if(is.null(tanc)) {
    tanc <- ape::vcv(tree)
  } else {
    tanc <- as.matrix(tanc)
  }
  
  N <- dim(tanc)[1]
  if(alpha > 0) {
    varanc <- 0.5 * (1 - exp(-2 * alpha * tanc)) / alpha
  } else {
    # limit case alpha==0
    varanc <- tanc
  }
  
  # time from the root to each tip
  ttips <- diag(tanc) 
  
  # distance-matrix between tips
  if(is.null(tauij)) {
    tauij <- sapply(ttips, function(.) . + ttips) - 2 * tanc  
  } else {
    tauij <- as.matrix(tauij)
  }
  
  covij <- exp(-alpha * tauij) * varanc * sigma^2 + diag(x = sigmae^2, nrow = N)
  if(corr) {
    sdi <- sqrt(diag(covij))
    sdi_sdj <- sapply(sdi, function(.) . * sdi)
    covij / sdi_sdj
  } else {
    covij
  }  
}

#' Distribution of the genotypic values under a POUMM fit
#'
#' @param tree an object of class phylo
#' @param z A numeric vector of size length(tree$tip.label) representing the trait
#'     values at the tip-nodes.
#' @param g0 A numeric value at the root of the tree, genotypic contribution.
#' @param alpha,theta,sigma Numeric values, parameters of the OU model.
#' @param sigmae Numeric non-negative value (default 0). Specifies the standard
#'   deviation of random environmental contribution (white noise) included in z.
#' @return a list with elements V.g, V.g_1, mu.g, V.e, V.e_1, mu.e, V.g.poumm, mu.g.poumm.
gPOUMM <- function(z, tree, g0, alpha, theta, sigma, sigmae) {
  N <- length(tree$tip.label)
  
  mu.g <- meanOU(g0, t = nodeTimes(tree, tipsOnly = TRUE), alpha = alpha, theta = theta)

  V.e <- diag(sigmae^2, nrow=N, ncol=N)
  V.e_1 <- V.e
  diag(V.e_1) <- 1/diag(V.e)
  mu.e <- z
  
  #V.g <- cov.poumm(tree, alpha, sigma, sigmae) ###???? shouldn't sigmae be 0 here?
  V.g <- covVTipsGivenTreePOUMM(tree, alpha, sigma)
  if(sigma > 0) {
    if(sigmae > 0) {
      V.g_1 <- chol2inv(chol(V.g))
      V.g.poumm <- try(chol2inv(chol(V.g_1 + V.e_1)), silent = TRUE)
      mu.g.poumm <- V.g.poumm %*% (V.g_1 %*% mu.g + V.e_1 %*% mu.e)  
    } else {
      V.g_1 <- chol2inv(chol(V.g))
      V.g.poumm <- try(chol2inv(chol(V.g_1)), silent = TRUE)
      mu.g.poumm <- z
      warning("For sigmae=0, gPOUMM returns z.")
    }
    
  } else {
    # sigma = 0
    warning("V.g.poumm is not defined for sigma == 0.")
    V.g.poumm <- V.g_1 <- matrix(as.double(NA), N, N)
    mu.g.poumm <- mu.g
  }

  if(inherits(V.g.poumm, 'try-error')) {
    warning(V.g.poumm)
  } 
  list(V.g = V.g, V.g_1 = V.g_1, mu.g = mu.g, V.e = V.e, V.e_1 = V.e_1,
       mu.e = mu.e, V.g.poumm = V.g.poumm, mu.g.poumm = mu.g.poumm)
}

r.squared <- function(obj) {
  fittedValues <- fitted(obj)
  grandMean <- mean(obj$pruneInfo$z)
  sum((fittedValues-grandMean)^2)/sum((obj$pruneInfo$z-grandMean)^2)
}

#' @importFrom stats nobs logLik
adj.r.squared <- function(obj, r.sq = r.squared(obj)) {
  1 - (1 - r.sq) * ((nobs(obj) - 1)/(nobs(obj) - attr(logLik(obj), 'df')))
}