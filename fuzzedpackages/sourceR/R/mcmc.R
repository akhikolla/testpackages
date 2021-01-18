##################################
# Name: mcmc.R                   #
# Created: 2016-06-17            #
# Author: Chris Jewell           #
# Modified by: Poppy Miller      #
# Copyright: Chris Jewell 2016   #
# Purpose: MCMC Updaters for ref #
#          class DAG model.      #
##################################

AdaptiveSingleSiteMRW <- R6::R6Class(
  "AdaptiveSingleSiteMRW",
  public = list(
    toupdate = NA,
    naccept = NA,
    ncalls = NA,
    acceptbatch = NA,
    batchsize = NA,
    node = NA,
    tune = NA,
    initialize = function(node, tune = 0.1,
                          batchsize = 50) {
      self$naccept <- 0
      self$ncalls <- 0
      self$acceptbatch <- 0
      self$batchsize <- batchsize
      self$tune <- tune
      self$node <- node
    },
    update = function() {
        self$ncalls <- self$ncalls + 1
        old_data <- self$node$getData()
        picur <- self$node$logPosterior()

        # Propose using MHRW
        self$node$data <- old_data * exp( rnorm(1, 0, self$tune) )
        pican <- self$node$logPosterior()

        alpha <- pican - picur + log(self$node$data/old_data)
        if (is.finite(alpha) &
            log(runif(1)) < alpha) {
          self$acceptbatch <- self$acceptbatch + 1
        }
        else {
          self$node$data <- old_data
        }
      private$adapt()
    },
    acceptance = function() {
      self$naccept / self$ncalls
    }
  ),
  private = list(
    adapt =  function() {
      if ((self$ncalls > 0) &
          (self$ncalls %% self$batchsize == 0)) {
        m <- ifelse(self$acceptbatch / self$batchsize > 0.44, 1,-1)
        self$tune <-
          exp(log(self$tune) + m * min(0.05, 1.0 / sqrt(self$ncalls)))
        if(!is.finite(self$tune)) self$tune <- 1e-8
        self$naccept <-
          self$naccept + self$acceptbatch
        self$acceptbatch <- 0
      }
    }
  )
)


#' AdaptiveMultiMRW
#'
#' This class implements an Adaptive Multi-site Metropolis random walk
#' algorithm.
#'
#' A multivariate Gaussian proposal is used, for which the proposal variance
#' is a scaled version of the evolving empirical posterior covariance matrix.
#' See Roberts and Rosenthal (2012) Examples of Adaptive MCMC. \emph{Journal of Computational
#' and Graphical Statistics}. \bold{18}:349--367.
#'
#' Please note that no checks are performed as to the suitability of this
#' algorithm for a particular \link{StochasticNode}.  It is up to the user
#' to use the correct update algorithm for the appropriate nodes.
#'
#' @docType class
#' @name AdaptiveMultiMRW
#' @importFrom R6 R6Class
#' @keywords DAG node MCMC
#' @return Object of \code{\link{AdaptiveMultiMRW}}
#' @format Object of \code{\link{R6Class}} with methods for updating a \link{Node} instance.
#' @field cov the current covariance
#' @field burnin the number of updates to burn in
#' @field tune the current tuning matrix
#' @field naccept the number of accepted proposals
#' @field ncalls the number of times \code{update} has been called
#' @field node the node to which the updater is attached
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(node, tune = rep(0.1, length(node$getData())), burning = 100)}}{constructor takes an instance of a \link{StochasticNode}
#'   node, initial tuning vector (diagonal of adaptive tuning matrix), and number of burnin calls.}
#'   \item{\code{update()}}{when called, updates \code{node}}
#'   \item{\code{acceptance()}}{return the acceptance rate}
#' }
AdaptiveMultiMRW <- R6::R6Class(
  "AdaptiveMultiMRW",
  public = list(
    cov = NA,
    burnin = NA,
    tune = NA,
    naccept = NA,
    ncalls = NA,
    node = NA,
    initialize = function(node, tune = rep(0.1, length(node$getData())), burnin = 100)
    {
      self$node <- node
      self$naccept <- 0
      self$ncalls <- 0
      self$burnin <- burnin
      self$cov <- makeOnlineCov(length(self$node$data))
      updateOnlineCov(self$cov, self$node$data)
      self$tune <- diag(tune)
    },
    update = function()
    {
      self$ncalls <- self$ncalls + 1
      pCov <- NULL
      if ((runif(1) < 0.95) & (self$ncalls > self$burnin)) {
        pCov <- getOnlineCov(self$cov)
      }
      else pCov <- self$tune

      old_data <- self$node$data
      picur <- self$node$logPosterior()
      self$node$data <- mvtnorm::rmvnorm(1, mean = self$node$data, sigma = 5.6644*pCov/length(self$node$data))
      pican <- self$node$logPosterior()
      if(log(runif(1)) < pican - picur) {
        self$naccept <- self$naccept + 1
      }
      else {
        self$node$data <- old_data
      }
    },
    acceptance = function() {
      self$naccept / self$ncalls
    }
  )
)



#' AdaptiveDirMRW
#'
#' This class implements an Adaptive Multi-site Metropolis random walk
#' algorithm, constrained so the parameter vector sums to 1.
#'
#' An adaptive multivariate Gaussian proposal is used for $d-1$ elements of a $d$-dimensional parameter
#' vector contained in \code{node}, with the $d$th element updated to ensure that the vector sums to 1.
#' This makes the updater useful for Dirichlet distributed random variables.
#'
#' For details of the adaptive scheme, see Roberts and Rosenthal (2012) Examples of Adaptive MCMC. \emph{Journal of Computational
#' and Graphical Statistics}. \bold{18}:349--367.
#'
#' Please note that no checks are performed as to the suitability of this
#' algorithm for a particular \link{StochasticNode}.  It is up to the user
#' to use the correct update algorithm for the appropriate nodes.
#'
#' @docType class
#' @name AdaptiveDirMRW
#' @importFrom R6 R6Class
#' @export
#' @keywords DAG node MCMC
#' @return Object of \code{\link{AdaptiveDirMRW}}
#' @format Object of \code{\link{R6Class}} with methods for updating a \link{DirichletNode} instance.
#' @field cov the current covariance
#' @field burnin the number of updates to burn in
#' @field tune the current tuning matrix
#' @field naccept the number of accepted proposals
#' @field ncalls the number of times \code{update} has been called
#' @field node the node to which the updater is attached
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(node, toupdate = function() 1:length(node$getData()), tune = rep(0.1, length(node$getData())), burning = 100)}}{constructor takes an instance of a \link{StochasticNode}
#'   node, function to choose the indices of the elements to update (by default all elements), initial tuning vector (diagonal of adaptive tuning matrix), and the number of calls between
#'   adaptations.}
#'   \item{\code{update()}}{when called, updates \code{node}}
#'   \item{\code{acceptance()}}{return the acceptance rate}
#' }
AdaptiveDirMRW <- R6::R6Class(
  "AdaptiveDirMRW",
  public = list(
    toupdate = NA,
    naccept = NA,
    ncalls = NA,
    acceptbatch = NA,
    batchsize = NA,
    node = NA,
    tune = NA,
    name = NA,
    initialize = function(node, toupdate = function() 1:length(node$getData()),
                          tune = rep(0.1, length(node$getData())),
                          batchsize = 50, name = "name") {
      self$toupdate = toupdate
      self$naccept <- rep(0, length(node$getData()))
      self$ncalls <- 0
      self$acceptbatch <-
        rep(0, length(node$getData()))
      self$batchsize <- batchsize
      self$tune <- tune
      self$node <- node
      self$name <- name
    },
    update = function() {
      self$ncalls <- self$ncalls + 1
      for (j in self$toupdate()) {
        old_data <- self$node$getData()
        picur <- self$node$logPosterior()

        # Propose using MHRW
        if(runif(1) < 0.95) {
        self$node$data[j] <-
          rnorm(1, old_data[j], self$tune[j])
        }
        else {
          self$node$data[j] <-
            rnorm(1, old_data[j], 0.01)
        }
        if (any(self$node$data < 0)) {
          self$node$data <- old_data
          next
        }

        self$node$data <-
          self$node$getData() / sum(self$node$getData())
        pican <- self$node$logPosterior()

        alpha <- pican - picur
        if (is.finite(alpha) &
            log(runif(1)) < alpha) {
          self$acceptbatch[j] <- self$acceptbatch[j] + 1
        }
        else {
          self$node$data <- old_data
        }
      }
      if (self$ncalls %% 50 == 0) private$adapt()
    },
    acceptance = function() {
      self$naccept / self$ncalls
    }
  ),
  private = list(
    adapt =  function() {
      if ((self$ncalls > 0) &
          (self$ncalls %% self$batchsize == 0)) {
        m <- ifelse(self$acceptbatch / self$batchsize > 0.44, 1,-1)
        self$tune <-
          exp(log(self$tune) + m * min(0.05, 1.0 / sqrt(self$ncalls)))
        self$naccept <-
          self$naccept + self$acceptbatch
        self$acceptbatch <- self$acceptbatch * 0
      }
    }
  )
)



#' AdaptiveLogDirMRW
#'
#' This class implements an Adaptive Multi-site logarithmic Metropolis-Hastings random walk
#' algorithm, constrained so the parameter vector sums to 1.
#'
#' An adaptive multivariate log-Gaussian proposal is used for $d-1$ elements of a $d$-dimensional parameter
#' vector contained in \code{node}, with the $d$th element updated to ensure that the vector sums to 1.
#' This makes the updater useful for Dirichlet distributed random variables, improving on \link{AdaptiveDirMRW} by
#' ensuring proposals do not go negative.
#'
#' For details of the adaptive scheme, see Roberts and Rosenthal (2012) Examples of Adaptive MCMC. \emph{Journal of Computational
#' and Graphical Statistics}. \bold{18}:349--367.
#'
#' Please note that no checks are performed as to the suitability of this
#' algorithm for a particular \link{StochasticNode}.  It is up to the user
#' to use the correct update algorithm for the appropriate nodes.
#'
#' @docType class
#' @name AdaptiveLogDirMRW
#' @importFrom R6 R6Class
#' @export
#' @keywords DAG node MCMC
#' @return Object of \code{\link{AdaptiveLogDirMRW}}
#' @format Object of \code{\link{R6Class}} with methods for updating a \link{DirichletNode} instance.
#' @field cov the current covariance
#' @field burnin the number of updates to burn in
#' @field tune the current tuning matrix
#' @field naccept the number of accepted proposals
#' @field ncalls the number of times \code{update} has been called
#' @field node the node to which the updater is attached
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(node, toupdate = function() 1:length(node$getData()), tune = rep(0.1, length(node$getData())), burning = 100)}}{constructor takes an instance of a \link{StochasticNode}
#'   node, function to choose the indices of the elements to update (by default all elements), initial tuning vector (diagonal of adaptive tuning matrix), and number of calls between
#'   adaptations.}
#'   \item{\code{update()}}{when called, updates \code{node}}
#'   \item{\code{acceptance()}}{return the acceptance rate}
#' }
AdaptiveLogDirMRW <- R6::R6Class(
  "AdaptiveDirMRW",
  public = list(
    toupdate = NA,
    naccept = NA,
    ncalls = NA,
    acceptbatch = NA,
    batchsize = NA,
    node = NA,
    tune = NA,
    name = NA,
    initialize = function(node, toupdate = function() 1:length(node$getData()),
                          tune = rep(0.1, length(node$getData())),
                          batchsize = 50, name = "name") {
      self$toupdate = toupdate
      self$naccept <- rep(0, length(node$getData()))
      self$ncalls <- rep(0, length(node$getData()))
      self$acceptbatch <-
        rep(0, length(node$getData()))
      self$batchsize <- batchsize
      self$tune <- tune
      self$node <- node
      self$name <- name
    },
    update = function() {
      updIdx <- self$toupdate()
      for (j in updIdx) {
        self$ncalls[j] <- self$ncalls[j] + 1
        old_data <- self$node$getData()
        picur <- self$node$logPosterior()

        # Propose using MHRW
        logq <- 0
        if(runif(1) < 0.95) {
          self$node$data[j] <- old_data[j] * exp( rnorm(1, 0, self$tune[j]) )
          logq <- log(self$node$data[j] / old_data[j])
        }
        else {
          self$node$data[j] <- rnorm(1, old_data[j], 0.1)
        }

        # Trap to avoid negative numbers
        if(any(self$node$data < 0)) {
          self$node$data <- old_data
          next
        }

        self$node$data <-
          self$node$getData() / sum(self$node$getData())
        pican <- self$node$logPosterior()

        alpha <- pican - picur + logq

        if (is.finite(alpha) &
            log(runif(1)) < alpha) {
          self$acceptbatch[j] <- self$acceptbatch[j] + 1
        }
        else {
          self$node$data <- old_data
        }
        private$adapt(j)
      }
    },
    acceptance = function() {
      self$naccept <- self$naccept + self$acceptbatch
      self$naccept / self$ncalls
    }
  ),
  private = list(
    adapt =  function(j) {
      if ((self$ncalls[j] > 0) &
          (self$ncalls[j] %% self$batchsize == 0)) {
        m <- ifelse(self$acceptbatch[j] / self$batchsize > 0.44, 1,-1)
        self$tune[j] <-
          exp(log(self$tune[j]) + m * min(0.05, 1.0 / sqrt(self$ncalls[j])))
        if(!is.finite(self$tune[j])) {
          warning('AdaptiveLogDirMRW tuning variance non-finite')
          self$tune[j] <- 1e-9
        }
        self$naccept[j] <-
          self$naccept[j] + self$acceptbatch[j]
        self$acceptbatch[j] <- 0
      }
    }
  )
)



AdaptiveLogDirMRW2 <- R6::R6Class(
  "AdaptiveDirMRW2",
  public = list(
    toupdate = NA,
    naccept = NA,
    ncalls = NA,
    acceptbatch = NA,
    batchsize = NA,
    node = NA,
    tune = NA,
    initialize = function(node, toupdate = function() 1:length(node$data),
                          tune = rep(0.1, length(node$data)),
                          batchsize = 50) {
      self$toupdate = toupdate
      self$naccept <- rep(0, length(node$data))
      self$ncalls <- rep(0, length(node$data))
      self$acceptbatch <-
        rep(0, length(node$data))
      self$batchsize <- batchsize
      self$tune <- tune
      self$node <- node
    },
    update = function() {
      updIdx <- self$toupdate()
      for (j in updIdx) {
        self$ncalls[j] <- self$ncalls[j] + 1
        old_data <- self$node$data
        picur <- self$node$logPosterior()

        # Propose using MHRW
        self$node$data[j] <- old_data[j] * exp( rnorm(1, 0, self$tune[j]) )
        pican <- self$node$logPosterior()

        alpha <- pican - picur + log(self$node$data[j]/old_data[j])
        if (is.finite(alpha) &
            log(runif(1)) < alpha) {
          self$acceptbatch[j] <- self$acceptbatch[j] + 1
        }
        else {
          self$node$data <- old_data
        }
        private$adapt(j)
      }
    },
    acceptance = function() {
      self$naccept <- self$naccept + self$acceptbatch
      self$naccept / self$ncalls
    }
  ),
  private = list(
    adapt =  function(j) {
      if ((self$ncalls[j] > 0) &
          (self$ncalls[j] %% self$batchsize == 0)) {
        m <- ifelse(self$acceptbatch[j] / self$batchsize > 0.44, 1,-1)
        self$tune[j] <-
          exp(log(self$tune[j]) + m * min(0.05, 1.0 / sqrt(self$ncalls[j])))
        self$naccept[j] <-
          self$naccept[j] + self$acceptbatch[j]
        self$acceptbatch[j] <- 0
      }
    }
  )
)


#' PoisGammaDPUpdate
#'
#' This class implements a marginal Gibbs sampler for a Dirichlet process
#' prior on the mean of a Poisson distributed random variable, with a Gamma-distributed
#' base function.
#'
#' The marginal Gibbs sampler is based on the description in Gelman et al. Bayesian Data
#' Analysis, 3rd edition, Chapter 23.
#'
#' Please note that no checks are performed as to the suitability of this
#' algorithm for a particular \link{StochasticNode}.  It is up to the user
#' to use the correct update algorithm for the appropriate nodes.
#'
#' @docType class
#' @name PoisGammaDPUpdate
#' @importFrom R6 R6Class
#' @export
#' @keywords DAG node Dirichlet-process
#' @return Object of \code{\link{PoisGammaDPUpdate}}
#' @format Object of \code{\link{R6Class}} with methods for updating a \link{DirichletProcessNode} instance.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(node)}}{constructor takes an instance of a \link{DirichletProcessNode} \code{node}}
#'   \item{\code{update()}}{when called, updates \code{node}}
#' }
PoisGammaDPUpdate <- R6::R6Class(
  "DirichletProcessUpdate",
  public = list(
    node = NA,

    # P0 params
    alpha = NULL,
    beta = NULL,

    initialize = function(node) {
      self$node <- node
      self$alpha <- node$baseShape
      self$beta <- node$baseRate
    },
    update =
      function() {

        # Get y data
        y <- sapply(self$node$children, function(child)
          child$getData()) %>% rowSums
        aX <- sapply(self$node$children, function(child)
          child$parents$offset$getData()) %>% rowSums

        # Update s
        nk <- table(self$node$s) # Expensive!
        for (i in 1:length(self$node$s)) {
          nk_local <- as.numeric(nk)
          names(nk_local) <- names(nk)
          nk_local[self$node$s[i]] <- nk_local[self$node$s[i]] - 1

          theta <- self$node$theta$find(names(nk))

          logpi <- log(nk_local) + y[i] * log(theta) - theta * aX[i]
          logpiplus <- log(self$node$conc) +
            log(self$beta ^ self$alpha / gamma(self$alpha)) +
            lgamma(self$alpha + y[i]) -
            (self$alpha + y[i]) * log(self$beta + aX[i])

          logpi <- c(logpi, logpiplus)
          logpi[is.nan(logpi)] <- -Inf
          logpi <- logpi - max(logpi)

          pi <- exp(logpi)

          s_old <- self$node$s[i]
          idx <- sample(length(pi),1,prob = pi)

          # Check if we're adding another cluster
          # n.b. could use has.key(model$s[i], model$theta)
          #      but almost certainly slower.
          if (idx == length(pi)) {
            self$node$s[i] <- dequeue(self$node$idBucket)
            self$node$theta$insert(self$node$s[i], rgamma(1, self$alpha + y[i], self$beta + aX[i]))
            nk[self$node$s[i]] <- 0 # Add new class, increment later
          }
          else {
            self$node$s[i] <- names(pi)[idx]
          }

          # Update table of s counts, remove a class from theta if necessary
          nk[s_old] <-
            nk[s_old] - 1
          if ((nk[s_old] == 0) &
              (s_old != self$node$s[i]))
            # Delete class
          {
            self$node$theta$erase(s_old)
            enqueue(self$node$idBucket, s_old)
            nk <- nk[names(nk) != s_old]
          }
          nk[self$node$s[i]] <- nk[self$node$s[i]] + 1
        }

        # Sample theta from the conditional posterior
        for (label in self$node$theta$keys())
        {
          sumy <- sum(y[self$node$s == label])
          aXs <- sum(aX[self$node$s == label])
          self$node$theta$insert(label, rgamma(1, shape = self$alpha + sumy, rate = self$beta + aXs))
        }
      }
  )
)
