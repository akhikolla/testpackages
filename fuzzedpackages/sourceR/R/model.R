#####################################################
# Name: model.R                                     #
# Author: Chris Jewell <c.jewell@lancaster.ac.uk>   #
# Modified by: Poppy Miller                         #
# Created: 2016-06-15                               #
# Copyright: Chris Jewell 2016                      #
# Purpose: Source attribution model creation        #
#####################################################


#' Builds the source attribution model. Is not intended to be used by a regular user.
#' Developers only here!
#' @return Object of \code{\link{R6Class}}.
#' @format \code{\link{R6Class}} object.
#' @field y 3D array of [type, time, location] of the number of human cases
#' @field X 3D array of the number of positive samples for each type, source and time
#' [type, source, time]
#' @field R 3D array of normalised relative prevalences for each timepoint
#' [type, source, time]
#' @field Time a character vector of timepoint ids matching time dimension in y and R
#' @field Location a character vector of location ids matching location dimension in y
#' @field Sources a character vector of source ids matching the source dimension in X
#' @field Type a character vector of type ids matching the type dimension in X
#' @field prev a 2D array (matrix) of [source, time].
#' @field a_q concentration parameter for the DP
#' @field a_theta shape parameter for the Gamma base distribution for the DP
#' @field b_theta rate parameter for the Gamma base distribution for the DP
#' @field a_r 3D array of [type, src, time] for the hyperprior on the relative prevalences R
#' @field a_alpha 3D array of [source, time, location] for the prior on the alpha parameters
#' @field s vector giving initial group allocation for each type for the DP
#' @field theta vector giving initial values for each group in the DP
#' @field alpha 3D array of [source, time, location] giving initial values for the alpha
#' parameters

DPModel_impl <- R6::R6Class(
  "DPModel_impl",
  public = list(
    # LambdaNode calculates \lambda^{\prime}_{itl} = \sum_{j=1}^{m} r_{ijt} k_{jt} \alpha_{jtl}
    #   N.B. q is multiplied on at a later stage, in the PoissonNode.
    LambdaNode = R6::R6Class("LambdaNode",
                              inherit = FormulaNode,
                              public = list(
                                initialize = function(k, rNodes, alpha, name) {
                                  super$initialize(name = name)
                                  self$addParent(k, 'k')
                                  self$parents$r <- rNodes  # N.B. attaches a list of nodes as a parent node
                                  # Todo: replace with R6 'NodesList' class
                                  sapply(rNodes, function(node) node$children[[self$name]] <- self)
                                  self$addParent(alpha, 'alpha')
                                },
                                getData = function() {
                                  r <- sapply(self$parents$r, function(x) x$getData())
                                  # identical(apply(campy[rownames(r), c(5, 2, 3, 4,7, 6)], 2, function(x) (x+0.000001)/sum(x+0.000001)), r)
                                  r %*% (self$parents$k$getData() * self$parents$alpha$getData())
                                }
                              )
    ),
    # Q
    qNodes = NULL,

    # Node lists
    alphaNodes = list(),
    yNodes = list(),
    rNodes = list(),

    initialize = function(y, X, R, Time, Location, Sources, Type, prev, a_q, a_theta,
                          b_theta, a_r, a_alpha, s, theta, alpha)
    {

      self$qNodes <- DirichletProcessNode$new(theta = theta, s = s, alpha = a_q,
                                              base = dgamma, shape = a_theta,
                                              rate = b_theta, name = 'q')
      # Node lists
      self$alphaNodes <- list()
      self$yNodes <- list()
      self$rNodes <- list()

      # Construct y, lambda, a, R, and k for location/time pairs
      for (time in unique(Time)) {
        # Prevalences
        k <- DataNode$new(data = setNames(prev[, time], Sources), name = paste('k',time,sep = '_'))

        # R_{tj}s
        self$rNodes[[time]] <- list()
        for (src in 1:length(Sources)) {
          # Dirichlet prior on r, as a result of Dirichlet/Multinomial conjugacy on R.
          a_r_full <- DataNode$new(data = a_r[, src, time] + X[, src, time], name = paste('a_r', time, src, sep = '_'))
          self$rNodes[[time]][[src]] <- DirichletNode$new(data = setNames(R[, src, time], Type),
                                                          alpha = a_r_full,
                                                          name = paste('r', time, src, sep = '_'))
        }
        names(self$rNodes[[time]]) <- Sources

        # Here we use y ~ Poisson(q \lambda), \lambda = r^T \alpha
        self$alphaNodes[[time]] <- list()
        self$yNodes[[time]] <- list()
        for (location in unique(Location)) {

          # a_t_l prior vector for alpha_t_l
          a_tl <- DataNode$new(data = a_alpha[, time, location],
                               name = paste('a_alpha', time, location, sep = '_'))

          # Location specific alpha
          alpha_tl <- DirichletNode$new(data = setNames(alpha[, time, location], Sources),
                                        alpha = a_tl,
                                        name = paste('alpha', time, location, sep = '_'))

          # Construct the lambda_i node
          lambdaPrime <- self$LambdaNode$new(k = k,
                                             rNodes = self$rNodes[[time]],
                                             alpha = alpha_tl,
                                             name = paste('lambda', time, location, sep = '_'))

          # From the point of view of q, lambda is fixed and known, therefore classed as the offset.
          y_tl <- PoissonNode$new(data = y[, time, location],
                                  lambda = self$qNodes,
                                  offset = lambdaPrime,
                                  name = paste('y',time, location,sep = '_'))

          self$alphaNodes[[time]][[location]] <- alpha_tl
          self$yNodes[[time]][[location]] <- y_tl
        }
      }
    }
  )
)
