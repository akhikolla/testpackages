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

library(testthat)
library(POUMM)
library(mvtnorm)


context("POUMM Likelihood")

set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")

N <- 100
tree <- ape::rtree(N)  

g0 <- 16
alpha <- 2
theta <- 4
sigma <- .2
sigmae <- .7
se = rep(0, N) #rexp(N, 1/.01)
z <- rVNodesGivenTreePOUMM(tree, g0, alpha, theta, sigma, sqrt(sigmae^2+se^2))

#' An "algebraic" implementation of the POUMM likelihood calculation based on
#' multivariate normal density function
#' This is an internal function used for tests only.
dVTipsGivenTreePOUMMg0Alg <- function(z, tree, alpha, theta, sigma, sigmae, se, g0) {
  tanc <- ape::vcv(tree)
  tipTimes <- diag(tanc)
  N <- length(tree$tip.label)
  
  D0 <- sapply(tipTimes, function(t) exp(-alpha*t))
  D1 <- 1-D0
  Vt <- covVTipsGivenTreePOUMM(tree, alpha, sigma, 0, tanc)
  Ve <- diag(sigmae^2+se^2, N, N)
  
  muz <- D0 * g0  +  D1 * theta
  
  dmvnorm(z, muz, Vt+Ve, log=TRUE)
}

pruneInfo <- POUMM:::pruneTree(tree, z[1:N], se)

test_that(
  "fastLik(R) vs algebraic lik", {
    expect_equivalent(likPOUMMGivenTreeVTips(z[1:N], tree, 0, theta, sigma, sqrt(sigmae^2+se^2), g0),
                      dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, 0, theta, sigma, sigmae, se, g0))
    expect_equivalent(likPOUMMGivenTreeVTips(z[1:N], tree, 0, theta, sigma, sqrt(0+se^2), g0),
                      dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, 0, theta, sigma, 0, se, g0))
    expect_equivalent(likPOUMMGivenTreeVTips(z[1:N], tree, alpha, theta, sigma, sqrt(sigmae^2+se^2), g0),
                      dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, alpha, theta, sigma, sigmae, se, g0))
    expect_equivalent(likPOUMMGivenTreeVTips(z[1:N], tree, alpha, theta, sigma, sqrt(0+se^2), g0),
                      dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, alpha, theta, sigma, 0, se, g0))
    expect_equivalent(likPOUMMGivenTreeVTips(z[1:N], tree, alpha, theta, sigma, sqrt(sigmae^2+0), g0),
                      dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, alpha, theta, sigma, sigmae, 0, g0))
  })

test_that(
  "fastLik(C++) vs algebraic lik", {
    expect_equivalent(likPOUMMGivenTreeVTipsC(pruneInfo$integrator, 0, theta, sigma, sigmae, g0),
                      dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, 0, theta, sigma, sigmae, se, g0))
    expect_equivalent(likPOUMMGivenTreeVTipsC(pruneInfo$integrator, 0, theta, sigma, 0, g0),
                      dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, 0, theta, sigma, 0, se, g0))
    expect_equivalent(likPOUMMGivenTreeVTipsC(pruneInfo$integrator, alpha, theta, sigma, sigmae, g0),
                      dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, alpha, theta, sigma, sigmae, se, g0))
    expect_equivalent(likPOUMMGivenTreeVTipsC(pruneInfo$integrator, alpha, theta, sigma, 0, g0),
                      dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, alpha, theta, sigma, 0, se, g0))
    expect_equivalent(likPOUMMGivenTreeVTips(z[1:N], tree, alpha, theta, sigma, sigmae, g0),
                      dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, alpha, theta, sigma, sigmae, 0, g0))
  })



