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

# test disabled in release-versions (CRAN) (takes too long)
if(POUMMIsADevRelease()) {
context("POUMM fit")
  
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
  
  test_that("fit passes without errors", {
    expect_is(fit1 <- POUMM(z = z[1:N], tree = tree, 
                            spec = specifyPOUMM(z[1:N], tree, nSamplesMCMC = 10000), 
                            verbose=FALSE), "POUMM")
    
    # end value not changed warnings from coda
    expect_warning(smmShort <- summary(fit1))
    expect_is(smmShort, "summary.POUMM")
    
    expect_warning(smmLong <- summary(fit1, mode="long"))
    expect_is(smmLong, "summary.POUMM")
    
    expect_warning(smmExpert <- summary(fit1, mode="expert"))
    expect_is(smmExpert, "summary.POUMM")
    
    expect_warning(plList <- plot(fit1, interactive = FALSE, doPlot = FALSE))
    expect_is(plList, "list")
  })
  
}