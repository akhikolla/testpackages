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

context("POUMM parameters") 

# At POUMM stationary state (equilibrium, t=Inf)
H2 <- POUMM::H2(alpha = 0.75, sigma = 1, sigmae = 1, t = Inf)     # 0.4
alpha <- POUMM::alpha(H2 = H2, sigma = 1, sigmae = 1, t = Inf)    # 0.75
sigma <- POUMM::sigmaOU(H2 = H2, alpha = 0.75, sigmae = 1, t = Inf) # 1
sigmae <- POUMM::sigmae(H2 = H2, alpha = 0.75, sigma = 1, t = Inf) # 1

test_that("stationary state (t = Inf)",
          expect_equal(c(H2, alpha, sigma, sigmae), c(0.4, 0.75, 1, 1)))

# At finite time t = 0.2
H2 <- POUMM::H2(alpha = 0.75, sigma = 1, sigmae = 1, t = 0.2)     # 0.1473309
alpha <- POUMM::alpha(H2 = H2, sigma = 1, sigmae = 1, t = 0.2)    # 0.75
sigma <- POUMM::sigmaOU(H2 = H2, alpha = 0.75, sigmae = 1, t = 0.2) # 1
sigmae <- POUMM::sigmae(H2  =  H2, alpha = 0.75, sigma = 1, t = 0.2) # 1

test_that("at finite time (t = 0.2)", 
          expect_equal(c(H2, alpha, sigma, sigmae), c(
            POUMM::H2(alpha, sigma, sigmae, 0.2), 0.75, 1, 1)))

