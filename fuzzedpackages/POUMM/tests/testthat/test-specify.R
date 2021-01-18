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

context("specifyPOUMM")
test_that("default calls no error", {
  expect_silent(specifyPOUMM())
  expect_silent(specifyPOUMM_ATH2tMeanSe())
  expect_silent(specifyPOUMM_ATH2tMeanSeG0())
  expect_silent(specifyPOUMM_ATS())
  expect_silent(specifyPOUMM_ATSG0())
  expect_silent(specifyPOUMM_ATSSeG0())
  expect_silent(specifyPMM())
  expect_silent(specifyPMM_H2tMeanSe())
  expect_silent(specifyPMM_H2tMeanSeG0())
  expect_silent(specifyPMM_SSeG0())
})