
library(testthat)
library(synthACS)

context("pull_acs_basetables")

test_that("returns results accurately - counties", {
  # create test geography and data
  ca_geo <- geo.make(state= 'CA', county= 'Los Angeles')
  ca_dat <- pull_acs_basetables(2014, 5, ca_geo, table_vec= c("B01001", "B01002", "B01003"))
  
  # test:
  synthACS:::confirm_macroACS_class(ca_dat)
})

test_that("returns results accurately - state", {
  # create test geography and data
  ca_geo <- geo.make(state= "CA")
  ca_dat <- pull_acs_basetables(2015, 5, ca_geo, c("B01001", "B01002", "B01003"))
  
  # test:
  synthACS:::confirm_macroACS_class(ca_dat)
})

context("pull_synth_data")

test_that("returns results accurately - counties", {
  # create test geography and data
  ca_geo <- geo.make(state= 'CA', county= 'Los Angeles')
  ca_dat <- pull_synth_data(2014, 5, ca_geo)
  
  # test:
  synthACS:::confirm_macroACS_class(ca_dat)
})

test_that("returns results accurately - state", {
  # create test geography and data
  ca_geo <- geo.make(state= "CA")
  ca_dat <- pull_synth_data(2016, 5, ca_geo)
  # test:
  synthACS:::confirm_macroACS_class(ca_dat)
  
  # create test geography and data
  ca_geo <- geo.make(state= "CA", county= '*')
  ca_dat <- pull_synth_data(2016, 5, ca_geo)
  # test:
  synthACS:::confirm_macroACS_class(ca_dat)
})