test_that("plot_RLumCarlo(): basic run", {
  testthat::skip_on_cran()
  local_edition(3)

  ## create data to plot
  results_par <- lapply(c(1.2, 1.4), function(r){
    run_MC_TL_TUN(
    s = 3.5e12,
    E = 1.45,
    rho = 0.015,
    r_c = r,
    times = 100:120,
    method = "seq"
  )})

  time <- run_MC_ISO_TUN(
      s = 3.5e12,
      E = 1.45,
      rho = 0.015,
      r_c = 0.5,
      times = 100:120,
      method = "seq")

  ## break plot function
  expect_error(plot_RLumCarlo(object = "error"), "'object' needs to be of class RLumCarlo_Model_Output!")
  expect_error(plot_RLumCarlo(object = list("error")), "At least one element in the list is not of class RLumCarlo_Model_Output!")
  expect_error(plot_RLumCarlo(object = results_par[[1]], FUN = function(){2}), regexp = "FUN has no argument!")

  ## simple run, one dataset
  expect_silent(plot_RLumCarlo(object = results_par[[1]]))
  expect_silent(plot_RLumCarlo(object = results_par[[1]], norm = TRUE))
  expect_silent(plot_RLumCarlo(object = results_par[[1]], plot_uncertainty = NULL))
  expect_silent(plot_RLumCarlo(object = results_par[[1]], plot_value = "sum"))
  expect_silent(plot_RLumCarlo(object = results_par[[1]], plot_value = NULL))
  expect_silent(plot_RLumCarlo(object = results_par[[1]], plot_uncertainty = "sd"))
  expect_silent(plot_RLumCarlo(object = results_par[[1]], plot_uncertainty = "var"))
  expect_silent(plot_RLumCarlo(object = results_par[[1]], plot_uncertainty = "range"))

  ## test FUN
  expect_silent(plot_RLumCarlo(object = results_par[[1]], FUN = function(x) x^2))

  ## plot the list
  expect_silent(plot_RLumCarlo(object = results_par, main = "Test 2"))
  expect_silent(plot_RLumCarlo(object = time, main = "Test 2"))
  expect_silent(plot_RLumCarlo(object = rep(results_par,4), main = "Rainbow"))

})

