test_that("basic run", {
  testthat::skip_on_cran()
  local_edition(3)

  ##break the function on purpose
  expect_error(run_MC_LM_OSL_TUN(method = "error"), "Allowed keywords for 'method' are either 'par' or 'seq'!")
  expect_error(run_MC_LM_OSL_TUN(output = "error"), "Allowed keywords for 'output' are either 'signal' or 'remaining_e'!")

  ## run function with basic setting
  ## sequential and multicore
  results_seq <- expect_silent(run_MC_LM_OSL_TUN(
    A = 1,
    rho = 1e-2,
    times = 0:1000,
    clusters = 15,
    N_e = 2,
    r_c = 0.01,
    delta.r = 1e-1,
    method = "seq"
  ))

  results_par <- expect_silent(run_MC_LM_OSL_TUN(
    A = 1,
    rho = 1e-2,
    times = 0:100,
    clusters = 15,
    N_e = 2,
    r_c = 0.01,
    delta.r = 1e-1,
    method = "par"
  ))

  ## create cluster system
  expect_silent(run_MC_LM_OSL_TUN(
    A = 1,
    rho = 1e-2,
    times = 0:1000,
    clusters = create_ClusterSystem(10),
    N_e = 20,
    r_c = 0.01,
    delta.r = 1e-1,
    method = "seq"
  ))

  ## check output
  expect_s3_class(results_par, class = "RLumCarlo_Model_Output")
  expect_length(results_par, 2)
  expect_s3_class(results_seq, class = "RLumCarlo_Model_Output")
  expect_length(results_seq, 2)

})
