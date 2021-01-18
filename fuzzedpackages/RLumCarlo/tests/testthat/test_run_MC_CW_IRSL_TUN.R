test_that("basic run", {
  testthat::skip_on_cran()
  local_edition(3)

  ##break the function on purpose
  expect_error(run_MC_CW_IRSL_TUN(method = "error"), "Allowed keywords for 'method' are either 'par' or 'seq'!")
  expect_error(run_MC_CW_IRSL_TUN(output = "error"), "Allowed keywords for 'output' are either 'signal' or 'remaining_e'!")

  ## run function with basic setting
  ## sequential and multi-core
  results_seq <- expect_silent(run_MC_CW_IRSL_TUN(
    A = 0.8,
    rho = 1e-4,
    times = 0:100,
    clusters = 10,
    N_e = 2,
    r_c = 0.05,
    delta.r = 0.1,
    method = "seq"
  ))

  results_par <- expect_silent(run_MC_CW_IRSL_TUN(
    A = 0.8,
    rho = 1e-4,
    times = 0:100,
    clusters = 10,
    N_e = 2,
    r_c = 0.05,
    delta.r = 0.1,
    method = "par"
  ))

  ##check cluster system
  expect_silent(run_MC_CW_IRSL_TUN(
    A = 0.8,
    rho = 1e-4,
    times = 0:100,
    clusters = create_ClusterSystem(10),
    N_e = 20,
    r_c = 0.05,
    delta.r = 0.1,
    method = "seq"
  ))

  ## check output
  expect_s3_class(results_par, class = "RLumCarlo_Model_Output")
  expect_length(results_par, 2)
  expect_s3_class(results_seq, class = "RLumCarlo_Model_Output")
  expect_length(results_seq, 2)

})
