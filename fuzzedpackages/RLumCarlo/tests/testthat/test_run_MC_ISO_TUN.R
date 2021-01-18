test_that("basic run", {
  testthat::skip_on_cran()
  local_edition(3)

  ##break the function on purpose
  expect_error(run_MC_ISO_TUN(method = "error"), "Allowed keywords for 'method' are either 'par' or 'seq'!")
  expect_error(run_MC_ISO_TUN(output = "error"), "Allowed keywords for 'output' are either 'signal' or 'remaining_e'!")

  ## run function with basic setting
  ## sequential and multicore
  results_seq <- expect_silent(run_MC_ISO_TUN(
    E = 0.5,
    s = 1e12,
    `T` = 150,
    rho = 1e-4,
    times = 0:100,
    clusters = 1e1,
    N_e = 2,
    r_c = 1e-4,
    delta.r = 0.5,
    method = "seq"
  ))

  expect_s3_class(run_MC_ISO_TUN(
     E = 1.2,
     s = 1e10,
     T = 200,
     N_e = 200,
     rho = 0.007,
     clusters = 10,
     times = seq(0, 5000)),
     "RLumCarlo_Model_Output")

  results_par <- expect_silent(run_MC_ISO_TUN(
    E = 0.5,
    s = 1e12,
    `T` = 150,
    rho = 1e-4,
    times = 0:100,
    clusters = 1e1,
    N_e = 2,
    r_c = 1e-4,
    delta.r = 0.5,
    method = "par"
  ))

  ##check cluster system
  expect_silent(run_MC_ISO_TUN(
    E = 0.5,
    s = 1e12,
    `T` = 150,
    rho = 1e-4,
    times = 0:100,
    clusters = create_ClusterSystem(10),
    N_e = 2,
    r_c = 1e-4,
    delta.r = 0.5,
    method = "seq"
  ))

  ## check output
  expect_s3_class(results_par, class = "RLumCarlo_Model_Output")
  expect_length(results_par, 2)
  expect_s3_class(results_seq, class = "RLumCarlo_Model_Output")
  expect_length(results_seq, 2)

})
