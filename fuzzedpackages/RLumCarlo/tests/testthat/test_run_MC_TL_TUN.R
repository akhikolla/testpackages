test_that("basic run", {
  testthat::skip_on_cran()
  local_edition(3)

  ##break the function on purpose
  expect_error(run_MC_TL_TUN(method = "error"), "Allowed keywords for 'method' are either 'par' or 'seq'!")
  expect_error(run_MC_TL_TUN(output = "error"), "Allowed keywords for 'output' are either 'signal' or 'remaining_e'!")

  ## run function with basic setting
  ## sequential and multicore
  results_seq <- expect_silent(run_MC_TL_TUN(
    s = 1e8,
    E = 0.5,
    rho = 1e-7,
    r_c = 0.1,
    times = seq(0,700,10),
    clusters = 100,
    N_e = 20,
    delta.r = 1e-3,
    method = "seq"
  ))

  results_par <- expect_silent(run_MC_TL_TUN(
    s = 1e8,
    E = 0.5,
    rho = 1e-7,
    r_c = 0.1,
    times = seq(0,700,10),
    clusters = 100,
    N_e = 20,
    method = "par",
    delta.r = 1e-3
  ))

  ## check cluster system
  expect_silent(run_MC_TL_TUN(
    s = 1e8,
    E = 0.5,
    rho = 1e-7,
    r_c = 0.1,
    times = seq(0,500,10),
    clusters = create_ClusterSystem(10),
    N_e = 100,
    method = "seq",
    delta.r = 1e-3
  ))

  ## check output
  expect_s3_class(results_par, class = "RLumCarlo_Model_Output")
  expect_length(results_par, 2)
  expect_s3_class(results_seq, class = "RLumCarlo_Model_Output")
  expect_length(results_seq, 2)

})
