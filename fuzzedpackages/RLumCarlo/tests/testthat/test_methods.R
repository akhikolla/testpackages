test_that("basic run", {
  testthat::skip_on_cran()
  local_edition(3)

  ## summary
  expect_s3_class(summary(run_MC_TL_TUN(
    s = 3.5e12,
    E = 1.45,
    rho = 0.015,
    r_c = 0.5,
    times = 100:500,
    method = "seq"
  )), "data.frame")

  expect_s3_class(summary(run_MC_TL_TUN(
    s = 3.5e12,
    E = 1.45,
    rho = 0.015,
    r_c = 0.5,
    times = 100:110,
    clusters = 2,
    method = "seq"
  )), "data.frame")


  ## check object dimension 2 (one cluster)
  expect_s3_class(summary(run_MC_TL_TUN(
    s = 3.5e12,
    E = 1.45,
    rho = 0.015,
    r_c = 0.5,
    times = 100:110,
    clusters = 1,
    method = "seq"
  )), "data.frame")

  ## check verbose (one cluster)
  expect_silent(summary(run_MC_TL_TUN(
    s = 3.5e12,
    E = 1.45,
    rho = 0.015,
    r_c = 0.5,
    times = 100:110,
    clusters = 1,
    method = "seq"
  ), verbose = FALSE))

})

test_that("check c()", {
  testthat::skip_on_cran()
  local_edition(3)

  ## create need output objects
  objectA <- run_MC_TL_LOC(
    s = 1e14,
    E = 0.9,
    times = 0:100,
    b = 1,
    method = "seq",
    clusters = 30,
    r = 1)

  objectB <- run_MC_TL_TUN(
    s = 1e12,
    E = 0.9,
    rho = 1,
    r_c = 0.1,
    times = 0:100,
    b = 1,
    clusters = 50,
    method = 'seq',
    delta.r = 1e-1)


  objectC <- objectB
  attr(objectC, which = "model") <- "run_MC_LM_OSL_TUN"

  objectD <- objectA
  objectD$time <- 1:1000

   ## crash function
   expect_error(c(objectA, objectD), regexp = "You cannot combine objects with different time vectors")

   ## warning
   expect_warning(c(objectA,objectC), regexp = "Stimulation modes differ")

   ## success
   expect_s3_class(c(objectA, objectB), class = "RLumCarlo_Model_Output")

})

