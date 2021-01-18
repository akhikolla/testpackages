test_that("utils.R: basic run", {
  testthat::skip_on_cran()
  local_edition(3)

  ## .registerClusters ... cores and verbose
  cl <- expect_type(.registerClusters(method = "par", cores = NULL, verbose = TRUE), "list")
  parallel::stopCluster(cl)

  ## .distribute_electrons
  expect_s3_class(
    .distribute_electrons(clusters = create_ClusterSystem(n = 10), N_system = 100),
    'data.frame'
  )

})

