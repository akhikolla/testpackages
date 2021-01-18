# example
#gp <- as.data.frame(fishflux::growth_params("Zebrasoma scopas"))

test_that("Simple corner cases", {
  needs_api()
  check_api()
  gp <- as.data.frame(fishflux::growth_params("Zebrasoma scopas"))
  expect_error(fishflux::gp("wrong name"))
  expect_s3_class(gp, "data.frame")
  expect_equal(gp, as.data.frame(fishflux::growth_params("Zebrasoma scopas")))
  expect_length(gp, 7)
  expect_equal(as.character(gp[1, 1]), "Zebrasoma scopas")
  expect_true(is.numeric(gp[1, 3]))
  expect_true(is.numeric(gp[1, 4]))
  expect_true(is.numeric(gp[1, 5]))
  expect_gt(gp[1, 3], 0)
  expect_gt(gp[1, 4], 0)
})
