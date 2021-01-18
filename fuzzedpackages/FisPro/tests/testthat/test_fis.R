context("fis")

test_that("fis infer", {
  fis_file <- system.file("extdata", "test.fis", package = "FisPro")
  fis <- new(fis, fis_file)

  expect_equal(fis$input_size, 2)
  expect_equal(fis$output_size, 2)
  
  expect_equal(fis$infer_output(c(0.5, 0.5), 0), 0.5)
  expect_equal(fis$infer_output(c(0.5, 0.5), 1), 0.5)
  expect_equal(fis$infer(c(0.5, 0.5)), c(0.5, 0.5))

  expect_equal(fis$infer_output(c(0.25, 0.75), 0), 0.333, tolerance = 1e-3)
  expect_equal(fis$infer_output(c(0.25, 0.75), 1), 0.361, tolerance = 1e-3)
  expect_equal(fis$infer(c(0.25, 0.75)), c(0.333, 0.361), tolerance = 1e-3)
})
