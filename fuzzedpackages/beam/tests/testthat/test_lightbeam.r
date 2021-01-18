data("TCPAprad")

test_that("NgreaterthanP", {

  # lightbeam
  res <- lightbeam(X = TCPAprad, verbose=TRUE)

  # Check results
  expect_equal(res@x[1], 0.01356278, tolerance=1e-5)
  expect_equal(res@x[10], 0.02367338, tolerance=1e-5)
  expect_equal(res@x[50], 0.003945791, tolerance=1e-5)
  expect_equal(res@x[90], 0.0001630148, tolerance=1e-5)

  # Check dim
  expect_equal(dim(res), c(189, 189))
  
  # Check class
  expect_equal(inherits(res, "dgCMatrix"), TRUE)
  
  # Check labs
  expect_equal(length(colnames(res)), 189)
  
})

test_that("NlowerthanP", {
  
  # lightbeam
  res <- lightbeam(X = t(TCPAprad), verbose=TRUE)
  
  # Check results
  expect_equal(res@x[1], 0.02458413, tolerance=1e-5)
  expect_equal(res@x[10], 0.004667081, tolerance=1e-5)
  expect_equal(res@x[50], 0.01377481, tolerance=1e-5)
  expect_equal(res@x[63], 0.09872211, tolerance=1e-5)
  
  # Check dim
  expect_equal(dim(res), c(164, 164))
  
  # Check class
  expect_equal(inherits(res, "dgCMatrix"), TRUE)
  
  # Check labs
  expect_equal(length(colnames(res)), 164)
  
})
