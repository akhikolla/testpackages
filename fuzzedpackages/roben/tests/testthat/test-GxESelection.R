test_that("check_selection_sparse", {
  iter = 5000
  fit=roben(X, Y, E, clin, iterations = iter)
  sel=GxESelection(fit)
  expect_equal(sel$method, "Median Probability Model (MPM)")
  expect_equal(dim(sel$indicator), c(ncol(E)+1, ncol(X)))
})


test_that("check_selection_nonsparse", {
  iter = 5000
  fit=roben(X, Y, E, clin, iterations = iter, sparse = FALSE)
  sel=GxESelection(fit, prob=0.9)
  expect_equal(sel$method, "90% credible interval")
  expect_equal(dim(sel$indicator), c(ncol(E)+1, ncol(X)))
})
