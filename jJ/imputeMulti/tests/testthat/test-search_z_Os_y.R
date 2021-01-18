
library(testthat)
library(imputeMulti)

context("int- search z_Os_y")


test_that("query text works for any number of columns",{
  data(tract2221)
  t2 <- tract2221[,1:4]
  x_p <- multinomial_stats(t2, "possible.obs")
  r1 <- data.frame(age= 4L, gender= 1L, mar= NA, edu= NA)
  r2 <- data.frame(age= 3L, gender= 1L, mar= NA, edu= 3L)

  expect_equal(imputeMulti:::create_search_query(x_p, r1, names(x_p)),
               "select rownames from x_possible where age= '30_34'and gender= 'Female'")
  expect_equal(imputeMulti:::create_search_query(x_p, r2, names(x_p)),
               "select rownames from x_possible where age= '25_29'and gender= 'Female'and edu_attain= 'grad_deg'")
  
})


test_that("works as designed", {
  # as designed means that x_possible, x_y, and z_Os_y have been specified as anticipated
  #  (variable names, factor variables, counts, etc)
  data(tract2221)
  t2 <- tract2221[,1:4]
  x_p <- multinomial_stats(t2, "possible.obs")
  x_y <- multinomial_stats(t2, "x_y")
  z_Os_y <- multinomial_stats(t2, "z_Os_y")
  
  incomp_ind <- imputeMulti:::search_z_Os_y(z_Os_y, x_p)
  
  expect_equal(length(incomp_ind), nrow(z_Os_y))
  expect_true(is.list(incomp_ind))
  expect_true(all(unlist(lapply(incomp_ind, is.numeric))))
  expect_true(all(unlist(lapply(incomp_ind, is.integer))))
  expect_true(all(unlist(lapply(incomp_ind, function(l) all(l %% 1 == 0)))))
  expect_true(all(unlist(lapply(incomp_ind, length)) > 0))
})

#context("int- count_sumStats")

