

library(testthat)
library(synthACS)

context("CRAN - derive synthetic microdata for geography set")

test_that("disaggregate_md works", {
  # load geography data, run disaggregation
  m_dat <- synthACS:::disaggregate_md(ca_dat$estimates)
  
  # test outputs
  expect_equal(length(m_dat), nrow(ca_dat$estimates[[1]]))
  expect_true(all(lapply(m_dat, length) == length(ca_dat$estimates)))
  expect_equal(names(m_dat[[1]]), names(ca_dat$estimates))
  expect_equal(sapply(m_dat[[10]], names), sapply(ca_dat$estimates, names))
  expect_true(all(unlist(lapply(m_dat[[30]], is.numeric))))
  
})


test_that("error checking", {
  diamonds <- data.frame(
    carat= rexp(100),
    cut= factor(sample(c("A", "B", "C"), size= 100, replace= TRUE)),
    x= runif(100, min= 0, max= 10),
    y= runif(100, min= 0, max= 10),
    x= runif(100, min= 0, max= 10)
  )
  expect_error(derive_synth_datasets(diamonds, parallel= FALSE))
  
  expect_error(derive_synth_datasets(ca_dat, parallel= TRUE, leave_cores= -1L))
  expect_error(derive_synth_datasets(ca_dat, parallel= TRUE, leave_cores= 2.5))
})
