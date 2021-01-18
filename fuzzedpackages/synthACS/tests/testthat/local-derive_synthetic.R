library(testthat)
library(synthACS)

context("local - derive synthetic microdata for geography set")

test_that("get correct results -- serial", {
  #-------------------------------
  testthat::skip_on_cran()
  testthat::skip_on_covr()
  testthat::skip_on_travis()
  #-------------------------------
  ## parallel == FALSE
  syn <- derive_synth_datasets(ca_dat, parallel= FALSE)
  
  # test class, structure, dimensions, etc
  expect_equal(class(syn), c("synthACS","list"))
  expect_true(is.synthACS(syn))
  expect_true(all(unlist(lapply(syn, is.macro_micro))))
  expect_equal(nrow(ca_dat$estimates[[1]]), length(syn))
  expect_true(all(unlist(lapply(syn, length)) == 2))
  expect_true(all(unlist(lapply(syn, function(l) ncol(l[[2]]))) == length(ca_dat$estimates) + 1))
  expect_true(all.equal(lapply(syn, function(l) names(l[[2]])),
                        replicate(length(syn), 
                                  c("age", "gender", "marital_status", "edu_attain", "emp_status", "nativity", 
                                    "pov_status", "geog_mobility", "ind_income", "race", "p"), simplify=FALSE), 
                        check.attributes= FALSE))
  
  # test total probabilities
  expect_true(all.equal(unlist(lapply(syn, function(l) sum(l[[2]]$p))), rep(1,58), 
                        tolerance= 1e-14, check.attributes= FALSE))
  
})

test_that("get correct results -- parallel", {
  #-------------------------------
  testthat::skip_on_cran()
  testthat::skip_on_covr()
  testthat::skip_on_travis()
  #-------------------------------
  ## parallel == TRUE
  syn <- derive_synth_datasets(ca_dat, parallel= TRUE)
  
  # test class, structure, dimensions, etc
  expect_equal(class(syn), c("synthACS","list"))
  expect_true(is.synthACS(syn))
  expect_true(all(unlist(lapply(syn, is.macro_micro))))
  expect_equal(nrow(ca_dat$estimates[[1]]), length(syn))
  expect_true(all(unlist(lapply(syn, length)) == 2))
  expect_true(all(unlist(lapply(syn, function(l) ncol(l[[2]]))) == length(ca_dat$estimates) + 1))
  expect_true(all.equal(lapply(syn, function(l) names(l[[2]])),
                        replicate(length(syn), 
                                  c("age", "gender", "marital_status", "edu_attain", "emp_status", "nativity", 
                                    "pov_status", "geog_mobility", "ind_income", "race", "p"), simplify=FALSE), 
                        check.attributes= FALSE))
  
  # test total probabilities
  expect_true(all.equal(unlist(lapply(syn, function(l) sum(l[[2]]$p))), rep(1,58), 
                        tolerance= 1e-14, check.attributes= FALSE))
  
})
