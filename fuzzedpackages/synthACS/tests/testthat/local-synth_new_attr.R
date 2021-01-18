
library(testthat)
library(synthACS)

#----------------------------------------------------------
context("LOCAL -- synthetic new attribute creation")
#----------------------------------------------------------
test_that("single synthetic dataset -- real df, fake ST (DF)", {
  # geography
  # load(...) 
  # set up symbol tables
  
  levels <- c("A", "B", "C")
  ST <- data.frame(marital_status= rep(levels(test_micro$marital_status), each= 7),
                   race= rep(levels(test_micro$race), 5))
  ST <- do.call("rbind", replicate(3, ST, simplify=FALSE))
  ST <- ST[order(ST$marital_status, ST$race),]
  ST$pct <- rep(c(0.33, 0.34, 0.33), 35)
  ST$levels <- rep(levels, 35)
  
  # run
  syn <- synthetic_new_attribute(df= test_micro, prob_name= "p", attr_name= "variable",
                                 conditional_vars= c("marital_status", "race"), sym_tbl= ST)
  
  # test output
  expect_equal(sum(syn$p), 1)
  expect_equal(tapply(syn$p, syn$gender, sum),
               tapply(test_micro$p, test_micro$gender, sum))
  expect_equal(tapply(syn$p, syn$marital_status, sum),
               tapply(test_micro$p, test_micro$marital_status, sum))
  expect_equal(tapply(syn$p, syn$race, sum),
               tapply(test_micro$p, test_micro$race, sum))
  expect_equal(nrow(test_micro) * length(levels), nrow(syn))
  expect_equal(ncol(test_micro), ncol(syn) - 1)
  expect_true(all.equal(as.vector(tapply(syn$p, syn$variable, sum)), c(0.33, 0.34, 0.33), check.attributes=FALSE))
})

test_that("single synthetic dataset -- real df, fake ST (DT)", {
  # set up symbol tables
  levels <- c("A", "B", "C")
  ST <- data.frame(marital_status= rep(levels(test_micro$marital_status), each= 7),
                   race= rep(levels(test_micro$race), 5))
  ST <- do.call("rbind", replicate(3, ST, simplify=FALSE))
  ST <- ST[order(ST$marital_status, ST$race),]
  ST$pct <- rep(c(0.33, 0.34, 0.33), 35)
  ST$levels <- rep(levels, 35)
  
  data.table::setDT(test_micro); data.table::setDT(ST)
  # run
  syn <- synthetic_new_attribute(df= test_micro, prob_name= "p", attr_name= "variable",
                                 conditional_vars= c("marital_status", "race"), sym_tbl= ST)
  
  # test output
  expect_equal(sum(syn$p), 1)
  expect_equal(tapply(syn$p, syn$gender, sum),
               tapply(test_micro$p, test_micro$gender, sum))
  expect_equal(tapply(syn$p, syn$marital_status, sum),
               tapply(test_micro$p, test_micro$marital_status, sum))
  expect_equal(tapply(syn$p, syn$race, sum),
               tapply(test_micro$p, test_micro$race, sum))
  expect_equal(nrow(test_micro) * length(levels), nrow(syn))
  expect_equal(ncol(test_micro), ncol(syn) - 1)
  expect_true(all.equal(as.vector(tapply(syn$p, syn$variable, sum)), c(0.33, 0.34, 0.33), check.attributes=FALSE))
})


test_that("single synthetic dataset -- real DF, ST (DF)", {
  # load data inputs
  test_micro <- syn[[1]][[2]]
  work_ST <- towork[[1]]
  class(work_ST) <- "data.frame"
  
  # run
  syn <- synthetic_new_attribute(df= test_micro, prob_name= "p", attr_name= "transit",
                                 conditional_vars= c("emp_status", "age"), sym_tbl= work_ST)
  
  # test output
  expect_equal(sum(syn$p), 1)
  expect_equal(tapply(syn$p, syn$emp_status, sum),
               tapply(test_micro$p, test_micro$emp_status, sum))
  expect_equal(tapply(syn$p, syn$age, sum),
               tapply(test_micro$p, test_micro$age, sum))
  expect_equal(tapply(syn$p, syn$race, sum),
               tapply(test_micro$p, test_micro$race, sum))
  expect_equal(ncol(test_micro), ncol(syn) - 1)
})

test_that("single synthetic dataset -- real DF, ST (DT)", {
  # load data inputs
  test_micro <- syn[[4]][[2]]
  work_ST <- towork[[4]]
  
  # run
  syn <- synthetic_new_attribute(df= test_micro, prob_name= "p", attr_name= "transit",
                                 conditional_vars= c("emp_status", "age"), sym_tbl= work_ST)
  
  # test output
  expect_equal(sum(syn$p), 1)
  expect_equal(tapply(syn$p, syn$emp_status, sum),
               tapply(test_micro$p, test_micro$emp_status, sum))
  expect_equal(tapply(syn$p, syn$age, sum),
               tapply(test_micro$p, test_micro$age, sum))
  expect_equal(tapply(syn$p, syn$race, sum),
               tapply(test_micro$p, test_micro$race, sum))
  expect_equal(ncol(test_micro), ncol(syn) - 1)
})

test_that("parallel - real DF, fake ST", {
  # create test inputs
  # load(...) 
  # set up symbol tables
  levels <- c("A", "B", "C")
  ST <- data.frame(marital_status= rep(levels(syn[[1]][[2]]$marital_status), each= 7),
                   race= rep(levels(syn[[2]][[2]]$race), 5))
  ST <- do.call("rbind", replicate(3, ST, simplify=FALSE))
  ST <- ST[order(ST$marital_status, ST$race),]
  ST$pct <- rep(c(0.33, 0.34, 0.33), 35)
  ST$levels <- rep(levels, 35)
  
  st_list <- replicate(4, ST, simplify= FALSE)
  
  # run
  syn2 <- all_geog_synthetic_new_attribute(syn, prob_name= "p", attr_name= "variable",
                                           conditional_vars= c("marital_status", "race"),
                                           st_list= st_list)
  # test 
  expect_equal(class(syn2), c("synthACS", "list"))
  expect_true(is.synthACS(syn2))
  expect_true(all(unlist(lapply(syn2, function(l) is.micro_synthetic(l[[2]])))))
  expect_true(all.equal(unlist(lapply(syn2, function(l) sum(l[[2]]$p))), rep(1, 4), check.attributes = FALSE)) 
  
  expect_equal(lapply(syn2, function(l) {tapply(l[[2]]$p, l[[2]]$marital_status, sum)}),
               lapply(syn,  function(l) {tapply(l[[2]]$p, l[[2]]$marital_status, sum)}) )
  expect_equal(lapply(syn2, function(l) {tapply(l[[2]]$p, l[[2]]$race, sum)}),
               lapply(syn,  function(l) {tapply(l[[2]]$p, l[[2]]$race, sum)}) )
  
  expect_true(all.equal(lapply(syn2, function(l) {as.vector(tapply(l[[2]]$p, l[[2]]$variable, sum))}), 
                        replicate(4, c(0.33, 0.34, 0.33), simplify = FALSE), check.attributes = FALSE))
  
})