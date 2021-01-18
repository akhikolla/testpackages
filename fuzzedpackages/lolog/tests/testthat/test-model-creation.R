context("Model Creation")
library(testthat)
#library(lolog)
library(ergm)


test_that("terms", {
  data(sampson)
  mod <- createCppModel(samplike ~ edges() + star(2:4) | 1:18)
  expect_true(class(mod) == "Rcpp_DirectedModel")
  f <- function() {
    b <- 1:18
    form <- a  ~ edges() + star(2:4) | 1:18
    form
  }
  terms <- lolog:::.prepModelTerms(f())
  expect_identical(terms, structure(
    list(
      stats = structure(list(
        star = list(2:4), edges = list()
      ), .Names = c("star",
                    "edges")),
      offsets = list(),
      vertexOrder = 1:18
    ),
    .Names = c("stats",
               "offsets", "vertexOrder")
  ))
  data(flo)
  flomarriage <- network(flo, directed = FALSE)
  mod <-
    createCppModel(flomarriage ~ edges() + star(2:4) + preferentialAttachment())
  expect_identical(mod$isIndependent(TRUE, TRUE),
                   c(TRUE, FALSE, FALSE, FALSE, FALSE))
  expect_identical(mod$isIndependent(FALSE, TRUE),
                   c(TRUE, TRUE, TRUE, TRUE, FALSE))
  
})
