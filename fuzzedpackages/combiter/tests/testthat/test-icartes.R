library(testthat)
library(combiter)

context("cartesian product iterator")

test_that("icartes goes through product of maxes", {

  for (n1 in 1:3) {
    for (n2 in 1:3) {
      for (n3 in 1:3)
      {
        nvec <- c(n1,n2,n3)
        x <- icartes(nvec)
        ct <- 0
        while (hasNext(x))
        {
          ct <- ct + 1
          nextElem(x)
        }
        expect_equal(ct, prod(nvec))

        # backward
        x <- icartes(nvec)
        ct <- 0
        while (hasPrev(x))
        {
          ct <- ct + 1
          prevElem(x)
        }
        expect_equal(ct, prod(nvec))
      }
    }
  }
})


test_that("icartes covers all cartesian product", {

  for (n1 in 1:3) {
    for (n2 in 1:3) {
      for (n3 in 1:3)
      {
        all_elems <- expand.grid(1:n1, 1:n2, 1:n3)
        all_elems <- lapply(1:nrow(all_elems), function(i) {
          as.integer(all_elems[i, ]) })
        nvec <- c(n1,n2,n3)
        x <- icartes(nvec)
        while (hasNext(x))
        {
          i <- nextElem(x)
          expect_false(is.na(match(list(i), all_elems)),
                       paste0(i, collapse=" "))
        }

      }
    }
  }
})


test_that("icartes elements are ordered lexicographically", {
  lexico_smaller <- function(a, b)
  {
    # check if a < b lexicograpically
    # assumes that a and b are vectors of the same length
    index <- c(which(a > b), which(a < b))
    if (length(index) == 0L) return(FALSE) # all elements are equal
    return(a[min(index)] < b[min(index)])
  }

  for (n1 in 1:3) {
    for (n2 in 1:3) {
      for (n3 in 1:3)
      {
        nvec <- c(n1,n2,n3)
        x <- icartes(nvec)
        i <- NULL
        while (hasNext(x))
        {
          j <- nextElem(x)
          if (!is.null(i)) expect_true(lexico_smaller(i, j),
                                       paste(paste0(i, collapse=" "),
                                             paste0(j, collase=" "), "?<"))
          i <- j
        }

        # backward
        x <- icartes(nvec)
        i <- NULL
        while (hasPrev(x))
        {
          j <- prevElem(x)
          if (!is.null(i)) expect_true(lexico_smaller(j, i),
                                       paste(paste0(j, collapse=" "),
                                             paste0(i, collase=" "), "?<"))
          i <- j
        }
      }
    }
  }
})


test_that("icartes rejects invalid elements", {
  expect_error(icartes(0))
  expect_error(icartes(c(-1, 3)))
  expect_error(icartes(c(4, -1, 3)))
  expect_error(icartes())
  expect_error(icartes(c(1.5, 1:10)))
  expect_error(icartes(3.0000000001))
  expect_error(icartes(c(4, 5, 1.2)))
})


