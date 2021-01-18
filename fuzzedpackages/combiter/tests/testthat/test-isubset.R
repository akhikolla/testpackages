library(testthat)
library(combiter)
context("subset iterator")

test_that("isubset goes through 2^n values", {
  for (n in 1:5)
  {
    x <- isubset(n)
    ct <- 0
    while (hasNext(x))
    {
      ct <- ct + 1
      nextElem(x)
    }
    expect_equal(ct, 2^n)

    # backward
    x <- isubset(n)
    ct <- 0
    while (hasPrev(x))
    {
      ct <- ct + 1
      prevElem(x)
    }
    expect_equal(ct, 2^n)
  }
})


test_that("all results of isubset are subsets of 1:n", {
  for (n in 1:5)
  {
    x <- isubset(n)
    while (hasNext(x))
    {
      i <- nextElem(x)
      expect_true(all(i %in% 1:n), paste0(i, collase=","))
    }

    # do the same for backward
    x <- isubset(n)
    while (hasPrev(x))
    {
      i <- prevElem(x)
      expect_true(all(i %in% 1:n), paste0(i, collase=","))
    }
  }
})


test_that("isubset elements are ordered by size first, and lexicographically if the sizes are equal", {
  lexico_smaller <- function(a, b)
  {
    # check if a < b lexicograpically
    # assumes that a and b are vectors of the same length
    index <- c(which(a > b), which(a < b))
    if (length(index) == 0L) return(FALSE) # all elements are equal
    return(a[min(index)] < b[min(index)])
  }

  subset_smaller <- function(a, b)
  {
    if (length(a) > length(b)) return(FALSE)
    if (length(a) < length(b)) return(TRUE)
    lexico_smaller(a, b)
  }

  for (n in 1:5)
  {
    x <- isubset(n)
    i <- NULL
    while (hasNext(x))
    {
      j <- nextElem(x)
      # requires i < j, but check only when i is not NULL
      if (!is.null(i)) {
        expect_true(subset_smaller(i, j))
      }
      i <- j
    }

    # backward
    x <- iperm(n)
    i <- NULL
    while (hasPrev(x))
    {
      j <- prevElem(x)
      # requires j < i
      if (!is.null(i)) {
        expect_true(subset_smaller(j, i))
      }
      i <- j
    }
  }
})


test_that("isubset rejects invalid elements", {
  expect_error(iperm(0))
  expect_error(iperm(-4))
  expect_error(iperm(1:2))
  expect_error(iperm(1.5))
  expect_error(iperm(3.0000000001))
})


