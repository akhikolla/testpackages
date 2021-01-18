library(testthat)
library(combiter)
context("combination iterator")

test_that("icomb goes through n C k values", {
  for (n in 1:5)
  {
    for (k in 1:n)
    {
      x <- icomb(n, k)
      ct <- 0
      while (hasNext(x))
      {
        ct <- ct + 1
        nextElem(x)
      }
      expect_equal(ct, choose(n, k))

      # backward
      x <- icomb(n, k)
      ct <- 0
      while (hasPrev(x))
      {
        ct <- ct + 1
        prevElem(x)
      }
      expect_equal(ct, choose(n, k))
    }
  }
})


test_that("icomb covers all combinations", {
  for (n in 1:5)
  {
    for (k in 1:n)
    {
      x <- icomb(n, k)
      allCombs <- combinat::combn(n, k, simplify = FALSE)
      while (hasNext(x))
      {
        i <- nextElem(x)
        expect_false(is.na(match(list(i), allCombs)))
      }

      # do the same for backward
      x <- icomb(n, k)
      while (hasPrev(x))
      {
        i <- prevElem(x)
        expect_false(is.na(match(list(i), allCombs)))
      }
    }
  }
})


test_that("icomb elements are ordered lexicographically", {
  lexico_smaller <- function(a, b)
  {
    # check if a < b lexicograpically
    # assumes that a and b are vectors of the same length
    index <- c(which(a > b), which(a < b))
    if (length(index) == 0L) return(FALSE) # all elements are equal
    return(a[min(index)] < b[min(index)])
  }

  for (n in 1:5)
  {
    for (k in 1:n)
    {
      x <- icomb(n, k)
      i <- NULL
      while (hasNext(x))
      {
        j <- nextElem(x)
        # requires i < j, but check only when i is not NULL
        if (!is.null(i)) {
          expect_true(lexico_smaller(i, j),
                      paste(paste0(i, collapse=" "),
                            paste0(j, collase=" "), "?<"))
        }
        i <- j
      }

      # backward
      x <- icomb(n, k)
      i <- NULL
      while (hasPrev(x))
      {
        j <- prevElem(x)
        # requires j < i
        if (!is.null(i)) {
          expect_true(lexico_smaller(j, i),
                      paste(paste0(j, collapse=" "),
                            paste0(i, collase=" "), "?<"))
        }
        i <- j
      }
    }
  }
})


test_that("icomb rejects invalid elements", {
  expect_error(icomb(0, 3))
  expect_error(icomb(-3, 3))
  expect_error(icomb(1:2, 3))
  expect_error(icomb(1.5, 3))
  expect_error(icomb(3.0000000001, 3))

  expect_error(icomb(10, 0))
  expect_error(icomb(10, -5))
  expect_error(icomb(10, 1:2))
  expect_error(icomb(10, 1.5))
  expect_error(icomb(10, 3.000000001))
})


