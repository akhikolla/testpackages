library(testthat)
library(combiter)
context("permutation iterator")

test_that("iperm goes through n P k values", {
  for (n in 1:4)
  {
    for (k in 1:n)
    {
      x <- iperm(n, k)
      ct <- 0
      while (hasNext(x))
      {
        ct <- ct + 1
        nextElem(x)
      }
      expect_equal(ct, choose(n,k)*factorial(k))


      # backward
      x <- iperm(n, k)
      ct <- 0
      while (hasPrev(x))
      {
        ct <- ct + 1
        prevElem(x)
      }
      expect_equal(ct, choose(n,k)*factorial(k))
    }
  }
})


test_that("iperm covers all permutations", {
  for (n in 1:4)
  {
    for (k in 1:n)
    {
      x <- iperm(n, k)
      allPerms <- combinat::permn(n) %>%
        lapply(`[`, 1:k) %>% unique()
      while (hasNext(x))
      {
        i <- nextElem(x)
        expect_false(is.na(match(list(i), allPerms)),
                     paste0(i, collapse=" "))
      }

      # do the same for backward
      x <- iperm(n, k)
      while (hasPrev(x))
      {
        i <- prevElem(x)
        expect_false(is.na(match(list(i), allPerms)),
                     paste0(i, collapse=" "))
      }
    }
  }
})


test_that("iperm elements are ordered lexicographically", {
  lexico_smaller <- function(a, b)
  {
    # check if a < b lexicograpically
    # assumes that a and b are vectors of the same length
    index <- c(which(a > b), which(a < b))
    if (length(index) == 0L) return(FALSE) # all elements are equal
    return(a[min(index)] < b[min(index)])
  }

  for (n in 1:4)
  {
    for (k in 1:n)
    {
      x <- iperm(n,k)
      i <- NULL
      while (hasNext(x))
      {
        j <- nextElem(x)
        # requires i < j, but check only when i is not NULL
        if (!is.null(i)) {
          expect_true(lexico_smaller(i, j))
        }
        i <- j
      }

      # backward
      x <- iperm(n,k)
      i <- NULL
      while (hasPrev(x))
      {
        j <- prevElem(x)
        # requires j < i
        if (!is.null(i)) {
          expect_true(lexico_smaller(j, i))
        }
        i <- j
      }
    }
  }
})


test_that("iperm rejects invalid elements", {
  expect_error(iperm(0))
  expect_error(iperm(-3))
  expect_error(iperm(1:2))
  expect_error(iperm(1.5))
  expect_error(iperm(3.0000000001))

  expect_error(iperm(3, 4))
  expect_error(iperm(3, -1))
  expect_error(iperm(3, 2.05))
  expect_error(iperm(3, 1.0000001))

})


