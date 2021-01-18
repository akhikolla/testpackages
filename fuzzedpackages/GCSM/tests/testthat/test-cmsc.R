test_that("cmsc() equals to cmsc_r()", {
  x <- rnorm(9)
  y <- rnorm(9)
  expect_equal(cmsc(x, y, rescale = FALSE), GCSM:::cmsc_r(x, y, rescale = FALSE))
  expect_equal(cmsc(x, y, rescale = TRUE),  GCSM:::cmsc_r(x, y, rescale = TRUE))
})

test_that("cmsc_sw() at the center equals to cmsc()", {
  x <- matrix(rnorm(9), nrow = 3)
  y <- matrix(rnorm(9), nrow = 3)
  expect_equal(cmsc(x, y, rescale = FALSE), cmsc_sw(x, y, rescale = FALSE, ksize=3)[2,2])
  expect_equal(cmsc(x, y, rescale = TRUE),  cmsc_sw(x, y, rescale = TRUE,  ksize=3)[2,2])
})

test_that("cmsc_tw() equals to for loop of cmsc()", {
  d <- c(3, 3, 9)
  x <- array(rnorm(prod(d)), d)
  y <- array(rnorm(prod(d)), d)
  xmin <- min(x, na.rm = T)
  xmax <- max(x, na.rm = T)
  ymin <- min(y, na.rm = T)
  ymax <- max(y, na.rm = T)
  s <- matrix(nrow = nrow(x), ncol = ncol(x))
  s_tw <- matrix(nrow = nrow(x), ncol = ncol(x))
  for (rn in 1:nrow(x)) {
    for (cn in 1:ncol(x)) {
      s[rn, cn] <- cmsc(x[rn, cn, ], y[rn, cn, ], xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
    }
  }
  s_tw <- cmsc_tw(x, y, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
  expect_equal(s, s_tw)
})
