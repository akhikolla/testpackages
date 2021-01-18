test_that("gcsm() equals to gcsm_r()", {
  x <- rnorm(9)
  y <- rnorm(9)
  expect_equal(gcsm(x, y, rescale = FALSE), GCSM:::gcsm_r(x, y, rescale = FALSE))
  expect_equal(gcsm(x, y, rescale = TRUE),  GCSM:::gcsm_r(x, y, rescale = TRUE))
})

test_that("gcsm_sw() at the center equals to gcsm()", {
  x <- matrix(rnorm(9), nrow = 3)
  y <- matrix(rnorm(9), nrow = 3)
  expect_equal(gcsm(x, y, rescale = FALSE), gcsm_sw(x, y, rescale = FALSE, ksize=3)[2,2])
  expect_equal(gcsm(x, y, rescale = TRUE),  gcsm_sw(x, y, rescale = TRUE,  ksize=3)[2,2])
})

test_that("gcsm_tw() equals to for loop of gcsm()", {
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
      s[rn, cn] <- gcsm(x[rn, cn, ], y[rn, cn, ], xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
    }
  }
  s_tw <- gcsm_tw(x, y, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
  expect_equal(s, s_tw)
})
