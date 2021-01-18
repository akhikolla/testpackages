## Hua created at Feb 27 2017
## Test the function "add.house.noise" in "house_noise_model.R".

library(testthat)
library(FunChisq)

context("Testing house and candle noise models")

test_that("Testing the function of adding house noise", {

  # Test house noise model
  t1 <- matrix(c(1,4,1,3,3,
                1,4,4,5,4,
                5,2,2,2,3,
                4,2,4,5,1,
                3,4,1,2,1),
              nrow = 5, byrow = TRUE)
  t1.XY <- add.house.noise(t1, 0.3, 0)
  t1.X <- add.house.noise(t1, 0.4, 2)
  t1.Y <- add.house.noise(t1, 0.5, 1)
  expect_identical(sum(t1), sum(t1.XY))
  expect_identical(colSums(t1), colSums(t1.X))
  expect_identical(rowSums(t1), rowSums(t1.Y))


  t2 <- matrix(c(1,4,2,3,
                 5,1,3,4,
                 2,2,3,4,
                 5,2,4,5),
               nrow = 4, byrow = TRUE)
  t2.XY <- add.house.noise(t2, 0.3, 0)
  t2.X <- add.house.noise(t2, 0.4, 2)
  t2.Y <- add.house.noise(t2, 0.5, 1)
  expect_identical(sum(t2), sum(t2.XY))
  expect_identical(colSums(t2), colSums(t2.X))
  expect_identical(rowSums(t2), rowSums(t2.Y))

  ts <- list(t1, 2*t1, 3*t1, t2, 2*t2, 3*t2)

  ts.XY <- add.house.noise(ts, 0.3, 0)
  ts.X <- add.house.noise(ts, 0.4, 2)
  ts.Y <- add.house.noise(ts, 0.5, 1)

  for(i in c(1:length(ts))){
    expect_identical(sum(ts[[i]]), sum(ts.XY[[i]]))
    expect_identical(colSums(ts[[i]]), colSums(ts.X[[i]]))
    expect_identical(rowSums(ts[[i]]), rowSums(ts.Y[[i]]))
  }


  # Test candle noise model
  t1 <- matrix(c(1,4,1,3,3,
                 1,4,4,5,4,
                 5,2,2,2,3,
                 4,2,4,5,1,
                 3,4,1,2,1),
               nrow = 5, byrow = TRUE)
  t1.XY <- add.candle.noise(t1, 0.3, 0)
  t1.X <- add.candle.noise(t1, 0.4, 2)
  t1.Y <- add.candle.noise(t1, 0.5, 1)
  expect_identical(sum(t1), sum(t1.XY))
  expect_identical(colSums(t1), colSums(t1.X))
  expect_identical(rowSums(t1), rowSums(t1.Y))


  t2 <- matrix(c(1,4,2,3,
                 5,1,3,4,
                 2,2,3,4,
                 5,2,4,5),
               nrow = 4, byrow = TRUE)
  t2.XY <- add.candle.noise(t2, 0.3, 0)
  t2.X <- add.candle.noise(t2, 0.4, 2)
  t2.Y <- add.candle.noise(t2, 0.5, 1)
  expect_identical(sum(t2), sum(t2.XY))
  expect_identical(colSums(t2), colSums(t2.X))
  expect_identical(rowSums(t2), rowSums(t2.Y))



  ts <- list(t1, 2*t1, 3*t1, t2, 2*t2, 3*t2)

  ts.XY <- add.candle.noise(ts, 0.3, 0)
  ts.X <- add.candle.noise(ts, 0.4, 2)
  ts.Y <- add.candle.noise(ts, 0.5, 1)

  for(i in c(1:length(ts))){
    expect_identical(sum(ts[[i]]), sum(ts.XY[[i]]))
    expect_identical(colSums(ts[[i]]), colSums(ts.X[[i]]))
    expect_identical(rowSums(ts[[i]]), rowSums(ts.Y[[i]]))
  }
})
