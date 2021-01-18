library(testthat)
library(OptCirClust)
library(Ckmeans.1d.dp)

context("Framed optimal clustering on simple and random data")

test_that("Framed optimal clustering", {
  x <- c(-1, 2, 4, 5, 6)
  result <- Ckmeans.1d.dp(x, 3, method = "linear")
  output <- FramedClust(X=x, K = 3, frame.size = 5,
                        first.frame = 1, last.frame = 1,
                        method = "linear.polylog")

  expect_equal(output$tot.withinss, result$tot.withinss)




  x <- c(-.9, 1, 1.1, 1.9, 2, 2.1)
  result <- Ckmeans.1d.dp(x, 3, method = "linear")
  output <- FramedClust(X=x, K = 3, frame.size = 6, first.frame = 1,
                        last.frame = 1, method = "linear.polylog")


  expect_equal(output$tot.withinss, result$tot.withinss)




  x <- c(-1, 2,-1, 2, 4, 5, 6,-1, 2,-1)

  result <- Ckmeans.1d.dp(x, 3, method = "linear")
  output <- FramedClust(X=x, K = 3, frame.size = 10, first.frame = 1,
                        last.frame = 1, method = "linear.polylog")

  expect_equal(output$tot.withinss, result$tot.withinss)




  x <- c(3, 2,-5.4, 0.1)
  result <- Ckmeans.1d.dp(x, 4, method = "linear")
  output <- FramedClust(X=x, K = 4, frame.size = 4, first.frame = 1,
                        last.frame = 1, method = "linear.polylog")

  expect_equal(output$tot.withinss, result$tot.withinss)





  x <- 1:10
  result <- Ckmeans.1d.dp(x, 2, method = "linear")
  output <- FramedClust(X=x, K = 2, frame.size = 10, first.frame = 1,
                        last.frame = 1, method = "linear.polylog")

  expect_equal(output$tot.withinss, result$tot.withinss)





  x <- c(-2.5,-2.5,-2.5,-2.5)
  result <- Ckmeans.1d.dp(x, 1, method = "linear")
  output <- FramedClust(X=x, K = 1, frame.size = 4, first.frame = 1,
                        last.frame = 1, method = "linear.polylog")

  expect_equal(output$tot.withinss, result$tot.withinss)




  x <- rep(1, 100)
  result <- Ckmeans.1d.dp(x, 1, method = "linear")
  output <- FramedClust(X=x, K = 1, frame.size = 100, first.frame = 1,
                        last.frame = 1, method = "linear.polylog")

  expect_equal(output$tot.withinss, result$tot.withinss)



  x <-
    c(-3, 2.2,-6, 7, 9, 11,-6.3, 75, 82.6, 32.3,-9.5, 62.5, 7, 95.2)
  result <- Ckmeans.1d.dp(x, 8, method = "linear")
  output <- FramedClust(X=x, K = 8, frame.size = 14, first.frame = 1,
                        last.frame = 1, method = "linear.polylog")

  expect_equal(output$tot.withinss, result$tot.withinss)


  x <- c(3.5, 3.6, 3.7, 3.1, 1.1, 0.9, 0.8, 2.2, 1.9, 2.1)
  result <- Ckmeans.1d.dp(x, 3, method = "linear")
  output <- FramedClust(X=x, K = 3, frame.size = 10, first.frame = 1,
                        last.frame = 1, method = "linear.polylog")

  expect_equal(output$tot.withinss, result$tot.withinss)





  x <- cos((-10:10))
  result <- Ckmeans.1d.dp(x, 2, method = "linear")
  output <- FramedClust(X=x, K = 2, frame.size = 21, first.frame = 1,
                        last.frame = 1, method = "linear.polylog")

  expect_equal(output$tot.withinss, result$tot.withinss)



  x <- dgamma(seq(1, 10, by = 0.5), shape = 2, rate = 1)

  result <- Ckmeans.1d.dp(x, 3, method = "linear")
  output <- FramedClust(X=x, K = 3, frame.size = 19, first.frame = 1,
                        last.frame = 1, method = "linear.polylog")

  expect_equal(output$tot.withinss, result$tot.withinss)











  x <- c(-1, 2, 4, 5, 6)
  result <- FramedClust(X=x, K = 3, frame.size = 5, first.frame = 1,
                        last.frame = 1, method = "Ckmeans.1d.dp")
  output <- FramedClust(X=x, K = 3, frame.size = 5, first.frame = 1,
                        last.frame = 1, method = "linear.polylog")

  expect_equal(output$tot.withinss, result$tot.withinss)
  # expect_equal(output$Border, result$Border)



  x <- c(-.9, 1, 1.1, 1.9, 2, 2.1)
  result <- FramedClust(X=x, K = 3, frame.size = 6, first.frame = 1,
                        last.frame = 1, method = "Ckmeans.1d.dp")
  output <- FramedClust(X=x, K = 3, frame.size = 6, first.frame = 1,
                        last.frame = 1, method = "linear.polylog")

  expect_equal(output$tot.withinss, result$tot.withinss)
  #  expect_equal(output$Border, result$Border)



  x <- c(-1, 2,-1, 2, 4, 5, 6,-1, 2,-1)

  result <- FramedClust(X=x, K = 3, frame.size = 10, first.frame = 1,
                        last.frame = 1, method = "Ckmeans.1d.dp")
  output <- FramedClust(X=x, K = 3, frame.size = 10, first.frame = 1,
                        last.frame = 1, method = "linear.polylog")

  expect_equal(output$tot.withinss, result$tot.withinss)
  #  expect_equal(output$Border, result$Border)



  x <- c(3, 2,-5.4, 0.1)
  result <- FramedClust(X=x, K = 4, frame.size = 4, first.frame = 1,
                        last.frame = 1, method = "Ckmeans.1d.dp")
  output <- FramedClust(X=x, K = 4, frame.size = 4, first.frame = 1,
                        last.frame = 1, method = "linear.polylog")

  expect_equal(output$tot.withinss, result$tot.withinss)
  #  expect_equal(output$Border, result$Border)




  x <- 1:10
  result <- FramedClust(X=x, K = 2, frame.size = 10, first.frame = 1,
                        last.frame = 1, method = "Ckmeans.1d.dp")
  output <- FramedClust(X=x, K = 2, frame.size = 10, first.frame = 1,
                        last.frame = 1, method = "linear.polylog")

  expect_equal(output$tot.withinss, result$tot.withinss)
  #  expect_equal(output$Border, result$Border)




  x <- c(-2.5,-2.5,-2.5,-2.5)
  result <- FramedClust(X=x, K = 1, frame.size = 4, first.frame = 1,
                        last.frame = 1, method = "Ckmeans.1d.dp")
  output <- FramedClust(X=x, K = 1, frame.size = 4, first.frame = 1,
                        last.frame = 1, method = "linear.polylog")

  expect_equal(output$tot.withinss, result$tot.withinss)
  #  expect_equal(output$Border, result$Border)



  x <- rep(1, 100)
  result <- FramedClust(X=x, K = 1, frame.size = 100, first.frame = 1,
                        last.frame = 1, method = "Ckmeans.1d.dp")
  output <- FramedClust(X=x, K = 1, frame.size = 100, first.frame = 1,
                        last.frame = 1, method = "linear.polylog")

  expect_equal(output$tot.withinss, result$tot.withinss)
  #  expect_equal(output$Border, result$Border)


  x <-
    c(-3, 2.2,-6, 7, 9, 11,-6.3, 75, 82.6, 32.3,-9.5, 62.5, 7, 95.2)
  result <- FramedClust(X=x, K = 8, frame.size = 14, first.frame = 1,
                        last.frame = 1, method = "Ckmeans.1d.dp")
  output <- FramedClust(X=x, K = 8, frame.size = 14, first.frame = 1,
                        last.frame = 1, method = "linear.polylog")

  expect_equal(output$tot.withinss, result$tot.withinss)
  #  expect_equal(output$Border, result$Border)

  x <- c(3.5, 3.6, 3.7, 3.1, 1.1, 0.9, 0.8, 2.2, 1.9, 2.1)
  result <- FramedClust(X=x, K = 3, frame.size = 10, first.frame = 1,
                        last.frame = 1, method = "Ckmeans.1d.dp")
  output <- FramedClust(X=x, K = 3, frame.size = 10, first.frame = 1,
                        last.frame = 1, method = "linear.polylog")

  expect_equal(output$tot.withinss, result$tot.withinss)
  #  expect_equal(output$Border, result$Border)




  x <- cos((-10:10))
  result <- FramedClust(X=x, K = 2, frame.size = 21, first.frame = 1,
                        last.frame = 1, method = "Ckmeans.1d.dp")
  output <- FramedClust(X=x, K = 2, frame.size = 21, first.frame = 1,
                        last.frame = 1, method = "linear.polylog")

  expect_equal(output$tot.withinss, result$tot.withinss)
  # expect_equal(output$Border, result$Border)


  x <- dgamma(seq(1, 10, by = 0.5), shape = 2, rate = 1)

  result <- FramedClust(X=x, K = 3, frame.size = 19, first.frame = 1,
                        last.frame = 1, method = "Ckmeans.1d.dp")
  output <- FramedClust(X=x, K = 3, frame.size = 19, first.frame = 1,
                        last.frame = 1, method = "linear.polylog")

  expect_equal(output$tot.withinss, result$tot.withinss)
  #  expect_equal(output$Border, result$Border)

  x <-
    c(rnorm(50, sd = 0.3),
      rnorm(50, mean = 10, sd = 0.3),
      rnorm(50, mean = 20, sd = 0.3))

  result <- FramedClust(X=x, K = 3, frame.size = 150, first.frame = 1,
                        last.frame = 1, method = "Ckmeans.1d.dp")
  output <- FramedClust(X=x, K = 3, frame.size = 150, first.frame = 1,
                        last.frame = 1, method = "linear.polylog")

  expect_equal(output$tot.withinss, result$tot.withinss)
  #  expect_equal(output$Border, result$Border)


  x <-
    c(rnorm(50, sd = 0.3),
      rnorm(50, mean = 1, sd = 0.3),
      rnorm(50, mean = 2, sd = 0.3))

  result <- FramedClust(X=x, K = 3, frame.size = 150, first.frame = 1,
                        last.frame = 1, method = "Ckmeans.1d.dp")
  output <- FramedClust(X=x, K = 3, frame.size = 150, first.frame = 1,
                        last.frame = 1, method = "linear.polylog")

  expect_equal(output$tot.withinss, result$tot.withinss)
  #  expect_equal(output$Border, result$Border)

})




context("Framed optimal clustering on more examples")


test_that("Brute force vs fast optimal framed clustering", {
  x <- c(1:10)

  result <- FramedClust(X=x, K = 3, frame.size = 5, first.frame = 1,
                        last.frame = 5, method = "Ckmeans.1d.dp")
  output <- FramedClust(X=x, K = 3, frame.size = 5, first.frame = 1,
                        last.frame = 5, method = "linear.polylog")

  expect_equal(output$tot.withinss, result$tot.withinss)


  x <- c(1:100)


  result <- FramedClust(X=x, K = 3, frame.size = 5, first.frame = 1,
                        last.frame = 95, method = "Ckmeans.1d.dp")
  output <- FramedClust(X=x, K = 3, frame.size = 5, first.frame = 1,
                        last.frame = 95, method = "linear.polylog")

  expect_equal(output$tot.withinss, result$tot.withinss)

  x <- c(1.5, 2.9, 3.4, 5.6, 6.7, 6.8, 7.9, 8.2, 9.5, 10.56)

  result <- FramedClust(X=x, K = 3, frame.size = 5, first.frame = 1,
                        last.frame = 5, method = "Ckmeans.1d.dp")
  output <- FramedClust(X=x, K = 3, frame.size = 5, first.frame = 1,
                        last.frame = 5, method = "linear.polylog")

  expect_equal(output$tot.withinss, result$tot.withinss)



  x <-
    c(
      1.500,
      5.600,
      -7.453,
      -6.000,
      -4.230,
      8.200,
      -1.240,
      6.700,
      2.900,
      -8.924,
      9.500,
      0.250,
      -3.450,
      7.900,
      -2.450,
      10.560,
      -5.650,
      3.400,
      -9.342,
      6.800
    )


  result <- FramedClust(X=x, K = 3, frame.size = 5, first.frame = 1,
                        last.frame = 15, method = "Ckmeans.1d.dp")
  output <- FramedClust(X=x, K = 3, frame.size = 5, first.frame = 1,
                        last.frame = 15, method = "linear.polylog")

  expect_equal(output$tot.withinss, result$tot.withinss)




  x <- seq(1, 100, by = 2)


  result <- FramedClust(X=x, K = 3, frame.size = 5, first.frame = 1,
                        last.frame = 45, method = "Ckmeans.1d.dp")
  output <- FramedClust(X=x, K = 3, frame.size = 5, first.frame = 1,
                        last.frame = 45, method = "linear.polylog")

  expect_equal(output$tot.withinss, result$tot.withinss)


  x <- seq(1, 100, by = 0.5)


  result <- FramedClust(X=x, K = 3, frame.size = 5, first.frame = 1,
                        last.frame = 194, method = "Ckmeans.1d.dp")
  output <- FramedClust(X=x, K = 3, frame.size = 5, first.frame = 1,
                        last.frame = 194, method = "linear.polylog")

  expect_equal(output$tot.withinss, result$tot.withinss)



  x <- seq(1, 100, by = 0.5)


  result <- FramedClust(X=x, K = 10, frame.size = length(x) / 2 , first.frame = 1,
                        last.frame = length(x) /2, method = "Ckmeans.1d.dp")
  output <- FramedClust(X=x, K = 10, frame.size = length(x) / 2 , first.frame = 1,
                        last.frame = length(x) /2, method = "linear.polylog")

  expect_equal(output$tot.withinss, result$tot.withinss)

  x <- seq(1, 1000, by = 0.5)


  result <- FramedClust(X=x, K = 10, frame.size = length(x) / 2 , first.frame = 1,
                        last.frame = length(x) - length(x) /2 - 1, method = "Ckmeans.1d.dp")
  output <- FramedClust(X=x, K = 10, frame.size = length(x) / 2 , first.frame = 1,
                        last.frame = length(x) - length(x) /2 -1, method = "linear.polylog")

  expect_equal(output$tot.withinss, result$tot.withinss)
})

context("Circular clustering")

test_that("Circular clustering", {
  x <-
    c(rnorm(50, sd = 0.3),
      rnorm(50, mean = 10, sd = 0.3),
      rnorm(50, mean = 20, sd = 0.3))

  x1 <- c(x, (x + 40))

  result <- FramedClust(X=x1, K = 3, frame.size = length(x) / 2, first.frame = 1,
                        last.frame =  length(x) -  length(x) /2, method = "Ckmeans.1d.dp")
  output <- FramedClust(X=x1, K = 3, frame.size = length(x) / 2, first.frame = 1,
                        last.frame =  length(x) -  length(x) /2, method = "linear.polylog")

  expect_equal(output$tot.withinss, result$tot.withinss)


  x <-
    c(rnorm(15, sd = 0.3),
      rnorm(45, mean = 10, sd = 0.3),
      rnorm(70, mean = 20, sd = 0.3))

  x1 <- c(x, (x + 40))


  result <- FramedClust(X=x1, K = 3, frame.size = length(x1) / 2, first.frame = 1,
                        last.frame = length(x1) - length(x1) /2, method = "Ckmeans.1d.dp")
  output <- FramedClust(X=x1, K = 3, frame.size = length(x1) / 2, first.frame = 1,
                        last.frame = length(x1) - length(x1) /2, method = "linear.polylog")

  expect_equal(output$tot.withinss, result$tot.withinss)


  x <-
    c(rnorm(50, sd = 0.3),
      rnorm(50, mean = 1, sd = 0.3),
      rnorm(50, mean = 2, sd = 0.3))

  x1 <- c(x, (x + 30))


  result <- FramedClust(X=x1, K = 3, frame.size = length(x1) / 2,
                        first.frame = 1, last.frame = length(x1) - length(x1) /2, method = "Ckmeans.1d.dp")
  output <- FramedClust(X=x1, K = 3, frame.size = length(x1) / 2,
                        first.frame = 1, last.frame = length(x1) - length(x1) /2, method = "linear.polylog")

  expect_equal(output$tot.withinss, result$tot.withinss)




  x <- c(1, 2, 10, 11, 12, 17, 18, 19, 30, 31)

  result <- CirClust(x, 3, Circumference = 32, method = "BOCC")
  output <- CirClust(x, 3, Circumference = 32)
  expect_equal(output$tot.withinss, result$tot.withinss)

  expect_equal(output$cluster[1], output$cluster[10])
  expect_equal(output$cluster[2], output$cluster[9])
  expect_equal(output$cluster[3], output$cluster[5])
  expect_equal(output$cluster[6], output$cluster[8])

  x <- c(1, 2, 10, 11, 12, 13, 17, 18, 19, 20, 30, 31)



  result <- CirClust(x, 3, Circumference = 31, method = "BOCC")
  output <- CirClust(x, 3, Circumference = 31)
  expect_equal(output$tot.withinss, result$tot.withinss)

  expect_equal(output$cluster[1], output$cluster[12])
  expect_equal(output$cluster[2], output$cluster[11])
  expect_equal(output$cluster[3], output$cluster[6])
  expect_equal(output$cluster[4], output$cluster[5])
  expect_equal(output$cluster[7], output$cluster[10])
  expect_equal(output$cluster[8], output$cluster[9])


  x <- c(1, 2, 3, 10, 11, 12, 13, 17, 18, 19, 20)


  result <- CirClust(x, 3, Circumference = 31, method = "BOCC")
  output <- CirClust(x, 3, Circumference = 31)
  expect_equal(output$tot.withinss, result$tot.withinss)

  expect_equal(output$cluster[1], output$cluster[3])
  expect_equal(output$cluster[2], output$cluster[3])
  expect_equal(output$cluster[4], output$cluster[7])
  expect_equal(output$cluster[5], output$cluster[6])
  expect_equal(output$cluster[8], output$cluster[11])
  expect_equal(output$cluster[9], output$cluster[10])


  x <-
    c(1, 2, 10, 11, 12, 13, 14, 15, 27, 28, 29, 30, 31, 32, 40, 41)


  result <- CirClust(x, 3, Circumference = 42, method = "BOCC")
  output <- CirClust(x, 3, Circumference = 42)
  expect_equal(output$tot.withinss, result$tot.withinss)

  expect_equal(output$cluster[1], output$cluster[16])
  expect_equal(output$cluster[2], output$cluster[15])
  expect_equal(output$cluster[3], output$cluster[8])
  expect_equal(output$cluster[9], output$cluster[14])



  x <-
    c(rnorm(25, sd = 0.3),
      rnorm(50, mean = 1, sd = 0.3),
      rnorm(100, mean = 2, sd = 0.3))


  result <- CirClust(x, 3, Circumference = 6, method = "BOCC")
  output <- CirClust(x, 3, Circumference = 6)
  expect_equal(output$tot.withinss, result$tot.withinss)


  x <-
    c(rnorm(25, sd = 0.3),
      rnorm(50, mean = 100, sd = 0.3),
      rnorm(100, mean = 200, sd = 0.3))


  result <- CirClust(x, 3, Circumference = 210, method = "BOCC")
  output <- CirClust(x, 3, Circumference = 210)
  expect_equal(output$tot.withinss, result$tot.withinss)


  x <-
    c(rnorm(25, sd = 0.3),
      rnorm(50, mean = 100, sd = 0.3),
      rnorm(100, mean = 200, sd = 0.3))


  result <- CirClust(x, 3, Circumference = 210, method = "BOCC")
  output <- CirClust(x, 3, Circumference = 210)
  expect_equal(output$tot.withinss, result$tot.withinss)


  x <-
    c(rnorm(25, sd = 3),
      rnorm(50, mean = 100, sd = 30),
      rnorm(100, mean = 200, sd = 60))


  result <- CirClust(x, 3, Circumference = 210, method = "BOCC")
  output <- CirClust(x, 3, Circumference = 210)
  expect_equal(output$tot.withinss, result$tot.withinss)


  x <-
    c(
      rnorm(25, sd = 3),
      rnorm(50, mean = 100, sd = 30),
      rnorm(100, mean = 200, sd = 60),
      rnorm(65, mean = 50, sd = 3),
      rnorm(110, mean = 150, sd = 30),
      rnorm(160, mean = 250, sd = 60)
    )


  result <- CirClust(x, 3, Circumference = 300, method = "BOCC")
  output <- CirClust(x, 3, Circumference = 300)
  expect_equal(output$tot.withinss, result$tot.withinss)

  x <-
    c(
      rnorm(25, sd = 3),
      rbinom(50, size = 5, prob = 0.5),
      rnorm(100, mean = 200, sd = 60),
      rbinom(65, size = 5, prob = 0.5),
      rnorm(110, mean = 150, sd = 30),
      rbinom(160, size = 5, prob = 0.5)
    )


  result <- CirClust(x, 3, Circumference = 300, method = "BOCC")
  output <- CirClust(x, 3, Circumference = 300)
  expect_equal(output$tot.withinss, result$tot.withinss)


  x <-
    c(
      rcauchy(25, location = 0, scale = 1),
      rbinom(50, size = 5, prob = 0.5),
      rcauchy(100, location = 0, scale = 1),
      rbinom(65, size = 5, prob = 0.5),
      rcauchy(110, location = 0, scale = 1),
      rbinom(160, size = 5, prob = 0.5)
    )


  result <- CirClust(x, 3, Circumference = 300, method = "BOCC")
  output <- CirClust(x, 3, Circumference = 300)
  expect_equal(output$tot.withinss, result$tot.withinss)


  x <-
    c(
      rnbinom(25, size = 10,  mu = 0),
      rnbinom(50, size = 10, mu = 10),
      rcauchy(100, location = 0, scale = 1),
      rnbinom(65, size = 10, mu = 20),
      rcauchy(110, location = 0, scale = 1),
      rnbinom(160, size = 10, mu = 30)
    )


  result <- CirClust(x, 3, Circumference = 300, method = "BOCC")
  output <- CirClust(x, 3, Circumference = 300)
  expect_equal(output$tot.withinss, result$tot.withinss)




  x <-
    c(
      rep(10, 100),
      rnbinom(50, size = 10, mu = 10),
      rcauchy(100, location = 0, scale = 1),
      rnbinom(65, size = 10, mu = 20),
      rcauchy(110, location = 0, scale = 1),
      rnbinom(160, size = 10, mu = 30)
    )


  result <- CirClust(x, 3, Circumference = 300, method = "BOCC")
  output <- CirClust(x, 3, Circumference = 300)
  expect_equal(output$tot.withinss, result$tot.withinss)




  x <- c(rep(1, 10), rep(10.45, 20), 30.5, 31.25, 32.75)

  result <- CirClust(x, 3, Circumference = 33, method = "BOCC")
  output <- CirClust(x, 3, Circumference = 33)
  expect_equal(output$tot.withinss, result$tot.withinss)



  x <-
    c(rep(1, 10),
      rep(10.45, 20),
      30.5,
      31.25,
      32.75,
      34.25,
      36.34,
      40,
      51,
      60)

  result <- CirClust(x, 3, Circumference = 33, method = "BOCC")
  output <- CirClust(x, 3, Circumference = 33)
  expect_equal(output$tot.withinss, result$tot.withinss)



  x <-
    c(
      rep(1, 10),
      rep(10.45, 20),
      30.5,
      31.25,
      32.75,
      34.25,
      36.34,
      40,
      51,
      60,
      1000,
      2000,
      3000,
      4000,
      5000
    )

  result <- CirClust(x, 3, Circumference = 33, method = "BOCC")
  output <- CirClust(x, 3, Circumference = 33)
  expect_equal(output$tot.withinss, result$tot.withinss)



  x <-
    c(
      -1000,
      -550.343,-450.54,-1.24,-2.45,-3.45,
      5.64,-6.34,
      rep(-10, 10),
      rep(1, 10),
      rep(10.45, 20),
      30.5,
      31.25,
      32.75,
      34.25,
      36.34,
      40,
      51,
      60,
      1000,
      2000,
      3000,
      4000,
      5000
    )

  result <- CirClust(x, 3, Circumference = 33, method = "BOCC")
  output <- CirClust(x, 3, Circumference = 33)
  expect_equal(output$tot.withinss, result$tot.withinss)


})

test_that("Framed clustering", {

  result <- FramedClust(1:10, 2, 4, 2, 7)
  output <- FramedClust(1:10, 2, 4, 2, 7,"Ckmeans.1d.dp")
  expect_equal(output$tot.withinss, result$tot.withinss)

  result <- FramedClust(1:10, 2, 4, 2, 6)
  output <- FramedClust(1:10, 2, 4, 2, 6,"Ckmeans.1d.dp")
  expect_equal(output$tot.withinss, result$tot.withinss)
})
