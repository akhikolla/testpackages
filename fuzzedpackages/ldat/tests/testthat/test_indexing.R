
context("Indexing")

test_that("indexing works", {
  x <- as_lvec(1:10)
  x_r <- as_rvec(x)

  x[4] <- 99
  x_r[4] <- 99
  expect_that(as_rvec(x), equals(x_r))
  expect_that(as_rvec(x[c(3,4)]), equals(x_r[c(3,4)]))
  expect_that(as_rvec(x[TRUE]), equals(x_r[TRUE]))

  x[1:4] <- 109
  x_r[1:4] <- 109
  expect_that(as_rvec(x), equals(x_r))
  expect_that(as_rvec(x[4]), equals(x_r[4]))
  expect_that(as_rvec(x[NA]), equals(x_r[NA]))

  x[range = c(3,6)] <- 299
  x_r[3:6] <- 299
  expect_that(as_rvec(x), equals(x_r))
  expect_that(as_rvec(x[range = c(2,3)]), equals(x_r[c(2,3)]))
})

test_that("indexing of ldat object", {
  x_r <- data.frame(a = 1:3, b = letters[1:3], stringsAsFactors = TRUE)
  x <- as_ldat(x_r)

  expect_equal(as.data.frame(x[]), x_r)
  expect_equal(as.data.frame(x[,]), x_r)
  
  expect_equal(as.data.frame(x[1:2, ]), x_r[1:2, ])  
  expect_equal(as.data.frame(x[1, ]), x_r[1, ])  
  
  expect_equal(as.data.frame(x[1:2]), x_r[1:2])  
  expect_equal(as.data.frame(x[1]), x_r[1])
  
  expect_equal(as.data.frame(x[1:2]), x_r[, 1:2, drop = FALSE])  
  expect_equal(as.data.frame(x[1]), x_r[, 1, drop = FALSE])  

  expect_equal(as.data.frame(x[1:2, drop = FALSE]), x_r[1:2])  
  expect_equal(as.data.frame(x[1, drop = FALSE]), x_r[1]) 
  
  # The following gave an error in version 0.2.0 with prerelease of R 4.0.0
  expect_equal(as.data.frame(x[1:2, , drop = FALSE]), x_r[1:2, , drop = FALSE])  
  expect_equal(as.data.frame(x[1, , drop = FALSE]), x_r[1, , drop = FALSE])  
  
  
  # Same as above with clone argument
  expect_equal(as.data.frame(x[clone = FALSE]), x_r)
  expect_equal(as.data.frame(x[, , clone = FALSE]), x_r)
  
  expect_equal(as.data.frame(x[1:2, , clone = FALSE]), x_r[1:2, ])  
  expect_equal(as.data.frame(x[1, , clone = FALSE]), x_r[1, ])  
  
  expect_equal(as.data.frame(x[1:2, clone = FALSE]), x_r[1:2])  
  expect_equal(as.data.frame(x[1, clone = FALSE]), x_r[1])
  
  expect_equal(as.data.frame(x[1:2, clone = FALSE]), x_r[, 1:2, drop = FALSE])  
  expect_equal(as.data.frame(x[1, clone = FALSE]), x_r[, 1, drop = FALSE])  
  
  expect_equal(as.data.frame(x[1:2, drop = FALSE, clone = FALSE]), x_r[1:2])  
  expect_equal(as.data.frame(x[1, drop = FALSE, clone = FALSE]), x_r[1]) 

  
  # drop = TRUE doesn't do anything
  expect_warning(x[1:2, , drop = TRUE])
  expect_warning(x[1, drop = TRUE])
  # Select non existent column
  expect_error(x[1:3])
})
