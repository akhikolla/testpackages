context("dat_transform")

# set up input variable
VChar <- c("high", "high", "median", "low", "zero")
VChar <- factor(VChar)
VInt <- c(1, 4, 5, 8, 1)
VNum <- c(1.2, 3.4, 5.6, 7.8, 9.0)
VStr <- c("TRUE", "FALSE", "TRUE", "TRUE", "TRUE")
VLessChar <- c("X", "X", "Y", "Y", "Z")
VLessChar <- factor(VLessChar, levels = c("W", "X", "Y", "Z"))
data <- data.frame(VChar, VInt, VNum, VStr, VLessChar, stringsAsFactors = TRUE)

level_list <- vector("list", length = 4)
level_list[[1]] <- c("zero", "low", "median", "high")
level_list[[2]] <- c(1, 2, 3, 4, 5, 6, 7, 8)
level_list[[3]] <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)
level_list[[4]] <- c("TRUE", "FALSE")
level_list[[5]] <- c("W", "X", "Y", "Z")

# test
test_that("Check dat_transform runs as expected", {

  databn <- sparsebnUtils::sparsebnData(data, type = "discrete")
  # dat_transform runs as expected
  expect_error(dat_transform(databn), NA)


  databn <- sparsebnUtils::sparsebnData(data, type = "discrete", levels = level_list)
  expect_error(dat_transform(databn), NA)

  level_list[[2]] <- 1:4
  databn <- sparsebnUtils::sparsebnData(data, type = "discrete", levels = level_list)
  expect_error(dat_transform(databn), NA)
})

test_that("Check dat_transform throw error as expected", {
  ### dat_transform throw an error if data.frame contain column other than type numeric and factor
  databn <- sparsebnUtils::sparsebnData(data, type = "discrete")
  databn$data[[4]] <- VStr
  expect_error(dat_transform(databn))

  ### dat_transform throw an error if the input level is incorrect
  level_list[[2]] <- 1:3
  databn <- sparsebnUtils::sparsebnData(data, type = "discrete", levels = level_list)
  expect_error(dat_transform(databn))

})

test_that("Check result for dat_transform", {
  # result if do not specify levels
  databn <- sparsebnUtils::sparsebnData(data, type = "discrete")
  result <- dat_transform(databn)
  expect_equal(result[[1]], c(0, 0, 2, 1, 3))
  expect_equal(result[[2]], c(0, 1, 2, 3, 0))
  expect_equal(result[[3]], c(0, 1, 2, 3, 4))
  expect_equal(result[[4]], c(1, 0, 1, 1, 1))
  expect_equal(result[[5]], c(0, 0, 1, 1, 2))

  # result if specify levels
  databn <- sparsebnUtils::sparsebnData(data, levels = level_list, type = "discrete")
  result <- dat_transform(databn)
  expect_equal(result[[1]], c(0, 0, 2, 1, 3))
  expect_equal(result[[2]], c(0, 3, 4, 7, 0))
  expect_equal(result[[3]], c(0, 2, 4, 6, 8))
  expect_equal(result[[4]], c(1, 0, 1, 1, 1))
  expect_equal(result[[5]], c(0, 0, 1, 1, 2))
})
