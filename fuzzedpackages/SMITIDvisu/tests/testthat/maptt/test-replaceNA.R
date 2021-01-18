context("Testing replaceNA")

test_that("Testing replaceNA with an empty vector", {
	vector <- c()
	result <- replaceNA(vector)
	expected <- c()
	expect_identical(result, expected)
})

test_that("Testing replaceNA with a simple vector", {
	vector <- c(1,2,3,NA,4)
	result <- replaceNA(vector)
	expected <- c(1,2,3,3,4)
	expect_identical(result, expected)
})

test_that("Testing replaceNA starting with NA", {
	vector <- c(NA,1,2,NA,NA,3,NA,4,NA)
	result <- replaceNA(vector)
	expected <- c(0,1,2,2,2,3,3,4,4)
	expect_identical(result, expected)
})
