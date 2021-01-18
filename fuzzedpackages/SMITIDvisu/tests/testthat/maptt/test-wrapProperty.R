context("Testing wrapProperty")

test_that("Testing wrapProperty with empty vector", {
	property <- c()
	result <- wrapProperty(property)
	expected <- '[]'
	expect_identical(result, expected)
})

test_that("Testing wrapProperty with NA vector", {
	property <- c(NA,NA)
	result <- wrapProperty(property)
	expected <- '[]'
	expect_identical(result, expected)
})

test_that("Testing wrapProperty with empty strings vector", {
	property <- c("","")
	result <- wrapProperty(property, wrap = '["%s"]', sep = '","')
	expected <- '["",""]'
	expect_identical(result, expected)
})

test_that("Testing wrapProperty with numeric vector", {
	property <- c(1,2,3)
	result <- wrapProperty(property)
	expected <- '[1,2,3]'
	expect_identical(result, expected)
})

test_that("Testing wrapProperty with character vector", {
	property <- c("a","b","c")
	result <- wrapProperty(property, wrap = '["%s"]', sep = '","')
	expected <- '["a","b","c"]'
	expect_identical(result, expected)
})

test_that("Testing wrapProperty with vector mixed with NAs", {
	property <- c(1,3,NA,NA,5)
	result <- wrapProperty(property)
	expected <- '[1,3,5]'
	expect_identical(result, expected)
})
