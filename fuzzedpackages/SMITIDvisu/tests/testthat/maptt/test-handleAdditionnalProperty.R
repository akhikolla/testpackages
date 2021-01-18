context("Testing handleAdditionnalProperty")

test_that("Testing handleAdditionnalProperty with empty vector", {
	feature <- data.frame(
		'prop1' = c(),
		'prop2' = c()
	)

	result <- handleAdditionnalProperty('prop1', feature)
	expected <- ''
	expect_identical(result, expected)

	result <- handleAdditionnalProperty('prop1', feature, c(''))
	expected <- ''
	expect_identical(result, expected)

	result <- handleAdditionnalProperty('prop1', feature, c('prop1'))
	expected <- '[]'
	expect_identical(result, expected)
})

test_that("Unique values by time", {
	propertyName <- 'prop4'
	feature <- data.frame(
		'prop3' = c(1,2,3),
		'prop4' = c("a","b","c")
	)

	# With numeric vector
	result <- handleAdditionnalProperty('prop3', feature)
	expected <- '1'
	expect_identical(result, expected)

	# With character vector (transformed into factor inside df)
	result <- handleAdditionnalProperty('prop4', feature)
	expected <- '"a"'
	expect_identical(result, expected)
})

test_that("Unique values by time starting with NA", {
	propertyName <- 'prop4'
	feature <- data.frame(
		'prop3' = c(NA,2,3),
		'prop4' = c(NA,"b","c")
	)

	# With numeric vector
	result <- handleAdditionnalProperty('prop3', feature)
	expected <- '2'
	expect_identical(result, expected)

	# With character vector (transformed into factor inside df)
	result <- handleAdditionnalProperty('prop4', feature)
	expected <- '"b"'
	expect_identical(result, expected)
})

test_that("Multiple values by time", {
	multipleValuesByTime <- c('prop1', 'prop2')
	feature <- data.frame(
		'prop1' = c(1,2,3),
		'prop2' = c("a","b","c")
	)

	# With numeric vector
	result <- handleAdditionnalProperty('prop1', feature, multipleValuesByTime)
	expected <- '[1,2,3]'
	expect_identical(result, expected)

	# With character vector (transformed into factor inside df)
	result <- handleAdditionnalProperty('prop2', feature, multipleValuesByTime)
	expected <- '["a","b","c"]'
	expect_identical(result, expected)
})

test_that("Multiple values by time mixed with NAs", {
	multipleValuesByTime <- c('prop1', 'prop2')
	feature <- data.frame(
		'prop1' = c(NA,1,2,NA,3,NA),
		'prop2' = c(NA,"a","b",NA,"c",NA)
	)

	# With numeric vector
	result <- handleAdditionnalProperty('prop1', feature, multipleValuesByTime)
	expected <- '[1,2,3]'
	expect_identical(result, expected)

	# With character vector (transformed into factor inside df)
	result <- handleAdditionnalProperty('prop2', feature, multipleValuesByTime)
	expected <- '["a","b","c"]'
	expect_identical(result, expected)
})
