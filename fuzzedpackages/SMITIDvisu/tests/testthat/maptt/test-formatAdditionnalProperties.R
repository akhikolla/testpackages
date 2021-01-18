context("Testing formatAdditionnalProperties")

test_that("Empty list", {
	list <- list()
	result <- formatAdditionnalProperties(list)
	expected <- ''

	expect_identical(result, expected)
})

test_that("Minimal list with no additionnal property", {
	list <- list(
		"time" = "[\"2019-01-02T00:00:01Z\"]",
		"coordinates" = "[[4.883519,43.91425]]"
	)
	result <- formatAdditionnalProperties(list)
	expected <- ''

	expect_identical(result, expected)
})

test_that("List with additionnal properties", {
	list <- list(
		"time" = "[\"2019-01-02T00:00:01Z\"]",
		"coordinates" = "[[4.883519,43.91425]]",
		"prop1" = "[[100,400]]",
		"prop2" = "[\"foo\"]"
	)
	result <- formatAdditionnalProperties(list)
	expected <- ',\"prop1\":[[100,400]],\"prop2\":[\"foo\"]'

	expect_identical(result, expected)
})
