context("Testing formatCoordinates")

test_that("Empty df", {

	df <- data.frame(
		"X"=c(4.8),
		"Y"=c(43.5)
	)

	result <- formatCoordinates(df)
	expected <- "[4.8,43.5]"
	expect_identical(result, expected)
})

test_that("Empty df", {

	df <- data.frame(
		"X"=c(4.8,5.5),
		"Y"=c(43.5,44.5)
	)

	result <- formatCoordinates(df)
	expected <- "[4.8,43.5]"
	expect_identical(result, expected)
})
