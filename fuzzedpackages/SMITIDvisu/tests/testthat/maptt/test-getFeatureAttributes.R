context("Testing getFeatureAttributes")

test_that("Minimal attributes", {
	df <- data.frame(
		"time" = c(1546300800,1546387200,1546473600),
		"X" = c(4.882602,4.883519,4.883259),
		"Y" = c(43.91555,43.91425,43.91428)
	)

	result <- getFeatureAttributes(df)
	expected <- list(
		"time" = "[1546300800,1546387200,1546473600]",
		"coordinates" = "[[4.882602,43.91555],[4.883519,43.91425],[4.883259,43.91428]]"
	)
	expect_identical(result, expected)
})

test_that("Additionnal properties (unique-value-by-time)", {
	df <- data.frame(
		"time" = c(1546300800,1546387200,1546473600,1546300800),
		"X" = c(4.882602,4.883519,4.883259,NA),
		"Y" = c(43.91555,43.91425,43.91428,NA),
		"prop1" = c(100,200,300,400),
		"prop2" = c("foo","bar","foo","bar")
	)

	result <- getFeatureAttributes(df)
	expected <- list(
		"time" = "[1546300800,1546387200,1546473600]",
		"coordinates" = "[[4.882602,43.91555],[4.883519,43.91425],[4.883259,43.91428]]",
		"prop1" = "[100,200,300]",
		"prop2" = "[\"foo\",\"bar\",\"foo\"]"
	)
	expect_identical(result, expected)
})

test_that("Additionnal properties (multiple-values-by-time)", {
	df <- data.frame(
		"time" = c(1546300800,1546387200,1546473600,1546300800),
		"X" = c(4.882602,4.883519,4.883259,NA),
		"Y" = c(43.91555,43.91425,43.91428,NA),
		"prop1" = c(100,200,300,400),
		"prop2" = c("foo","bar","foo","bar")
	)

	result <- getFeatureAttributes(df, c('prop1'))
	expected <- list(
		"time" = "[1546300800,1546387200,1546473600]",
		"coordinates" = "[[4.882602,43.91555],[4.883519,43.91425],[4.883259,43.91428]]",
		"prop1" = "[[100,400],[200],[300]]",
		"prop2" = "[\"foo\",\"bar\",\"foo\"]"
	)
	expect_identical(result, expected)

	result <- getFeatureAttributes(df, c('prop2'))
	expected <- list(
		"time" = "[1546300800,1546387200,1546473600]",
		"coordinates" = "[[4.882602,43.91555],[4.883519,43.91425],[4.883259,43.91428]]",
		"prop1" = "[100,200,300]",
		"prop2" = "[[\"foo\",\"bar\"],[\"bar\"],[\"foo\"]]"
	)
	expect_identical(result, expected)

	result <- getFeatureAttributes(df, c('prop1','prop2'))
	expected <- list(
		"time" = "[1546300800,1546387200,1546473600]",
		"coordinates" = "[[4.882602,43.91555],[4.883519,43.91425],[4.883259,43.91428]]",
		"prop1" = "[[100,400],[200],[300]]",
		"prop2" = "[[\"foo\",\"bar\"],[\"bar\"],[\"foo\"]]"
	)
	expect_identical(result, expected)
})
