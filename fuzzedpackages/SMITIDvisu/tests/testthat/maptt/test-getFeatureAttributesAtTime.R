context("Testing getFeatureAttributesAtTime")

test_that("Empty data frame", {
	df <- data.frame()

	result <- getFeatureAttributesAtTime(df, 0)
	expected <- data.frame()
	expect_identical(result, expected)
})

test_that("Minimal data frame (time, X, Y)", {
	df <- data.frame(
		"time" = c(1546300800,1546387200,1546473600,1546300800),
		"X" = c(4.882602,4.883519,4.883259,NA),
		"Y" = c(43.91555,43.91425,43.91428,NA)
	)

	result <- getFeatureAttributesAtTime(df, 1546300800)
	expected <- list(
		"time" = "1546300800",
		"coordinates" = "[4.882602,43.91555]"
	)
	expect_identical(result, expected)

	result <- getFeatureAttributesAtTime(df, 1546387200)
	expected <- list(
		"time" = "1546387200",
		"coordinates" = "[4.883519,43.91425]"
	)
	expect_identical(result, expected)
})

test_that("Data-frame with unique-value-by-time columns", {
	df <- data.frame(
		"time" = c(1546300800,1546387200,1546473600,1546300800),
		"X" = c(4.882602,4.883519,4.883259,NA),
		"Y" = c(43.91555,43.91425,43.91428,NA),
		"prop1" = c(100,200,300,400),
		"prop2" = c("foo","bar","foo","bar")
	)

	result <- getFeatureAttributesAtTime(df, 1546300800)
	expected <- list(
		"time" = "1546300800",
		"coordinates" = "[4.882602,43.91555]",
		"prop1" = "100",
		"prop2" = "\"foo\""
	)
	expect_true(all(result == expected))

	result <- getFeatureAttributesAtTime(df, 1546387200)
	expected <- list(
		"time" = "1546387200",
		"coordinates" = "[4.883519,43.91425]",
		"prop1" = "200",
		"prop2" = "\"bar\""
	)
	expect_true(all(result == expected))
})

test_that("Data-frame with multiple-values-by-time columns", {
	df <- data.frame(
		"time" = c(1546300800,1546387200,1546473600,1546300800),
		"X" = c(4.882602,4.883519,4.883259,NA),
		"Y" = c(43.91555,43.91425,43.91428,NA),
		"prop1" = c(100,200,300,400),
		"prop2" = c("foo","bar","foo","bar")
	)

	result <- getFeatureAttributesAtTime(df, 1546300800, c('prop1'))
	expected <- list(
		"time" = "1546300800",
		"coordinates" = "[4.882602,43.91555]",
		"prop1" = "[100,400]",
		"prop2" = "\"foo\""
	)
	expect_true(all(result == expected))

	result <- getFeatureAttributesAtTime(df, 1546387200, c('prop2'))
	expected <- list(
		"time" = "1546387200",
		"coordinates" = "[4.883519,43.91425]",
		"prop1" = "200",
		"prop2" = "[\"bar\"]"
	)
	expect_true(all(result == expected))
})
