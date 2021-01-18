context("Testing getFeature")

test_that("Minimal attributes", {
	df <- data.frame(
		"id" = c("100", "100", "100"),
		"time" = c(1546300800,1546387200,1546473600),
		"X" = c(4.882602,4.883519,4.883259),
		"Y" = c(43.91555,43.91425,43.91428)
	)

	result <- getFeature(df, "100")
	expected <- "{\"type\":\"Feature\",\"properties\":{\"id\":\"100\",\"linestringTimestamps\":[1546300800,1546387200,1546473600]},\"geometry\":{\"type\":\"LineString\",\"coordinates\":[[4.882602,43.91555],[4.883519,43.91425],[4.883259,43.91428]]}}"
	expect_identical(result, expected)

	df <- data.frame(
		"id" = c(100, 100, 100),
		"time" = c(1546300800,1546387200,1546473600),
		"X" = c(4.882602,4.883519,4.883259),
		"Y" = c(43.91555,43.91425,43.91428)
	)

	result <- getFeature(df, 100)
	expected <- "{\"type\":\"Feature\",\"properties\":{\"id\":100,\"linestringTimestamps\":[1546300800,1546387200,1546473600]},\"geometry\":{\"type\":\"LineString\",\"coordinates\":[[4.882602,43.91555],[4.883519,43.91425],[4.883259,43.91428]]}}"
	expect_identical(result, expected)
})

test_that("Additionnal properties (unique-value-by-time)", {
	df <- data.frame(
		"id" = c(100,100,100,100),
		"time" = c(1546300800,1546387200,1546473600,1546300800),
		"X" = c(4.882602,4.883519,4.883259,NA),
		"Y" = c(43.91555,43.91425,43.91428,NA),
		"prop1" = c(100,200,300,400),
		"prop2" = c("foo","bar","foo","bar")
	)

	result <- getFeature(df, 100)
	expected <- "{\"type\":\"Feature\",\"properties\":{\"id\":100,\"linestringTimestamps\":[1546300800,1546387200,1546473600],\"prop1\":[100,200,300],\"prop2\":[\"foo\",\"bar\",\"foo\"]},\"geometry\":{\"type\":\"LineString\",\"coordinates\":[[4.882602,43.91555],[4.883519,43.91425],[4.883259,43.91428]]}}"
	expect_identical(result, expected)
})

test_that("Additionnal properties with multiple-values-by-time", {
	df <- data.frame(
		"id" = c(100,100,100,100),
		"time" = c(1546300800,1546387200,1546473600,1546300800),
		"X" = c(4.882602,4.883519,4.883259,NA),
		"Y" = c(43.91555,43.91425,43.91428,NA),
		"prop1" = c(100,200,300,400),
		"prop2" = c("foo","bar","foo","bar")
	)

	result <- getFeature(df, 100, c('prop1'))
	expected <- "{\"type\":\"Feature\",\"properties\":{\"id\":100,\"linestringTimestamps\":[1546300800,1546387200,1546473600],\"prop1\":[[100,400],[200],[300]],\"prop2\":[\"foo\",\"bar\",\"foo\"]},\"geometry\":{\"type\":\"LineString\",\"coordinates\":[[4.882602,43.91555],[4.883519,43.91425],[4.883259,43.91428]]}}"
	expect_identical(result, expected)

	result <- getFeature(df, 100, c('prop2'))
	expected <- "{\"type\":\"Feature\",\"properties\":{\"id\":100,\"linestringTimestamps\":[1546300800,1546387200,1546473600],\"prop1\":[100,200,300],\"prop2\":[[\"foo\",\"bar\"],[\"bar\"],[\"foo\"]]},\"geometry\":{\"type\":\"LineString\",\"coordinates\":[[4.882602,43.91555],[4.883519,43.91425],[4.883259,43.91428]]}}"
	expect_identical(result, expected)

	result <- getFeature(df, 100, c('prop1','prop2'))
	expected <- "{\"type\":\"Feature\",\"properties\":{\"id\":100,\"linestringTimestamps\":[1546300800,1546387200,1546473600],\"prop1\":[[100,400],[200],[300]],\"prop2\":[[\"foo\",\"bar\"],[\"bar\"],[\"foo\"]]},\"geometry\":{\"type\":\"LineString\",\"coordinates\":[[4.882602,43.91555],[4.883519,43.91425],[4.883259,43.91428]]}}"
	expect_identical(result, expected)
})
