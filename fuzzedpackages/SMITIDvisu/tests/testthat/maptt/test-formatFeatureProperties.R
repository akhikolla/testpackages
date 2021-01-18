context("Testing formatFeatureProperties")

test_that("Minimal attributes", {
	list <- list(
		"time" = "[1548979200]",
		"coordinates" = "[[4.883519,43.91425]]"
	)

	result <- formatFeatureProperties("\"ID1\"", list)
	expected <- "\"id\":\"ID1\",\"linestringTimestamps\":[1548979200]"
	expect_identical(result, expected)

	result <- formatFeatureProperties("100", list)
	expected <- "\"id\":100,\"linestringTimestamps\":[1548979200]"
	expect_identical(result, expected)
})

test_that("Additionnal properties", {
	list <- list(
		"time" = "[1548979200]",
		"coordinates" = "[[4.883519,43.91425]]",
		"prop1" = "[[100,400]]",
		"prop2" = "[\"foo\"]"
	)

	result <- formatFeatureProperties("\"ID1\"", list)
	expected <- "\"id\":\"ID1\",\"linestringTimestamps\":[1548979200],\"prop1\":[[100,400]],\"prop2\":[\"foo\"]"
	expect_identical(result, expected)
})
