context("Testing formatFeatureGeometry")

test_that("Minimal attributes", {
	list <- list(
		"time" = "[\"2019-01-02T00:00:01Z\"]",
		"coordinates" = "[[4.883519,43.91425]]"
	)

	result <- formatFeatureGeometry(list)
	expected <- "\"type\":\"LineString\",\"coordinates\":[[4.883519,43.91425]]"
	expect_identical(result, expected)
})
