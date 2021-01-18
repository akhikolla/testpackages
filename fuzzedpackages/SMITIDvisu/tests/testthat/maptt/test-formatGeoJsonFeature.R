context("Testing formatGeoJsonFeature")

test_that("Minimal attributes", {
	properties <- "\"id\":100,\"linestringTimestamps\":[\"2019-01-01T00:00:01Z\",\"2019-01-02T00:00:01Z\",\"2019-01-03T00:00:01Z\"]"
	geometry <- "\"type\":\"LineString\",\"coordinates\":[[4.882602,43.91555],[4.883519,43.91425],[4.883259,43.91428]]"
	result <- formatGeoJsonFeature(properties, geometry)

	expected <- "{\"type\":\"Feature\",\"properties\":{\"id\":100,\"linestringTimestamps\":[\"2019-01-01T00:00:01Z\",\"2019-01-02T00:00:01Z\",\"2019-01-03T00:00:01Z\"]},\"geometry\":{\"type\":\"LineString\",\"coordinates\":[[4.882602,43.91555],[4.883519,43.91425],[4.883259,43.91428]]}}"
	expect_identical(result, expected)
})
