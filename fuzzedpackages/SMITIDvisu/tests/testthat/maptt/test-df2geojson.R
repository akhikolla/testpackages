context("Testing df2geojson")

test_that("Testing empty df", {
	df <- data.frame(
		"id" = c(),
		"time" = c(),
		"X" = c(),
		"Y" = c()
	) # same as : data.frame()

	result <- df2geojson(df)
	expected <- '{\"type\":\"FeatureCollection\",\"features\":[]}'
	expect_identical(result, expected)
})

test_that("Testing df with minimal columns", {
	df <- data.frame(
		"id"=c(113,113),
		"time"=c("2019-01-02","2019-01-03"),
		"X"=c(4.882602,4.883519),
		"Y"=c(43.91555,43.91425)
	)

	result <- df2geojson(df)
	expected <- "{\"type\":\"FeatureCollection\",\"features\":[{\"type\":\"Feature\",\"properties\":{\"id\":113,\"linestringTimestamps\":[1546387200000,1546473600000]},\"geometry\":{\"type\":\"LineString\",\"coordinates\":[[4.882602,43.91555],[4.883519,43.91425]]}}]}"
	expect_identical(result, expected)
})

test_that("Testing df with unsorted values & NA", {
	df <- data.frame(
		"id" = c(100,100,100,100),
		"time" = c("2019-01-01","2019-01-02","2019-01-03","2019-01-01"),
		"X" = c(4.882602,4.883519,4.883259,NA),
		"Y" = c(43.91555,43.91425,43.91428,NA),
		"prop1" = c(100,200,300,400),
		"prop2" = c("foo","bar","foo","bar")
	)

	result <- df2geojson(df)
	expected <- "{\"type\":\"FeatureCollection\",\"features\":[{\"type\":\"Feature\",\"properties\":{\"id\":100,\"linestringTimestamps\":[1546300800000,1546387200000,1546473600000],\"prop1\":[100,200,300],\"prop2\":[\"foo\",\"bar\",\"foo\"]},\"geometry\":{\"type\":\"LineString\",\"coordinates\":[[4.882602,43.91555],[4.883519,43.91425],[4.883259,43.91428]]}}]}"
	expect_identical(result, expected)
})

test_that("Testing df with unsorted values & NA & string ids", {
	df <- data.frame(
		"id" = c("113","113","113","113","116","116","116"),
		"time" = c("2019-01-01","2019-01-02","2019-01-03","2019-01-01","2019-01-01","2019-01-09","2019-01-05"),
		"status" = c("Contact","","Removal","","","","Removal"),
		"infectedby" = c("0",NA,NA,"116","0",NA,NA),
		"X" = c(4.882602,4.883519,4.883259,NA,4.884697,4.883721,NA),
		"Y" = c(43.91555,43.91425,43.91428,NA,43.91454,43.91466,NA)
	)

	result <- df2geojson(df, c('infectedby'))
	expected <- "{\"type\":\"FeatureCollection\",\"features\":[{\"type\":\"Feature\",\"properties\":{\"id\":\"113\",\"linestringTimestamps\":[1546300800000,1546387200000,1546473600000],\"infectedby\":[[\"0\",\"116\"],[],[]],\"status\":[\"Contact\",\"\",\"Removal\"]},\"geometry\":{\"type\":\"LineString\",\"coordinates\":[[4.882602,43.91555],[4.883519,43.91425],[4.883259,43.91428]]}},{\"type\":\"Feature\",\"properties\":{\"id\":\"116\",\"linestringTimestamps\":[1546300800000,1546646400000,1546992000000],\"infectedby\":[[\"0\"],[],[]],\"status\":[\"\",\"Removal\",\"\"]},\"geometry\":{\"type\":\"LineString\",\"coordinates\":[[4.884697,43.91454],[4.884697,43.91454],[4.883721,43.91466]]}}]}"
	expect_identical(result, expected)
})
