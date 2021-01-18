context("Testing prepareData")

test_that("Empty df", {
	df <- data.frame()

	result <- prepareData(df)
	expected <- data.frame()
	expect_identical(result, expected)
})

test_that("df with minimal attributes & julien date", {
	df <- data.frame(
		"id" = c(1,2,3),
		"time" = c("2019-01-01T00:00:00","2019-01-02T00:00:00","2019-01-03T00:00:00")
	)

	result <- prepareData(df)
	expected <- data.frame(
	  "id" = c(1,2,3),
		"time" = c("1546300800000","1546387200000","1546473600000"),
		stringsAsFactors = FALSE
	)
	expect_identical(result, expected)
})

test_that("df with unsorted minimal attributes & julien date", {
	df <- data.frame(
		"id" = c(4,4,1,2,1,3),
		"time" = c(1559347200,1546300800,1546300800,1548979200,1559260800,1559260800)
	)

	result <- prepareData(df)

	# Reset row names for the test comparison
	rownames(result) <- NULL
	# Drop unused factor levels for the test comparison
	# levels(result$status) <- droplevels(result$status)

	expected <- data.frame(
	  "id" = c(1,1,2,3,4,4),
		"time" = c("1546300800000","1559260800000","1548979200000","1559260800000","1546300800000","1559347200000"),
		stringsAsFactors = FALSE
	)
	expect_identical(result, expected)
})

test_that("df with unsorted minimal attributes & julien date", {
	df <- data.frame(
		"id" = c(4,4,1,2,1,3),
		"time" = c(1559347200,1546300800,1546300800,1548979200,1559260800,1559260800),
		stringsAsFactors = FALSE
	)

	result <- prepareData(df)

	# Reset row names for the test comparison
	rownames(result) <- NULL
	# Drop unused factor levels for the test comparison
	# levels(result$status) <- droplevels(result$status)
	expected <- data.frame(
	  "id" = c(1,1,2,3,4,4),
		"time" = c("1546300800000","1559260800000","1548979200000","1559260800000","1546300800000","1559347200000"),
		stringsAsFactors = FALSE
	)
	expect_identical(result, expected)
})
