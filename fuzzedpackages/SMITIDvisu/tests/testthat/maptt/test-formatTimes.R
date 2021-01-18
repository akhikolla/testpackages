context("Testing formatTimes")

test_that("Date", {
	timeColumn <- c(
		as.Date("2019-01-01"),
		as.Date("2019-02-01"),
		as.Date("2020-01-01")
	)
	result <- formatTimes(timeColumn)
	expected <- c("1546300800000","1548979200000","1577836800000")
	expect_identical(result, expected)
})

test_that("julian day", {
	timeColumn <- c(0,1,30)

	result <- formatTimes(timeColumn)

	today <- Sys.Date()
	expected <- c(today, today+1, today+30)
	expected <- as.numeric(as.POSIXct(expected, origin = "1970-01-01", tz = "GMT"))
	expected <- sapply(expected, FUN = "*", 1000)
  expected <- format(expected, scientific = FALSE)

	expect_identical(result, expected)
})

test_that("timestamp", {
	timeColumn <- c(1546300800,1548979200,1577836800)
	result <- formatTimes(timeColumn)
	expected <- c("1546300800000","1548979200000","1577836800000")
	expect_identical(result, expected)
})

test_that("string julian days", {
	timeColumn <- c("0","1","30")

	result <- formatTimes(timeColumn)

	today <- Sys.Date()
	expected <- c(today, today+1, today+30)
	expected <- as.numeric(as.POSIXct(expected, origin = "1970-01-01", tz = "GMT"))
	expected <- sapply(expected, FUN = "*", 1000)
  expected <- format(expected, scientific = FALSE)

	expect_identical(result, expected)
})

test_that("string timestamp", {
	timeColumn <- c("1546300800","1548979200","1577836800")
	result <- formatTimes(timeColumn)
	expected <- c("1546300800000","1548979200000","1577836800000")
	expect_identical(result, expected)
})

test_that("Date string", {
	times <- c(
		"2019-01-01T00:00:00",
		"2019-01-02T00:00:00",
		"2019-01-03T00:00:00",
		"2019-01-08T00:00:00",
		"2019-01-09T00:00:00"
	)

	result <- formatTimes(times)
	expected <- c("1546300800000","1546387200000","1546473600000","1546905600000","1546992000000")
	expect_identical(result, expected)
})
