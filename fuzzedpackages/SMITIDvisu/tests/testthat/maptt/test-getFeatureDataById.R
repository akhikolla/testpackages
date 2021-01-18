context("Testing getFeatureDataById")

test_that("Empty df", {
	df <- data.frame()
	result <- getFeatureDataById(df, 1)
	expected <- data.frame()
	expect_identical(result, expected)
})

test_that("df containing a unique id", {
	df <- data.frame(
		"id"=c(1,1),
		"time"=c(1,2),
		"X"=c(4.882602,4.883519),
		"Y"=c(43.91555,43.91425)
	)
	result <- getFeatureDataById(df,1)
	expected <- data.frame(
		"time"=c(1,2),
		"X"=c(4.882602,4.883519),
		"Y"=c(43.91555,43.91425)
	)
	expect_identical(result, expected)
})

test_that("df containing several ids", {
	df <- data.frame(
		"id"=c(1,2),
		"time"=c(1,2),
		"X"=c(4.882602,4.883519),
		"Y"=c(43.91555,43.91425)
	)
	result <- getFeatureDataById(df,1)
	expected <- data.frame(
		"time"=c(1),
		"X"=c(4.882602),
		"Y"=c(43.91555)
	)
	expect_identical(result, expected)
})

test_that("df with one different column name", {
	df <- data.frame(
		"id"=c(1,2,1),
		"time"=c(1,2,2),
		"X"=c(4.882602,4.883519,5.5),
		"Y"=c(43.91555,43.91425,44.5)
	)
	result <- getFeatureDataById(df,1)
	expected <- data.frame(
		"otherName"=c(1,2),
		"X"=c(4.882602,5.5),
		"Y"=c(43.91555,44.5)
	)
	expect_false(identical(result, expected))
})

test_that("df with a different value", {
	df <- data.frame(
		"id"=c(1,2,1),
		"time"=c(1,2,2),
		"X"=c(4.882602,4.883519,5.5),
		"Y"=c(43.91555,43.91425,44.5)
	)
	result <- getFeatureDataById(df,1)
	expected <- data.frame(
		"time"=c(1,3),
		"X"=c(4.882602,5.5),
		"Y"=c(43.91555,44.5)
	)
	expect_false(identical(result, expected))
})

test_that("Additionnal columns", {
	df <- data.frame(
		"id" = c(113,113,113,113,116,116,116),
		"time" = c(0,1,2,0,0,8,4),
		"status" = c("Contact","","Removal","","","","Removal"),
		"infectedby" = c(0,NA,NA,116,0,NA,NA),
		"X" = c(4.882602,4.883519,4.883259,NA,4.884697,4.883721,NA),
		"Y" = c(43.91555,43.91425,43.91428,NA,43.91454,43.91466,NA)
	)
	result <- getFeatureDataById(df,116)

	# Reset row names for the test comparison
	rownames(result) <- NULL
	# Drop unused factor levels for the test comparison
	#levels(result$status) <- droplevels(result$status)
	
	expected <- data.frame(
		"time" = c(0,8,4),
		"status" = c("","","Removal"),
		"infectedby" = c(0,NA,NA),
		"X" = c(4.884697,4.883721,NA),
		"Y" = c(43.91454,43.91466,NA)
	)

	expect_identical(result, expected)
})
