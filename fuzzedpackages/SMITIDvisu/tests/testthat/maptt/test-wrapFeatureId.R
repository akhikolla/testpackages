context("Testing wrapFeatureId")

test_that("Using a NA id", {
	id <- NA
	result <- wrapFeatureId(id)
	expected <- NA
	expect_identical(result, expected)
})

test_that("Using an empty string", {
	id <- ''
	result <- wrapFeatureId(id)
	expected <- '""'
	expect_identical(result, expected)
})

test_that("Using a numeric id", {
	id <- 15
	result <- wrapFeatureId(id)
	expected <- 15
	expect_identical(result, expected)
})

test_that("Using a string id", {
	id <- "ID"
	result <- wrapFeatureId(id)
	expected <- '"ID"'
	expect_identical(result, expected)
	expected <- '"ID2"'
	expect_false(identical(result, expected))
})

test_that("Using a factor level id", {
	id <- factor("Level2", levels = c("Level1", "Level2"))
	result <- wrapFeatureId(id)
	expected <- '"Level2"'
	expect_identical(result, expected)
})
