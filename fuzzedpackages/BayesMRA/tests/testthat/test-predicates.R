context("testing predicates")

test_that("is_numeric function", {
    expect_false(is_numeric(c(2, 2), 1))
    expect_false(is_numeric(NA, 1))
    expect_false(is_numeric("NA", 1))
    expect_false(is_numeric(NULL, 1))
    expect_true(is_numeric(matrix(1:6, 3, 2), 6))
})

test_that("is_positive_numeric function", {
    expect_true(is_positive_numeric(1, 1))
    expect_error(is_positive_numeric(-3))
    expect_false(is_positive_numeric(c(2, 2, 2), 2))
    expect_false(is_positive_numeric(NA))
    expect_false(is_positive_numeric("NA"))
    expect_false(is_positive_numeric(NULL))
    expect_true(is_positive_numeric(1:6, 6))
    expect_true(is_positive_numeric(matrix(1:6, 3, 2), 6))
})

test_that("is_numeric_vector function", {
    expect_error(is_numeric_vector(c(2, 2)))
    expect_error(is_numeric_vector(matrix(1:6, 3, 2)))
    expect_false(is_numeric_vector(c(2, 2), 3))
    expect_false(is_numeric_vector(c(NA, 2), 2))
    expect_false(is_numeric_vector(NA))
    expect_false(is_numeric_vector("NA"))
    expect_false(is_numeric_vector(NULL))
    expect_false(is_numeric_vector(matrix(1:6, 3, 2), 6))
    expect_true(is_numeric_vector(c(2, 2), 2))

})

test_that("is_numeric_matrix function", {
    expect_error(is_numeric_matrix(c(2, 2), 2))
    ## add in testing on the input dimensions
    expect_error(is_numeric_matrix(matrix(c(2, 2), NA, 2)))
    expect_false(is_numeric_matrix(matrix(1:6, 3, 2), 2, 2))
    expect_false(is_numeric_matrix(matrix(c(1:5, NA), 3, 2), 2, 2))
    expect_false(is_numeric_matrix(matrix(c(1:5, "NA"), 3, 2), 2, 2))
    expect_false(is_numeric_matrix(NULL))
    expect_true(is_numeric_matrix(matrix(1:6, 3, 2), 3, 2))
})

test_that("is_sympd_matrix function", {
    ## add in testing on the input dimensions
    expect_error(is_numeric_matrix(matrix(c(2, 2), NA, 2)))
    expect_false(is_sympd_matrix(matrix(1:4, 2, 2), 2))
    expect_false(is_sympd_matrix(matrix(rep(NA, 4), 2, 2), 2))
    expect_false(is_sympd_matrix(NULL, 2))
    expect_false(is_sympd_matrix(c(2, 2), 2))
    expect_false(is_sympd_matrix(matrix(c(1:5, NA), 3, 2), 2))
    expect_false(is_sympd_matrix(matrix(1:6, 3, 2), 3))
    expect_false(is_sympd_matrix(matrix(1:6, 3, 2), 2))
    expect_error(is_sympd_matrix(NULL, 2, NA))
    expect_error(is_sympd_matrix(matrix(1:6, 3, 2), 3, 2))
    expect_true(is_sympd_matrix(diag(4), 4))
    expect_true(is_sympd_matrix(exp(- as.matrix(dist(1:4))), 4))
})


test_that("is_integer", {
    expect_error(is_integer(1L))
    expect_false(is_integer(2.5, 1))
    expect_false(is_integer(NA, 1))
    expect_false(is_integer("NA", 1))
    expect_false(is_integer(NULL, 1))
    expect_false(is_integer(1L:6L, 4))
    expect_true(is_integer(1, 1))
    expect_true(is_integer(1L, 1))
    expect_true(is_integer(1L:6L, 6))
    expect_true(is_integer(matrix(1:6, 3, 2), 6))
    expect_true(is_integer(2.0, 1))
})

test_that("is_integer_matrix", {
    expect_error(is_integer_matrix(1L))
    expect_error(is_integer_matrix(2.5, 1))
    expect_false(is_integer_matrix(NA, 1, 1))
    expect_false(is_integer_matrix("NA", 1))
    expect_false(is_integer_matrix(NULL, 1))
    expect_false(is_integer_matrix(1L:6L, 4, 1))
    expect_false(is_integer_matrix(1, 1, 1))
    expect_false(is_integer_matrix(matrix(1:6 + 0.5, 3, 2), 3, 2))
    expect_false(is_integer_matrix(matrix(1:6, 2, 3), 3, 2))
    expect_true(is_integer_matrix(matrix(1:6, 2, 3), 2, 3))
    expect_true(is_integer_matrix(matrix(1L:6L, 2, 3), 2, 3))
    expect_true(is_integer_matrix(matrix(2.0, 2, 2), 2, 2))
})


test_that("is_positive_integer function", {
    expect_true(is_positive_integer(1, 1))
    expect_false(is_positive_integer(1.1, 1))
    expect_error(is_positive_integer(-3))
    expect_false(is_positive_integer(-3, 1))
    expect_false(is_positive_integer(c(2, 2, 2), 2))
    expect_false(is_positive_integer(NA))
    expect_false(is_positive_integer("NA"))
    expect_false(is_positive_integer(NULL))
    expect_true(is_positive_integer(1:6, 6))
    expect_false(is_positive_integer(1:6 + 0.5, 6))
    expect_true(is_positive_integer(matrix(1:6, 3, 2), 6))
})
