# Check that the functions disjoint_maxima(), all_disjoint_maxima(),
# sliding_maxima() and all_maxima() agree and give the correct results
# on some simple examples and on the newlyn data

my_tol <- 1e-10

context("disjoint_maxima")

x <- 1:9
index <- c(3, 6, 9)
temp <- all_disjoint_maxima(x, b = 3)
y_mat <- matrix(x[index], 3, 1)
x_mat <- cbind(1:9)
test_that("x = 1:9, b = 3, all_disjoint_maxima", {
  testthat::expect_equal(temp$y_mat, y_mat, tol = my_tol)
})
test_that("x = 1:9, b = 3, all_disjoint_maxima input values", {
  testthat::expect_equal(temp$x_mat, x_mat, tol = my_tol)
})
temp2 <- disjoint_maxima(x, b = 3)
test_that("x = 1:9, b = 3, disjoint_maxima", {
  testthat::expect_equal(temp$y_mat, as.matrix(temp2$y), tol = my_tol)
})
test_that("x = 1:9, b = 3, disjoint_maxima input values", {
  testthat::expect_equal(temp$x_mat, as.matrix(temp2$x), tol = my_tol)
})

x <- 1:11
temp <- all_disjoint_maxima(x, b = 3)
index <- c(3, 6, 9, 4, 7, 10, 5, 8, 11)
y_mat <- matrix(x[index], 3, 3)
x_mat <- cbind(1:9, 2:10, 3:11)
test_that("x = 1:11, b = 3, all_disjoint_maxima", {
  testthat::expect_equal(temp$y_mat, y_mat, tol = my_tol)
})
test_that("x = 1:11, b = 3, all_disjoint_maxima input values", {
  testthat::expect_equal(temp$x_mat, x_mat, tol = my_tol)
})
temp2 <- disjoint_maxima(x, b = 3)
test_that("x = 1:11, b = 3, all_disjoint_maxima", {
  testthat::expect_equal(temp$y_mat[, 1], temp2$y, tol = my_tol)
})
test_that("x = 1:11, b = 3, all_disjoint_maxima input values", {
  testthat::expect_equal(temp$x_mat[, 1], temp2$x, tol = my_tol)
})

context("sliding_maxima")

# all_maxima() vs truth and all_maxima() vs sliding_maxima()

x <- 1:9
temp <- all_maxima(x, b = 3)
index <- 3:9
y_vec <- x[index]
x_vec <- 1:9
test_that("x = 1:9, b = 3, all_maxima", {
  testthat::expect_equal(temp$ys, y_vec, tol = my_tol)
})
test_that("x = 1:9, b = 3, all_maxima input values", {
  testthat::expect_equal(temp$xs, x_vec, tol = my_tol)
})
temp2 <- sliding_maxima(x, b = 3)
test_that("x = 1:9, b = 3, all_maxima", {
  testthat::expect_equal(temp2$y, y_vec, tol = my_tol)
})
test_that("x = 1:9, b = 3, all_maxima input values", {
  testthat::expect_equal(temp2$x, x_vec, tol = my_tol)
})

x <- 1:11
temp <- all_maxima(x, b = 3)
index <- c(3, 6, 9, 4, 7, 10, 5, 8, 11)
y_mat <- matrix(x[index], 3, 3)
x_mat <- cbind(1:9, 2:10, 3:11)
test_that("x = 1:11, b = 3, all_maxima", {
  testthat::expect_equal(temp$yd, y_mat, tol = my_tol)
})
test_that("x = 1:11, b = 3, all_maxima input values", {
  testthat::expect_equal(temp$xd, x_mat, tol = my_tol)
})
temp2 <- disjoint_maxima(x, b = 3)
test_that("x = 1:11, b = 3, all_disjoint_maxima", {
  testthat::expect_equal(temp$yd[, 1], temp2$y, tol = my_tol)
})
test_that("x = 1:11, b = 3, all_disjoint_maxima input values", {
  testthat::expect_equal(temp$xd[, 1], temp2$x, tol = my_tol)
})
temp3 <- all_maxima(x, b = 3, which_dj = "first")
temp4 <- all_maxima(x, b = 3, which_dj = "last")
test_that("x = 1:11, b = 3, all_maxima, which_dj = first", {
  testthat::expect_equal(temp$yd[, 1, drop = FALSE], temp3$yd, tol = my_tol)
})
test_that("x = 1:11, b = 3, all_maxima, which_dj = first, input values", {
  testthat::expect_equal(temp$xd[, 1, drop = FALSE], temp3$xd, tol = my_tol)
})
test_that("x = 1:11, b = 3, all_maxima, which_dj = first", {
  testthat::expect_equal(temp$yd[, 3, drop = FALSE], temp4$yd, tol = my_tol)
})
test_that("x = 1:11, b = 3, all_maxima, which_dj = first, input values", {
  testthat::expect_equal(temp$xd[, 3, drop = FALSE], temp4$xd, tol = my_tol)
})


context("all_maxima, newlyn")

# Check that the function all_maxima() gives the same results as the
# functions all_disjoint_maxima(), disjoint_maxima() and sliding_maxima(),
# using the newlyn data

# Sliding maxima and disjoint maxima
a_res <- all_maxima(newlyn, 100)
# All disjoint maxima
d_res <- all_disjoint_maxima(newlyn, 100)
# Sliding maxima
s_res <- sliding_maxima(newlyn, 100)
# Disjoint maxima, starting only from the first observation
d1_res <- disjoint_maxima(newlyn, 100)
# Disjoint maxima, ending only from the lastst observation
d2_res <- disjoint_maxima(newlyn, 100, which_dj = "last")

test_that("newlyn: disjoint contributing values", {
  testthat::expect_identical(a_res$xd, d_res$x_mat)
})
test_that("newlyn: disjoint maxima", {
  testthat::expect_identical(a_res$yd, d_res$y_mat)
})
test_that("newlyn: sliding contributing values", {
  testthat::expect_identical(a_res$xs, s_res$x)
})
test_that("newlyn: sliding maxima", {
  testthat::expect_identical(a_res$ys, s_res$y)
})
test_that("newlyn: disjoint contributing values vs disjoint_maxima()", {
  testthat::expect_identical(a_res$xd[, 1], d1_res$x)
})
test_that("newlyn: disjoint maxima vs disjoint_maxima()", {
  testthat::expect_identical(a_res$yd[, 1], d1_res$y)
})
test_that("newlyn: disjoint contributing values vs disjoint_maxima(), last", {
  testthat::expect_identical(a_res$xd[, ncol(a_res$xd)], d2_res$x)
})
test_that("newlyn: disjoint maxima vs disjoint_maxima(), last", {
  testthat::expect_identical(a_res$yd[, ncol(a_res$yd)], d2_res$y)
})

# Sliding maxima and (first) disjoint maxima
a_first_res <- all_maxima(newlyn, 100, which_dj = "first")
test_that("newlyn: which_dj = first", {
  testthat::expect_identical(a_res$xd[, 1, drop = FALSE], a_first_res$xd)
})
test_that("newlyn: which_dj = first, input values", {
  testthat::expect_identical(a_res$yd[, 1, drop = FALSE], a_first_res$yd)
})
# Sliding maxima and (first) disjoint maxima
a_last_res <- all_maxima(newlyn, 100, which_dj = "last")
test_that("newlyn: which_dj = last", {
  testthat::expect_identical(a_res$xd[, ncol(a_res$xd), drop = FALSE],
                             a_last_res$xd)
})
test_that("newlyn: which_dj = last, input values", {
  testthat::expect_identical(a_res$yd[, ncol(a_res$yd), drop = FALSE],
                             a_last_res$yd)
})
