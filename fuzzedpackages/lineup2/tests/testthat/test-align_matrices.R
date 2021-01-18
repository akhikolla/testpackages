context("Align rows or columns of two matrices")

test_that("align_matrix_rows and align_matrix_cols work", {

    set.seed(20171202)

    n <- p <- q <-8
    m <- 10

    x <- matrix(rnorm(n*p), ncol=p)
    y <- matrix(rnorm(m*q), ncol=q)

    rownames(x) <- sample(letters, n)
    rownames(y) <- sample(letters, m)

    colnames(x) <- sample(LETTERS, p)
    colnames(y) <- sample(LETTERS, q)

    expect_equal(align_matrix_rows(x,y), list(x=x[c(2,3,7,8),], y=y[c(4,7,5,9),]))

    expect_equal(align_matrix_cols(x,y), list(x=x[,c(7,8)], y=y[,c(1,5)]))


})
