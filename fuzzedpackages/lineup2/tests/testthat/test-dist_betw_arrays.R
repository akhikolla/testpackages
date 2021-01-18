context("Distance between arrays")

test_that("dist_betw_arrays() works", {

    set.seed(20171130)
    d <- c(5, 2, 8)
    x <- array(stats::rnorm(prod(d)), dim=d)
    y <- array(stats::rnorm(prod(d)), dim=d)
    dimnames(x) <- dimnames(y) <- list(paste0("ind", 1:d[1]),
                                       paste0("c", 1:d[2]),
                                       paste0("v", 1:d[3]))

    X <- cbind(x[,1,], x[,2,])
    Y <- cbind(y[,1,], y[,2,])
    colnames(X) <- colnames(Y) <- paste0("col", 1:prod(d[2:3]))

    expect_equal(dist_betw_arrays(x,y), dist_betw_matrices(X,Y))


    # different numbers of rows in the two arrays
    d2 <- c(8, 2, 8)
    y2 <- array(stats::rnorm(prod(d2)), dim=d2)
    dimnames(y2) <- list(paste0("ind", 1:d2[1]),
                         paste0("c", 1:d2[2]),
                         paste0("v", 1:d2[3]))

    Y2 <- cbind(y2[,1,], y2[,2,])
    colnames(Y2) <- paste0("col", 1:prod(d2[2:3]))

    expect_equal(dist_betw_arrays(x,y2), dist_betw_matrices(X,Y2))

    # give error if different shape
    d3 <- c(5, 3, 8)
    y3 <- array(stats::rnorm(prod(d3)), dim=d3)
    dimnames(y3) <- list(paste0("ind", 1:d3[1]),
                         paste0("c", 1:d3[2]),
                         paste0("v", 1:d3[3]))
    expect_error(dist_betw_arrays(x,y3))

    d4 <- c(5, 2, 8, 2)
    y4 <- array(stats::rnorm(prod(d4)), dim=d4)
    dimnames(y4) <- list(paste0("ind", 1:d4[1]),
                         paste0("c", 1:d4[2]),
                         paste0("v", 1:d4[3]),
                         paste0("vv", 1:d4[4]))
    expect_error(dist_betw_arrays(x,y4))

})
