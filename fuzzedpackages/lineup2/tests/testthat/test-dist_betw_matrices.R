context("Distance between matrices")

test_that("dist_betw_matrices() works", {

    set.seed(20171130)
    n <- 5
    p <- 10
    x <- matrix(stats::rnorm(n*p), nrow=n)
    expect_equivalent(dist_betw_matrices(x, x), as.matrix(stats::dist(x))/sqrt(p))

    rownames(x) <- LETTERS[1:n]
    expect_equal(dist_betw_matrices(x, x), as.matrix(stats::dist(x))/sqrt(p))
    expect_equal(dist_betw_matrices(x, x, "mad"), as.matrix(stats::dist(x, method="manhattan"))/p)

    m <- 3
    y <- matrix(stats::rnorm(m*p), nrow=m)
    rownames(y) <- letters[1:m]

    expect_error(dist_betw_matrices(x, y[,1:5]))
    expect_error(dist_betw_matrices(x, y[,1:5], "mad"))

    expect_equal(dist_betw_matrices(x, y), as.matrix(stats::dist(rbind(x,y)))[1:n, n + 1:m]/sqrt(p))
    expect_equal(dist_betw_matrices(x, y, "mad"),
                 as.matrix(stats::dist(rbind(x,y), method="manhattan"))[1:n, n + 1:m]/p)

    zx <- x
    zy <- y
    zx[1,1] <- zx[2,4] <- zx[3,4] <- NA
    zy[1,1] <- zy[2,3] <- zy[3,5] <- NA

    expect_equal(dist_betw_matrices(zx, zy),
                 as.matrix(stats::dist(rbind(zx,zy)))[1:n, n + 1:m]/sqrt(p))

    expect_equal(dist_betw_matrices(zx, zy, "mad"),
                 as.matrix(stats::dist(rbind(zx,zy), method="manhattan"))[1:n, n + 1:m]/p)

})


test_that("dist_betw_matrices() works with multi-core", {

    if(isnt_karl()) skip("These tests only run locally")

    set.seed(20171130)
    n <- 5
    p <- 10
    x <- matrix(stats::rnorm(n*p), nrow=n)
    expect_equivalent(dist_betw_matrices(x, x, cores=0), as.matrix(stats::dist(x))/sqrt(p))

    rownames(x) <- LETTERS[1:n]
    expect_equal(dist_betw_matrices(x, x, cores=0), as.matrix(stats::dist(x))/sqrt(p))
    expect_equal(dist_betw_matrices(x, x, "mad", cores=0), as.matrix(stats::dist(x, method="manhattan"))/p)

    m <- 3
    y <- matrix(stats::rnorm(m*p), nrow=m)
    rownames(y) <- letters[1:m]

    expect_error(dist_betw_matrices(x, y[,1:5], cores=0))
    expect_error(dist_betw_matrices(x, y[,1:5], "mad", cores=0))

    expect_equal(dist_betw_matrices(x, y, cores=0), as.matrix(stats::dist(rbind(x,y)))[1:n, n + 1:m]/sqrt(p))
    expect_equal(dist_betw_matrices(x, y, "mad", cores=0),
                 as.matrix(stats::dist(rbind(x,y), method="manhattan"))[1:n, n + 1:m]/p)

    zx <- x
    zy <- y
    zx[1,1] <- zx[2,4] <- zx[3,4] <- NA
    zy[1,1] <- zy[2,3] <- zy[3,5] <- NA

    expect_equal(dist_betw_matrices(zx, zy, cores=0),
                 as.matrix(stats::dist(rbind(zx,zy)))[1:n, n + 1:m]/sqrt(p))

    expect_equal(dist_betw_matrices(zx, zy, "mad", cores=0),
                 as.matrix(stats::dist(rbind(zx,zy), method="manhattan"))[1:n, n + 1:m]/p)

})
