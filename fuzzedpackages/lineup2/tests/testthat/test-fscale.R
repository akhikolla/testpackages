context("Standardize a vector")

test_that("fscale (in c++) works", {

    set.seed(20171201)

    x <- rnorm(24, 20, 5)
    y <- rnorm(24, 20, 5)
    expect_equivalent(fscale(x, y), as.numeric(scale(x)))

    # a couple of NAs in x
    x[2] <- x[24] <- NA
    expect_equivalent(fscale(x, y), as.numeric(scale(x)))

    # a couple of NAs in y as well
    y[3] <- y[12] <- NA
    z <- x
    z[3] <- z[12] <- NA
    expect_equivalent(fscale(x, y), as.numeric(scale(z)))

})
