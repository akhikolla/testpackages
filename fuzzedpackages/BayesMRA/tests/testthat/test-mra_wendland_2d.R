context("testing wendland MRA functions")

test_that("wendland_basis", {
    expect_error(wendland_basis(1:10, 1:10), "radius must be a single positive numeric value")
    expect_error(wendland_basis(-10:10, 5), "d must be nonnegative")
    expect_error(wendland_basis(c(1:10, NA), 2), "d must not contain missing values")
    expect_equal(wendland_basis(1:5, 3), c(0.37717827566936, 0.013971447441955, 0, 0, 0))
})

test_that("illegal input for mra_wendland_2d", {
    locs <- NULL
    expect_error(mra_wendland_2d(locs), "locs must be a numeric matrix with N rows and 2 columns")
    locs <- matrix(NA, 10, 3)
    expect_error(mra_wendland_2d(locs), "locs must be a numeric matrix with N rows and 2 columns")
    locs <- matrix(NA, 10, 2)
    expect_error(mra_wendland_2d(locs), "locs must be a numeric matrix with N rows and 2 columns")
    locs <- matrix("aaa", 10, 2)
    expect_error(mra_wendland_2d(locs), "locs must be a numeric matrix with N rows and 2 columns")

    locs <- matrix(1:20, 10, 2)
    expect_error(mra_wendland_2d(locs, M = 3.5), "the number of resolutions M must be a positive integer")
    expect_error(mra_wendland_2d(locs, M = NA), "the number of resolutions M must be a positive integer")
    expect_error(mra_wendland_2d(locs, M = "aaa"), "the number of resolutions M must be a positive integer")
    expect_error(mra_wendland_2d(locs, n_neighbors = 3.5), "n_neighbors must be a positive integer")
    expect_error(mra_wendland_2d(locs, n_neighbors = NA), "n_neighbors must be a positive integer")
    expect_error(mra_wendland_2d(locs, n_neighbors = "aaa"), "n_neighbors must be a positive integer")
    expect_error(mra_wendland_2d(locs, n_coarse_grid = 3.5), "n_coarse_grid must be a positive integer")
    expect_error(mra_wendland_2d(locs, n_coarse_grid = NA), "n_coarse_grid must be a positive integer")
    expect_error(mra_wendland_2d(locs, n_coarse_grid = "aaa"), "n_coarse_grid must be a positive integer")
    expect_error(mra_wendland_2d(locs, n_padding = 3.5), "n_padding must be a positive integer")
    expect_error(mra_wendland_2d(locs, n_padding = NA), "n_padding must be a positive integer")
    expect_error(mra_wendland_2d(locs, n_padding = "aaa"), "n_padding must be a positive integer")
    expect_error(mra_wendland_2d(locs, use_spam = 3.5), "use_spam must be either TRUE or FALSE")
    expect_error(mra_wendland_2d(locs, use_spam = NA), "use_spam must be either TRUE or FALSE")
    expect_error(mra_wendland_2d(locs, use_spam = "aaa"), "use_spam must be either TRUE or FALSE")
    expect_error(mra_wendland_2d(locs, M = 4, n_coarse_grid = 3), "There are too many resolutions to form a reliable grid. Reduce M and try again.")

    locs <- matrix(1:6, 3, 2)
    expect_s3_class(mra_wendland_2d(locs, M = 2), "mra_wendland_2d")
    expect_s3_class(mra_wendland_2d(locs, M = 2, use_spam = FALSE), "mra_wendland_2d")

})


test_that("illegal input for mra_wendland_2d_pred", {
    locs <- matrix(1:20, 10, 2)
    MRA <- mra_wendland_2d(locs)
    locs_pred <- matrix(NA, 20, 2)
    expect_error(mra_wendland_2d_pred(locs, locs_pred, MRA), "locs_pred must be a numeric matrix with N rows and 2 columns")

    locs_pred <- matrix(1:30, 10, 3)
    expect_error(mra_wendland_2d_pred(locs, locs_pred, MRA), "locs_pred must be a numeric matrix with N rows and 2 columns")

    locs_pred <- matrix("11", 10, 2)
    expect_error(mra_wendland_2d_pred(locs, locs_pred, MRA), "locs_pred must be a numeric matrix with N rows and 2 columns")

    locs_pred <- 1:10
    expect_error(mra_wendland_2d_pred(locs, locs_pred, MRA), "locs_pred must be a numeric matrix with N rows and 2 columns")

    locs <- matrix(1:30, 10, 3)
    expect_error(mra_wendland_2d_pred(locs, locs_pred, MRA), "locs must be a numeric matrix with N rows and 2 columns")

    locs <- matrix(1:20, 10, 2)
    locs_pred <- matrix(1:20, 10, 2)
    expect_error(mra_wendland_2d_pred(locs, locs_pred, MRA, use_spam = 3.5), "use_spam must be either TRUE or FALSE")
    expect_error(mra_wendland_2d_pred(locs, locs_pred, MRA, use_spam = NA), "use_spam must be either TRUE or FALSE")
    expect_error(mra_wendland_2d_pred(locs, locs_pred, MRA, use_spam = "aaa"), "use_spam must be either TRUE or FALSE")

    class(MRA) <- NULL
    expect_error(mra_wendland_2d_pred(locs, locs_pred, MRA), 'MRA must be of class "mra_wendland_2d"')

    class(MRA) <- "XXX"
    expect_error(mra_wendland_2d_pred(locs, locs_pred, MRA), 'MRA must be of class "mra_wendland_2d"')

    MRA <- mra_wendland_2d(locs)
    expect_s3_class(mra_wendland_2d_pred(locs, locs_pred, MRA), "mra_wendland_2d_pred")
    expect_s3_class(mra_wendland_2d_pred(locs, locs_pred, MRA, use_spam = FALSE), "mra_wendland_2d_pred")
})
