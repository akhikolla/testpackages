context("testing precision matrices")

# test_that("make_Q", {
#
#     expect_error(make_Q(4, "aaa"), "phi must be between -1 and 1")
#     expect_error(make_Q(4, NA), "phi must be between -1 and 1")
#     expect_error(make_Q(4, -2), "phi must be between -1 and 1")
#     expect_error(make_Q(4, 2), "phi must be between -1 and 1")
#     make_Q(4, 0.9)
#     # expect_error(wendland_basis(1:10, 1:10), "radius must be a single positive numeric value")
#     # expect_error(wendland_basis(-10:10, 5), "d must be nonnegative")
#     # expect_error(wendland_basis(c(1:10, NA), 2), "d must not contain missing values")
#     # expect_equal(wendland_basis(1:5, 3), c(0.37717827566936, 0.013971447441955, 0, 0, 0))
# })

test_that("make_Q_alpha_2d", {

    expect_identical(
        make_Q_alpha_2d(4, 0.5),
        structure(list(new("spam",
                           entries = c(2, -0.5, -0.5, -0.5, 3, -0.5, -0.5,
                                       -0.5, 3, -0.5, -0.5, -0.5, 2, -0.5, -0.5, 3, -0.5, -0.5, -0.5,
                                       -0.5, 4, -0.5, -0.5, -0.5, -0.5, 4, -0.5, -0.5, -0.5, -0.5, 3,
                                       -0.5, -0.5, 3, -0.5, -0.5, -0.5, -0.5, 4, -0.5, -0.5, -0.5, -0.5,
                                       4, -0.5, -0.5, -0.5, -0.5, 3, -0.5, -0.5, 2, -0.5, -0.5, -0.5,
                                       3, -0.5, -0.5, -0.5, 3, -0.5, -0.5, -0.5, 2),
                           colindices = c(1L,
                                          2L, 5L, 1L, 2L, 3L, 6L, 2L, 3L, 4L, 7L, 3L, 4L, 8L, 1L, 5L, 6L,
                                          9L, 2L, 5L, 6L, 7L, 10L, 3L, 6L, 7L, 8L, 11L, 4L, 7L, 8L, 12L,
                                          5L, 9L, 10L, 13L, 6L, 9L, 10L, 11L, 14L, 7L, 10L, 11L, 12L, 15L,
                                          8L, 11L, 12L, 16L, 9L, 13L, 14L, 10L, 13L, 14L, 15L, 11L, 14L,
                                          15L, 16L, 12L, 15L, 16L),
                           rowpointers = c(1L, 4L, 8L, 12L, 15L,
                                           19L, 24L, 29L, 33L, 37L, 42L, 47L, 51L, 54L, 58L, 62L, 65L),
                           dimension = c(16L, 16L))),
                  class = "Q_alpha"
        )
    )

    n_dims <- c(4, 16)
    phi <- rep(1, 3)
    expect_error(make_Q_alpha_2d(n_dims, phi), "n_dims and phi must both be vectors of length M.")

    n_dims <- c(4, 16, 32)
    phi <- rep(1, 2)
    expect_error(make_Q_alpha_2d(n_dims, phi), "n_dims and phi must both be vectors of length M.")

    n_dims <- c(4, 16)
    phi <- rep(NA, 2)
    expect_error(make_Q_alpha_2d(n_dims, phi), "phi must be a numeric vector of length M with values between -1 and 1.")
    phi <- rep(2, 2)
    expect_error(make_Q_alpha_2d(n_dims, phi), "phi must be a numeric vector of length M with values between -1 and 1.")
    phi <- rep(-2, 2)
    expect_error(make_Q_alpha_2d(n_dims, phi), "phi must be a numeric vector of length M with values between -1 and 1.")

    expect_error(make_Q_alpha_2d(4, 1.5), "phi must be a numeric vector of length M with values between -1 and 1.")
    expect_error(make_Q_alpha_2d(4, -1.5), "phi must be a numeric vector of length M with values between -1 and 1.")
    expect_error(make_Q_alpha_2d(c(2, 4), c(0.5, 1.5)), "phi must be a numeric vector of length M with values between -1 and 1.")
    expect_error(make_Q_alpha_2d(c(2, 4), c(0.5, NA)), "phi must be a numeric vector of length M with values between -1 and 1.")
    expect_error(make_Q_alpha_2d(c(2, 4), c(0.5, "aaa")), "phi must be a numeric vector of length M with values between -1 and 1.")
    expect_error(make_Q_alpha_2d(c(2, NA), c(0.5, 0.5)), "n_dims must be a vector of integers of length M.")
    expect_error(make_Q_alpha_2d(c(2, "aaa"), c(0.5, 0.5)), "n_dims must be a vector of integers of length M.")
    expect_error(make_Q_alpha_2d(c(2, 3.5), c(0.5, 0.5)), "n_dims must be a vector of integers of length M.")

    n_dims <- c(4, 16)
    phi <- rep(1, 2)

    expect_error(make_Q_alpha_2d(n_dims, phi, use_spam = "TRUE"), "use_spam must be either TRUE or FALSE.")
    expect_error(make_Q_alpha_2d(n_dims, phi, use_spam = 3.5), "use_spam must be either TRUE or FALSE.")
    expect_error(make_Q_alpha_2d(n_dims, phi, use_spam = NA), "use_spam must be either TRUE or FALSE.")

    expect_error(make_Q_alpha_2d(n_dims, phi, prec_model = "AAA"), 'The only valid options for prec_model are \"CAR\" and \"SAR\".')
    expect_error(make_Q_alpha_2d(n_dims, phi, prec_model = NA), 'The only valid options for prec_model are \"CAR\" and \"SAR\".')
    expect_error(make_Q_alpha_2d(n_dims, phi, prec_model = 32), 'The only valid options for prec_model are \"CAR\" and \"SAR\".')

})


test_that("make_Q_alpha_tau2", {

    n_dims <- c(4, 16)
    phi    <- c(1, 1)
    tau2 <- c(2, 2)
    Q_alpha <- make_Q_alpha_2d(n_dims, phi)
    class(Q_alpha) <- "aaa"
    expect_error(make_Q_alpha_tau2(Q_alpha, tau2), 'Q_alpha must by of class "Q_alpha" which is the output of make_Q_alpha_2d\\(\\)')
    class(Q_alpha) <- "Q_alpha"
    tau2 <- 1:3
    expect_error(make_Q_alpha_tau2(Q_alpha, tau2), "Q_alpha must be a list of length M and tau2 must be a positive numeric vector of length M.")
    tau2 <- c(1, NA)
    expect_error(make_Q_alpha_tau2(Q_alpha, tau2), "tau2 must be a positive numeric vector of length M.")
    tau2 <- c(1, "aaa")
    expect_error(make_Q_alpha_tau2(Q_alpha, tau2), "tau2 must be a positive numeric vector of length M.")
    tau2 <- c(1, -1)
    expect_error(make_Q_alpha_tau2(Q_alpha, tau2), "tau2 must be a positive numeric vector of length M.")

    tau2 <- 1:2
    expect_error(make_Q_alpha_tau2(Q_alpha, tau2, use_spam = "TRUE"), "use_spam must be either TRUE or FALSE.")
    expect_error(make_Q_alpha_tau2(Q_alpha, tau2, use_spam = 3.5), "use_spam must be either TRUE or FALSE.")
    expect_error(make_Q_alpha_tau2(Q_alpha, tau2, use_spam = NA), "use_spam must be either TRUE or FALSE.")

    expect_equal({
        n_dims <- c(2, 4)
        phi    <- c(1, 1)
        Q_alpha <- make_Q_alpha_2d(n_dims, phi, prec_model = "SAR")
        tau2 <- c(3, 5)
        make_Q_alpha_tau2(Q_alpha, tau2)
    },
    {
        new("spam", entries = c(4.5, -3, -3, 1.5, -3, 4.5, 1.5, -3, -3,
                                1.5, 4.5, -3, 1.5, -3, -3, 4.5, 7.5, -4.16666666666667, 0.833333333333333,
                                -4.16666666666667, 1.25, 0.833333333333333, -4.16666666666667,
                                6.66666666666667, -3.33333333333333, 0.833333333333333, 1.11111111111111,
                                -2.91666666666667, 0.833333333333333, 0.416666666666667, 0.833333333333333,
                                -3.33333333333333, 6.66666666666667, -4.16666666666667, 0.833333333333333,
                                -2.91666666666667, 1.11111111111111, 0.416666666666667, 0.833333333333333,
                                -4.16666666666667, 7.5, 1.25, -4.16666666666667, 0.833333333333333,
                                -4.16666666666667, 1.11111111111111, 6.66666666666667, -2.91666666666667,
                                0.416666666666667, -3.33333333333333, 0.833333333333333, 0.833333333333333,
                                1.25, -2.91666666666667, 0.833333333333333, -2.91666666666667,
                                6.25, -2.5, 0.416666666666667, 0.833333333333333, -2.5, 0.625,
                                0.416666666666667, 0.833333333333333, -2.91666666666667, 1.25,
                                0.416666666666667, -2.5, 6.25, -2.91666666666667, 0.625, -2.5,
                                0.833333333333333, 0.416666666666667, 1.11111111111111, -4.16666666666667,
                                0.416666666666667, -2.91666666666667, 6.66666666666667, 0.833333333333333,
                                -3.33333333333333, 0.833333333333333, 0.833333333333333, -3.33333333333333,
                                0.833333333333333, 6.66666666666667, -2.91666666666667, 0.416666666666667,
                                -4.16666666666667, 1.11111111111111, 0.416666666666667, 0.833333333333333,
                                -2.5, 0.625, -2.91666666666667, 6.25, -2.5, 0.416666666666667,
                                1.25, -2.91666666666667, 0.833333333333333, 0.416666666666667,
                                0.625, -2.5, 0.833333333333333, 0.416666666666667, -2.5, 6.25,
                                -2.91666666666667, 0.833333333333333, -2.91666666666667, 1.25,
                                0.833333333333333, 0.833333333333333, -3.33333333333333, 0.416666666666667,
                                -2.91666666666667, 6.66666666666667, 1.11111111111111, -4.16666666666667,
                                0.833333333333333, -4.16666666666667, 1.25, 7.5, -4.16666666666667,
                                0.833333333333333, 0.416666666666667, 1.11111111111111, -2.91666666666667,
                                0.833333333333333, -4.16666666666667, 6.66666666666667, -3.33333333333333,
                                0.833333333333333, 0.416666666666667, 0.833333333333333, -2.91666666666667,
                                1.11111111111111, 0.833333333333333, -3.33333333333333, 6.66666666666667,
                                -4.16666666666667, 0.833333333333333, 1.25, -4.16666666666667,
                                0.833333333333333, -4.16666666666667, 7.5),
            colindices = c(1L,
                           2L, 3L, 4L, 1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L, 5L,
                           6L, 7L, 9L, 10L, 13L, 5L, 6L, 7L, 8L, 9L, 10L, 11L, 14L, 5L,
                           6L, 7L, 8L, 10L, 11L, 12L, 15L, 6L, 7L, 8L, 11L, 12L, 16L, 5L,
                           6L, 9L, 10L, 11L, 13L, 14L, 17L, 5L, 6L, 7L, 9L, 10L, 11L, 12L,
                           13L, 14L, 15L, 18L, 6L, 7L, 8L, 9L, 10L, 11L, 12L, 14L, 15L,
                           16L, 19L, 7L, 8L, 10L, 11L, 12L, 15L, 16L, 20L, 5L, 9L, 10L,
                           13L, 14L, 15L, 17L, 18L, 6L, 9L, 10L, 11L, 13L, 14L, 15L, 16L,
                           17L, 18L, 19L, 7L, 10L, 11L, 12L, 13L, 14L, 15L, 16L, 18L, 19L,
                           20L, 8L, 11L, 12L, 14L, 15L, 16L, 19L, 20L, 9L, 13L, 14L, 17L,
                           18L, 19L, 10L, 13L, 14L, 15L, 17L, 18L, 19L, 20L, 11L, 14L, 15L,
                           16L, 17L, 18L, 19L, 20L, 12L, 15L, 16L, 18L, 19L, 20L),
            rowpointers = c(1L,
                            5L, 9L, 13L, 17L, 23L, 31L, 39L, 45L, 53L, 64L, 75L, 83L, 91L,
                            102L, 113L, 121L, 127L, 135L, 143L, 149L),
            dimension = c(20L,
                          20L)
            )
    })

    expect_identical({
        n_dims <- c(2, 4)
        phi    <- c(1, 1)
        Q_alpha <- make_Q_alpha_2d(n_dims, phi)
        tau2 <- c(3, 5)
        make_Q_alpha_tau2(Q_alpha, tau2)
    },
    {
        new("spam",
            entries = c(6, -3, -3, -3, 6, -3, -3, 6, -3, -3,
                        -3, 6, 10, -5, -5, -5, 15, -5, -5, -5, 15, -5, -5, -5, 10, -5,
                        -5, 15, -5, -5, -5, -5, 20, -5, -5, -5, -5, 20, -5, -5, -5, -5,
                        15, -5, -5, 15, -5, -5, -5, -5, 20, -5, -5, -5, -5, 20, -5, -5,
                        -5, -5, 15, -5, -5, 10, -5, -5, -5, 15, -5, -5, -5, 15, -5, -5,
                        -5, 10),
            colindices = c(1L, 2L, 3L, 1L, 2L, 4L, 1L, 3L, 4L, 2L,
                           3L, 4L, 5L, 6L, 9L, 5L, 6L, 7L, 10L, 6L, 7L, 8L, 11L, 7L, 8L,
                           12L, 5L, 9L, 10L, 13L, 6L, 9L, 10L, 11L, 14L, 7L, 10L, 11L, 12L,
                           15L, 8L, 11L, 12L, 16L, 9L, 13L, 14L, 17L, 10L, 13L, 14L, 15L,
                           18L, 11L, 14L, 15L, 16L, 19L, 12L, 15L, 16L, 20L, 13L, 17L, 18L,
                           14L, 17L, 18L, 19L, 15L, 18L, 19L, 20L, 16L, 19L, 20L),
            rowpointers = c(1L,
                            4L, 7L, 10L, 13L, 16L, 20L, 24L, 27L, 31L, 36L, 41L, 45L, 49L,
                            54L, 59L, 63L, 66L, 70L, 74L, 77L),
            dimension = c(20L, 20L))
    })

    expect_equal({
        n_dims <- c(2, 4)
        phi    <- c(1, 1)
        Q_alpha <- make_Q_alpha_2d(n_dims, phi, prec_model = "SAR", use_spam = FALSE)
        tau2 <- c(3, 5)
        make_Q_alpha_tau2(Q_alpha, tau2, use_spam = FALSE)
    },
    {
        new("dgCMatrix",
            i = c(0L, 1L, 2L, 3L, 0L, 1L, 2L, 3L, 0L, 1L,
                  2L, 3L, 0L, 1L, 2L, 3L, 4L, 5L, 6L, 8L, 9L, 12L, 4L, 5L, 6L,
                  7L, 8L, 9L, 10L, 13L, 4L, 5L, 6L, 7L, 9L, 10L, 11L, 14L, 5L,
                  6L, 7L, 10L, 11L, 15L, 4L, 5L, 8L, 9L, 10L, 12L, 13L, 16L, 4L,
                  5L, 6L, 8L, 9L, 10L, 11L, 12L, 13L, 14L, 17L, 5L, 6L, 7L, 8L,
                  9L, 10L, 11L, 13L, 14L, 15L, 18L, 6L, 7L, 9L, 10L, 11L, 14L,
                  15L, 19L, 4L, 8L, 9L, 12L, 13L, 14L, 16L, 17L, 5L, 8L, 9L, 10L,
                  12L, 13L, 14L, 15L, 16L, 17L, 18L, 6L, 9L, 10L, 11L, 12L, 13L,
                  14L, 15L, 17L, 18L, 19L, 7L, 10L, 11L, 13L, 14L, 15L, 18L, 19L,
                  8L, 12L, 13L, 16L, 17L, 18L, 9L, 12L, 13L, 14L, 16L, 17L, 18L,
                  19L, 10L, 13L, 14L, 15L, 16L, 17L, 18L, 19L, 11L, 14L, 15L, 17L,
                  18L, 19L),
            p = c(0L, 4L, 8L, 12L, 16L, 22L, 30L, 38L, 44L, 52L,
                  63L, 74L, 82L, 90L, 101L, 112L, 120L, 126L, 134L, 142L, 148L),
            Dim = c(20L, 20L),
            Dimnames = list(NULL, NULL),
            x = c(4.5,
                  -3, -3, 1.5, -3, 4.5, 1.5, -3, -3, 1.5, 4.5, -3, 1.5, -3,
                  -3, 4.5, 7.5, -4.16666666666667, 0.833333333333333, -4.16666666666667,
                  1.25, 0.833333333333333, -4.16666666666667, 6.66666666666667,
                  -3.33333333333333, 0.833333333333333, 1.11111111111111, -2.91666666666667,
                  0.833333333333333, 0.416666666666667, 0.833333333333333,
                  -3.33333333333333, 6.66666666666667, -4.16666666666667, 0.833333333333333,
                  -2.91666666666667, 1.11111111111111, 0.416666666666667, 0.833333333333333,
                  -4.16666666666667, 7.5, 1.25, -4.16666666666667, 0.833333333333333,
                  -4.16666666666667, 1.11111111111111, 6.66666666666667, -2.91666666666667,
                  0.416666666666667, -3.33333333333333, 0.833333333333333,
                  0.833333333333333, 1.25, -2.91666666666667, 0.833333333333333,
                  -2.91666666666667, 6.25, -2.5, 0.416666666666667, 0.833333333333333,
                  -2.5, 0.625, 0.416666666666667, 0.833333333333333, -2.91666666666667,
                  1.25, 0.416666666666667, -2.5, 6.25, -2.91666666666667, 0.625,
                  -2.5, 0.833333333333333, 0.416666666666667, 1.11111111111111,
                  -4.16666666666667, 0.416666666666667, -2.91666666666667,
                  6.66666666666667, 0.833333333333333, -3.33333333333333, 0.833333333333333,
                  0.833333333333333, -3.33333333333333, 0.833333333333333,
                  6.66666666666667, -2.91666666666667, 0.416666666666667, -4.16666666666667,
                  1.11111111111111, 0.416666666666667, 0.833333333333333, -2.5,
                  0.625, -2.91666666666667, 6.25, -2.5, 0.416666666666667,
                  1.25, -2.91666666666667, 0.833333333333333, 0.416666666666667,
                  0.625, -2.5, 0.833333333333333, 0.416666666666667, -2.5,
                  6.25, -2.91666666666667, 0.833333333333333, -2.91666666666667,
                  1.25, 0.833333333333333, 0.833333333333333, -3.33333333333333,
                  0.416666666666667, -2.91666666666667, 6.66666666666667, 1.11111111111111,
                  -4.16666666666667, 0.833333333333333, -4.16666666666667,
                  1.25, 7.5, -4.16666666666667, 0.833333333333333, 0.416666666666667,
                  1.11111111111111, -2.91666666666667, 0.833333333333333, -4.16666666666667,
                  6.66666666666667, -3.33333333333333, 0.833333333333333, 0.416666666666667,
                  0.833333333333333, -2.91666666666667, 1.11111111111111, 0.833333333333333,
                  -3.33333333333333, 6.66666666666667, -4.16666666666667, 0.833333333333333,
                  1.25, -4.16666666666667, 0.833333333333333, -4.16666666666667,
                  7.5),
            factors = list())
    })

    expect_identical({
        n_dims <- c(2, 4)
        phi    <- c(1, 1)
        Q_alpha <- make_Q_alpha_2d(n_dims, phi, use_spam = FALSE)
        tau2 <- c(3, 5)
        make_Q_alpha_tau2(Q_alpha, tau2, use_spam = FALSE)
    },
    {
        new("dgCMatrix",
            i = c(0L, 1L, 2L, 0L, 1L, 3L, 0L, 2L, 3L, 1L,
                  2L, 3L, 4L, 5L, 8L, 4L, 5L, 6L, 9L, 5L, 6L, 7L, 10L, 6L, 7L,
                  11L, 4L, 8L, 9L, 12L, 5L, 8L, 9L, 10L, 13L, 6L, 9L, 10L, 11L,
                  14L, 7L, 10L, 11L, 15L, 8L, 12L, 13L, 16L, 9L, 12L, 13L, 14L,
                  17L, 10L, 13L, 14L, 15L, 18L, 11L, 14L, 15L, 19L, 12L, 16L, 17L,
                  13L, 16L, 17L, 18L, 14L, 17L, 18L, 19L, 15L, 18L, 19L),
            p = c(0L,
                  3L, 6L, 9L, 12L, 15L, 19L, 23L, 26L, 30L, 35L, 40L, 44L, 48L,
                  53L, 58L, 62L, 65L, 69L, 73L, 76L),
            Dim = c(20L, 20L),
            Dimnames = list(NULL, NULL),
            x = c(6, -3, -3, -3, 6, -3, -3, 6, -3, -3, -3,
                  6, 10, -5, -5, -5, 15, -5, -5, -5, 15, -5, -5, -5, 10, -5, -5,
                  15, -5, -5, -5, -5, 20, -5, -5, -5, -5, 20, -5, -5, -5, -5, 15,
                  -5, -5, 15, -5, -5, -5, -5, 20, -5, -5, -5, -5, 20, -5, -5, -5,
                  -5, 15, -5, -5, 10, -5, -5, -5, 15, -5, -5, -5, 15, -5, -5, -5,
                  10),
            factors = list()
        )
    })

    # locs <- matrix(1:20, 10, 2)
    # MRA <- mra_wendland_2d(locs)
    # locs_pred <- matrix(NA, 20, 2)
    # expect_error(mra_wendland_2d_pred(locs, locs_pred, MRA), "locs_pred must be a numeric matrix with N rows and 2 columns")
    #
    # locs_pred <- matrix(1:30, 10, 3)
    # expect_error(mra_wendland_2d_pred(locs, locs_pred, MRA), "locs_pred must be a numeric matrix with N rows and 2 columns")
    #
    # locs_pred <- matrix("11", 10, 2)
    # expect_error(mra_wendland_2d_pred(locs, locs_pred, MRA), "locs_pred must be a numeric matrix with N rows and 2 columns")
    #
    # locs_pred <- 1:10
    # expect_error(mra_wendland_2d_pred(locs, locs_pred, MRA), "locs_pred must be a numeric matrix with N rows and 2 columns")
    #
    # locs <- matrix(1:30, 10, 3)
    # expect_error(mra_wendland_2d_pred(locs, locs_pred, MRA), "locs must be a numeric matrix with N rows and 2 columns")
    #
    # locs <- matrix(1:20, 10, 2)
    # locs_pred <- matrix(1:20, 10, 2)
    # expect_error(mra_wendland_2d_pred(locs, locs_pred, MRA, use_spam = 3.5), "use_spam must be either TRUE or FALSE")
    # expect_error(mra_wendland_2d_pred(locs, locs_pred, MRA, use_spam = NA), "use_spam must be either TRUE or FALSE")
    # expect_error(mra_wendland_2d_pred(locs, locs_pred, MRA, use_spam = "aaa"), "use_spam must be either TRUE or FALSE")
    #
    # class(MRA) <- NULL
    # expect_error(mra_wendland_2d_pred(locs, locs_pred, MRA), 'MRA must be of class "mra_wendland_2d"')
    #
    # class(MRA) <- "XXX"
    # expect_error(mra_wendland_2d_pred(locs, locs_pred, MRA), 'MRA must be of class "mra_wendland_2d"')
})
