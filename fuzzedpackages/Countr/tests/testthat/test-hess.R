test_that("test .checkHess()", {

    ## NULL
    expect_false(.checkHess(NULL))
    ## any NA
    expect_false(.checkHess(NA))
    expect_false(.checkHess(c(1, NA)))

    ## matrix
    expect_false(.checkHess(c(1.1, 2.2), 2))
    expect_error(.checkHess(matrix(0, nrow = 2, ncol = 3)),
                 'argument "nPars" is missing, with no default')
    expect_false(.checkHess(matrix(0, nrow = 2, ncol = 3), 2))
    expect_false(.checkHess(diag(c(1.1, 2.2)), 3))
                              
    expect_false(.checkHess(diag(c(1.1, 2.2)), 3))

    expect_true( .checkHess(diag(c(1.1, 2.2)), 2))
    expect_true( .checkHess(matrix(0, nrow = 3, ncol = 3), 3))

                 
})
