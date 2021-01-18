test_that("Test suite aaa.R",{  # tests of signature

    expect_silent(signature())
    expect_silent(signature(Inf))

    
    signature(0)
    expect_true(basis(1)*basis(1) == as.scalar(1))
    expect_true(basis(2)*basis(2) == as.scalar(1))
    expect_true(basis(1)*basis(2) ==  basis(1:2))
    expect_true(basis(2)*basis(1) == -basis(1:2))

    signature(1)
    expect_true(basis(1)*basis(1) ==  as.scalar(1))
    expect_true(basis(2)*basis(2) == -as.scalar(1))
    expect_true(basis(1)*basis(2) ==  basis(1:2))
    expect_true(basis(2)*basis(1) == -basis(1:2))
    expect_true(1/basis(1) ==  basis(1))
    expect_true(1/basis(2) == -basis(2))

    signature(-1)  # reduces to Grassman algebra
    expect_true(basis(1)*basis(1) == as.clifford(0))
    expect_true(basis(2)*basis(2) == as.clifford(0))
    expect_true(basis(1)*basis(2) ==  basis(1:2))
    expect_true(basis(2)*basis(1) == -basis(1:2))
    expect_error(1/basis(1))
    expect_error(1/basis(2))

    signature(0)
    expect_true(1/basis(1) == basis(1))
    expect_true(1/basis(2) == basis(2))

    signature(1)
    expect_true(1/basis(1) == +basis(1))
    expect_true(1/basis(2) == -basis(2))

    expect_true(basis(1) %^% basis(2) + basis(2) %^% basis(1) == as.clifford(0))

    signature(0)
    options("show_signature" = TRUE)
    signature(2)

    signature(0)

    x <- clifford(list(6,8,9,22),1:4)
    expect_true(as.1vector(as.vector(x)) == x)

})
