test_that("test .process_formula() works correctly", {
    ## ============================== weibull =================================
    getParNames("weibull") # "scale" "shape"

    .process_formula(y ~ a + b | c + d, "weibull")

    ## =============================== gamma ===================================
    getParNames("gamma") # "rate"  "shape"

    expect_identical(.process_formula(y ~ a + b, "gamma"), NULL)

    abc_gamma <- list(formula = y ~ a + b, anc = list(shape = ~ c))
    expect_identical(.process_formula(y ~ a + b | c, "gamma"), abc_gamma)


    expect_identical(.process_formula(y | shape ~ a + b | c, "gamma"),
                     .process_formula(y ~ a + b | c, "gamma") )

    expect_identical(.process_formula(y | shape ~ a + b, "gamma"),
                     .process_formula(y | shape ~ a + b | a + b, "gamma") )
})
