test_that("test for issue 1", {
    ## the example from issue#1 (and Stackoverflow)
    D <- matrix(c(1:9))
    Numb <- matrix(c(15, 4, 5, 3, 1, 2, 1, 1, 1))
    Number <- data.frame(Numb)
    ## this was giving error on Windows with optimx_2018-7.10
    a_weib <- renewalCount(formula = D ~ 1, data = Number, dist = "weibull",
                           computeHessian = FALSE,  control = renewal.control(trace = 0))
    expect_true(TRUE)
})
