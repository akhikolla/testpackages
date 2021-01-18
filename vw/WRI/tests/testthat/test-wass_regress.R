test_that("function wass_regress works for a trivial example", {

        p1 = pnorm(2.5)
        p0 = pnorm(-2.5)
        t_vec = unique(c(seq(0, 0.05, 0.001), seq(0.05, 0.95, 0.05), seq(0.95, 1, 0.001)))

        # Q0 is the standard normal distribution shifted 5 units to left, truncated (-2.5, 2.5) - 5
        quantile_curve_0 = qnorm((p1-p0)*t_vec + p0) - 5
        # when the single predictor x = 0
        quantile_curve_0_plus1 = quantile_curve_0 + 1
        quantile_curve_0_minus1 = quantile_curve_0 - 1

        # Q0 is the standard normal distribution shifted 5 units to right, truncated (-2.5, 2.5) + 5
        quantile_curve_1 = qnorm((p1-p0)*t_vec + p0) + 5
        # when the single predictor x = 1
        quantile_curve_1_plus1 = quantile_curve_1 + 1
        quantile_curve_1_minus1 = quantile_curve_1 - 1

        Ymat = matrix(c(quantile_curve_0_plus1,
                        quantile_curve_0,
                        quantile_curve_0_minus1,
                        quantile_curve_1_plus1,
                        quantile_curve_1,
                        quantile_curve_1_minus1), nrow = 6, byrow = TRUE)
        Xfit_df = data.frame(X1 = rep(c(0, 1), each = 3))

        res = wass_regress(rightside_formula = ~., Xfit_df = Xfit_df,
                            Ytype = 'quantile', Ymat = Ymat, Sup = t_vec)

        expect_equal(wass_R2(res), 0.974, tolerance = 1e-2)
        expect_equal(quantile_curve_0, res$Qfit[1, ], tolerance = 1e-5)
        expect_equal(quantile_curve_1, res$Qfit[4, ], tolerance = 1e-5)

})
