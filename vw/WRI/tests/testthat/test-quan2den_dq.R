test_that("function quan2den_dq works for standard normal distribution", {

        p1 = pnorm(2.5)
        p0 = pnorm(-2.5)
        t_vec = unique(c(seq(0, 0.05, 0.001), seq(0.05, 0.95, 0.05), seq(0.95, 1, 0.001)))

        # closed form curves
        quantile_curve = qnorm((p1-p0)*t_vec + p0)
        quantile_density_curve = (p1 - p0)/dnorm(quantile_curve)
        quantile_density_prime_curve = quantile_density_curve^3 * quantile_curve * dnorm(quantile_curve)/(p1-p0)
        density_curve = dnorm(quantile_curve)/(p1-p0)

        # give an equally quantile density curve n-by-m (2-by-101) matrix
        quantileCurves = matrix(rep(qnorm((p1-p0)*t_vec + p0), 2), nrow = 2, byrow = TRUE)
        res = quan2den_qd(quantileCurves, t_vec)

        # mean(abs(quantile_curve)) is 1.729919
        expect_equal(quantile_curve,res$Qobs[1, ], tolerance = 0.001)
        # mean(quantile_density_curve) is 16.27296
        expect_equal(quantile_density_curve, res$qobs[1, ], tolerance = 0.05)
        # mean(abs(quantile_density_prime_curve)) is 871.4039
        expect_equal(quantile_density_prime_curve, res$qobs_prime[1, ], tolerance = 0.05)
        # mean(density_curve) is 0.105482
        expect_equal(density_curve,res$fobs[1, ], tolerance = 1e-3) # mean is 16.27296

})
