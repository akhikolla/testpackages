# example
mr <- fishflux::metabolic_rate(temp = 27, m_max = 600, m = 300, asp = 3,
                         troph = 2, f = 2, growth_g_day = 0.05, B0 = 0.2, a = 0.6 )

test_that("Simple corner cases", {
  expect_gt(min(mr), 0)
  expect_length(mr, 4)
  expect_equal(nrow(mr), 1)
  expect_true(is.numeric(mr$Total_metabolic_rate_C_g_d))
  expect_true(is.numeric(mr$Total_metabolic_rate_j_d))
  expect_true(is.numeric(mr$Resting_metabolic_rate_j_d))
  expect_true(is.numeric(mr$Cost_growth_j_g))
})
