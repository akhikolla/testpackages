# example
sens <- suppressMessages(suppressWarnings(fishflux::sensitivity(TL = 10, param = list(k_sd = 0.2, Dn_sd = 0.2, Dc_sd = 0.1),
                                                                par = c("k_sd","Dn_sd","Dc_sd"), out = c("Ic", "In", "Ip", "Gc"))))

test_that("Simple corner cases", {
  expect_gt(min(sens), 0)
  expect_length(sens, 4)
  expect_equal(nrow(sens), 3)
  expect_true(is.numeric(sens$Ic_CI))
  expect_true(is.numeric(sens$In_CI))
  expect_true(is.numeric(sens$Ip_CI))
  expect_true(is.numeric(sens$Gc_CI))
})
