test_that("S4 class glmdisc available", {
  expect_s4_class(methods::new(
    Class = "glmdisc",
    parameters = list(),
    best.disc = list(),
    performance = list(),
    disc.data = data.frame(),
    cont.data = data.frame()
  ), "glmdisc")
})
