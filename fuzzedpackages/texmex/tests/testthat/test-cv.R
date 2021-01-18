context("cv")

test_that("cv appears to work ok", {
  skip_on_cran()
  skip_on_travis()

  set.seed(1234)

  x <- rgev(1000, xi = .2, mu = 0, sigma = 1)

  g <- evm(x, family = gev)

  set.seed(4321)
  cg <- cv(g, range = seq(.1, 64, length.out = 25))

  expect_true(all(diff(cg$cv$estimate) < 0),
              label = "Parameter estimates shrink with increasing penalty")

  set.seed(4321)
  cvL1 <- cv(g, range = seq(.1, 64, length.out = 25), penalty = "lasso")

  expect_true(all(cg$cv$estimate != cvL1$cv$estimatae),
              label = "Specifying penalty = 'lasso' gives different output")

  x <- rgev(50, xi = .2, mu = 0, sigma = 1)
  g50 <- evm(x, family = gev)
  g50num <-

  cg50 <- cv(g50, range = seq(.1, 64, length.out = 25), folds = 5)

  expect_equal(range(cg50$data$fold), c(1, 5),
               label = "Specifying number of cv folds passes through")

  cg50all <- cv(g50, range = seq(.1, 10, length.out = 25), folds = 50)
  cg50all_2 <- cv(g50, range = seq(.1, 10, length.out = 25), folds = 50)


  expect_equal(cg50all$cv, cg50all_2$cv,
               label = "Leave-one-out CV produces the same results from one run to another")

  cg50_num <- cv(g50, range = seq(.1, 64, length.out = 25), folds = 5)


  # Switch to gpd
  m <- evm(log(ALT.M / ALT.B), data = liver, qu = .7)
  cm <- cv(m, range = seq(.1, 9, length.out = 25), folds = 182)

})
