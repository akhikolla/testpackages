context("test-getpdfs.R")

tt <- seq(0, 5, .01)
pars <- structure(c(0.8, 2, 0.2, 0.5, 0.5, 0.5, 0, 0.8, 3, 0.2, 0.5, 0.5, 0.5, 0, 0.8, 4, 0.2, 0.5, 0.5, 0.5, 0),
                  .Dim = c(7L, 3L))

ncondition <- NCOL(pars)
mm <- matrix(0, ncondition * 2, ncondition)
mm[1:dim(mm)[1L] + dim(mm)[1L] * rep(1:dim(mm)[2L] - 1L, each = 2)] <- 1

test_that("getPdf == getPdfC (1)", {

  pars.list <- unlist(apply(pars, 2, list), recursive = FALSE)
  m1 <- DstarM:::getPdf(pars.list = pars.list, tt = tt, DstarM = FALSE, mm = mm,
                        oscPdf = FALSE, fun.density = Voss.density, args.density = list())

  m2 <- DstarM:::getPdfC(tt = tt, pars = pars, mm = mm, DstarM = FALSE, oscPdf = FALSE, precision = 3)

  expect_equal(m1, m2, label = "test")
})

pars2 <- pars + runif(21, -.1, .1)

test_that("getPdf == getPdfC (2)", {

  pars.list <- unlist(apply(pars, 2, list), recursive = FALSE)
  m1 <- DstarM:::getPdf(pars.list = pars.list, tt = tt, DstarM = FALSE, mm = mm,
                        oscPdf = FALSE, fun.density = Voss.density, args.density = list())

  m2 <- DstarM:::getPdfC(tt = tt, pars = pars, mm = mm, DstarM = FALSE, oscPdf = FALSE, precision = 3)

  expect_equal(m1, m2, label = "test")

})
