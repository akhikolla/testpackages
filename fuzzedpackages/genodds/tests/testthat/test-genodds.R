context("Input Sanitation for genodds()")

test_that("Error if response is empty", {
  expect_error(
    genodds(NULL,alteplase$treat)
  )
})

test_that("Error if group is empty", {
  expect_error(
    genodds(alteplase$mRS,NULL)
  )
})


test_that("Error if response and group are of different length", {
  expect_error(
    genodds(alteplase$mRS[1:100],alteplase$treat)
  )
})


test_that("Error if response and strata are of different length", {
  expect_error(
    genodds(alteplase$mRS,alteplase$treat,alteplase$time[100])
  )
})

test_that("Error if treat contains more than two groups", {
  df <- alteplase
  df$treat <- as.character(df$treat)
  df$treat[c(100,600,651)] <- "Other"

  expect_error(
    genodds(df$mRS,df$treat,df$time)
  )
})

test_that("Warning if response contains NA", {
  df <- alteplase
  df$mRS[c(100,600,651)] <- NA

  expect_warning(
    genodds(df$mRS,df$treat,df$time)
  )
})

test_that("Warning if group contains NA", {
  df <- alteplase
  df$treat[c(100,600,651)] <- NA

  expect_warning(
    genodds(df$mRS,df$treat,df$time)
  )
})


test_that("Warning if strata contains NA", {
  df <- alteplase
  df$time[c(100,600,651)] <- NA

  expect_warning(
    genodds(df$mRS,df$treat,df$time)
  )
})


context("Calculation tests for genodds()")

test_that("pooled odds equal to stratum odds when stratum=NULL", {

  obj <- genodds(alteplase$mRS,alteplase$treat)

  expect_equal(exp(obj$pooled_lnodds), obj$results[[1]]$odds)
})

test_that("pooled p value equal to stratum p when stratum=NULL", {

  obj <- genodds(alteplase$mRS,alteplase$treat)

  expect_equal(obj$pooled_p, obj$results[[1]]$p)
})


context("Input Sanitation for genodds.power()")

test_that("Error if neither N nor power are specified",{

          p0 <- c(0.224,0.191,0.082,0.133,0.136,0.043,0.191)
          p1 <- c(0.109,0.199,0.109,0.120,0.194,0.070,0.200)

          expect_error(
            genodds.power(p0,p1)
          )
          })


test_that("Error if power isn't in [0,1]",{

  p0 <- c(0.224,0.191,0.082,0.133,0.136,0.043,0.191)
  p1 <- c(0.109,0.199,0.109,0.120,0.194,0.070,0.200)

  expect_error(
    genodds.power(p0,p1,power=-1)
  )

  expect_error(
    genodds.power(p0,p1,power=2)
  )

})

test_that("Error if N isn't >=2",{

  p0 <- c(0.224,0.191,0.082,0.133,0.136,0.043,0.191)
  p1 <- c(0.109,0.199,0.109,0.120,0.194,0.070,0.200)

  expect_error(
    genodds.power(p0,p1,N=1)
  )
})

test_that("Error p0 and p1 are different lengths",{

  p0 <- c(0.224,0.191,0.082,0.133,0.136,0.043,0.191)
  p1 <- c(0.109,0.199,0.109,0.120,0.194,0.070)

  expect_error(
    genodds.power(p0,p1,N=1)
  )
})


context("Calculation Tests for genodds.power()")

test_that("Required sample size is conservative", {

  p0 <- c(0.224,0.191,0.082,0.133,0.136,0.043,0.191)
  p1 <- c(0.109,0.199,0.109,0.120,0.194,0.070,0.200)
  power_in <- seq(0.1,0.9,0.1)
  N <- genodds.power(p0,p1,power=power_in)

  power_out <- genodds.power(p0,p1,N=N)

  expect_true(prod(power_out>power_in)==1)
})

test_that("Power at sample size within 5e-2 of input power", {

  p0 <- c(0.224,0.191,0.082,0.133,0.136,0.043,0.191)
  p1 <- c(0.109,0.199,0.109,0.120,0.194,0.070,0.200)
  power_in <- seq(0.1,0.9,0.1)
  N <- genodds.power(p0,p1,power=power_in)

  power_out <- genodds.power(p0,p1,N=N)

  expect_true(prod((power_out-power_in)<5e-2)==1)
})

test_that("Theoretical power within 1% of empirical power", {

  set.seed(1908141806)

  p0 <- c(0.224,0.191,0.082,0.133,0.136,0.043,0.191)
  p1 <- c(0.109,0.199,0.109,0.120,0.194,0.070,0.200)

  power_in <- 0.8
  alpha=0.05
  N <- genodds.power(p0,p1,power=power_in)

  empiricalPower <- sapply(1:1e4,function(x){
    response <- c(
                  sample(0:6,size = ceiling(N/2),prob = p0,replace=TRUE),
                  sample(0:6,size = ceiling(N/2),prob = p1,replace=TRUE)
    )

    group <- c(
      rep(0,ceiling(N/2)),
      rep(1,ceiling(N/2))
    )

    obj <- genodds(response,group,alpha=alpha)$pooled_p<=alpha
  })

  empiricalPower <- mean(empiricalPower)

  expect_lte(abs(power_in - empiricalPower),1e-2)
})

