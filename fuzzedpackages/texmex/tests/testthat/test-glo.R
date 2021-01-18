context("glo")

test_that("glo behaves as it should",{
  skip_on_cran()

  set.seed(20101111)

  # large sample estimation, testing Normality of estimators with and without covariates in model
  myTest <- function(Mu, Phi, Xi, ..., label){
    ests <- lapply(1:nreps, function(i) evm(data=data.frame(y=x[,i],x=covariate), y=y,family=glo,...))
    par <- sapply(1:nreps,function(i)ests[[i]]$par)
    se <- sapply(1:nreps,function(i)ests[[i]]$se)
    true <- cbind(Mu,Phi,Xi)
    p.vals <- apply((t(par)-true)/t(se),2,function(x)shapiro.test(x)$p.value)

    expect_gte(min(p.vals), alpha, label=label)
  }

  set.seed(20101111)

  #*************************************************************
  # Set up simulation test info

  ntests <- 21 #total number of parameters estimated in tests below:
  nreps <- 20
  nsim <- 10000
  alpha <- 0.05/ntests # bonferroni
  p <- matrix(runif(6*nreps, -1, 1),ncol=6)
  p[,3] <- p[,3]
  p[,5] <- p[,5]/3
  p[,6] <- p[,6]/4
  covariate <- seq(0.1,0.9,length=nsim)

  #*************************************************************
  # Test stationary model

  x <- sapply(1:nreps,function(i)rglo(nsim,mu=p[i,1],sigma=exp(p[i,3]), xi=p[i,5]))

  myTest(Mu=p[,1], Phi=p[,3], Xi=p[,5], label="glo: no covariates, p-value")

  #*************************************************************
  # Test covariate in mu only

  x <- sapply(1:nreps,function(i)rglo(nsim,mu=p[i,1] + p[i,2]*covariate,sigma=exp(p[i,3]), xi=p[i,5]))

  myTest(Mu=p[,1:2],Phi=p[,3], Xi=p[,5], mu=~x, label="glo: covariates in mu only, p-value")

  #*************************************************************
  # Test covariate in sigma only

  x <- sapply(1:nreps,function(i)rglo(nsim,mu=p[i,1],sigma=exp(p[i,3]+p[i,4]*covariate), xi=p[i,5]))

  myTest(Mu=p[,1],Phi=p[,3:4], Xi=p[,5], phi=~x,label="glo: covariates in phi only, p-value")

  #*************************************************************
  # Test covariate in xi only

  x <- sapply(1:nreps,function(i)rglo(nsim,mu=p[i,1],sigma=exp(p[i,3]), xi=p[i,5]+p[i,6]*covariate))

  myTest(Mu=p[,1],Phi=p[,3], Xi=p[,5:6], xi=~x, label="glo: covariates in xi only, p-value")

  #*************************************************************
  # Test covariate in all variables

  x <- sapply(1:nreps,function(i)rglo(nsim,mu=p[i,1]+ p[i,2]*covariate,sigma=exp(p[i,3]+p[i,4]*covariate), xi=p[i,5]+p[i,6]*covariate))

  myTest(Mu=p[,1:2],Phi=p[,3:4], Xi=p[,5:6], mu=~x,phi=~x,xi=~x, label="glo: covariates in all parameters, p-value")

})

test_that("pglo behaves as it should", {
  skip_on_cran()

  set.seed(20101111)

  probabilities <- runif(10)

  xi.values <- c(0, seq(-0.3, 0.3, length.out=10))

  core.sanity.test <- function(xi) {
    randoms <- sort(rglo(10, 0, 1, xi))
    expect_true(all(diff(pglo(randoms, 0, 1, xi)) >= 0),
                label="pglo: ascending")
    expect_true(all(diff(pglo(randoms, 0,1,xi,lower.tail=FALSE)) <= 0),
                label="pglo: descending")

    expect_equal(log(pglo(randoms, 0, 1,xi)),
                 pglo(randoms, 0, 1, xi, log.p=TRUE),
                 label="pglo: log.p (1)")
    expect_equal(log(pglo(randoms, 0, 1, xi, lower.tail=FALSE)),
                 pglo(randoms, 0, 1, xi, lower.tail=FALSE, log.p=TRUE),
                 label="pglo: log.p (2)")

    mu <- runif(1, -5, 5)
    sigma <- rexp(1)

    expect_equal(pglo(randoms, 0, 1,xi),
                 pglo(mu + sigma * randoms, mu, sigma, xi),
                 label="pglo: shift and scale")
  } # Close core.sanity.tests

  qglo.comparison <- function(xi) {

    quantiles <- qglo(probabilities, 0, 1, xi)
    my.probs  <- pglo(quantiles, 0, 1, xi)
    expect_equal(my.probs, probabilities, label="pglo:straighttest")

    quantiles <- qglo(probabilities, 0, 1, xi, lower.tail=FALSE)
    my.probs  <- pglo(quantiles, 0, 1, xi, lower.tail=FALSE)
    expect_equal(my.probs, probabilities, label="pglo:lowertail")

    my.probs.2 <- pglo(quantiles, 0, 1, xi, lower.tail=TRUE)
    expect_equal(probabilities+my.probs.2,
                 rep(1, length(probabilities)),
                 label="pglo: tail flip")
  } # Close qglo.comparison

  lapply(xi.values, core.sanity.test)
  lapply(xi.values, qglo.comparison)
})

test_that("qglo behaves as it should", {
  skip_on_cran()

  ## get the probabilities that we'll use and sort them
  ## into ascending order for safekeeping
  set.seed(20101111)
  probabilities <- sort(runif(10))

  core.sanity.test <- function(xi) {
    base.quantiles <- qglo(probabilities, 0, 1, xi)
    ## check that the values are ascending
    expect_true(all(diff(base.quantiles)>=0), "qglo:ascendingquantiles")
    ## and check that we're descending correctly for upper tail
    bq2 <- qglo(probabilities, 0, 1, xi, lower.tail=FALSE)
    expect_true(all(diff(bq2)<=0), "qglo:descendingquantiles")
    ## does lower.tail work
    expect_equal(base.quantiles,
                 qglo(1 - probabilities, 0, 1, xi, lower.tail=FALSE),
                 label="qglo: lower.tail works correctly")
    ## does log.p work?
    expect_equal(base.quantiles,
                 qglo(log(probabilities), 0, 1, xi, log.p=TRUE),
                 label="qglo: log.p works")
    ## check shift and scale property
    sigma <- rexp(1)
    mu    <- runif(1, -5, 5)
    shifted <- mu + sigma * base.quantiles
    expect_equal(shifted, qglo(probabilities, mu, sigma, xi),
                 label="qglo: shift and scale")
  } # Close core.sanity.test

  lapply(c(0, seq(-5, 5, length.out=10)), core.sanity.test)

  ## known values
  expect_equal(0, qglo(0.5, 0, 1, 0),
               label="qglo: median match at zero xi")
  xi <- seq(-2, 2, length=10)
  expect_equal(qglo(0.5, 0, 1, xi),
               rep(0,10),
               label="qglo: median match at nonzero xi")
})

test_that("dglo behaves as it should", {
  skip_on_cran()

  local.dglo <- function(x,loc,scale,shape){ # function is not vectorised
    if(shape == 0){
      Ex <- exp(-(x-loc)/scale)
      out <- Ex/(scale* (1+Ex)^2)
    } else {
      Ox <- 1+shape/scale * (x-loc)
      out <- Ox^{-1/shape - 1} / (scale * (1+ Ox^(-1/shape))^2)
    }
    out
  }

  myTest <- function(mu, sig, xi, label){
    myd <- sapply(1:nreps, function(i) dglo(x[,i], mu[i],sig[i], xi[i]))
    ed <- sapply(1:nreps, function(i) local.dglo(x[,i], loc=mu[i], scale=sig[i], shape=xi[i]))
    expect_equal(ed, myd, label=label)
  }

  set.seed(20101111)

  #*************************************************************
  # Test dglo. Note that local.dglo is NOT vectorized.

  nreps <- 100
  nsim <- 1000
  p <- matrix(runif(2*nreps, -1, 1),ncol=2)
  p[, 1] <- p[, 1] + 1
  mu <- rep(runif(nreps,0,1))

  x <- sapply(1:nreps,
              function(i)rglo(nsim,mu=mu[i],sigma=p[i,1], xi=p[i,2]))

  myTest(mu=mu,sig=p[,1], xi=p[,2], label="dglo: random xi")

  #*************************************************************
  # Test dglo when some or all of xi == 0

  p[sample(1:nreps,nreps/2),2] <- 0
  x <- sapply(1:nreps,function(i)rglo(nsim,mu=mu[i],sigma=p[i,1],xi=p[i,2]))
  myTest(mu=mu,sig=p[,1], xi=p[,2], label="dglo: some zero xi")

  p[,2] <-  0
  x <- sapply(1:nreps,function(i)rglo(nsim,mu=mu[i],sigma=p[i,1],xi=p[i,2]))
  myTest(mu=mu,sig=p[,1], xi=p[,2], label="dglo: all zero xi")

  #*************************************************************
  # Test vectorization of dglo.

  mu <- rnorm(nsim)
  sig <- runif(nsim, 0, 2)
  xi <- runif(nsim)

  x <- rglo(nsim, mu, sig, xi)
  myd <- dglo(x, mu, sig, xi)

  ed <- sapply(1:nsim, function(i) local.dglo(x[i], loc=mu[i], scale=sig[i], shape=xi[i]))
  expect_equal(ed, myd, label="dglo:vectorisation")

  #*************************************************************
  # test log.d argument

  ld <- dglo(x,mu,sig,xi,log.d=TRUE)
  expect_equal(myd, exp(ld), label="dglo:logdensity")
}
)

test_that("rglo behaves as it should", {
  skip_on_cran()
  ## so, how do we test an RNG...
  num.simple <- 1000
  num.quantile <- 1e6

  xi.values  <- c(0, seq(-5, 5, length.out=10))
  test.quantiles <- c(0.25, 0.5, 0.75)

  core.sanity.test <- function(xi) {
    seed <- as.integer(runif(1, -1, 1)*(2**30))
    set.seed(seed)
    samples <- rglo(num.simple, 0, 1, xi)
    expect_that(length(samples), equals(num.simple),
                "rglo: output of correct length")
    if (xi > 0) {
      expect_true(all(samples>=-1/xi), "rglo:lowerboundcheck")
    } else if (xi < 0) {
      expect_true(all(samples<=-1/xi), "rglo:upperboundcheck")
    }
    ## scale and shift property
    sigma <- rexp(1)
    mu    <- runif(1, -5, 5)
    shifted <- mu + sigma * samples
    set.seed(seed)
    expect_that(shifted, equals(rglo(num.simple, mu, sigma, xi)), "rglo: scale and shift")
  } # Close core.sanity.test

  quantile.test <- function(xi) {
    ## here are the sampled quantiles
    quantiles <- quantile(pglo(rglo(num.quantile, 0, 1, xi),
                               0, 1, xi),
                          probs=test.quantiles,
                          names=FALSE,na.rm=TRUE)
    ## this is a bit crude, but hey...
    expect_equal(test.quantiles, quantiles, tolerance=0.02, label="rglo: quantile test")
  } # Close quantile.test
  lapply(xi.values, core.sanity.test)
  lapply(xi.values, quantile.test)
}
)
