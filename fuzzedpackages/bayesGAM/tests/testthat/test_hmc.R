test_that("hmc testing", {
  require(stats); require(graphics); require(SemiPar)
  # Linear regression example
  # linear: identity
  f1 <- bayesGAM(weight ~ np(height), data = women, iter=500,
                   family = gaussian(link="identity"), 
                 chains=1, seed=521)
  cf1 <- coef(f1)
  expect_equal(length(cf1), 8)
  expect_true(sum(abs(cf1)) > 0)
  
  # linear: log
  f2 <- bayesGAM(weight ~ np(height), data = women, iter=500,
                 family = gaussian(link="log"), 
                 chains=1, beta=st(c(35, 0, 4)), seed=521)
  cf2 <- coef(f2)
  expect_equal(length(cf2), 8)
  expect_true(sum(abs(cf2)) > 0)

  # binomial models
  set.seed(651)
  X <- matrix(rnorm(100*5), ncol=5)
  y <- sample(x=c(0,1), size=100, replace=TRUE, prob=c(0.5, 0.5))
  
  bdat <- data.frame(y, X)
  
  # binomial: logit
  f4 <- bayesGAM(y ~ . - 1, data=bdat, chains=1, iter=500,
                 family=binomial(link="logit"),
                    spcontrol = list(qr = TRUE), 
                 seed=521)

  cf4 <- coef(f4)
  expect_equal(length(cf4), 5)
  expect_true(sum(abs(cf4)) > 0)

  # binomial: logit2 qr=FALSE
  f5 <- bayesGAM(y ~ . - 1, data=bdat, chains=1, iter=500,
                 family=binomial(link="logit"),
                 spcontrol = list(qr = FALSE), 
                 seed=521)

  cf5 <- coef(f5)
  expect_equal(length(cf5), 5)
  expect_true(sum(abs(cf5)) > 0)

  # binomial: probit
  f6 <- bayesGAM(y ~ . - 1, data=bdat, chains=1, iter=500,
                 family=binomial(link="probit"),
                 spcontrol = list(qr = TRUE), 
                 seed=521)

  cf6 <- coef(f6)
  expect_equal(length(cf6), 5)
  expect_true(sum(abs(cf6)) > 0)

  # binomial: logit3 with . in formula
  dat<- data.frame(abs(X),
                   y=y)
  f8 <- bayesGAM(y ~ ., data=dat, chains=1, iter=500,
                 family=binomial(link="logit"),
                 spcontrol = list(qr = TRUE), 
                 seed=521)

  cf8 <- coef(f8)
  expect_equal(length(cf8), 6)
  expect_true(sum(abs(cf8)) > 0)

  # binomial: cloglog
  data("PimaIndiansDiabetes2", package = "mlbench")
  PimaIndiansDiabetes2$y <- ifelse(PimaIndiansDiabetes2$diabetes=="pos",1,0)
  f9 <- bayesGAM(y ~ pregnant+pedigree+age, chains=1, 
                 data=PimaIndiansDiabetes2, iter=500,
                 family=binomial(link="cloglog"),
                 spcontrol = list(qr = TRUE), seed=521)

  cf9 <- coef(f9)
  expect_equal(length(cf9), 4)
  expect_true(sum(abs(cf9)) > 0)

  # Poisson: log
  ## Dobson (1990) Page 93: Randomized Controlled Trial :
  counts <- c(18,17,15,20,10,20,25,13,12)
  outcome <- gl(3,1,9)
  treatment <- gl(3,3)
  pdata <- data.frame(counts, outcome, treatment)

  f10 <- bayesGAM(counts ~ outcome + treatment, 
                  data=pdata, chains=1, iter=500,
                  family = poisson(link="log"),
                       spcontrol = list(qr = TRUE), 
                  seed=521)

  cf10 <- coef(f10)
  expect_equal(length(cf10), 5)
  expect_true(sum(abs(cf10)) > 0)

  # bivariate smoothing
  data(scallop)
  set.seed(982)
  f13 <- bayesGAM(log(tot.catch+1) ~ np(longitude, latitude),
                 data=scallop, cores=1, chains=1, iter=500,
                 seed=521)
  
  cf13 <- coef(f13)
  expect_equal(length(cf13), 42)
  expect_true(sum(abs(cf13)) > 0)
  
  # autoregressive
  f15 <- bayesGAM(lh ~ L(lh), family=gaussian, cores=1, chains=1, 
                  seed=521)
  cf15 <- coef(f15)
  expect_equal(length(cf15), 3)
  expect_true(sum(abs(cf15)) > 0)
  
  # posterior predict
  set.seed(981)
  pp <- posterior_predict(f13, draws=50)
  ppdim <- dim(pp@pp[[1]])
  expect_equal(ppdim, c(50, 148))
  expect_equal(sum(is.na(pp@pp)), 0)
  
  # predict
  set.seed(432)
  f <- bayesGAM(weight ~ np(height), data = women,
                family = gaussian, iter=500, chains=1)
  newheights <- with(women, rnorm(10, mean = mean(height)), sd=sd(height))
  women2 <- data.frame(height=newheights)
  
  pred <- predict(f, women2, draws=100)
  
  ppdim <- dim(pred@pp[[1]])
  expect_equal(ppdim, c(100, 10))
  expect_equal(sum(is.na(pred@pp)), 0)
  
  # random intercept
  dat <- data.frame(id = rep(1:5, each=10),
                    y1 = rnorm(50),
                    y2 = rexp(50),
                    x = runif(50))
  
  f16 <- bayesGAM(cbind(y1, y2) ~ np(x), random = ~factor(id), data=dat, 
                chains=1, iter=500)
  
  cf16 <- coef(f16)
  expect_equal(length(cf16), 55)
  expect_true(sum(abs(cf16)) > 0)
})



