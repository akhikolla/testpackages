## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(HTLR)
library(bayesplot)

## -----------------------------------------------------------------------------
SEED <- 1234

n <- 510
p <- 2000

means <- rbind(
  c(0, 1, 0),
  c(0, 0, 0),
  c(0, 0, 1),
  c(0, 0, 1),
  c(0, 0, 1),
  c(0, 0, 1),
  c(0, 0, 1),
  c(0, 0, 1),
  c(0, 0, 1),
  c(0, 0, 1)
) * 2

means <- rbind(means, matrix(0, p - 10, 3))

A <- diag(1, p)

A[1:10, 1:3] <-
  rbind(
    c(1, 0, 0),
    c(2, 1, 0),
    c(0, 0, 1),
    c(0, 0, 1),
    c(0, 0, 1),
    c(0, 0, 1),
    c(0, 0, 1),
    c(0, 0, 1),
    c(0, 0, 1),
    c(0, 0, 1)
  )

set.seed(SEED)
dat <- gendata_FAM(n, means, A, sd_g = 0.5, stdx = TRUE)
str(dat)

## -----------------------------------------------------------------------------
# require(corrplot)
cor(dat$X[ , 1:11]) %>% corrplot::corrplot(tl.pos = "n")

## -----------------------------------------------------------------------------
set.seed(SEED)
dat <- split_data(dat$X, dat$y, n.train = 500)
str(dat)

## -----------------------------------------------------------------------------
set.seed(SEED)
system.time(
  fit.t <- htlr(dat$x.tr, dat$y.tr)
)
print(fit.t)

## -----------------------------------------------------------------------------
set.seed(SEED)
system.time(
  fit.t2 <- htlr(X = dat$x.tr, y = dat$y.tr, 
                 prior = htlr_prior("t", df = 1, logw = -20, sigmab0 = 1500), 
                 iter = 4000, init = "bcbc", keep.warmup.hist = T)
)
print(fit.t2)

## -----------------------------------------------------------------------------
summary(fit.t2, features = c(1:10, 100, 200, 1000, 2000), method = median)

## -----------------------------------------------------------------------------
post.t <- as.matrix(fit.t2, k = 2)
## signal parameters
mcmc_intervals(post.t, pars = c("Intercept", "V1", "V2", "V3", "V1000"))

## -----------------------------------------------------------------------------
as.matrix(fit.t2, k = 2, include.warmup = T) %>%
  mcmc_trace(c("V1", "V1000"), facet_args = list("nrow" = 2), n_warmup = 2000)

## -----------------------------------------------------------------------------
y.class <- predict(fit.t, dat$x.te, type = "class")
y.class
print(paste0("prediction accuracy of model 1 = ", 
             sum(y.class == dat$y.te) / length(y.class)))

y.class2 <- predict(fit.t2, dat$x.te, type = "class")
print(paste0("prediction accuracy of model 2 = ", 
             sum(y.class2 == dat$y.te) / length(y.class)))


## -----------------------------------------------------------------------------
predict(fit.t, dat$x.te, type = "response") %>%
  evaluate_pred(y.true = dat$y.te)

