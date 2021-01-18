## ----setup, include=FALSE, cache=FALSE----------------------------------------
knitr::render_sweave()
knitr::opts_chunk$set(prompt = TRUE,
                      fig.show = "hide",
                      warning = FALSE,
                      error = FALSE,
                      message = FALSE,
                      echo = FALSE,
                      cache = TRUE,
                      fig.path = "Figures/article-")
base::options(continue = "+  ", prompt = "R> ")

## ----presvrunmodel, eval=TRUE, results='hide'---------------------------------
set.seed(1)
library("stochvol")
data("exrates")
ind <- which(exrates$date >= as.Date("2008-03-01") &
  exrates$date <= as.Date("2012-03-01"))
CHF_price <- exrates$CHF[ind]

## ----svrunmodel, eval=FALSE, echo=TRUE----------------------------------------
#  set.seed(1)
#  library("stochvol")
#  data("exrates")
#  ind <- which(exrates$date >= as.Date("2008-03-01") &
#    exrates$date <= as.Date("2012-03-01"))
#  CHF_price <- exrates$CHF[ind]
#  res_sv <- svsample(CHF_price, designmatrix = "ar1")

## ----presvtrunmodel, results='hide', eval=TRUE, dependson="presvrunmodel"-----
set.seed(2)
CHF_logret <- 100 * logret(CHF_price)

## ----svtrunmodel, echo=TRUE, eval=FALSE---------------------------------------
#  set.seed(2)
#  CHF_logret <- 100 * logret(CHF_price)
#  res_svt <- svtsample(CHF_logret, designmatrix = "ar0")

## ----svlrunmodel, echo=TRUE, eval=TRUE, dependson="presvtrunmodel"------------
set.seed(3)
X <- cbind(constant = 1,
  100 * logret(exrates$USD[ind]),
  100 * logret(exrates$JPY[ind]))
res_svl <- svlsample(CHF_logret, designmatrix = X, thin = 10)

## ----svlplot, echo=TRUE, dependson="svlrunmodel", fig.height=4----------------
plot(res_svl, showobs = FALSE,
     dates = exrates$date[ind[-1]])

## ----svlbetaplot, echo=2, dependson="svlrunmodel", fig.width=6.7, fig.height=3.5----
opar <- par(mar = c(2.5, 1.5, 0.5, 0.5), mfrow = c(3, 2), mgp = c(1.7, 0.5, 0))
for (i in seq_len(3)) {
  coda::traceplot(svbeta(res_svl)[, i])
  coda::densplot(svbeta(res_svl)[, i], show.obs = FALSE)
}
par(opar)

## ----printsummary, echo=TRUE, eval=TRUE, results="markup"---------------------
summary(res_svl, showlatent = FALSE)

## ----svlpredict, echo=TRUE, eval=TRUE-----------------------------------------
set.seed(4)
pred_ind <- seq(tail(ind, 1), length.out = 25)
pred_X <- cbind(constant = 1,
  100 * logret(exrates$USD[pred_ind]),
  100 * logret(exrates$JPY[pred_ind]))
pred_svl <- predict(res_svl, 24, newdata = pred_X)

## ----plotsvlpred, echo=TRUE, eval=TRUE, fig.height=3.5, fig.width=10----------
opar <- par(mgp = c(1.7, 0.5, 0))
obs_CHF <- 100 * logret(exrates$CHF[pred_ind])
ts.plot(cbind(t(apply(predy(pred_svl), 2, quantile, c(0.05, 0.5, 0.95))),
  obs_CHF), xlab = "Periods ahead", lty = c(rep(1, 3), 2),
  col = c("gray80", "black", "gray80", "red"))
par(opar)

## ----svroll, echo=TRUE, eval=FALSE--------------------------------------------
#  set.seed(5)
#  res <- svsample_roll(CHF_logret, n_ahead = 1,
#    forecast_length = 30,
#    refit_window = "moving",
#    calculate_quantile = c(0.01, 0.05),
#    calculate_predictive_likelihood = TRUE)

## ----printpriordefault, echo=TRUE, eval=FALSE---------------------------------
#  svsample(CHF_logret, priormu = c(0, 100), priorphi = c(5, 1.5),
#    priorsigma = 1, priorbeta = c(0, 10000))
#  svtsample(CHF_logret, priormu = c(0, 100), priorphi = c(5, 1.5),
#    priorsigma = 1, priorbeta = c(0, 10000), priornu = 0.1)
#  svlsample(CHF_logret, priormu = c(0, 100), priorphi = c(5, 1.5),
#    priorsigma = 1, priorbeta = c(0, 10000), priorrho = c(4, 4))
#  svtlsample(CHF_logret, priormu = c(0, 100), priorphi = c(5, 1.5),
#    priorsigma = 1, priorbeta = c(0, 10000), priornu = 0.1,
#    priorrho = c(4, 4))

## ----printpriorspecdefault, echo=TRUE, eval=FALSE-----------------------------
#  ps <- specify_priors(
#    mu = sv_normal(mean = 0, sd = 100),
#    phi = sv_beta(shape1 = 5, shape2 = 1.5),
#    sigma2 = sv_gamma(shape = 0.5, rate = 0.5),
#    nu = sv_infinity(),
#    rho = sv_constant(0),
#    latent0_variance = "stationary",
#    beta = sv_multinormal(mean = 0, sd = 10000, dim = 1))
#  svsample(CHF_logret, priorspec = ps)

## ----eval=FALSE---------------------------------------------------------------
#  y <- svsim(50)$y
#  svsample(y, expert = list(correct_model_misspecification = TRUE))

