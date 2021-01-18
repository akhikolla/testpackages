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
#used_packages <- c("LSD", "RColorBrewer")
#for (p in used_packages) {
#  if (!require(p, character.only = TRUE)) {
#    install.packages(p)
#  }
#}

## ----presvrunmodel, eval=TRUE, echo=FALSE-------------------------------------
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

## ----presvtrunmodel, echo=FALSE, eval=TRUE, dependson="presvrunmodel"---------
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
res_svl <- svlsample(CHF_logret, designmatrix = X)

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

## ----fsvprepdata, echo=TRUE, fig.width=10, fig.height=5-----------------------
library("factorstochvol")
library("zoo")
data("exrates", package = "stochvol")
m <- 6
n <- 1000
y <- 100 * logret(tail(exrates[, seq_len(m)], n + 1))
y <- zoo(y, order.by = tail(exrates$date, n))
plot(y, main = "", xlab = "Time")

## ----preorder, echo=TRUE------------------------------------------------------
preorder(y, factors = 2)

## ----findrestrict, echo=TRUE--------------------------------------------------
findrestrict(y, factors = 2)

## ----runmodel, echo=TRUE------------------------------------------------------
set.seed(1)
res <- fsvsample(y, factors = 2, draws = 10000,
  zeromean = FALSE, thin = 10, quiet = TRUE)

## ----printrres, echo = TRUE---------------------------------------------------
res

## ----covn, echo = TRUE--------------------------------------------------------
dim(cov_n <- covmat(res))

## ----logdetcovn, echo = 2:5, fig.width = 10, fig.height=3.5-------------------
opar <- par(mfrow = c(1, 2), mgp = c(1.7, 0.5, 0), mar = c(3, 3, 1, 1))
logdet <- function (x) log(det(x))
logdet_n <- apply(cov_n[,,,1], 3, logdet)
ts.plot(logdet_n)
acf(logdet_n, main = "")
par(opar)

## ----covess, echo = TRUE------------------------------------------------------
round(apply(cov_n, 1:2, coda::effectiveSize))

## ----corimageplot, echo=2-----------------------------------------------------
opar <- par(mfrow = c(1, 3), xpd = TRUE)
corimageplot(res, these = seq(1, n, length.out = 3), plotCI = "circle",
  plotdatedist = 2, date.cex = 1.1)
par(opar)

## ----voltimeplot, echo=2:3, fig.width = 10, fig.height = 3, cache.rebuild = TRUE----
opar <- par(mgp = c(1.7, 0.5, 0), mar = c(2, 1.5, 1, 0.5))
palette(RColorBrewer::brewer.pal(7, "Dark2")[-5])
voltimeplot(res, legend = "top")
par(opar)

## ----cortimeplot, echo=2:4, fig.width = 10, fig.height = 5, cache.rebuile = TRUE----
opar <- par(mfrow = c(2, 1), mgp = c(1.7, 0.5, 0), mar = c(2, 1.5, 1, 0.5))
palette(RColorBrewer::brewer.pal(6, "Dark2"))
cortimeplot(res, 1)
cortimeplot(res, 2)
par(opar)

## ----comtimeplot, echo=2, fig.height = 6.5------------------------------------
opar <- par(mgp = c(1.7, 0.5, 0), mar = c(3, 3, 1, 1))
comtimeplot(res, maxrows = 6)
par(opar)

## ----loadplot2, fig.width=4.5, fig.height=4.5---------------------------------
opar <- par(mgp = c(1.7, 0.5, 0), mar = c(2.7, 2.7, 2, 0.5))
facloadpairplot(res)
facloadcredplot(res)
par(opar)

## ----varplot, fig.width=10, fig.height=4--------------------------------------
opar <- par(mgp = c(1.7, 0.5, 0), mar = c(2.7, 2.7, 2, 0.5))
logvartimeplot(res, show = "fac")
par(opar)

## ----varplot2, fig.width=7, fig.height=6.5------------------------------------
opar <- par(mgp = c(1.7, 0.5, 0), mar = c(2.7, 2.7, 2, 0.5))
logvartimeplot(res, show = "idi", maxrows = 6)
par(opar)

## ----evdiag, fig.width=10, fig.height=4, echo=2:4, results = 'hide'-----------
opar <- par(mgp = c(1.7, 0.5, 0), mar = c(2.7, 2.7, 2, 0.5))
set.seed(6)
largemodel <- fsvsample(y, factors = 6)
evdiag(largemodel)
par(opar)

## ----predcov1, echo=TRUE------------------------------------------------------
set.seed(4)
predcor1 <- predcor(res)
round(apply(predcor1[,,,1], 1:2, mean), 2)
round(apply(predcor1[,,,1], 1:2, sd), 2)

## ----preddist, fig.height = 6, fig.width = 9, echo = TRUE---------------------
set.seed(5)
predcov_1 <- predcov(res)
effectivedraws <- res$config$draws/res$config$thin
preddraws <- matrix(NA_real_, effectivedraws, m)
for (i in seq_len(effectivedraws))
  preddraws[i,] <- chol(predcov_1[,,i,1]) %*% rnorm(m)
plotlims <- quantile(preddraws, c(0.01, 0.99))
LSD::heatpairs(preddraws, labels = colnames(y),
  cor.cex = 1.5, gap = 0.3, xlim = plotlims, ylim = plotlims)

## ----echo = TRUE--------------------------------------------------------------
set.seed(6)
predloglik(res, matrix(0, nrow = 2, ncol = m), ahead = 1:2, each = 10)

