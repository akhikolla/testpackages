## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, tidy.opts=list(width.cutoff=80), tidy=TRUE, comment=NA)

## ----message = FALSE----------------------------------------------------------
require("regmed")

## ---- loaddat-----------------------------------------------------------------
data(regmed_example)

y <- regmed_example$y
x <- regmed_example$x
med <- regmed_example[, -c(1,2)]

## ---- gridfit-----------------------------------------------------------------
fit.grid <- regmed.grid(x, med, y, lambda.vec= c(seq(from=1, to=0, by = -.1)), frac.lasso=.8)

## ---- methodsgrid-------------------------------------------------------------
plot.regmed.grid(fit.grid)

## ---- fitbest-----------------------------------------------------------------
fit.trim <- trim.best(fit.grid)
summary.regmed(fit.trim)

## ---- plotregmed--------------------------------------------------------------
plot(fit.trim, cex=.6)

## ----med.subset---------------------------------------------------------------
## choose subset of mediators that are in fit.trim
which.med <- colnames(med) %in% dimnames(fit.trim$alpha)[[1]]
med.selected <- med[, which.med]
trauma=x
cortisol=y
fit.lam2 <- regmed.fit(trauma, med.selected, cortisol, lambda = 0.2, frac.lasso=.8)
summary(fit.lam2)

fit.lam0 <- regmed.fit(trauma, med.selected, cortisol, lambda = 0.0, frac.lasso=.8)
summary(fit.lam0)

## ---- plotfit-----------------------------------------------------------------
plot(fit.lam2, lty=2, lwd=3, cex=.7)

## ---- lavaan------------------------------------------------------------------
library(lavaan)
## choose subset of mediators that are in fit.trim
which.med <- colnames(med) %in% dimnames(fit.trim$alpha)[[1]]
med.selected <- med[, which.med]
mediator.names <- dimnames(med.selected)[[2]]

## setup lavaan model
med.model <- lavaan.model(y.name="y", x.name="x",
                          med.name=mediator.names, medcov=fit.trim$MedCov)

## setup data for lavaan
dat <- data.frame( cbind(scale(x, center=TRUE, scale=TRUE),scale(y, center=TRUE, scale=TRUE),scale(med.selected, center=TRUE, scale=TRUE)) )
names(dat) <- c("x","y",dimnames(med.selected)[[2]])

## fit sem with lavaan

fit.lavaan <- sem(model=med.model, data=dat)
summary(fit.lavaan)

