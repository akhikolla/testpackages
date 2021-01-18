## ----setup, include=FALSE-----------------------------------------------------
library(knitr)
opts_chunk$set(cache = TRUE)
options(digits = 3)


## ----eu-stock-markets---------------------------------------------------------
ret <- diff(log(EuStockMarkets)) * 100
plot(ret)

## ----pairs-metrics------------------------------------------------------------
library(LambertW)  # this will load the `moments` package as well
data_metrics <- function(x) {
  c(mean = mean(x), sd = sd(x), min = min(x), max = max(x), 
    skewness = skewness(x), kurtosis = kurtosis(x))
}
ret.metrics <- t(apply(ret, 2, data_metrics))
ret.metrics

## ----gaussianize-returns------------------------------------------------------
library(LambertW)
ret.gauss <- Gaussianize(ret, type = "h", method = "IGMM")
colnames(ret.gauss) <- gsub("\\.X", "", colnames(ret.gauss))

plot(ret.gauss)

## ----metrics-gaussianized-----------------------------------------------------
ret.gauss.metrics <- t(apply(ret.gauss, 2, data_metrics))
ret.gauss.metrics

## ----DAX-FTSE, fig.width = 16, fig.height = 8---------------------------------
layout(matrix(1:2, ncol = 2, byrow = TRUE))
plot(ret[, "FTSE"], ret[, "DAX"])
grid()
plot(ret.gauss[, "DAX"], ret.gauss[, "FTSE"])
grid()

## ----fit-models, include = TRUE, eval = FALSE---------------------------------
#  # try these models on your own
#  mod <- lm(FTSE ~ DAX + SMI + CAC, data = ret)
#  mod.robust <- rlm(FTSE ~ DAX + SMI + CAC, data = ret)
#  mod.gauss <- lm(FTSE ~ DAX + SMI + CAC, data = ret.gauss)
#  
#  summary(mod)
#  summary(mod.robust)
#  summary(mod.gauss)

## ----vars-ret-----------------------------------------------------------------
library(vars)
mod.vars <- vars::VAR(ret[, c("DAX", "CAC")], p = 6)
causality(mod.vars, "DAX")$Granger
causality(mod.vars, "CAC")$Granger


## ----vars-ret-gauss-----------------------------------------------------------
mod.vars.gauss <- vars::VAR(ret.gauss[, c("DAX", "CAC")], p = 6)
causality(mod.vars.gauss, "DAX")$Granger
causality(mod.vars.gauss, "CAC")$Granger

