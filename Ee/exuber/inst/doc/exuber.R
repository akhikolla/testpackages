## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message=FALSE
)
options(cli.width = 80)

## ----echo=FALSE, message=FALSE------------------------------------------------
library(exuber)

## ----dataset------------------------------------------------------------------
stocks <- aggregate(EuStockMarkets, nfrequency = 52, mean)

## ----estimation---------------------------------------------------------------
est_stocks <- radf(stocks, lag = 1)

## ----summary------------------------------------------------------------------
summary(est_stocks)

## ----diagnostics--------------------------------------------------------------
diagnostics(est_stocks)

## ----datestamp----------------------------------------------------------------
# Minimum duration of an explosive period 
rot = psy_ds(stocks) # log(n) ~ rule of thumb

dstamp_stocks <- datestamp(est_stocks, min_duration = rot)
dstamp_stocks

## ----datestamp-dummy----------------------------------------------------------
dummy <- attr(dstamp_stocks, "dummy")
tail(dummy)

## ----plot-radf, fig.width = 9-------------------------------------------------
autoplot(est_stocks)

## ----plot-datestaemp, fig.width = 7, fig.height=2-----------------------------
datestamp(est_stocks) %>% 
  autoplot()

