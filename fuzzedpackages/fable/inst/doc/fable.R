## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.height = 4,
  fig.width = 7
)

## ----setup, message=FALSE-----------------------------------------------------
library(fable)
library(tsibble)
library(dplyr)

## ----data---------------------------------------------------------------------
tourism_melb <- tourism %>%
  filter(Region == "Melbourne")
tourism_melb %>%
  group_by(Purpose) %>%
  slice(1)

## ----plot---------------------------------------------------------------------
tourism_melb %>%
  autoplot(Trips)

## ----mdl----------------------------------------------------------------------
fit <- tourism_melb %>%
  model(
    ets = ETS(Trips ~ trend("A")),
    arima = ARIMA(Trips)
  )
fit

## ----coef---------------------------------------------------------------------
fit %>%
  select(Region, State, Purpose, arima) %>%
  coef()

## ----glance-------------------------------------------------------------------
fit %>%
  glance()

## ----report-------------------------------------------------------------------
fit %>%
  filter(Purpose == "Holiday") %>%
  select(ets) %>%
  report()

## ----augment------------------------------------------------------------------
fit %>%
  augment()

## ----accuracy-----------------------------------------------------------------
fit %>%
  accuracy() %>%
  arrange(MASE)

## ----fc-----------------------------------------------------------------------
fc <- fit %>%
  forecast(h = "5 years")
fc

## ----fc-hilo------------------------------------------------------------------
fc %>%
  hilo(level = c(80, 95))

## ---- eval = FALSE------------------------------------------------------------
#  library(tidyr)
#  fc %>%
#    mutate(interval = hilo(.distribution, 80)) %>%
#    unpack_hilo(interval)

## ----fc-plot, fig.height=10---------------------------------------------------
fc %>%
  autoplot(tourism_melb)

