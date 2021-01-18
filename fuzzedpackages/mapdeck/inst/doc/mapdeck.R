## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "# ",
  eval = F
)

## ----packages, eval = TRUE, message = FALSE-----------------------------------
library(mapdeck)

## ---- fig.width=6-------------------------------------------------------------
#  key <- 'abc'    ## put your own token here
#  mapdeck(token = key)

## ---- eval = T----------------------------------------------------------------
set_token('abc')
mapdeck_tokens()

## ---- fig.width=6-------------------------------------------------------------
#  mapdeck(token = key, style = 'mapbox://styles/mapbox/dark-v9')

## -----------------------------------------------------------------------------
#  mapdeck_style(style = 'dark')

