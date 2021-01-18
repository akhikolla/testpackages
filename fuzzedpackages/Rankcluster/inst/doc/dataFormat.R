## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library("Rankcluster")

## -----------------------------------------------------------------------------
data(words)
head(words$data)

## -----------------------------------------------------------------------------
data(big4)
head(big4$data)
big4$m

