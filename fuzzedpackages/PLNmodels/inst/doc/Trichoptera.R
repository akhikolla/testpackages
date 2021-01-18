## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  screenshot.force = FALSE, 
  echo = TRUE,
  rows.print = 5,
  message = FALSE,
  warning = FALSE)

## ----packages-----------------------------------------------------------------
library(PLNmodels)

## ----trichoptera, echo = FALSE, fig.align='center', fig.cap = "Macronema Zebratum captured by Y. Dubuc at Donacona (Québec), 06-20-2001."----
knitr::include_graphics("figures/macronema_zebratum.jpg")

## ----data_load----------------------------------------------------------------
data(trichoptera)

## ----prepare_data-------------------------------------------------------------
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)

## ----data_str-----------------------------------------------------------------
str(trichoptera)

## ----responses_trichoptera----------------------------------------------------
trichoptera$Abundance %>% head() %>% knitr::kable()

## ----count diplay, fig.width = 7, fig.cap = "log-counts in the trichoptera data set"----
corrplot::corrplot(
  t(log(1 + trichoptera$Abundance)),
  is.corr = FALSE,
  addgrid.col = NA
)

## ----covariates_trichoptera---------------------------------------------------
dplyr::select(trichoptera, -Offset, -Abundance) %>% head() %>% knitr::kable()

## ----offset_trichopera--------------------------------------------------------
trichoptera$Offset

