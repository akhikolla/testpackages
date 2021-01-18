## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(PandemicLP)

## ---- fig.asp=2/(sqrt(5)+1), fig.width=8, fig.align='center', warnings = FALSE----
MGdata = load_covid("Brazil","MG")
plot(MGdata)$new

## -----------------------------------------------------------------------------
# Commenting to reduce vignette runtime
# MGestimated = pandemic_model(MGdata, case_type = "deaths", covidLPconfig = TRUE)
# Load precomputed MCMC
download.file("http://github.com/CovidLP/PandemicLP/raw/master/temp/PandemicLP.rda","./PandemicLP.rda")
load("./PandemicLP.rda")

## -----------------------------------------------------------------------------
MGestimated

## ---- fig.asp=2/(sqrt(5)+1), fig.align='center', fig.width=4------------------
traceplot(MGestimated)+theme(legend.position = "")
density(MGestimated)

## -----------------------------------------------------------------------------
MGpredicted = posterior_predict(MGestimated,horizonLong=200)
MGpredicted

## -----------------------------------------------------------------------------
pandemic_stats(MGpredicted)

## ---- fig.asp=2/(sqrt(5)+1), fig.width=8, fig.align='center', warnings = FALSE----
MGplots = plot(MGpredicted,term="both")
MGplots$long
MGplots$short

