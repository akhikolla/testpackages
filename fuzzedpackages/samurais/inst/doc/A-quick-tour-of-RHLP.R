## ----setup, include=FALSE------------------------------------------------
library(knitr)
knitr::opts_chunk$set(
	fig.align = "center",
	fig.height = 5.5,
	fig.width = 6,
	warning = FALSE,
	collapse = TRUE,
	dev.args = list(pointsize = 10),
	out.width = "90%",
	par = TRUE
)
knit_hooks$set(par = function(before, options, envir)
  { if (before && options$fig.show != "none") 
       par(family = "sans", mar = c(4.1,4.1,1.1,1.1), mgp = c(3,1,0), tcl = -0.5)
})

## ---- message = FALSE, echo = FALSE--------------------------------------
library(samurais)

## ------------------------------------------------------------------------
data("univtoydataset")

## ------------------------------------------------------------------------
K <- 5 # Number of regimes (mixture components)
p <- 3 # Dimension of beta (order of the polynomial regressors)
q <- 1 # Dimension of w (order of the logistic regression: to be set to 1 for segmentation)
variance_type <- "heteroskedastic" # "heteroskedastic" or "homoskedastic" model

## ------------------------------------------------------------------------
n_tries <- 1
max_iter = 1500
threshold <- 1e-6
verbose <- TRUE
verbose_IRLS <- FALSE

## ------------------------------------------------------------------------
rhlp <- emRHLP(univtoydataset$x, univtoydataset$y, K, p, q, 
               variance_type, n_tries, max_iter, threshold, verbose, verbose_IRLS)

## ------------------------------------------------------------------------
rhlp$summary()

## ------------------------------------------------------------------------
rhlp$plot(what = "regressors")

## ------------------------------------------------------------------------
rhlp$plot(what = "estimatedsignal")

## ------------------------------------------------------------------------
rhlp$plot(what = "loglikelihood")

