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
data("multivtoydataset")

## ------------------------------------------------------------------------
K <- 5 # Number of regimes (states)
p <- 3 # Dimension of beta (order of the polynomial regressors)
variance_type <- "heteroskedastic" # "heteroskedastic" or "homoskedastic" model

## ------------------------------------------------------------------------
n_tries <- 1
max_iter <- 1500
threshold <- 1e-6
verbose <- TRUE

## ------------------------------------------------------------------------
mhmmr <- emMHMMR(multivtoydataset$x, multivtoydataset[,c("y1", "y2", "y3")], 
                 K, p, variance_type, n_tries, max_iter, threshold, verbose)

## ------------------------------------------------------------------------
mhmmr$summary()

## ------------------------------------------------------------------------
mhmmr$plot(what = "predicted")

## ------------------------------------------------------------------------
mhmmr$plot(what = "filtered")

## ------------------------------------------------------------------------
mhmmr$plot(what = "regressors")

## ------------------------------------------------------------------------
mhmmr$plot(what = "smoothed")

## ------------------------------------------------------------------------
mhmmr$plot(what = "loglikelihood")

