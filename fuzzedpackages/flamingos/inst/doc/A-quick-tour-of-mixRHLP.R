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
library(flamingos)

## ------------------------------------------------------------------------
data("toydataset")
x <- toydataset$x
Y <- t(toydataset[,2:ncol(toydataset)])

## ------------------------------------------------------------------------
K <- 3 # Number of clusters
R <- 3 # Number of regimes (polynomial regression components)
p <- 1 # Degree of the polynomials
q <- 1 # Order of the logistic regression (by default 1 for contiguous segmentation)
variance_type <- "heteroskedastic" # "heteroskedastic" or "homoskedastic" model

## ------------------------------------------------------------------------
n_tries <- 1
max_iter <- 1000
threshold <- 1e-5
verbose <- TRUE
verbose_IRLS <- FALSE
init_kmeans <- TRUE

## ---- echo=TRUE----------------------------------------------------------
mixrhlp <- emMixRHLP(X = x, Y = Y, K, R, p, q, variance_type, init_kmeans, 
                     n_tries, max_iter, threshold, verbose, verbose_IRLS)

## ------------------------------------------------------------------------
mixrhlp$summary()

## ------------------------------------------------------------------------
mixrhlp$plot(what = "estimatedsignal")

## ------------------------------------------------------------------------
mixrhlp$plot(what = "regressors")

## ------------------------------------------------------------------------
mixrhlp$plot(what = "loglikelihood")

