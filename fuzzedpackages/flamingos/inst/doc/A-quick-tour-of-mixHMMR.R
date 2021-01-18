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
R <- 3 # Number of regimes/states
p <- 1 # Degree of the polynomial regression
variance_type <- "heteroskedastic" # "heteroskedastic" or "homoskedastic" model

## ------------------------------------------------------------------------
ordered_states <- TRUE
n_tries <- 1
max_iter <- 1000
init_kmeans <- TRUE
threshold <- 1e-6
verbose <- TRUE

## ---- echo=TRUE----------------------------------------------------------
mixhmmr <- emMixHMMR(X = x, Y = Y, K, R, p, variance_type, ordered_states, 
                     init_kmeans, n_tries, max_iter, threshold, verbose)

## ------------------------------------------------------------------------
mixhmmr$summary()

## ------------------------------------------------------------------------
mixhmmr$plot()

