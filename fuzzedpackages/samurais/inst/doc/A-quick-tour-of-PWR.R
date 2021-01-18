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
K <- 5 # Number of segments
p <- 3 # Polynomial degree

## ------------------------------------------------------------------------
pwr <- fitPWRFisher(univtoydataset$x, univtoydataset$y, K, p)

## ------------------------------------------------------------------------
pwr$summary()

## ------------------------------------------------------------------------
pwr$plot(what = "regressors")

## ------------------------------------------------------------------------
pwr$plot(what = "segmentation")

