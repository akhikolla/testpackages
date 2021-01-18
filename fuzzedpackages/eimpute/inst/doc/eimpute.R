## ---- include = FALSE, eval=TRUE----------------------------------------------
knitr::opts_chunk$set(comment = "#>", warning = FALSE, eval = TRUE, message = FALSE, collapse = TRUE)
library(eimpute)

## ---- eval=FALSE--------------------------------------------------------------
#  install.packages("eimpute")

## ---- eval=FALSE--------------------------------------------------------------
#  library(devtools)
#  install_github("Mamba413/eimpute", build_vignettes = TRUE)

## -----------------------------------------------------------------------------
m <- 6
n <- 5
r <- 3
x_na <- incomplete.generator(m, n, r)
x_na

## -----------------------------------------------------------------------------
x_impute <- eimpute(x_na, r)
x_impute[["x.imp"]]

