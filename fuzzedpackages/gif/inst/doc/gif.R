## ---- echo=FALSE, message=FALSE, warning=FALSE--------------------------------
knitr::opts_chunk$set(comment = "#", warning = FALSE, eval = TRUE, message = FALSE)
set.seed(1)
library(gif)

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("gif")

## ---- eval=FALSE--------------------------------------------------------------
#  library(devtools)
#  install_github("Mamba413/gif/R-package", build_vignettes = TRUE)

## ---- echo=FALSE--------------------------------------------------------------
library(gif)

## -----------------------------------------------------------------------------
set.seed(1)
n <- 200
p <- 100
Omega <- diag(1, p, p)
for(i in 1:(p - 1)) {
  Omega[i, i + 1] <- 0.5
  Omega[i + 1, i] <- 0.5
}
x <- ggm.generator(n, Omega)

## -----------------------------------------------------------------------------
non_zero_num <- sum(Omega != 0) - p
res <- hgt(x, size = non_zero_num / 2)
Omega_hat <- as.matrix(res[["Omega"]])
head(Omega_hat[, 1:6])
active.entry <- res[["active.entry"]]
head(active.entry)

## ---- eval=FALSE--------------------------------------------------------------
#  non_zero_index <- which(as.matrix(Omega) != 0, arr.ind = TRUE)
#  active.entry <- non_zero_index[which(non_zero_index[,1] < non_zero_index[,2]),]
#  res <- hgt(x, active.entry = active.entry)

## -----------------------------------------------------------------------------
res <- sgt(x, lambda = 0.01)
res[["is.acyclic"]]

