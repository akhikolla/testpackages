## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, tidy.opts=list(width.cutoff=80), tidy=TRUE, comment=NA)

## ----message = FALSE----------------------------------------------------------
require(vcpen)

## ---- loaddat-----------------------------------------------------------------

data(vcexample)
ls()
head(dose)
head(doseinfo)
response[1:10]

## ---- kerns-------------------------------------------------------------------
nvc <- 1+length(unique(doseinfo[,2]))
id <- 1:nrow(dose)

## vcs for genetic kernel matrices
Kerns <- vector("list", length=nvc)
for(i in 1:(nvc-1)){
  ## below uses kernel_linear, but users can replace this with their choice of function to 
  ## create other types of kernel matrices.
  Kerns[[i]] <- kernel_linear(dose[,grep(i, doseinfo[,2])])
  rownames(Kerns[[i]]) <- id
  colnames(Kerns[[i]]) <- id  
}
## vc for residual variance requires identity matrix
Kerns[[nvc]] <- diag(nrow(dose))
rownames(Kerns[[nvc]]) <- id
colnames(Kerns[[nvc]]) <- id


## ---- runvcpen6---------------------------------------------------------------
fit  <- vcpen(response, covmat, Kerns)
summary(fit)

## ---- runvcpen1---------------------------------------------------------------
fit.frac1  <- vcpen(response, covmat, Kerns, frac1 = .1)
summary(fit.frac1)

## ---- vcinit------------------------------------------------------------------
vcinit <- minque(response, covmat, Kerns, n.iter=2)
names(vcinit)
vcinit$beta
vcinit$vc

