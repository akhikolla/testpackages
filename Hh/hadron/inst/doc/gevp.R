## ----echo=FALSE---------------------------------------------------------------
library(knitr)
library(hadron)

## ---- eval=FALSE--------------------------------------------------------------
#  Time <- 48
#  correlatormatrix <- cf()
#  for(i in c(1:4)) {
#    tmp <- readbinarycf(files=paste0("corr", i, ".dat"), T=Time)
#    correlatormatrix <- c(correlatormatrix, tmp)
#  }
#  rm(tmp)

## ---- eval=FALSE--------------------------------------------------------------
#  getorderedfilelist <- function(path="./", basename="onlinemeas",
#                                 last.digits=4, ending="")
#  getconfignumbers <- function(ofiles, basename="onlinemeas",
#                               last.digits=4, ending="")
#  getorderedconfigindices <- function(path="./", basename="onlinemeas",
#                                      last.digits=4, ending="")

## -----------------------------------------------------------------------------
data(correlatormatrix)

## ---- cache=TRUE--------------------------------------------------------------
boot.R <- 150
boot.l <- 1
seed <- 1433567
correlatormatrix <- bootstrap.cf(cf=correlatormatrix,
                                 boot.R=boot.R,
                                 boot.l=boot.l,
                                 seed=seed)

## ---- warning=FALSE-----------------------------------------------------------
plot(correlatormatrix, log="y",
     xlab=c("t/a"), ylab="C(t)")

## ---- cache=TRUE--------------------------------------------------------------
t0 <- 4
correlatormatrix.gevp <- bootstrap.gevp(cf=correlatormatrix, t0=t0,
                                        element.order=c(1,2,3,4),
                                        sort.type="values")

## ---- warning=FALSE-----------------------------------------------------------
pc1 <- gevp2cf(gevp=correlatormatrix.gevp, id=1)
pc2 <- gevp2cf(gevp=correlatormatrix.gevp, id=2)
plot(pc1, col="red", pch=21, log="y", xlab="t", ylab="C(t)")
plot(pc2, rep=TRUE, col="blue", pch=22)

## ---- warning=FALSE-----------------------------------------------------------
pc1.matrixfit <- matrixfit(cf=pc1, t1=6, t2=21, useCov=TRUE,
                           parlist=array(c(1,1), dim=c(2,1)),
                           sym.vec=c("cosh"), fit.method="lm")
plot(pc1.matrixfit, do.qqplot=FALSE,
     xlab="t", ylab="C(t)")

## -----------------------------------------------------------------------------
summary(pc1.matrixfit)

## ---- warning=FALSE-----------------------------------------------------------
pc1.matrixfit <- matrixfit(cf=pc1, t1=3, t2=20, useCov=TRUE,
                           parlist=array(c(1,1), dim=c(2,1)),
                           sym.vec=c("cosh"), fit.method="lm",
                           model="pc")
plot(pc1.matrixfit, do.qqplot=FALSE,
     xlab="t", ylab="C(t)")

## ---- warning=FALSE-----------------------------------------------------------
plot(pc1.matrixfit, do.qqplot=FALSE,
     xlab="t", ylab="C(t)", plot.raw=FALSE)
abline(h=1, lty=2)

## -----------------------------------------------------------------------------
pc1.effectivemass <- fit.effectivemass(cf=bootstrap.effectivemass(cf=pc1),
                                       t1=5, t2=20)
plot(pc1.effectivemass, col="red", pch=21, ylim=c(0,1),
     xlab="t", ylab="Meff")

