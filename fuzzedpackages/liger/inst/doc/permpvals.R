## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- echo=TRUE---------------------------------------------------------------
library(liger)
# load gene set
data("org.Hs.GO2Symbol.list")  
# get universe
universe <- unique(unlist(org.Hs.GO2Symbol.list))
# get a gene set
gs <- org.Hs.GO2Symbol.list[[1]]
# fake dummy example where everything in gene set is perfectly enriched
vals <- rnorm(length(universe), 0, 10)
names(vals) <- universe
set.seed(0)
vals[gs] <- rnorm(length(gs), 10, 10)

## ----test1--------------------------------------------------------------------
set.seed(0)
gsea(values=vals, geneset=gs, plot=TRUE, n.rand=10)

## ----test2--------------------------------------------------------------------
set.seed(0)
gsea(values=vals, geneset=gs, plot=TRUE, n.rand=1e3)

## ----test3--------------------------------------------------------------------
set.seed(0)
gsea(values=vals, geneset=gs, plot=TRUE, n.rand=1e5)

## ---- echo=TRUE---------------------------------------------------------------
sessionInfo()

