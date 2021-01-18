## ----echo=FALSE, include=FALSE----------------------------------------------------------
require(knitr)
options(width = 90, prompt="> ")
knit_hooks$set(fig.mar=function() par(mar=c(5.1,4.1,3.1,2.1)))
opts_chunk$set(out.width='0.4\\textwidth', fig.align='center', highlight=TRUE, comment=NA, fig.height=6, fig.width=6)
opts_knit$set(unnamed.chunk.label='gaston')

## ----prompton, echo=FALSE---------------------------------------------------------------
opts_chunk$set(prompt=TRUE, continue = " ");

## ----promptoff, echo=FALSE--------------------------------------------------------------
opts_chunk$set(prompt=FALSE, continue=" ");

## ----echo=FALSE-------------------------------------------------------------------------
opts_chunk$set(prompt=TRUE, continue = " ");

## ----desc, include=FALSE, echo=FALSE----------------------------------------------------
require(gaston)
desc <- packageDescription("gaston")

## ----fig.mar=TRUE, fig.width=14, fig.height = 7, out.width='0.55\\textwidth'------------
x <- as.bed.matrix(AGT.gen, AGT.fam, AGT.bim)
standardize(x) <- "p"  # needed for matrix product below

K <- GRM(x)
set.seed(17); 
# SNP effects, drown in a normal distribution
u <- rnorm( ncol(x), sd = sqrt(1/ncol(x)) ); 
# Simulated phenotype 
y <- (x %*% u) + rnorm( nrow(x) , sd = 0.7)
# fiting the linear model (note: above simulation is 
# done with tau = sigma2 = 1)
fit <- lmm.diago(y, eigenK = eigen(K), verbose=FALSE )
str(fit)
# retrieving BLUPs for u
BLUP_u <- fit$tau * as.vector(fit$Py %*% x) / (ncol(x) - 1)
# comparison with true effect values
par(mfrow = c(1,2))
plot(u, BLUP_u)
abline(0, 1, col = "red", lty = 3)
# these values allow to recompute the BLUP of omega
plot(x %*% BLUP_u, fit$BLUP_omega)
abline(0, 1, col = "red", lty = 3)

## ---------------------------------------------------------------------------------------
x <- as.bed.matrix(matrix( rbinom(1000*5, 2, 0.5), ncol = 5))
x@snps$chr <- 1  # needed to compute the GRM

## ----fig.mar=TRUE, fig.width=7, fig.height = 7, out.width='0.55\\textwidth'-------------
x <- as.bed.matrix(AGT.gen, AGT.fam, AGT.bim)
standardize(x) <- "p" 
p <-x@p # save the frequency of alleles A2

set.seed(17);
u <- rnorm( ncol(x), sd = sqrt(1/ncol(x)) );
y <- (x %*% u) + rnorm( nrow(x) , sd = 0.7)
# training set : 403 first individuals
I.tr <- 1:403
x.tr <- x[I.tr, ]
x.tr@p <- p # use allele frequencies computed on whole sample
K.tr <- GRM(x.tr)
y.tr <- y[I.tr]
fit <- lmm.diago(y.tr, eigenK = eigen(K.tr), verbose = FALSE)
BLUP_u <- fit$tau * as.vector(fit$Py %*% x.tr) / (ncol(x.tr) - 1)
# prediction on remaining individuals
x1 <- x[-I.tr,]
x1@p <- p # use same allele frequenciesa
# predicted values for y
BLUP_y <- x1 %*% BLUP_u
# compare with simulated value
plot(BLUP_y, y[-I.tr])
abline(0, 1, col = "red", lty = 3)

## ----fig.mar=TRUE, fig.width=7, fig.height = 7, out.width='0.55\\textwidth'-------------
K <- GRM(x)
BLUP_y1 <- fit$tau * K[ -I.tr, I.tr ] %*% fit$Py
plot(BLUP_y, y[-I.tr])
abline(0, 1, col = "red", lty = 3)

## ---------------------------------------------------------------------------------------
x <- as.bed.matrix(AGT.gen, AGT.fam, AGT.bim)
standardize(x) <- "p"  # needed for matrix product below

K <- GRM(x)
set.seed(17); 
# SNP effects, drown in a normal distribution
u <- rnorm( ncol(x), sd = sqrt(1/ncol(x)) ); 
# Simulated phenotype 
y <- (x %*% u) + rnorm( nrow(x) , sd = 0.7)

## ---------------------------------------------------------------------------------------
# fiting the linear model (note: above simulation is 
# done with tau = sigma2 = 1)
fit <- lmm.diago(y, eigenK = eigen(K), verbose=FALSE )

