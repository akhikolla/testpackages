## ---- echo = FALSE------------------------------------------------------------
knitr::opts_chunk$set(
 fig.width  = 5 ,
 fig.height = 3.5,
 fig.align  = 'center'
)

## -----------------------------------------------------------------------------
library(cellWise)
library(robustHD) # for the TopGear data


## -----------------------------------------------------------------------------
set.seed(1)
X = rnorm(100) 
X[50] = NA
qqnorm(X)

## -----------------------------------------------------------------------------
ML.out <- transfo(X, type = "bestObj", robust=FALSE)
ML.out$lambdahat

## -----------------------------------------------------------------------------
ML.out$objective
ML.out$Xt[45:55] 
ML.out$muhat 
ML.out$sigmahat 
ML.out$Zt[45:55]
qqnorm(ML.out$Zt); abline(0,1)

## -----------------------------------------------------------------------------
plot(ML.out$weights)
ML.out$ttypes 

## -----------------------------------------------------------------------------
RewML.out <- transfo(X, type = "bestObj", robust=TRUE)
RewML.out$lambdahat 
RewML.out$objective 
RewML.out$Xt[45:55]
RewML.out$muhat 
RewML.out$sigmahat 
RewML.out$Zt[45:55] 
qqnorm(RewML.out$Zt); abline(0,1)

## -----------------------------------------------------------------------------
plot(RewML.out$weights) 
X[61] 
RewML.out$ttypes 

## -----------------------------------------------------------------------------
X = exp(X) 


## -----------------------------------------------------------------------------
ML.out <- transfo(X, type = "BC", robust=FALSE)
ML.out$lambdahat 
ML.out$objective 
ML.out$Xt[45:55] 
ML.out$Zt[45:55] 
qqnorm(ML.out$Zt); abline(0,1)

## -----------------------------------------------------------------------------
plot(ML.out$weights) 
ML.out$ttypes 


## -----------------------------------------------------------------------------
RewML.out <- transfo(X, type = "bestObj", robust=TRUE)
RewML.out$lambdahat 
RewML.out$objective 
RewML.out$Xt[45:55]
RewML.out$Zt[45:55] 
qqnorm(RewML.out$Zt); abline(0,1)
RewML.out$ttypes 

## -----------------------------------------------------------------------------
data(TopGear) 

## -----------------------------------------------------------------------------
CD.out   <- checkDataSet(TopGear)
colnames(CD.out$remX)
# remove the subjective variable `Verdict':
X        <- CD.out$remX[,-12] 
carnames <- TopGear[CD.out$rowInAnalysis, 1:2]
X        <- pmax(X, 1e-8) # avoids zeroes

## -----------------------------------------------------------------------------
lambdarange <- c(-4, 6) # for illustration, is the default anyway
ML.out    <- transfo(X, lambdarange = lambdarange, type = "BC", robust=F)
RewML.out <- transfo(X, lambdarange = lambdarange, type = "BC", robust=T)

## -----------------------------------------------------------------------------
ML.out$lambdahat[7] 
RewML.out$lambdahat[7] 

## -----------------------------------------------------------------------------
qqnorm(X[, 7], main = "", cex.lab = 1.5, cex.axis = 1.5)

qqnorm(ML.out$Zt[, 7], main = "Classical transform",
       cex.lab = 1.5, cex.axis = 1.5)
abline(0, 1)

qqnorm(RewML.out$Zt[, 7], main = "Robustly transformed",
       cex.lab = 1.5, cex.axis = 1.5)
abline(0, 1)

## -----------------------------------------------------------------------------
ML.out$lambdahat[8] 
RewML.out$lambdahat[8] 

## -----------------------------------------------------------------------------
qqnorm(X[, 8], main = "Original variable", 
       cex.lab = 1.5, cex.axis = 1.5)

qqnorm(ML.out$Zt[, 8], main = "Classical transform", 
       cex.lab = 1.5, cex.axis = 1.5)
abline(0, 1)

qqnorm(RewML.out$Zt[, 8], main = "Robustly transformed", 
       cex.lab = 1.5, cex.axis = 1.5)
abline(0, 1)


## -----------------------------------------------------------------------------
data("data_glass")

## -----------------------------------------------------------------------------
X <- as.matrix(data_glass[, -c(1:13)])
X <- X[, 1:500]
Z <- scale(X, center=FALSE, robustbase::colMedians(X))
dim(Z)

## -----------------------------------------------------------------------------
ML.out    <- transfo(Z, type = "YJ", robust=F)
RewML.out <- transfo(Z, type = "YJ", robust=T)

## -----------------------------------------------------------------------------
indcells_clas = which(abs(ML.out$Zt) > sqrt(qchisq(0.99, 1)))
indcells_rob  = which(abs(RewML.out$Zt) > sqrt(qchisq(0.99, 1)))
n = dim(ML.out$Zt)[1]
d = dim(ML.out$Zt)[2]; d
nrowsinblock = 5
rowlabels = rep("", floor(n/nrowsinblock));
rowtitle = "" 
ncolumnsinblock = 5
columnlabels = rep("",floor(d/ncolumnsinblock));
columnlabels[3] = "1";
columnlabels[62] = "wavelengths";
columnlabels[floor(d/ncolumnsinblock)] = "500"
columntitle = "" 

CM_clas = cellMap(D = ML.out$Zt, R = ML.out$Zt, 
                   indcells = indcells_clas,
                   indrows = NULL, 
                   rowlabels = rowlabels,
                   columnlabels = columnlabels,
                   mTitle = "YJ transformed variables by ML",
                   rowtitle = rowtitle,
                   columntitle = columntitle,
                   nrowsinblock = nrowsinblock, 
                   ncolumnsinblock = ncolumnsinblock,
                   columnangle = 0, 
                   autolabel = FALSE,
                   colContrast = 1.0,
                   drawCircles = F)
plot(CM_clas)

CM_rob = cellMap(D = RewML.out$Zt, R = RewML.out$Zt, 
                  indcells = indcells_rob,
                  indrows = NULL, 
                  rowlabels = rowlabels,
                  mTitle = "YJ transformed variables by RewML",
                  columnlabels = columnlabels,
                  rowtitle = rowtitle,
                  columntitle = columntitle,
                  nrowsinblock = nrowsinblock, 
                  ncolumnsinblock = ncolumnsinblock,
                  columnangle = 0, 
                  autolabel = FALSE,
                  colContrast = 1.0,
                  drawCircles = F)
plot(CM_rob)

# pdf("Glass_YJ_ML_RewML.pdf",width=10,height=6)
# gridExtra::grid.arrange(CM_clas, CM_rob,ncol=1)
# dev.off()

## -----------------------------------------------------------------------------
data("data_dposs") # in package cellWise
n = nrow(data_dposs); n 
ncol(data_dposs) 

## -----------------------------------------------------------------------------
missmat = is.na(data_dposs)
sizemat = nrow(missmat)*ncol(missmat); sizemat 
100*sum(as.vector(missmat))/sizemat 

missrow = length(which(rowSums(missmat) > 0))
100*missrow/nrow(missmat) 

# Missingness by band:

# F band:
300*sum(as.vector(missmat[,1:7]))/sizemat           
100*length(which(rowSums(missmat[,1:7]) > 0))/20000 

# J band:
300*sum(as.vector(missmat[,8:14]))/sizemat           
100*length(which(rowSums(missmat[,8:14]) > 0))/20000  

# N band:
300*sum(as.vector(missmat[,15:21]))/sizemat           
100*length(which(rowSums(missmat[,15:21]) > 0))/20000  

# So typically the whole band is missing or not.
# We focus on the J band which has the most available rows.

indx = which(rowSums(missmat[,8:14]) ==0)
dpossJ = data_dposs[indx,8:14]
dim(dpossJ) 

## -----------------------------------------------------------------------------
par(mfrow = c(1, 1))
boxplot(scale(dpossJ)) 

## -----------------------------------------------------------------------------
transfoJ_YJ <- transfo(dpossJ, type = "YJ", robust = T)
plot(transfoJ_YJ$lambdahat) 
dpossJ_YJ = transfoJ_YJ$Xt

## -----------------------------------------------------------------------------
DDCPars = list(fastDDC=F,fracNA=0.5)
MacroPCAPars = list(DDCpars=DDCPars,scale=TRUE,silent=T)
MacroPCAdpossJ = MacroPCA(dpossJ,k=4,MacroPCApars=MacroPCAPars) 
MacroPCAdpossJ_YJ = MacroPCA(dpossJ_YJ,k=4,MacroPCApars=MacroPCAPars) 

## -----------------------------------------------------------------------------
MacroPCAdpossJ_YJ$scores[, 1] <- -MacroPCAdpossJ_YJ$scores[, 1] 
cols <- rep("black", dim(MacroPCAdpossJ$scores)[1])
cols[c(98, 10894)] <- "darkorange"
cols[which(MacroPCAdpossJ_YJ$SD > 14)] <- "skyblue3"
cols[which(MacroPCAdpossJ_YJ$SD > 25)] <- "firebrick"

# pdf("dpossJ_scores_YJ.pdf")
pairs(MacroPCAdpossJ_YJ$scores,gap=0,main="",pch=19,col=cols)
# dev.off()

## -----------------------------------------------------------------------------
# pdf("dpossJ_outliermap_YJ.pdf",height=5,width=5)
outlierMap(MacroPCAdpossJ_YJ, title="Outlier map of transformed data", col=cols, labelOut=FALSE)
# dev.off()

# pdf("dpossJ_outliermap_rawdata.pdf",height=5,width=5)
outlierMap(MacroPCAdpossJ, title="Outlier map of raw data", col=cols, labelOut=FALSE)
# dev.off()

