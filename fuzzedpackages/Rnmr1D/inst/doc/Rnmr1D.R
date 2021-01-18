## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning=FALSE, fig.align=TRUE)

## ----init-object0, echo=TRUE, include=TRUE-------------------------------
library(Rnmr1D)
data_dir <- system.file("extra", package = "Rnmr1D")
RAWDIR <- file.path(data_dir, "CD_BBI_16P02")
CMDFILE <- file.path(data_dir, "NP_macro_cmd.txt")
SAMPLEFILE <- file.path(data_dir, "Samples.txt")

## ----init-object1, echo=TRUE, include=TRUE-------------------------------
samples <- read.table(SAMPLEFILE, sep="\t", header=T,stringsAsFactors=FALSE)
samples

## ----init-object2, echo=TRUE, include=TRUE-------------------------------
CMDTEXT <- readLines(CMDFILE)
CMDTEXT[grep("^#$", CMDTEXT, invert=TRUE)]

## ----eval0, echo=TRUE, eval=TRUE-----------------------------------------
out <- Rnmr1D::doProcessing(RAWDIR, cmdfile=CMDFILE, samplefile=SAMPLEFILE, ncpu=2)

## ----proc2, echo=TRUE, eval=TRUE-----------------------------------------
ls(out)

out$infos

## ----plot1, echo=TRUE, fig.align='center', fig.width=12, fig.height=6----
plotSpecMat(out$specMat, ppm_lim=c(0.5,5), K=0.33)

## ----plot2, echo=TRUE, fig.align='center', fig.width=12, fig.height=6----
plotSpecMat(out$specMat, ppm_lim=c(0.5,5), K=0, pY=0.1)

## ----plot3, echo=TRUE, fig.align='center', fig.width=12, fig.height=6----
plotSpecMat(out$specMat, ppm_lim=c(0.5,5), K=0.33, asym=0)
cols<- c( rep("blue",length(out$samples$Treatment)));
cols[out$samples$Treatment=="stress"] <- "red"
plotSpecMat(out$specMat, ppm_lim=c(0.5,5), K=0.67, dppm_max=0, cols=cols)

## ----proc3, echo=TRUE, eval=TRUE-----------------------------------------
specMat.new <- Rnmr1D::doProcCmd(out, 
     c( "bucket aibin 10.2 10.5 0.2 3 0", "9.5 4.9", "4.8 0.5", "EOL" ), ncpu=2, debug=TRUE)
out$specMat <- specMat.new

## ----proc4, echo=TRUE, eval=TRUE-----------------------------------------
outMat <- Rnmr1D::getBucketsDataset(out, norm_meth='CSN')
outMat[, 1:10]

## ----proc5a, echo=TRUE, eval=TRUE----------------------------------------
options(warn=-1)
clustcor <- Rnmr1D::getClusters(outMat, method='corr', cval=0, dC=0.003, ncpu=2)

## ----proc5b, echo=TRUE, eval=TRUE----------------------------------------
options(warn=-1)
clusthca <- Rnmr1D::getClusters(outMat, method='hca', vcutusr=0.11)

## ----proc5c, echo=TRUE, eval=TRUE----------------------------------------
clustcor$clustertab[1:20, ]

clustcor$clusters$C3      # same as clustcor$clusters[['C3']]

## ----plot5a, echo=TRUE, fig.align='center', fig.width=12, fig.height=6----
g1 <- ggplotCriterion(clustcor)
#ggplotPlotly(g1, width=820, height=400)
g1
g2 <- ggplotCriterion(clusthca)
#ggplotPlotly(g2, width=820, height=400)
g2

## ----plot5b, echo=TRUE, fig.align='center', fig.width=12, fig.height=6----
layout(matrix(1:2, 1, 2,byrow = TRUE))

hist(simplify2array(lapply(clustcor$clusters, length)), 
     breaks=20, main="CORR", xlab="size", col="darkcyan")
mtext("clusters size distribution", side = 3)

hist(simplify2array(lapply(clusthca$clusters, length)), 
     breaks=20, main="HCA", xlab="size", col="darkcyan")
mtext("clusters size distribution", side = 3)

## ----plot5c, echo=TRUE, fig.align='center', fig.width=12, fig.height=6----
g3 <- ggplotClusters(outMat,clustcor)
#ggplotPlotly(g3, width=820, height=400)
g3
g4 <- ggplotClusters(outMat,clusthca)
#ggplotPlotly(g4, width=820, height=400)
g4


## ----proc6, echo=TRUE, eval=TRUE-----------------------------------------
pca <- prcomp(outMat,retx=TRUE,scale=T, rank=2)
sd <- pca$sdev
eigenvalues <- sd^2
evnorm <- (100*eigenvalues/sum(eigenvalues))[1:10]

## ----plot6a, echo=TRUE, fig.align='center', fig.width=12, fig.height=6----
g5 <- ggplotScores(pca$x, 1, 2, groups=out$samples$Treatment, EV=evnorm , gcontour="polygon")
#ggplotPlotly(g5, width=820, height=450)
g5

## ----plot6b, echo=TRUE, fig.align='center', fig.width=12, fig.height=10----
g6 <- ggplotLoadings(pca$rotation, 1, 2, associations=clusthca$clustertab, EV=evnorm, main=sprintf("Loadings - Crit=%s",clusthca$vcrit), gcontour="ellipse" )
#ggplotPlotly(g6, width=820, height=650)
g6

## ----plot6c, echo=TRUE, fig.align='center', fig.width=12, fig.height=10----
outMat.merged <- Rnmr1D::getMergedDataset(outMat, clusthca, onlycluster=TRUE)
pca.merged <- prcomp(outMat.merged,retx=TRUE,scale=T, rank=2)
g7 <- ggplotLoadings(pca.merged$rotation, 1, 2, associations=NULL, EV=evnorm)
#ggplotPlotly(g7, width=820, height=650)
g7

## ----proc101, echo=TRUE, eval=TRUE---------------------------------------
data_dir <- system.file("extra", package = "Rnmr1D")
RAWDIR <- file.path(data_dir, "CD_BBI_16P02")

## ----proc102, echo=TRUE, eval=TRUE---------------------------------------
procParams <- Spec1rProcpar
procParams$LOGFILE <- ""
procParams$VENDOR <- 'bruker'
procParams$INPUT_SIGNAL <- 'fid'
procParams$LB <- 0.3
procParams$ZEROFILLING <- TRUE
procParams$ZFFAC <- 4
procParams$OPTPHC1 <- TRUE
procParams$TSP <- TRUE

## ----proc103, echo=TRUE, eval=TRUE---------------------------------------
metadata <- generateMetadata(RAWDIR, procParams)
metadata

## ----proc104, echo=TRUE, eval=TRUE---------------------------------------
ID <- 1
ACQDIR <- metadata$rawids[ID,1]
spec <- Spec1rDoProc(Input=ACQDIR,param=procParams)

## ----proc105, echo=TRUE, eval=TRUE---------------------------------------
ls(spec)

## ----plot100, echo=TRUE, fig.align='center', fig.width=9, fig.height=6----
plot( spec$ppm, spec$int, type="l", col="blue", 
                xlab="ppm", ylab="Intensities", 
                xlim=c( spec$pmax, spec$pmin ), ylim=c(0, max(spec$int/100)) )
legend("topleft", legend=metadata$samples[ID,1])

