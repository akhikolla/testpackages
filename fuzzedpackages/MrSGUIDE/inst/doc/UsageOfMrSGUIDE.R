## ----setup, include=FALSE-----------------------------------------------------
library(knitr)
knitr::opts_chunk$set(echo = TRUE, fig.align='center', fig.width=18, fig.heigh=9)

## ---- eval=FALSE--------------------------------------------------------------
#  install.packages("devtools")

## ---- eval=FALSE--------------------------------------------------------------
#  library(devtools)
#  install_github("BaconZhou/MrSGUIDE")

## ----quick-example-data-------------------------------------------------------
set.seed(1234)

N = 400
np = 3

numX <- matrix(rnorm(N * np), N, np) ## numerical features
gender <- sample(c('Male', 'Female'), N, replace = TRUE)
country <- sample(c('US', 'UK', 'China', 'Japan'), N, replace = TRUE)

z <- sample(c(0, 1), N, replace = TRUE) # Binary treatment assignment

y1 <- numX[, 1] + 1 * z * (gender == 'Female') + rnorm(N)
y2 <- numX[, 2] + 2 * z * (gender == 'Female') + rnorm(N)

train <- data.frame(numX, gender, country, z, y1, y2)
role <- c(rep('n', 3), 'c', 'c', 'r', 'd', 'd')

## ----quick-example-fit--------------------------------------------------------
library(MrSGUIDE)
mrsobj <- MrSFit(dataframe = train, role = role)

## ----quick-example-print------------------------------------------------------
printTree(mrsobj = mrsobj)

## ----quick-example-print-detailF----------------------------------------------
printTree(mrsobj = mrsobj, details = FALSE)

## ----check-packages, eval=FALSE-----------------------------------------------
#  for (pack in c('visNetwork', 'ggplot2')) {
#    if(pack %in% rownames(installed.packages()) == FALSE) {install.packages(pack)}
#  }

## ----quick-example-plot-------------------------------------------------------
plotObj <- plotTree(mrsobj = mrsobj)

## ----quick-example-plot-treeplot----------------------------------------------
plotObj$treeplot

## ----quick-example-plot-nodeTreat---------------------------------------------
plotObj$nodeTreat

## ----quick-example-plot-trtPlot, fig.align='center', fig.width=20-------------
plotObj$trtPlot

## ----predictTree--------------------------------------------------------------
newx <- train[1,]
predictNode <- predictTree(mrsobj = mrsobj, newx, type='node')
predictY <- predictTree(mrsobj = mrsobj, newx, type='outcome')
predictY

## -----------------------------------------------------------------------------
writeTex(mrsobj, file = 'test.tex')

## -----------------------------------------------------------------------------
NCOL(train) == length(role)

## ----bestK-0------------------------------------------------------------------
mrsobj <- MrSFit(dataframe = train, role = role, bestK = 0)
plotObj <- plotTree(mrsobj)
plotObj$treeplot

## ----bestK-1------------------------------------------------------------------
mrsobj <- MrSFit(dataframe = train, role = role, bestK = 1)
plotObj <- plotTree(mrsobj)
plotObj$treeplot

## ----boot---------------------------------------------------------------------
mrsobj <- MrSFit(dataframe = train, role = role, bestK = 1, bootNum = 50, alpha = 0.05)
plotObj <- plotTree(mrsobj)
plotObj$treeplot

## ----boot-nodetreat-----------------------------------------------------------
plotObj$nodeTreat

## ----plot-boot----------------------------------------------------------------
plotObj$trtPlot

## ----maxDepth-1---------------------------------------------------------------
mrsobj <- MrSFit(dataframe = train, role = role, bestK = 1, maxDepth = 1)
plotObj <- plotTree(mrsobj)
plotObj$treeplot

## ----minTrt-30----------------------------------------------------------------
mrsobj <- MrSFit(dataframe = train, role = role, bestK = 0, minTrt = 30)
plotObj <- plotTree(mrsobj)
plotObj$treeplot

## ----minData-50---------------------------------------------------------------
mrsobj <- MrSFit(dataframe = train, role = role, bestK = 0, minData = 50)
plotObj <- plotTree(mrsobj)
plotObj$treeplot

## ----CVFolds-0----------------------------------------------------------------
mrsobj <- MrSFit(dataframe = train, role = role, 
                 bestK = 0, maxDepth = 5, minTrt = 1, 
                 minData = 2, 
                 CVFolds = 0)
plotObj <- plotTree(mrsobj)
plotObj$treeplot

## ----CVSE---------------------------------------------------------------------
mrsobj <- MrSFit(dataframe = train, role = role, 
                 bestK = 0, maxDepth = 5, minTrt = 1, 
                 minData = 2, CVSE = 0.5,
                 CVFolds = 10)
plotObj <- plotTree(mrsobj)
plotObj$treeplot

