## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=6, 
  fig.height=4,
  fig.align = "center"
)

## ----setup---------------------------------------------------------------
library(MIDASwrappeR)
data("ArtificialDistributionChange")
head(ArtificialDistributionChange)
ArtificialDistributionChange$row <- 1:nrow(ArtificialDistributionChange)

## ----distribution--------------------------------------------------------
hist(subset(ArtificialDistributionChange,ArtificialDistributionChange$row<90000)$dst,freq=F)
hist(subset(ArtificialDistributionChange,ArtificialDistributionChange$row>=90000)$dst,freq=F)

## ----scoring-------------------------------------------------------------
ArtificialDistributionChange$score <- getMIDASScore(ArtificialDistributionChange)

plot(x=tail(ArtificialDistributionChange, 20000)$row, y=tail(ArtificialDistributionChange, 20000)$score, pch=20)

## ----check---------------------------------------------------------------
head(subset(ArtificialDistributionChange, ArtificialDistributionChange$score > 8))

## ----undirected----------------------------------------------------------
ArtificialDistributionChange$score <- getMIDASScore(ArtificialDistributionChange, undirected=T)[c( F, T )]
plot(x=tail(ArtificialDistributionChange, 20000)$row, y=tail(ArtificialDistributionChange, 20000)$score, pch=20)

## ----norelations---------------------------------------------------------
ArtificialDistributionChange$score <- getMIDASScore(ArtificialDistributionChange, norelations=T)
plot(x=tail(ArtificialDistributionChange, 20000)$row, y=tail(ArtificialDistributionChange, 20000)$score, pch=20)

## ----alphaphigh----------------------------------------------------------
ArtificialDistributionChange$score <- getMIDASScore(ArtificialDistributionChange, alpha = .9)
plot(x=tail(ArtificialDistributionChange, 20000)$row, y=tail(ArtificialDistributionChange, 20000)$score, pch=20)

## ----alphahighedges------------------------------------------------------
aggregate(score ~ dst, data=subset(ArtificialDistributionChange,ArtificialDistributionChange$times==907) , max)

## ----alphalow------------------------------------------------------------
ArtificialDistributionChange$score <- getMIDASScore(ArtificialDistributionChange, alpha = .1)
plot(x=tail(ArtificialDistributionChange, 20000)$row, y=tail(ArtificialDistributionChange, 20000)$score, pch=20)

## ----bucket--------------------------------------------------------------
ArtificialDistributionChange$score <- getMIDASScore(ArtificialDistributionChange, alpha = .9, buckets=10)
plot(x=tail(ArtificialDistributionChange, 20000)$row, y=tail(ArtificialDistributionChange, 20000)$score, pch=20)

