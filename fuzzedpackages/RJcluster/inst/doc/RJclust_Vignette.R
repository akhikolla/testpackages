## ---- include = FALSE---------------------------------------------------------
# knitr::opts_chunk$set(
#   collapse = TRUE
# )

## ----setup, echo = FALSE------------------------------------------------------
library(RJcluster)

## -----------------------------------------------------------------------------
data = generateSimulationData()
dim(data$X)

## -----------------------------------------------------------------------------
num_cut_example = sqrt(ncol(data$X))
num_cut_example = floor(num_cut_example)
print(num_cut_example)
clust = RJclust(data$X, num_cut_example)

## -----------------------------------------------------------------------------
clust$G
table(clust$classification)

