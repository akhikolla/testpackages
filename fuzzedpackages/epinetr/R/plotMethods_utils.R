# Split an incidence matrix into matrices per number of vertices per
# edge
splitMatrix <- function(mm) {
  mm[mm > 0] <- 1
  mm[mm < 0] <- 0

  sums <- colSums(mm)
  numNodes <- unique(sort(sums))
  numNodes <- numNodes[numNodes > 1]
  retmat <- list()

  for (i in 1:length(numNodes)) retmat[[i]] <- as.matrix(mm[, which(sums ==
      numNodes[i])])

  return(retmat)
}
