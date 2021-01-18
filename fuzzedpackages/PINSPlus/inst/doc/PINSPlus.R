## ----message=FALSE------------------------------------------------------------
library(PINSPlus)
data(AML2004)
data <- as.matrix(AML2004$Gene)

## -----------------------------------------------------------------------------
system.time(result <- PerturbationClustering(data = data, verbose = FALSE))

## ---- eval=FALSE--------------------------------------------------------------
#  result <- PerturbationClustering(data = data)

## -----------------------------------------------------------------------------
result$k

## -----------------------------------------------------------------------------
result$cluster

## -----------------------------------------------------------------------------
condition <- seq(unique(AML2004$Group[, 2]))
names(condition) = unique(AML2004$Group[, 2])
plot(prcomp(AML2004$Gene)$x, col = result$cluster, 
     pch = condition[AML2004$Group[, 2]], main = "AML2004")
legend("bottomright", legend = paste("Cluster ", sort(unique(result$cluster)), sep = ""),
        fill = sort(unique(result$cluster)))
legend("bottomleft", legend = names(condition), pch = condition)

## ----eval=FALSE---------------------------------------------------------------
#  result <- PerturbationClustering(data = data, kMax = 5,
#                                   clusteringMethod = "kmeans")

## ----eval=FALSE---------------------------------------------------------------
#  result <- PerturbationClustering(data = data, kMax = 5,
#                                   clusteringMethod = "pam")

## ----eval=FALSE---------------------------------------------------------------
#  result <- PerturbationClustering(data = data, kMax = 5,
#                                   clusteringMethod = "hclust")

## ----eval=FALSE---------------------------------------------------------------
#  result <- PerturbationClustering(
#      data = data,
#      clusteringMethod = "kmeans",
#      clusteringOptions = list(nstart = 100, iter.max = 500),
#      verbose = FALSE
#  )

## ----eval=FALSE---------------------------------------------------------------
#  result <- PerturbationClustering(data = data,
#      clusteringFunction = function(data, k){
#      # this function must return a vector of cluster
#      kmeans(x = data, centers = k, nstart = k*10, iter.max = 2000)$cluster
#  })

## ----eval=FALSE---------------------------------------------------------------
#  result <- PerturbationClustering(data = data,
#                                   perturbMethod = "noise",
#                                   perturbOptions = list(noise = 1.23))

## ----eval=FALSE---------------------------------------------------------------
#  result <- PerturbationClustering(data = data,
#                                   perturbMethod = "noise",
#                                   perturbOptions = list(noisePercent = 10))

## ----eval=FALSE---------------------------------------------------------------
#  result <- PerturbationClustering(data = data,
#                                   perturbMethod = "subsampling",
#                                   perturbOptions = list(percent = 80))

## ----eval=FALSE---------------------------------------------------------------
#  result <- PerturbationClustering(data = data, perturbFunction = function(data){
#      rowNum <- nrow(data)
#      colNum <- ncol(data)
#      epsilon <-
#          matrix(
#              data = rnorm(rowNum * colNum, mean = 0, sd = 1.23456),
#              nrow = rowNum, ncol = colNum
#          )
#  
#      list(
#          data = data + epsilon,
#          ConnectivityMatrixHandler = function(connectivityMatrix, iter, k) {
#              connectivityMatrix
#          }
#      )
#  })

## ----eval=FALSE---------------------------------------------------------------
#  sampleNum <- 50000 # Number of samples
#  geneNum <- 5000 # Number of genes
#  subtypeNum <- 3 # Number of subtypes
#  
#  # Generate expression matrix
#  exprs <- matrix(rnorm(sampleNum*geneNum, 0, 1), nrow = sampleNum, ncol = geneNum)
#  rownames(exprs) <- paste0("S", 1:sampleNum) # Assign unique names for samples
#  
#  # Generate subtypes
#  group <- sort(rep(1:subtypeNum, sampleNum/subtypeNum + 1)[1:sampleNum])
#  names(group) <- rownames(exprs)
#  
#  # Make subtypes separate
#  for (i in 1:subtypeNum) {
#    exprs[group == i, 1:100 + 100*(i-1)] <- exprs[group == i, 1:100 + 100*(i-1)] + 2
#  }
#  
#  # Plot the data
#  library(irlba)
#  exprs.pca <- irlba::prcomp_irlba(exprs, n = 2)$x
#  plot(exprs.pca, main = "PCA")

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(1)
#  t1 <- Sys.time()
#  result <- PerturbationClustering(data = exprs.pca, ncore = 1)
#  t2 <- Sys.time()

## ----eval=FALSE---------------------------------------------------------------
#  t2-t1

## ----eval=FALSE---------------------------------------------------------------
#  result$k

## ----eval=FALSE---------------------------------------------------------------
#  subtype <- result$cluster

## ----eval=FALSE---------------------------------------------------------------
#  if (!require("mclust")) install.packages("mclust")
#  library(mclust)
#  ari <- mclust::adjustedRandIndex(subtype, group)
#  

## ----eval=FALSE---------------------------------------------------------------
#  
#  colors <- as.numeric(as.character(factor(subtype)))
#  
#  plot(exprs.pca, col = colors, main = "Cluster assigments for simulation data")
#  
#  legend("topright", legend = paste("ARI:", ari))
#  
#  legend("bottomright", fill = unique(colors),
#      legend = paste("Group ",
#                     levels(factor(subtype)), ": ",
#                     table(subtype)[levels(factor(subtype))], sep = "" )
#      )

## ----eval=FALSE---------------------------------------------------------------
#  # Load the kidney cancer carcinoma data
#  data(KIRC)
#  # SubtypingOmicsData`'s input data must be a list of
#  # numeric matrices that have the same number of rows:
#  dataList <- list (as.matrix(KIRC$GE), as.matrix(KIRC$ME), as.matrix(KIRC$MI))
#  names(dataList) <- c("GE", "ME", "MI")
#  # Run `SubtypingOmicsData`:
#  result <- SubtypingOmicsData(dataList = dataList)

## ----eval=FALSE---------------------------------------------------------------
#  result <- SubtypingOmicsData(
#      dataList = dataList,
#      clusteringMethod = "kmeans",
#      clusteringOptions = list(nstart = 50)
#  )

## ----eval=FALSE---------------------------------------------------------------
#  library(survival)
#  cluster1=result$cluster1;cluster2=result$cluster2
#  a <- intersect(unique(cluster2), unique(cluster1))
#  names(a) <- intersect(unique(cluster2), unique(cluster1))
#  a[setdiff(unique(cluster2), unique(cluster1))] <-
#      seq(setdiff(unique(cluster2), unique(cluster1))) + max(cluster1)
#  colors <- a[levels(factor(cluster2))]
#  coxFit <- coxph(
#       Surv(time = Survival, event = Death) ~ as.factor(cluster2),
#       data = KIRC$survival,
#       ties = "exact"
#  )
#  mfit <- survfit(Surv(Survival, Death == 1) ~ as.factor(cluster2), data = KIRC$survival)
#  plot(
#       mfit, col = colors, main = "Survival curves for KIRC, level 2",
#       xlab = "Days", ylab = "Survival",lwd = 2
#  )
#  legend("bottomright",
#      legend = paste(
#          "Cox p-value:", round(summary(coxFit)$sctest[3], digits = 5), sep = ""
#      )
#  )
#  legend(
#      "bottomleft",
#      fill = colors,
#      legend = paste("Group ", levels(factor(cluster2)), ": ",
#          table(cluster2)[levels(factor(cluster2))], sep =""
#      )
#  )

