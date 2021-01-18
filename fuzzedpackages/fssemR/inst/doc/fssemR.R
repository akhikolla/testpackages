## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = ">"
)

## ------------------------------------------------------------------------
library(fssemR)
library(network)
library(ggnetwork)
library(Matrix)

## ------------------------------------------------------------------------
n = c(100, 100)    # number of observations in two conditions
p = 20             # number of genes in our simulation
k = 3              # each gene has nonzero 3 cis-eQTL effect
sigma2 = 0.01      # simulated noise variance
prob = 3           # average number of edges connected to each gene
type = "DG"        # `fssemR` also offers simulated ER and directed graph (DG) network
dag  = TRUE        # if DG is simulated, user can select to simulate DAG or DCG
## seed = as.numeric(Sys.time())  # any seed acceptable
seed = 1234        # set.seed(100)
set.seed(seed)
data = randomFSSEMdata2(n = n, p = p, k = p * k, sparse = prob / 2, df = 0.3, 
                        sigma2 = sigma2, type = type, dag = T)

## ---- echo=FALSE, eval=TRUE----------------------------------------------
# genes 1 to 20 are named as g1, g2, ..., g20
rownames(data$Vars$B[[1]]) = colnames(data$Vars$B[[1]]) = paste("g", seq(1, p), sep = "")
rownames(data$Vars$B[[2]]) = colnames(data$Vars$B[[2]]) = paste("g", seq(1, p), sep = "")
rownames(data$Data$Y[[1]]) = rownames(data$Data$Y[[2]]) = paste("g", seq(1, p), sep = "")
names(data$Data$Sk) = paste("g", seq(1, p), sep = "")
# qtl 1 to qtl 60 are named as rs1, rs2, ..., rs60
rownames(data$Vars$F) = paste("g", seq(1, p), sep = "")
colnames(data$Vars$F) = paste("rs", seq(1, p * k), sep = "")
rownames(data$Data$X[[1]]) = rownames(data$Data$X[[2]]) = paste("rs", seq(1, p * k), sep = "")

## ---- fig.align="center", fig.cap = "Simulated GRN under condition 1"----
# data$Vars$B[[1]]    ## simulated GRN under condition 1
GRN_1 = network(t(data$Vars$B[[1]]) != 0, matrix.type = "adjacency", directed = TRUE)
plot(GRN_1, displaylabels = TRUE, label = network.vertex.names(GRN_1), label.cex = 0.5)

## ---- fig.align="center", fig.cap = "Simulated GRN under condition 2"----
# data$Vars$B[[2]]    ## simulated GRN under condition 2
GRN_2 = network(t(data$Vars$B[[2]]) != 0, matrix.type = "adjacency", directed = TRUE)
plot(GRN_2, displaylabels = TRUE, label = network.vertex.names(GRN_2), label.cex = 0.5)

## ---- fig.align="center", fig.cap = "Simulated differential GRN (GRN2 - GRN1), up-regulated are red and down-regulated are blue"----
# data$Vars$B[[2]]    ## simulated GRN under condition 2
diffGRN = network(t(data$Vars$B[[2]] - data$Vars$B[[1]]) != 0, matrix.type = "adjacency", directed = TRUE)
ecol = 3 - sign(t(data$Vars$B[[2]] - data$Vars$B[[1]]))
plot(diffGRN, displaylabels = TRUE, label = network.vertex.names(GRN_2), label.cex = 0.5, edge.col = ecol)

## ------------------------------------------------------------------------
library(Matrix)
print(Matrix(data$Vars$F, sparse = TRUE))

## ------------------------------------------------------------------------
head(data$Data$Y[[1]])
head(data$Data$Y[[2]])

## ------------------------------------------------------------------------
head(data$Data$X[[1]] - 1)
head(data$Data$X[[2]] - 1)

## ------------------------------------------------------------------------
head(data$Data$Sk)

## ------------------------------------------------------------------------
Xs  = data$Data$X     ## eQTL's genotype data
Ys  = data$Data$Y     ## gene expression data
Sk  = data$Data$Sk    ## cis-eQTL indices
gamma = cv.multiRegression(Xs, Ys, Sk, ngamma = 50, nfold = 5, n = data$Vars$n, 
                           p = data$Vars$p, k = data$Vars$k)
fit0   = multiRegression(data$Data$X, data$Data$Y, data$Data$Sk, gamma, trans = FALSE,
                         n = data$Vars$n, p = data$Vars$p, k = data$Vars$k)

## ------------------------------------------------------------------------
fitOpt <- opt.multiFSSEMiPALM2(Xs = Xs, Ys = Ys, Bs = fit0$Bs, Fs = fit0$Fs, Sk = Sk,
                               sigma2 = fit0$sigma2, nlambda = 10, nrho = 10,
                               p = data$Vars$p, q = data$Vars$k, wt = TRUE)

fit <- fitOpt$fit

## ------------------------------------------------------------------------
cat("Power of two estimated GRNs = ", 
    (TPR(fit$Bs[[1]], data$Vars$B[[1]]) + TPR(fit$Bs[[2]], data$Vars$B[[2]])) / 2)
cat("FDR of two estimated GRNs = ", 
    (FDR(fit$Bs[[1]], data$Vars$B[[1]]) + FDR(fit$Bs[[2]], data$Vars$B[[2]])) / 2)
cat("Power of estimated differential GRN = ", 
    TPR(fit$Bs[[1]] - fit$Bs[[2]], data$Vars$B[[1]] - data$Vars$B[[2]]))
cat("FDR of estimated differential GRN = ", 
    FDR(fit$Bs[[1]] - fit$Bs[[2]], data$Vars$B[[1]] - data$Vars$B[[2]]))

## ---- fig.align="center", fig.cap = "estimated differential GRN by fssemR"----
# data$Vars$B[[2]]    ## simulated GRN under condition 2
diffGRN = network(t(fit$Bs[[2]] - fit$Bs[[1]]) != 0, matrix.type = "adjacency", directed = TRUE)
# up-regulated edges are colored by `red` and down-regulated edges are colored by `blue`
ecol = 3 - sign(t(fit$Bs[[2]] - fit$Bs[[1]]))
plot(diffGRN, displaylabels = TRUE, label = network.vertex.names(GRN_2), label.cex = 0.5, edge.col = ecol)

## ------------------------------------------------------------------------
diffGRN = Matrix::Matrix(fit$Bs[[1]] - fit$Bs[[2]], sparse = TRUE)
diffGRN

## ------------------------------------------------------------------------
sessionInfo()

