## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  fig.width=4.5, 
  fig.height=3.5,
  collapse = TRUE,
  comment = "#>"
)

## ----packages------------------------------------------------------------
#nstall.packages("robustbase")
#install.packages("rrcov")
library(ltsspca)
library(rrcov)
library(robustbase)

## ----simulation----------------------------------------------------------
# Clean data
# dataM <- dataSim(n = 200, p = 20, eps = 0)
# Simulation setting 1: outliers which are outlying 
# in the first two variables in the second block
# dataM <- dataSim(n = 200, p = 20, eps = 0.2, setting = "1")
# Simulation setting 2: score outliers
# dataM <- dataSim(n = 200, p = 20, bLength = 4, eps = 0.2, setting = "2")
# Simulation setting 2: orthogonal outliers in Hubert, et al. (2016)
dataM <- dataSim(n = 200, p = 20, bLength = 4, 
                 eps = 0.2, setting = "3", eta = 25)
# get the data
x <- dataM$data
# get the true loading vector
v <- svd(dataM$R)$v[,1:2]

#initial LTS-SPCA
ltsspcaMI <- ltsspca(x = x, kmax = 5, alpha = 0.5)
#reweighted LTS-SPCA
ltsspcaMR <- ltsspcaRw(x = x , obj = ltsspcaMI)
# there are 2 PCs
# compute the angle value
print(Angle(v,ltsspcaMI$loading[,1:2]))
print(Angle(v,ltsspcaMR$loading[,1:2]))
# visualize the loading matrix
matplot(ltsspcaMR$loading[,1:2],type="b",ylab="Loadings",xlab="Variables")
# detect outliers
diagM <-  mydiagPlot(x = x, obj = ltsspcaMR, k =2)

## ----realdata------------------------------------------------------------
data("Glass")
x <- data.matrix(Glass)
#initial LTS-SPCA
ltsspcaMI <- ltsspca(x = x, kmax = 10, alpha = 0.5) ## it takes about 1 minutes
#reweighted LTS-SPCA
ltsspcaMR <- ltsspcaRw(x = x , obj = ltsspcaMI)
# scree plot
plot(ltsspcaMR$eigenvalues,type="b",ylab="Explained Variance",lwd=2)
# select 4 PCs
# visualize the loading matrix
matplot(ltsspcaMR$loadings[,1:4],type="l",lwd=2,ylab="Loadings",xlab="Variables")
# detect outliers
diagM <-  mydiagPlot(x = x, obj = ltsspcaMR, k =4)
#print(diagM$ws.od)
#print(diagM$ws.sd)

