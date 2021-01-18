% RMKL Readme
# RMKL
Integrating multiple heterogeneous high throughput data sources is an emerging topic of interest in cancer research. Making decisions based upon metabolomic, genomic, etc. data sources can lead to better prognosis or diagnosis than simply using clinical data alone. Support vector machines (SVM) are not suitable for analyzing multiple data sources in a single analysis. SVM employs the kernel trick, thus it is able to construct nonlinear classification boundaries. However, obtaining an optimal kernel type, hyperparameter, and regularization parameter is a challenging task. Cross validation can be employed to make these selections. Ultimately, there may not be one single optimal kernel, but rather a convex combination of several kernel representations of the data. Methods that produce a classification rule based on a convex combination of candidate kernels are referred to as multiple kernel learning (MKL) methods

You can install RMKL from GitHub. If you already have a previous version of RMKL installed, you can use that to install the new version:

```{r}
install.packages("devtools")
library(devtools)
devtools::install_github("cwilso6/RMKL")
library(RMKL)
```
# Requirements
In order for RMKL to work properly, the following packages are required:
* caret
* kernlab

# Examples 
For our benchmark example, we have two groups which are drawn from a bivariate normal distribution where the mean of one groupis fixed and the group means shift to provide different amounts of overlap of the two groups.
This example is meant to illustrate that each MKL implentation selects the correct kernel. For instance, as the hyperparameter gets smaller (using kernlab's parameterization) the resulting decision boundary is almost linear. On the other hand, as the hyperparameter gets larger then decision boundary becomes more jagged and circular. Further details of this are highlighted "Multiple-kernel learning for genomic data mining and prediction" on bioRvix doi: https://doi.org/10.1101/415950.

## Loading data
```{r}
data(benchmark.data)
# The data sets are organized in a a list. Each entry of the list is a 100x3 matrix with each row consisting of a x- and y- coordinate, and a group label (-1,1).
#Below is a summary of the mean of each group for each mean structure.
lapply(1:length(data), function(a) aggregate(x = data[[a]][,1:2], by=list(data[[a]][,3]), mean))
```
## Using RMKL
```{r}
 data.mkl=data[[1]]
 kernels=rep('radial',2)
 sigma=c(2,1/20)
 train.samples=sample(1:nrow(data.mkl),floor(0.7*nrow(data.mkl)),replace=FALSE)
 degree=sapply(1:length(kernels), function(a) ifelse(kernels[a]=='p',2,0))
 #Kernels.gen splts the data into a training and test set, and generates the desired kernel matrices.
 #Here we generate two gaussisan kernel matrices with sigma hyperparameter 2 and 0.05
 K=kernels.gen(data=data.mkl[,1:2],train.samples=train.samples,kernels=kernels,sigma=sigma,degree=degree,scale=rep(0,length(kernels)))
 C=0.05 #Cost parameter for DALMKL
 K.train=K$K.train
 K.test=K$K.test
  
  # parameters set up
  cri_outer = 0.01 # criterion for outer cycle, 0.01 is default by author
  cri_inner = cri_outer/10000  #criterion for inner cycle, this ratio is default by author
  calpha = 10 ### Lagrangian duality constraint parameter, must be positive, 10 is default by author
  max_iter_outer = 500 # maximum number of iterations in outer cycle
  max_iter_inner = 500 # maximum number of iterations in inner cycle
  ytr=data.mkl[train.samples,3]
  k.train=simplify2array(K.train) #Converts list of kernel matrices in to an array with is appropriate for C++ code
  k.test=simplify2array(K.test)
  
  #Implement DALMKL with the hinge loss function
  spicy_svmb1n=SpicySVM(k.train, ytr, C, cri_outer, cri_inner, max_iter_outer, max_iter_inner, calpha)
  spicysvmb1n_results=predictspicy(spicy_svmb1n$alpha,spicy_svmb1n$b, k = k.test)
  cm.DALMKL.svm=confusionMatrix(factor(sign(spicysvmb1n_results),levels=c(-1,1)), factor(data.mkl[-train.samples,3],levels=c(-1,1)))
  cm.DALMKL.svm
  
  #Implement DALMKL with a logistic loss function
  spicy_logib1n=SpicyLogit(k.train, ytr, C, cri_outer, cri_inner, max_iter_outer, max_iter_inner, calpha)
  spicylogib1n_results=predictspicy(spicy_logib1n$alpha,spicy_logib1n$b, k = k.test)
  cm.DALMKL.logi=confusionMatrix(factor(sign(spicylogib1n_results),levels=c(-1,1)), factor(data.mkl[-train.samples,3],levels=c(-1,1)))
  cm.DALMKL.logi
 
  #Convert C parameter from DALMKL implenetation to SimpleMKL and SEMKL implementation to make the four implementations comparible.
  C_SEMKL=C.convert(K.train,spicy_logib1n$model,C)
  
  #Implement SimpleMKL
  SimpleMKL.model=SimpleMKL.classification(k=K.train,data.mkl[train.samples,3], penalty=C_SEMKL)
  cm.SimpleMKL=confusionMatrix(factor(prediction.Classification(SimpleMKL.model,ktest=K.test,data.mkl[train.samples,3])$predict,       levels=c(-1,1)),factor(data.mkl[-train.samples,3],levels=c(-1,1)))
  cm.SimpleMKL
  
  #Implement SEMKL
  SEMKL.model=SEMKL.classification(k=K.train,data.mkl[train.samples,3], penalty=C_SEMKL)
  cm.SEMKL=confusionMatrix(factor(prediction.Classification(SEMKL.model,ktest=K.test,data.mkl[train.samples,3])$predict,
  levels=c(-1,1)),factor(data.mkl[-train.samples,3],levels=c(-1,1)))
  cm.SEMKL
  
  #Selecting a plot in the middle to show the benefit of MKL over SVM
plot(data[[4]][,-3],col=data[[4]][,3]+3,main='Benchmark Data',pch=19,xlab='X1', ylab='X2')

#Using the radial kernel with both hyperparameter individually and then in a combined analysis
C=100
K=kernels.gen(data=data[[4]][,1:2],train.samples=train.samples,kernels=kernels,
              sigma=sigma,degree=degree,scale=rep(0,length(kernels)))
K.train=K$K.train
K.test=K$K.test
#MKL with only one candidate kernel is equivalent to SVM
#SVM with radial hyperparameter 2
rbf2=SEMKL.classification(k = list(K.train[[1]]),outcome = data[[4]][train.samples,3],penalty = C)

#SVM with radial hyperparameter 1/20
rbf.05=SEMKL.classification(k=list(K.train[[2]]),outcome = data[[4]][train.samples,3],penalty = C)
kernels.gen(data=data[[4]][,1:2],train.samples=train.samples,kernels=kernels,sigma=sigma,degree=degree,scale=rep(0,length(kernels)))
domain=seq(1,8,0.1)
grid=cbind(c(replicate(length(domain), domain)),c(t(replicate(length(domain), domain))))
predict.data=rbind(data[[4]][train.samples,1:2],grid)
kernels.predict=kernels.gen(data=predict.data,train.samples=1:length(train.samples),kernels=kernels,
            sigma=sigma,degree=degree,scale=rep(0,length(kernels)))

predict2=prediction.Classification(rbf2, ktest = list(kernels.predict$K.test[[1]]),
                          train.outcome = data[[4]][train.samples,3])$

predict.05=prediction.Classification(rbf.05, ktest = list(kernels.predict$K.test[[2]]),
                                   train.outcome = data[[4]][train.samples,3])

#Contour plot of the predicted values using the model where a single kernel was used
filled.contour(domain,domain, matrix(predict2$predict,length(domain),length(domain)),
               col = colorRampPalette(c('indianred1','lightskyblue'))(2),
               main='Classication Rule Hyperparameter=2', 
               plot.axes={points(data[[4]][,-3],col=data[[4]][,3]+3,pch=18,cex=1.5)})


filled.contour(domain,domain, matrix(predict.05$predict,length(domain),length(domain)),
               col = colorRampPalette(c('indianred1','lightskyblue'))(2),
               main='Classication Rule Hyperparameter=0.05',
               plot.axes={points(data[[4]][,-3],col=data[[4]][,3]+3,pch=18,cex=1.5)})
###################################################################################################
#Use the optimal model with the combination of kernels

predict.combined=prediction.Classification(SEMKL.model[[4]], ktest = kernels.predict$K.test,
                                   train.outcome = data[[4]][train.samples,3])
filled.contour(domain,domain, matrix(predict.combined$predict,length(domain),length(domain)),
               col = colorRampPalette(c('indianred1','lightskyblue'))(2),
               main='Classication Rule MKL', 
               plot.axes={points(data[[4]][,-3],col=data[[4]][,3]+3,pch=18,cex=1.5)})
```
Realizations that fall in the light blue region will be classified as 1, while the points that fall in the light red region will be classified as -1. The points are the original observations. Notice that the two groups do overlap, and that a radial kernel with a large hyperparameter is able to classify in areas with overlap, while a radial kernel with a small hyperparameter can not. The kernel wieghts for this example are 0.9997 for a radial kernel 2 as a hyperparameter, and 0.0002 for radial kernel with 1/20 as a hyper parameter. 

