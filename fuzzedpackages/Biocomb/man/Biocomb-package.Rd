\name{Biocomb-package}
\alias{Biocomb-package}
\alias{Biocomb}
\docType{package}
\title{
Tools for Data Mining
}
\description{
  Functions to make the data analysis with the emphasis on biological data.
They can deal with both numerical and nominal features. Biocomb includes functions for several feature ranking, feature selection algorithms. The feature ranking is based on several criteria: information gain, symmetrical uncertainty, chi-squared statistic etc. There are a number of features selection algorithms: Chi2 algorithm, based on chi-squared test,  fast correlation-based filter algorithm, feature weighting algorithm (RelieF), sequential forward search algorithm (CorrSF), Correlation-based feature selection algorithm (CFS). Package includes several classification algorithms with embedded feature selection and validation schemes. It includes also the functions for calculation of feature AUC (Area Under the ROC Curve) values with statistical significance analysis, calculation of Area Above the RCC (AAC) values. For two- and multi-class problems it is possible to use functions for HUM (hypervolume under manifold) calculation and construction 2D- and 3D- ROC curves. Relative Cost Curves (RCC) are provided to estimate the classifier performance under unequal misclassification costs.\cr
  Biocomb has a special function to deal with missing values, including different imputing schemes.
}
\details{
\tabular{ll}{
Package: \tab Biocomb\cr
Type: \tab Package\cr
Version: \tab 0.3\cr
Date: \tab 2016-08-14\cr
License: \tab GPL (>= 3)\cr
}
Biocomb package presents the functions for two stages of data mining process: feature selection and classification. One of the main functions of Biocomb is the \code{\link{select.process}} function. It presents the infrostructure to perform the feature ranking or feature selection for the data set with two or more class labels. Functions \code{\link{compute.aucs}}, \code{\link{select.inf.gain}}, \code{\link{select.inf.symm}} and \code{\link{select.inf.chi2}} calculate the different criterion measure for each feature in the dataset. Function \code{\link{select.fast.filter}} realizes the fast correlation-based filter method. Function \code{\link{chi2.algorithm}} performes Chi2 discretization algorithms with feature selection. Function \code{\link{select.forward.Corr}} is designed for the sequential forward features search according to the correlation measure. Function \code{\link{select.forward.wrapper}} is the realization of the wrapper feature selection method with sequential forward search strategy.
The auxiliary function  \code{\link{ProcessData}} performs the discretization of the numerical features and is called from the several functions for feature selection.
The second main function of the Biocomb is \code{\link{classifier.loop}} which presents the infrastructure for the classifier construction with the embedded feature selection and using the different schemes for the performance validation.
The functions \code{\link{compute.aucs}}, \code{\link{compute.auc.permutation}}, \code{\link{pauc}}, \code{\link{pauclog}}, \code{\link{compute.auc.random}} are the functions for calculation of feature AUC (Area Under  the ROC Curve) values with statistical significance analysis.
The functions \code{\link{plotRoc.curves}} is assigned for the construction of the ROC curve in 2D-space.
The functions \code{\link{cost.curve}} plots the RCC and calculates the corresponding AAC to estimate the classifier performance under unequal misclassification costs problem.
The function \code{\link{input_miss}} deals with missing value problem and realizes the two methods of missing value imputing.
The function \code{\link{generate.data.miss}} allows to generate the dataset with missing values from the input dataset in order to test the algorithms, which are designed to deal with missing values problem.
The functions \code{\link{CalculateHUM_seq}}, \code{\link{CalculateHUM_ROC}}, \code{\link{CalculateHUM_Plot}} are for HUM calculation and construction 2D- and 3D- ROC curves.
}
\section{Function}{
\tabular{ll}{
\code{\link{select.process}} \tab Perform the features ranking or features selection\cr
\code{\link{compute.aucs}} \tab Calculate the AUC values \cr
\code{\link{select.inf.gain}} \tab Calculate the Information Gain criterion \cr
\code{\link{select.inf.symm}} \tab Calculate the Symmetrical uncertainty criterion\cr
\code{\link{select.inf.chi2}} \tab Calculate the chi-squared statistic\cr
\code{\link{select.fast.filter}} \tab Select the feature subset with fast correlation-based filter method\cr
\code{\link{chi2.algorithm}} \tab Select the feature subset with Chi2 discretization algorithm.\cr
\code{\link{select.forward.Corr}} \tab Select the feature subset with forward search strategy and correlation measure\cr
\code{\link{select.forward.wrapper}} \tab Select the feature subset with a wrapper method\cr
\code{\link{ProcessData}} \tab Perform the discretization of the numerical features\cr
\code{\link{classifier.loop}} \tab Perform the classification with the embedded feature selection\cr
\code{\link{pauc}} \tab Calculate the p-values of the statistical significance of the two-class difference\cr
\code{\link{pauclog}} \tab Calculate the logarithm of p-values of the statistical significance\cr
\code{\link{compute.auc.permutation}} \tab Compute the p-value of the significance of the AUC using the permutation test\cr
\code{\link{compute.auc.random}} \tab Compute the p-value of the significance of the AUC using random sample generation\cr
\code{\link{plotRoc.curves}} \tab Plot the ROC curve in 2D-space \cr
\code{\link{CalculateHUM_seq}} \tab Calculate a maximal HUM value and the corresponding permutation of class labels\cr
\code{\link{CalculateHUM_Ex}} \tab Calculate the HUM values with exaustive serach for specified number of class labels \cr
\code{\link{CalculateHUM_ROC}} \tab Function to construct and plot the 2D- or 3d-ROC curve \cr
\code{\link{CalcGene}} \tab Compute the HUM value for one feature \cr
\code{\link{CalcROC}} \tab Compute the point coordinates to plot the 2D- or 3D-ROC curve \cr
\code{\link{CalculateHUM_Plot}} \tab Plot the 2D-ROC curve \cr
\code{\link{Calculate3D}} \tab Plot the 3D-ROC curve \cr
\code{\link{cost.curve}} \tab Plot the RCC and calculate the AAC for unequal misclassification costs\cr
\code{\link{input_miss}} \tab Perform the missing values imputation \cr
\code{\link{generate.data.miss}} \tab Generate the dataset with missing values \cr
  }
}

\section{Dataset}{
  This package comes with two simulated datasets and a real dataset of leukemia patients with 72 cases and 101 features. The last feature is the class (disease labels).
}

\section{Installing and using}{
  To install this package, make sure you are connected to the internet and issue the following command in the R prompt:

  \preformatted{
    install.packages("Biocomb")
  }

  To load the package in R:

  \preformatted{
    library(Biocomb)
  }
}

\author{
Natalia Novoselova, Junxi Wang,Frank Pessler,Frank Klawonn

Maintainer: Natalia Novoselova <novos65@mail.ru>
}
\references{
  H. Liu and L. Yu. "Toward Integrating Feature Selection Algorithms for Classification and Clustering", IEEE Trans. on Knowledge and Data Engineering, pdf, 17(4), 491-502, 2005.\cr
  L. Yu and H. Liu. "Feature Selection for High-Dimensional Data: A Fast Correlation-Based Filter Solution". In Proceedings of The Twentieth International Conference on Machine Leaning (ICML-03), Washington, D.C. pp. 856-863. August 21-24, 2003.\cr
  Y. Wang, I.V. Tetko, M.A. Hall, E. Frank, A. Facius, K.F.X. Mayer, and H.W. Mewes, "Gene Selection from Microarray Data for Cancer Classification?A Machine Learning Approach," Computational Biology and Chemistry, vol. 29, no. 1, pp. 37-46, 2005.\cr
  Olga Montvida and Frank Klawonn Relative cost curves: An alternative to AUC and an extension to 3-class problems,Kybernetika 50 no. 5, 647-660, 2014
}

\keyword{package}
\keyword{feature selection}
\keyword{classification}
\keyword{information gain}
\keyword{Relative cost curves}
\keyword{HUM}
\keyword{roc}
\keyword{auc}

\seealso{
CRAN packages \pkg{arules} or \pkg{discretization}
for feature discretization.
CRAN packages \pkg{pROC} for ROC curves.
CRAN packages \pkg{FSelector} for chi-squared test, forward search strategy.
CRAN packages \pkg{pamr} for nearest shrunken centroid classifier,
CRAN packages \pkg{MASS}, \pkg{e1071}, \pkg{randomForest},\pkg{class},
\pkg{nnet}, \pkg{rpart} are used in this package.
}
\examples{
data(data_test)
# class label must be factor
data_test[,ncol(data_test)]<-as.factor(data_test[,ncol(data_test)])

# Perform the feature selection using the fast correlation-based filter algorithm
disc="MDL"
threshold=0.2
attrs.nominal=numeric()
out=select.fast.filter(data_test,disc.method=disc,threshold=threshold,
attrs.nominal=attrs.nominal)

# Perform the classification with cross-validation of results
out=classifier.loop(data_test,classifiers=c("svm","lda","rf"),
 feature.selection="auc", flag.feature=FALSE,method.cross="fold-crossval")

# Calculate the coordinates for 2D- or 3D- ROC curve and the optimal threshold point
\dontrun{data(data_test)
xllim<--4
xulim<-4
yllim<-30
yulim<-110

attrs.no=1
pos.Class<-levels(data_test[,ncol(data_test)])[1]
add.legend<-TRUE

aacs<-rep(0,length(attrs.no))
color<-c(1:length(attrs.no))

out <- cost.curve(data_test,attrs.no, pos.Class,col=color[1],add=F,
 xlim=c(xllim,xulim),ylim=c(yllim,yulim))
}
}
