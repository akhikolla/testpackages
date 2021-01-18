\name{MAINT.Data-package}
\alias{MAINT.Data-package}
\alias{MAINT.Data}
\docType{package}
\title{Modelling and Analizing Interval Data}

\description{MAINT.Data implements methodologies for modelling Interval Data by Normal and Skew-Normal distributions,
considering four different possible configurations structures for the variance-covariance matrix.
It introduces a data class for representing interval data and includes functions and methods for parametric modelling and analysing of interval data.
It performs maximum likelihood and trimmed maximum likelihood estimation, statistical tests, as well as (M)ANOVA, Discriminant Analysis and Gaussian Model Based Clustering.
}

\details{In the classical model of multivariate data analysis, data is represented in a data-array where n ``individuals" (usually in rows) take exactly one value for each variable (usually in columns).
Symbolic Data Analysis (see, e.g., Noirhomme-Fraiture and Brito (2011)) provides a framework where new variable types allow to take directly into account variability and/or uncertainty associated to each single ``individual",
by allowing multiple, possibly weighted, values for each variable. 
New variable types - interval, categorical multi-valued and modal variables - have been introduced.\cr
We focus on the analysis of interval data, i.e., where elements are described by variables whose values are intervals.
Parametric inference methodologies based on probabilistic models for interval variables are developed in Brito and Duarte Silva (2011) where each interval is represented by its midpoint and log-range,for which Normal and Skew-Normal (Azzalini and Dalla Valle (1996)) distributions are assumed.
The intrinsic nature of the interval variables leads to special structures of the variance-covariance matrix, which are represented by four different possible configurations.\cr
MAINT.Data implements the proposed methodologies in R, introducing a data class for representing interval data; it
includes functions for modelling and analysing interval data, in particular maximum likelihood and trimmed maximum likelihood (Duarte Silva, Filzmoser and Brito (2017)) estimation, and statistical tests for the different considered configurations. 
Methods for (M)ANOVA, Discriminant Analysis (Duarte Silva and Brito (2015)) and model based clustering (Brito, Duarte Silva and Dias (2015))  of this data class are also provided.

\tabular{ll}{
Package: \tab MAINT.Data\cr
Type: \tab Package\cr
Version: \tab 2.5.0\cr
Date: \tab 2020-09-24\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
LazyData: \tab yes\cr
}
}

\author{Pedro Duarte Silva <psilva@porto.ucp.pt>\cr 
        Paula Brito <mpbrito.fep.up.pt>

Maintainer: Pedro Duarte Silva <psilva@porto.ucp.pt>\cr 
}

\references{
Azzalini, A. and Dalla Valle, A. (1996), The multivariate skew-normal distribution. \emph{Biometrika} \bold{83}(4), 715--726.\cr

Brito, P. and Duarte Silva, A. P. (2012), Modelling Interval Data with Normal and Skew-Normal Distributions. \emph{Journal of Applied Statistics} \bold{39}(1), 3--20.\cr

Brito, P., Duarte Silva, A. P. and Dias, J. G. (2015), Probabilistic Clustering of Interval Data. \emph{Intelligent Data Analysis} \bold{19}(2), 293--313.\cr

Duarte Silva, A.P. and Brito, P. (2015), Discriminant analysis of interval data: An assessment of parametric and distance-based   approaches. \emph{Journal of Classification} \bold{39}(3), 516--541.\cr

Duarte Silva, A.P., Filzmoser, P. and Brito, P. (2017), Outlier detection in interval data. \emph{Advances in Data Analysis and Classification}, 1--38.\cr

Noirhomme-Fraiture, M. and Brito, P. (2011), Far Beyond the Classical Data Models: Symbolic Data Analysis. \emph{Statistical Analysis and Data Mining} \bold{4}(2), 157--170.
}

\keyword{MAINT.Data}
\keyword{Interval Data}
\keyword{Symbolic Data Analysis}
\keyword{MANOVA for Interval Data}
\keyword{Parametric modelling of Interval Data}
\keyword{Statistical tests for Interval Data}

\examples{
# Create an Interval-Data object containing the intervals for 899 observations 
# on the temperatures by quarter in 60 Chinese meteorological stations.

ChinaT <- IData(ChinaTemp[1:8],VarNames=c("T1","T2","T3","T4"))

#Display the first and last observations

head(ChinaT)
tail(ChinaT)

#Print summary statistics

summary(ChinaT)

#Create a new data set considering only the Winter (1st and 4th) quarter intervals

ChinaWT <- ChinaT[,c(1,4)]

# Estimate normal distribution parameters by maximum likelihood, assuming 
# the classical (unrestricted) covariance configuration Case 1

ChinaWTE.C1 <- mle(ChinaWT,CovCase=1)
cat("Winter temperatures of China -- normal maximum likelhiood estimation results:\n")
print(ChinaWTE.C1)
cat("Standard Errors of Estimators:\n") ; print(stdEr(ChinaWTE.C1))

# Estimate normal distribution parameters by maximum likelihood, 
# assuming that one of the C2, C3 or C4 restricted covariance configuration cases hold

ChinaWTE.C234 <- mle(ChinaWT,CovCase=2:4)
cat("Winter temperatures of China -- normal maximum likelihood estimation results:\n")
print(ChinaWTE.C234)
cat("Standard Errors of Estimators:\n") ; print(stdEr(ChinaWTE.C234))

# Estimate normal distribution  parameters robustly by fast maximun trimmed likelihood, 
# assuming that one of the C2, C3 or C4 restricted covariance configuration cases hold

\dontrun{
ChinaWTE.C234 <- fasttle(ChinaWT,CovCase=2:4)
cat("Winter temperatures of China -- normal maximum trimmed likelhiood estimation results:\n")
print(ChinaWTE.C234)

# Estimate skew-normal distribution  parameters 

ChinaWTE.SkN <- mle(ChinaWT,Model="SKNormal")
cat("Winter temperatures of China -- Skew-Normal maximum likelhiood estimation results:\n")
print(ChinaWTE.SkN)
cat("Standard Errors of Estimators:\n") ; print(stdEr(ChinaWTE.SkN))
}

#MANOVA tests assuming that configuration case 1 (unrestricted covariance) 
# or 3 (MidPoints independent of Log-Ranges) holds.  

ManvChinaWT.C13 <- MANOVA(ChinaWT,ChinaTemp$GeoReg,CovCase=c(1,3))
cat("Winter temperatures of China -- MANOVA by geografical regions results:\n")
print(ManvChinaWT.C13)

#Linear Discriminant Analysis

ChinaWT.lda <- lda(ManvChinaWT.C13)
cat("Winter temperatures of China -- linear discriminant analysis results:\n")
print(ChinaWT.lda)
cat("lda Prediction results:\n")
print(predict(ChinaWT.lda,ChinaWT)$class)

\dontrun{
#Estimate error rates by ten-fold cross-validation 

CVlda <- DACrossVal(ChinaWT,ChinaTemp$GeoReg,TrainAlg=lda,
CovCase=BestModel(H1res(ManvChinaWT.C13)),CVrep=1)

#Robust Quadratic Discriminant Analysis

ChinaWT.rqda <- Robqda(ChinaWT,ChinaTemp$GeoReg)
cat("Winter temperatures of China -- robust quadratic discriminant analysis results:\n")
print(ChinaWT.rqda)
cat("robust qda prediction results:\n")
print(predict(ChinaWT.rqda,ChinaWT)$class)
}

# Create an Interval-Data object containing the intervals of loan data
# (from the Kaggle Data Science platform) aggregated by loan purpose

LbyPIdt <- IData(LoansbyPurpose_minmaxDt,
  VarNames=c("ln-inc","ln-revolbal","open-acc","total-acc")) 

print(LbyPIdt)

\dontrun{

#Fit homoscedastic Gaussian mixtures with up to six components

mclustres <- Idtmclust(LbyPIdt,G=1:6)
plotInfCrt(mclustres,legpos="bottomright")
print(mclustres)

#Display the results of the best mixture according to the BIC

summary(mclustres,parameters=TRUE,classification=TRUE)
pcoordplot(mclustres)

}


}
