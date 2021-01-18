\name{lda-methods}
\docType{methods}
\alias{lda}
\alias{lda-methods}
\alias{lda,IData-method}
\alias{lda,IdtMxtNDE-method}
\alias{lda,IdtClMANOVA-method}
\alias{lda,IdtLocNSNMANOVA-method}

\title{Linear Discriminant Analysis of Interval Data}

\description{lda performs linear discriminant analysis of Interval Data based on classic estimates of a mixture of  Gaussian models.}

\usage{
\S4method{lda}{IData}(x, grouping, prior="proportions", CVtol=1.0e-5, egvtol=1.0e-10, 
  subset=1:nrow(x), CovCase=1:4, SelCrit=c("BIC","AIC"), silent=FALSE, k2max=1e6, \dots )

\S4method{lda}{IdtMxtNDE}(x, prior="proportions", selmodel=BestModel(x), egvtol=1.0e-10,
  silent=FALSE, k2max=1e6, \dots )

\S4method{lda}{IdtClMANOVA}( x, prior="proportions", selmodel=BestModel(H1res(x)),
  egvtol=1.0e-10, silent=FALSE, k2max=1e6, \dots )

\S4method{lda}{IdtLocNSNMANOVA}( x, prior="proportions", 
  selmodel=BestModel(H1res(x)@NMod), egvtol=1.0e-10, silent=FALSE, k2max=1e6, \dots )

}

\arguments{

  \item{x}{An object of class \code{\linkS4class{IData}}, \code{\linkS4class{IdtMxtNDE}}, \code{\linkS4class{IdtClMANOVA}} or \code{\linkS4class{IdtLocNSNMANOVA}} with either the original Interval Data, an estimate of a mixture of gaussian models for Interval Data, or the results of an Interval Data MANOVA, from which the discriminant analysis will be based.}

  \item{grouping}{Factor specifying the class for each observation.}

  \item{prior}{The prior probabilities of class membership.  If unspecified, the class proportions for the training set are used.  If
present, the probabilities should be specified in the order of the factor levels.}

  \item{CVtol}{Tolerance level for absolute value of the coefficient of variation of non-constant variables. When a MidPoint or LogRange has an absolute value within-groups coefficient of variation below CVtol, it is considered to be a constant.}

  \item{egvtol}{Tolerance level for the eigenvalues of the product of the inverse within by the between covariance matrices. When a  eigenvalue has an absolute value below egvtol, it is considered to be zero.}

  \item{subset}{An index vector specifying the cases to be used in the analysis.}

  \item{CovCase}{Configuration of the variance-covariance matrix: a set of integers between 1 and 4.}

  \item{SelCrit}{The model selection criterion.}

  \item{silent}{A boolean flag indicating whether a warning message should be printed if the method fails.}
 
  \item{selmodel}{Selected model from a list of candidate models saved in object x.}

  \item{k2max}{Maximal allowed l2-norm condition number for correlation matrices. Correlation matrices with condition number above k2max are considered to be numerically singular, leading to degenerate results.}

  \item{\dots}{Other named arguments.}
}

\references{
Brito, P., Duarte Silva, A. P. (2012), Modelling Interval Data with Normal and Skew-Normal Distributions. \emph{Journal of Applied Statistics} \bold{39}(1), 3--20.\cr

Duarte Silva, A.P. and Brito, P. (2015), Discriminant analysis of interval data: An assessment of parametric and distance-based   approaches. \emph{Journal of Classification} \bold{39}(3), 516--541.
}

\seealso{ \code{\link{qda}}, \code{\link{snda}}, \code{\link{Roblda}}, \code{\link{Robqda}}, \code{\linkS4class{IData}}, \code{\linkS4class{IdtMxtNDE}}, \code{\linkS4class{IdtClMANOVA}}, 
 \code{\linkS4class{IdtLocNSNMANOVA}} }

\keyword{methods}
\keyword{Interval Data}
\keyword{Discriminant Analysis Interval Data}

\examples{

# Create an Interval-Data object containing the intervals for 899 observations 
# on the temperatures by quarter in 60 Chinese meteorological stations.

ChinaT <- IData(ChinaTemp[1:8],VarNames=c("T1","T2","T3","T4"))

#Linear Discriminant Analysis

ChinaT.lda <- lda(ChinaT,ChinaTemp$GeoReg)
cat("Temperatures of China -- linear discriminant analysis results:\n")
print(ChinaT.lda)
cat("lda Prediction results:\n")
print(predict(ChinaT.lda,ChinaT)$class)

\dontrun{
#Estimate error rates by ten-fold cross-validation replicated 20 times  

CVlda <- DACrossVal(ChinaT,ChinaTemp$GeoReg,TrainAlg=lda,CovCase=CovCase(ChinaT.lda))
summary(CVlda[,,"Clerr"])
}

}
