\name{snda-methods}
\docType{methods}
\alias{snda}
\alias{snda-methods}
\alias{snda,IData-method}
\alias{snda,IdtLocSNMANOVA-method}
\alias{snda,IdtLocNSNMANOVA-method}
\alias{snda,IdtGenSNMANOVA-method}
\alias{snda,IdtGenNSNMANOVA-method}

\title{Skew-Normal Discriminant Analysis of Interval Data}

\description{snda performs discriminant analysis of Interval Data based on estimates of mixtures of Skew-Normal models}

\usage{

\S4method{snda}{IData}(x, grouping, prior="proportions", CVtol=1.0e-5, subset=1:nrow(x),
  CovCase=1:4, SelCrit=c("BIC","AIC"), Mxt=c("Loc","Gen"), k2max=1e6, \dots )

\S4method{snda}{IdtLocSNMANOVA}( x, prior="proportions", selmodel=BestModel(H1res(x)),
  egvtol=1.0e-10, silent=FALSE, k2max=1e6, \dots )

\S4method{snda}{IdtLocNSNMANOVA}( x, prior="proportions",
  selmodel=BestModel(H1res(x)@SNMod), egvtol=1.0e-10, silent=FALSE, k2max=1e6, \dots )

\S4method{snda}{IdtGenSNMANOVA}( x, prior="proportions", selmodel=BestModel(H1res(x)),
  silent=FALSE, k2max=1e6, \dots )

\S4method{snda}{IdtGenNSNMANOVA}( x, prior="proportions",
  selmodel=BestModel(H1res(x)@SNMod), silent=FALSE, k2max=1e6, \dots )

}

\arguments{

  \item{x}{An object of class \code{\linkS4class{IData}}, \code{\linkS4class{IdtLocSNMANOVA}}, \code{\linkS4class{IdtLocNSNMANOVA}},\code{\linkS4class{IdtGenSNMANOVA}} or \code{\linkS4class{IdtGenNSNMANOVA}} with either the original Interval Data, or the results of a Interval Data Skew-Normal MANOVA, from which the discriminant analysis will be based.}

  \item{grouping}{Factor specifying the class for each observation.}

  \item{prior}{The prior probabilities of class membership.  If unspecified, the class proportions for the training set are used.  If
present, the probabilities should be specified in the order of the factor levels.}

  \item{CVtol}{Tolerance level for absolute value of the coefficient of variation of non-constant variables. When a MidPoint or LogRange has an absolute value within-groups coefficient of variation below CVtol, it is considered to be a constant.}

  \item{subset}{An index vector specifying the cases to be used in the analysis.}

  \item{CovCase}{Configuration of the variance-covariance matrix: a set of integers between 1 and 4.}

  \item{SelCrit}{The model selection criterion.}

 \item{Mxt}{Indicates the type of mixing distributions to be considered. Current alternatives are \dQuote{Loc} (location model -- groups differ only on the location parameters of a Skew-Normal model) and \dQuote{Gen} (general model -- groups differ on all parameters of a Skew-Normal models).}

  \item{silent}{A boolean flag indicating whether a warning message should be printed if the method fails.} 
 
  \item{selmodel}{Selected model from a list of candidate models saved in object x.}

  \item{egvtol}{Tolerance level for the eigenvalues of the product of the inverse within by the between covariance matrices. When a  eigenvalue has an absolute value below egvtol, it is considered to be zero.}

  \item{k2max}{Maximal allowed l2-norm condition number for correlation matrices. Correlation matrices with condition number above k2max are considered to be numerically singular, leading to degenerate results.}

  \item{\dots}{Other named arguments.}
}

\references{
Azzalini, A. and Dalla Valle, A. (1996), The multivariate skew-normal distribution. \emph{Biometrika} \bold{83}(4), 715--726.\cr

Brito, P., Duarte Silva, A. P. (2012), Modelling Interval Data with Normal and Skew-Normal Distributions. \emph{Journal of Applied Statistics} \bold{39}(1), 3--20.

Duarte Silva, A.P. and Brito, P. (2015), Discriminant analysis of interval data: An assessment of parametric and distance-based   approaches. \emph{Journal of Classification} \bold{39}(3), 516--541.\cr
}

\seealso{\code{\link{lda}}, \code{\link{qda}}, \code{\link{Roblda}}, \code{\link{Robqda}}, \code{\linkS4class{IData}}, \code{\linkS4class{IdtLocSNMANOVA}}, \code{\linkS4class{IdtLocNSNMANOVA}}, \code{\linkS4class{IdtGenSNMANOVA}},\code{\linkS4class{IdtGenSNMANOVA}} }

\keyword{methods}
\keyword{Interval Data}
\keyword{Discriminant Analysis Interval Data}

\examples{

\dontrun{

# Create an Interval-Data object containing the intervals for 899 observations 
# on the temperatures by quarter in 60 Chinese meteorological stations.

ChinaT <- IData(ChinaTemp[1:8],VarNames=c("T1","T2","T3","T4"))

# Skew-Normal based discriminant analysis, asssuming that the different regions differ
# only in location parameters

ChinaT.locsnda <- snda(ChinaT,ChinaTemp$GeoReg,Mxt="Loc")

cat("Temperatures of China -- SkewNormal location model discriminant analysis results:\n")
print(ChinaT.locsnda)

#Estimate error rates by three-fold cross-validation without replication  

CVlocsnda <- DACrossVal(ChinaT,ChinaTemp$GeoReg,TrainAlg=snda,Mxt="Loc",
  CovCase=CovCase(ChinaT.locsnda),kfold=3,CVrep=1)

summary(CVlocsnda[,,"Clerr"])

# Skew-Normal based discriminant analysis, asssuming that the different regions may differ
# in all SkewNormal parameters

ChinaT.gensnda <- snda(ChinaT,ChinaTemp$GeoReg,Mxt="Gen")

cat("Temperatures of China -- SkewNormal general model discriminant analysis results:\n")
print(ChinaT.gensnda)

#Estimate error rates by three-fold cross-validation without replication  

CVgensnda <- DACrossVal(ChinaT,ChinaTemp$GeoReg,TrainAlg=snda,Mxt="Gen",
  CovCase=CovCase(ChinaT.gensnda),kfold=3,CVrep=1)

summary(CVgensnda[,,"Clerr"])

}

}
