\name{fulltle-methods}
\docType{methods}
\alias{fulltle-methods}
\alias{fulltle,IData-method}
\alias{fulltle}   % Check if this sintax is really correct !!!
\title{ Methods for Function fulltle in Package \sQuote{MAINT.Data}}

\description{Performs maximum trimmed likelihood estimation by an exact algorithm (full enumeratiom of all k-trimmed subsets)}

\usage{fulltle(Idt, CovCase = 1:4, SelCrit = c("BIC", "AIC"), alpha =
                 0.75, use.correction = TRUE, getalpha = "TwoStep",
                 rawMD2Dist = c("ChiSq", "HardRockeAsF",
                 "HardRockeAdjF"), MD2Dist = c("ChiSq",
                 "CerioliBetaF"), eta = 0.025, multiCmpCor = c("never",
                 "always", "iterstep"), outlin = c("MidPandLogR",
                 "MidP", "LogR"), reweighted = TRUE, k2max=1e6, 
                 force = FALSE, \dots)
}

\arguments{
    \item{Idt}{An IData object representing interval-valued entities.}

    \item{CovCase}{Configuration of the variance-covariance matrix: a set of integers between 1 and 4.}

    \item{SelCrit}{The model selection criterion.}

    \item{alpha}{Numeric parameter controlling the size of the subsets over which the trimmed likelihood is maximized; roughly alpha*nrow(Idt) observations are used for computing the trimmed likelihood. Note that when argument \sQuote{getalpha} is set to \dQuote{TwoStep} the final value of \sQuote{alpha} is estimated by a two-step procedure and the value of argument \sQuote{alpha} is only used to specify the size of the samples used in the first step. Allowed values are between 0.5 and 1.}

    \item{use.correction}{whether to use finite sample correction factors; defaults to \code{TRUE}.}

    \item{getalpha}{Argument specifying if the \sQuote{alpha} parameter (roughly the percentage of the sample used for computing the trimmed likelihood) should be estimated from the data, or if the value of the argument \sQuote{alpha} should be used instead. When set to \dQuote{TwoStep}, \sQuote{alpha} is estimated by a two-step procedure with the value of argument \sQuote{alpha} specifying the size of the samples used in the first step. Otherwise, with the value of argument \sQuote{alpha} is used directly.}   

    \item{rawMD2Dist}{The assumed reference distribution of the raw MCD squared distances, which is used to find to cutoffs defining the observations kept in one-step reweighted MCD estimates. Alternatives are \sQuote{ChiSq}, \sQuote{HardRockeAsF} and \sQuote{HardRockeAdjF}, respectivelly for the usual Chi-squared, and the asymptotic and adjusted scaled F distributions proposed by Hardin and Rocke (2005).}    

    \item{MD2Dist}{The assumed reference distributions used to find cutoffs defining the observations assumed as outliers. Alternatives are \dQuote{ChiSq} and \dQuote{CerioliBetaF} respectivelly for the usual Chi-squared, and the Beta and F distributions proposed by Cerioli (2010).}    

    \item{eta}{Nominal size of the null hypothesis that a given observation is not an outlier. Defines the raw MCD Mahalanobis distances cutoff used to choose the observations kept in the reweightening step.}

    \item{multiCmpCor}{Whether a multicomparison correction of the nominal size (eta) for the outliers tests should be performed. Alternatives are: \sQuote{never} -- ignoring the multicomparisons and testing all entities at the \sQuote{eta} nominal level. \sQuote{always} -- testing all n entitites at 1.- (1.-\sQuote{eta}^(1/n)); and \sQuote{iterstep} -- use the iterated rule proposed by  Cerioli (2010), \emph{i.e.}, make an initial set of tests using the nominal size 1.- (1-\sQuote{eta}^(1/n)), and if no outliers are detected stop. Otherwise, make a second step testing for outliers at the \sQuote{eta} nominal level.}   

    \item{outlin}{The type of outliers to be consideres. \dQuote{MidPandLogR} if outliers may be present in both MidPpoints and LogRanges, \dQuote{MidP} if outliers are only present in MidPpoints, or \dQuote{LogR} if outliers are only present in LogRanges.}

    \item{reweighted}{should a (Re)weighted estimate of the covariance matrix be used in the computation of the trimmed likelihood or just a \dQuote{raw} covariance estimate; default is (Re)weighting.}

  \item{k2max}{Maximal allowed l2-norm condition number for correlation matrices. Correlation matrices with condition number above k2max are considered to be numerically singular, leading to degenerate results.}

   \item{force}{A boolean flag indicating whether, for moderate or large data sets the algorithm should proceed anyway, regardless of an expected long excution time, due to exponential explosions in the number of different subsets that need to be avaluated by fulltle}

    \item{\dots}{Further arguments to be passed to internal functions of \sQuote{fulltle}.} 
}

\value{An object of class \code{\linkS4class{IdtE}} with the fulltle estimates, the value of the comparison criterion used to select the covariance configurations and the robust squared Mahalanobis distances.} 

\references{
Brito, P., Duarte Silva, A. P. (2012), Modelling Interval Data with Normal and Skew-Normal Distributions. \emph{Journal of Applied Statistics} \bold{39}(1), 3--20.\cr

Cerioli, A. (2010), Multivariate Outlier Detection with High-Breakdown Estimators.
\emph{Journal of the American Statistical Association} \bold{105} (489), 147--156.\cr 

Duarte Silva, A.P., Filzmoser, P. and Brito, P. (2017), Outlier detection in interval data. \emph{Advances in Data Analysis and Classification}, 1--38.\cr

Hadi, A. S. and Luceno, A. (1997), Maximum trimmed likelihood estimators: a unified approach, examples, and algorithms.
\emph{Computational Statistics and Data Analysis} \bold{25}(3), 251--272.\cr

Hardin, J. and Rocke, A. (2005), The Distribution of Robust Distances.
\emph{Journal of Computational and Graphical Statistics} \bold{14}, 910--927.\cr
}

\keyword{methods}
\keyword{trimmed maximum likelihood estimator}
\keyword{interval data}

\seealso{ \code{\link{fasttle}}, \code{\link{getIdtOutl}} }


\examples{


\dontrun{

# Create an Interval-Data object containing the intervals for characteristics 
# of 27 cars models.

CarsIdt <- IData(Cars[1:8],VarNames=c("Price","EngineCapacity","TopSpeed","Acceleration"))

# Estimate parameters by the full trimmed maximum likelihood estimator, 
# using a two-step procedure to select the trimming parameter, a reweighed 
# MCD estimate, and the classical 97.5\% chi-squared quantile cut-offs.

CarsTE1 <- fulltle(CarsIdt)
cat("Cars data -- normal maximum trimmed likelihood estimation results:\n")
print(CarsTE1)
		
# Estimate parameters by the full trimmed maximum likelihood estimator, using
# the triming parameter that maximizes breakdown, and a reweighed MCD estimate
# based on the 97.5\% quantiles of Hardin and Rocke adjusted F distributions.
		
CarsTE2 <- fulltle(CarsIdt,alpha=0.5,getalpha=FALSE,rawMD2Dist="HardRockeAdjF")
cat("Cars data -- normal maximum trimmed likelihood estimation results:\n")
print(CarsTE2)
		
# Estimate parameters by the full trimmed maximum likelihood estimator, using 
# a two-step procedure to select the trimming parameter, and a reweighed MCD estimate 
# based on Hardin and Rocke adjusted F distributions, 95\% quantiles, and 
# the Cerioli Beta and F distributions together with his iterated procedure 
# to identify outliers in the first step.
		
CarsTE3 <- fulltle(CarsIdt,rawMD2Dist="HardRockeAdjF",eta=0.05,MD2Dist="CerioliBetaF",
multiCmpCor="iterstep")
cat("Cars data -- normal maximum trimmed likelihood estimation results:\n")
print(CarsTE3)

}

}
