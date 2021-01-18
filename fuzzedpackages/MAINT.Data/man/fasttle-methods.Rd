\name{fasttle-methods}
\docType{methods}
\alias{fasttle-methods}
\alias{fasttle,IData-method}
\alias{fasttle}
\title{ Methods for Function fasttle in Package \sQuote{MAINT.Data}}

\description{Performs maximum trimmed likelihood estimation by the fasttle algorithm}

\usage{fasttle(Idt,
    CovCase=1:4,
    SelCrit=c("BIC","AIC"),
    alpha=control@alpha,
    nsamp = control@nsamp,
    seed=control@seed,
    trace=control@trace,
    use.correction=control@use.correction,
    ncsteps=control@ncsteps,
    getalpha=control@getalpha,
    rawMD2Dist=control@rawMD2Dist,				
    MD2Dist=control@MD2Dist,
    eta=control@eta,
    multiCmpCor=control@multiCmpCor,				
    getkdblstar=control@getkdblstar,
    outlin=control@outlin,
    trialmethod=control@trialmethod,
    m=control@m,
    reweighted = control@reweighted,
    k2max = control@k2max, 
    otpType=control@otpType,
    control=RobEstControl(), \dots)
}

\arguments{
    \item{Idt}{An IData object representing interval-valued entities.}

    \item{CovCase}{Configuration of the variance-covariance matrix: a set of integers between 1 and 4.}

    \item{SelCrit}{The model selection criterion.}

    \item{alpha}{Numeric parameter controlling the size of the subsets over which the trimmed likelihood is maximized; roughly alpha*nrow(Idt) observations are used for computing the trimmed likelihood.  Note that when argument \sQuote{getalpha} is set to \dQuote{TwoStep} the final value of \sQuote{alpha} is estimated by a two-step procedure and the value of argument \sQuote{alpha} is only used to specify the size of the samples used in the first step. Allowed values are between 0.5 and 1.}

    \item{nsamp}{Number of subsets used for initial estimates.}  

    \item{seed}{Initial seed for random generator, like \code{\link{.Random.seed}}, see \code{\link{rrcov.control}}.}

    \item{trace}{Logical (or integer) indicating if intermediate results should be printed; defaults to \code{FALSE}.}

    \item{use.correction}{ whether to use finite sample correction factors; defaults to \code{TRUE}.}

    \item{ncsteps}{The maximum number of concentration steps used each iteration of the fasttle algorithm.}

    \item{getalpha}{Argument specifying if the \sQuote{alpha} parameter (roughly the percentage of the sample used for computing the trimmed likelihood) should be estimated from the data, or if the value of the argument \sQuote{alpha} should be used instead. When set to \dQuote{TwoStep}, \sQuote{alpha} is estimated by a two-step procedure with the value of argument \sQuote{alpha} specifying the size of the samples used in the first step. Otherwise, with the value of argument \sQuote{alpha} is used directly.}   

    \item{rawMD2Dist}{The assumed reference distribution of the raw MCD squared distances, which is used to find to cutoffs defining the observations kept in one-step reweighted MCD estimates. Alternatives are \sQuote{ChiSq},\sQuote{HardRockeAsF} and \sQuote{HardRockeAdjF}, respectivelly for the usual Chi-squared, and the asymptotic and adjusted scaled F distributions proposed by Hardin and Rocke (2005).}

    \item{MD2Dist}{The assumed reference distributions used to find cutoffs defining the observations assumed as outliers. Alternatives are \dQuote{ChiSq} and \dQuote{CerioliBetaF} respectivelly for the usual Chi-squared, or the Beta and F distributions proposed by Cerioli (2010).}    

    \item{eta}{Nominal size for the null hypothesis that a given observation is not an outlier. Defines the raw MCD Mahalanobis distances cutoff used to choose the observations kept in the reweightening step.}

    \item{multiCmpCor}{Whether a multicomparison correction of the nominal size (eta) for the outliers tests should be performed. Alternatives are: \sQuote{never} -- ignoring the multicomparisons and testing all entities at \sQuote{eta} nominal level. \sQuote{always} -- testing all n entitites at 1.- (1.-\sQuote{eta}^(1/n)); and \sQuote{iterstep} -- use the iterated rule proposed by  Cerioli (2010), \emph{i.e.}, make an initial set of tests using the nominal size 1.- (1-\sQuote{eta}^(1/n)), and if no outliers are detected stop. Otherwise, make a second step testing for outliers at the \sQuote{eta} nominal level. }   

    \item{getkdblstar}{Argument specifying the size of the initial small (in order to minimize the probability of outliers) subsets. If set to the string \dQuote{Twopplusone} (default) the initial sets have twice the number of interval-value variables plus one (i.e., they are the smaller samples that lead to a non-singular covariance estimate). Otherwise, an integer with the size of the initial sets.}

    \item{outlin}{The type of outliers to be considered. \dQuote{MidPandLogR} if outliers may be present in both MidPpoints and LogRanges, \dQuote{MidP} if outliers are only present in MidPpoints, or \dQuote{LogR} if outliers are only present in LogRanges.}

    \item{trialmethod}{The method to find a trial subset used to initialize each replication of the fasttle algorithm. The current options are \dQuote{simple} (default) that simply selects \sQuote{kdblstar} observations at random, and \dQuote{Poolm} that divides the original sample into \sQuote{m} non-overlaping subsets, applies the \sQuote{simple trial} and the refinement methods to each one of them, and merges the results into a trial subset.}

    \item{m}{Number of non-overlaping subsets used by the trial method when the argument of \sQuote{trialmethod} is set to 'Poolm'.}

    \item{reweighted}{Should a (Re)weighted estimate of the covariance matrix be used in the computation of the trimmed likelihood or just a \dQuote{raw} covariance estimate; default is (Re)weighting.}

    \item{k2max}{Maximal allowed l2-norm condition number for correlation matrices. Correlation matrices with condition number above k2max are considered to be numerically singular, leading to degenerate results.}

    \item{otpType}{The amount of output returned by fasttle. Current options are \dQuote{SetMD2andEst} (default) which returns an \sQuote{IdtSngNDRE} object with the fasttle estimates, \cr
a vector with the final trimmed subset elements used to compute these estimates and the corresponding robust squared Mahalanobis distances, and \cr
\dQuote{SetMD2EstandPrfSt} wich returns an \sQuote{IdtSngNDRE} object with the previous slots plust a list of some performance statistics concerning the algorithm execution.}

  \item{control}{a list with estimation options - this includes those above provided in the function specification. See
\code{\link{RobEstControl}} for the defaults.  If \code{control} is supplied, the parameters from it will be used.
If parameters are passed also in the invocation statement, they will override the corresponding elements of the control object.}

  \item{\dots}{Further arguments to be passed to internal functions of \code{fasttle}.} 
}

\value{An object of class \code{\linkS4class{IdtE}} with the fasttle estimates, the value of the comparison criterion used to select the covariance configurations, the robust squared Mahalanobis distances, and optionally (if argument \sQuote{otpType} is set to true) performance statistics concerning the algorithm execution.} 

\references{
Brito, P., Duarte Silva, A. P. (2012), Modelling Interval Data with Normal and Skew-Normal Distributions. \emph{Journal of Applied Statistics} \bold{39}(1), 3--20.\cr

Cerioli, A. (2010), Multivariate Outlier Detection with High-Breakdown Estimators.
\emph{Journal of the American Statistical Association} \bold{105} (489), 147--156.\cr 

Duarte Silva, A.P., Filzmoser, P. and Brito, P. (2017), Outlier detection in interval data. \emph{Advances in Data Analysis and Classification}, 1--38.\cr

Hadi, A. S. and Luceno, A. (1997), Maximum trimmed likelihood estimators: a unified approach, examples, and algorithms.
\emph{Computational Statistics and Data Analysis} \bold{25}(3), 251--272.\cr

Hardin, J. and Rocke, A. (2005), The Distribution of Robust Distances.
\emph{Journal of Computational and Graphical Statistics} \bold{14}, 910--927.\cr

Todorov V. and Filzmoser P. (2009), An Object Oriented Framework for Robust Multivariate Analysis. \emph{Journal of Statistical Software} \bold{32}(3), 1--47.
}

\keyword{methods}
\keyword{trimmed maximum likelihood estimator}
\keyword{interval data}


\seealso{ \code{\link{fulltle}}, \code{\linkS4class{RobEstControl}}, \code{\link{getIdtOutl}}, \code{\linkS4class{IdtSngNDRE}}  }


\examples{

\dontrun{

# Create an Interval-Data object containing the intervals of temperatures by quarter 
# for 899 Chinese meteorological stations.

ChinaT <- IData(ChinaTemp[1:8])

# Estimate parameters by the fast trimmed maximum likelihood estimator, 
# using a two-step procedure to select the trimming parameter, a reweighted 
# MCD estimate, and the classical 97.5\% chi-squared quantile cut-offs.

Chinafasttle1 <- fasttle(ChinaT)
cat("China maximum trimmed likelihood estimation results =\n")
print(Chinafasttle1)

# Estimate parameters by the fast trimmed maximum likelihood estimator, using 
# the triming parameter that maximizes breakdown, and a reweighted MCD estimate 
# based on the 97.5\% quantiles of Hardin and Rocke adjusted F distributions.

Chinafasttle2 <- fasttle(ChinaT,alpha=0.5,getalpha=FALSE,rawMD2Dist="HardRockeAdjF")
cat("China maximum trimmed likelihood estimation results =\n")
print(Chinafasttle2)

# Estimate parameters by the fast trimmed maximum likelihood estimator, using a two-step procedure
# to select the triming parameter, a reweighed MCD estimate based on Hardin and Rocke adjusted 
# F distributions, and 95\% quantiles, and the Cerioli Beta and F distributions together
# with Cerioli iterated procedure to identify outliers in the first step.

Chinafasttle3 <- fasttle(ChinaT,rawMD2Dist="HardRockeAdjF",eta=0.05,MD2Dist="CerioliBetaF",
multiCmpCor="iterstep")
cat("China maximum trimmed likelihood estimation results =\n")
print(Chinafasttle3)

}

}
