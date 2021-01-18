\name{RobEstControl-class}
\docType{class}
\alias{RobEstControl-class}

\title{Class 'RobEstControl' - contains control parameters for the robust estimation of parametric interval data models.}
\description{
    This class extends the \code{CovControlMcd} class
    and contains control parameters for the robust estimation of parametric interval data models.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("RobEstControl", ...)}
or by calling the constructor-function \code{RobEstControl}.
}
\section{Slots}{
  \describe{
    \item{\code{alpha}:}{Inherited from class \code{"CovControlMcd"}. Numeric parameter controlling the size of the subsets over which the trimmed likelihood is maximized; roughly alpha*nrow(Idt) observations are used for computing the trimmed likelihood. Allowed values are between 0.5 and 1.
Note that when argument \sQuote{getalpha} is set to \dQuote{TwoStep} the final value of \sQuote{alpha} is estimated by a two-step procedure and the value of argument \sQuote{alpha} is only used to specify the size of the samples used in the first step.}
    \item{\code{nsamp}:}{Inherited from class \code{"CovControlMcd"}. Number of subsets used for initial estimates.}
    \item{\code{scalefn}:}{Inherited from class \code{"CovControlMcd"} and not used in the package \sQuote{Maint.Data.}}
    \item{\code{maxcsteps}:}{Inherited from class \code{"CovControlMcd"} and not used in the package \sQuote{Maint.Data.}}
    \item{\code{seed}:}{Inherited from class \code{"CovControlMcd"}. Starting value for random generator. Default is \code{seed = NULL.}}
    \item{\code{use.correction}:}{Inherited from class \code{"CovControlMcd"}. Whether to use finite sample correction factors.
Default is \code{use.correction=TRUE}.}
    \item{\code{trace}, \code{tolSolve}:}{Inherited from class \code{"CovControl"}.} 
    \item{\code{ncsteps}:}{The maximum number of concentration steps used each iteration of the fasttle algorithm.}
    \item{\code{getalpha}:}{Argument specifying if the \sQuote{alpha} parameter (roughly the percentage of the sample used for computing the trimmed likelihood) should be estimated from the data, or if the value of the argument \sQuote{alpha} should be used instead. When set to \dQuote{TwoStep}, \sQuote{alpha} is estimated by a two-step procedure with the value of argument \sQuote{alpha} specifying the size of the samples used in the first step. Otherwise, with the value of argument \sQuote{alpha} is used directly.}   
    \item{\code{rawMD2Dist}:}{The assumed reference distribution of the raw MCD squared distances, which is used to find to cutoffs defining the observations kept in one-step reweighted MCD estimates. Alternatives are \sQuote{ChiSq},\sQuote{HardRockeAsF} and \sQuote{HardRockeAdjF}, respectivelly for the usual Chi-squared, and the asymptotic and adjusted scaled F distributions proposed by Hardin and Rocke (2005).}    
    \item{\code{MD2Dist}:}{The assumed reference distributions used to find cutoffs defining the observations assumed as outliers. Alternatives are \dQuote{ChiSq} and \dQuote{CerioliBetaF} respectivelly for the usual Chi-squared, and the Beta and F distributions proposed by Cerioli (2010).}    
    \item{\code{eta}:}{Nominal size of the null hypothesis that a given observation is not an outlier. Defines the raw MCD Mahalanobis distances cutoff used to choose the observations kept in the reweightening step.}
    \item{\code{multiCmpCor}:}{Whether a multicomparison correction of the nominal size (eta) for the outliers tests should be performed. Alternatives are: \sQuote{never} -- ignoring the multicomparisons and testing all entities at \sQuote{eta}. \sQuote{always} -- testing all n entitites at 1.- (1.-\sQuote{eta}^(1/n)); and \sQuote{iterstep} -- as suggested by Cerioli (2010), make an initial set of tests using the nominal size 1.- (1-\sQuote{eta}^(1/n)), and if no outliers were detected stop. Otherwise, make a second step testing for outliers at \sQuote{eta}. }   
    \item{\code{getkdblstar}:}{Argument specifying the size of the initial small (in order to minimize the probability of outliers) subsets. If set to the string \dQuote{Twopplusone} (default) the initial sets have twice the number of interval-value variables plus one (i.e., they are the smaller samples that lead to a non-singular covariance estimate). Otherwise, an integer with the size of the initial sets.}
    \item{\code{k2max}:}{Maximal allowed l2-norm condition number for correlation matrices. Correlation matrices with condition number above k2max are considered to be numerically singular, leading to degenerate results.}
    \item{\code{outlin}:}{The type of outliers to be consideres. \dQuote{MidPandLogR} if outliers may be present in both MidPpoints and LogRanges, \dQuote{MidP} if outliers are only present in MidPpoints, or \dQuote{LogR} if outliers are only present in LogRanges.}
    \item{\code{trialmethod}:}{The method to find a trial subset used to initialize each replication of the fasttle algorithm. The current options are \dQuote{simple} (default) that simply selects \sQuote{kdblstar} observations at random, and \dQuote{Poolm} that divides the original sample into \sQuote{m} non-overlaping subsets, applies the \sQuote{simple trial} and the refinement methods to each one of them, and merges the results into a trial subset.}
    \item{\code{m}:}{Number of non-overlaping subsets used by the trial method when the argument of \sQuote{trialmethod} is set to 'Poolm'.}
    \item{\code{reweighted}:}{Should a (Re)weighted estimate of the covariance matrix be used in the computation of the trimmed likelihood or just a \dQuote{raw} covariance estimate; default is (Re)weighting.}
    \item{\code{otpType}:}{The amount of output returned by fasttle. Current options are \dQuote{OnlyEst} (default) where only an \sQuote{IdtE} object with the fasttle estimates is returned, \dQuote{SetMD2andEst} which returns a list with an \sQuote{IdtE} object of fasttle estimates, a vector with the final trimmed subset elements used to compute these estimates and the corresponding robust squared Mahalanobis distances, and \dQuote{SetMD2EstandPrfSt} wich returns a list with the previous three components plust a list of some performance statistics concerning the algorithm execution.}
  }
}

\section{Extends}{
    Class \code{\linkS4class{CovControlMcd}}, directly.
    Class \code{\linkS4class{CovControl}} by CovControlMcd, distance 2.
}

\section{Methods}{
No methods defined with class "RobEstControl" in the signature.
}

\seealso{
  \code{\link{RobEstControl}}, \code{\link{fasttle}}, \code{\link{RobMxtDEst}}, \code{\link{Roblda}}, \code{\link{Robqda}}
}

\references{
Cerioli, A. (2010), Multivariate Outlier Detection with High-Breakdown Estimators.
\emph{Journal of the American Statistical Association} \bold{105} (489), 147--156. 

Duarte Silva, A.P., Filzmoser, P. and Brito, P. (2017), Outlier detection in interval data. \emph{Advances in Data Analysis and Classification}, 1--38.\cr

Hardin, J. and Rocke, A. (2005), The Distribution of Robust Distances.
\emph{Journal of Computational and Graphical Statistics} \bold{14}, 910--927.\cr

Todorov V. and Filzmoser P. (2009), An Object Oriented Framework for Robust Multivariate Analysis. \emph{Journal of Statistical Software} \bold{32}(3), 1--47.
}


\keyword{classes}
\keyword{robust}
\keyword{interval data}
