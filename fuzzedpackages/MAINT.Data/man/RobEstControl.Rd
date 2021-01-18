\name{RobEstControl}
\alias{RobEstControl}

\title{ Constructor function for objects of class RobEstControl}
\description{
  This function will create a control object of class \code{\linkS4class{RobEstControl}}
  containing the control parameters for the robust estimation functions \code{\link{fasttle}}, 
  \code{\link{RobMxtDEst}}, \code{\link{Roblda}} and \code{\link{Robqda}}.  
}

\usage{
RobEstControl(alpha=0.75, nsamp=500,  seed=NULL, trace=FALSE, use.correction=TRUE,
  ncsteps=200, getalpha="TwoStep", rawMD2Dist="ChiSq", MD2Dist="ChiSq", eta=0.025,
  multiCmpCor="never",  getkdblstar="Twopplusone", outlin="MidPandLogR", 
  trialmethod="simple", m=1, reweighted=TRUE, k2max=1e6, otpType="SetMD2andEst")
}

\arguments{

  \item{alpha}{Numeric parameter controlling the size of the subsets over which the trimmed likelihood is maximized; roughly alpha*nrow(Idt) observations are used for computing the trimmed likelihood. Allowed values are between 0.5 and 1.
Note that when argument \sQuote{getalpha} is set to \dQuote{TwoStep} the final value of \sQuote{alpha} is estimated by a two-step procedure and the value of argument \sQuote{alpha} is only used to specify the size of the samples used in the first step.}
  \item{nsamp}{Number of subsets used for initial estimates.}
  \item{seed}{Starting value for random generator.}
  \item{trace}{Whether to print intermediate results.}
  \item{use.correction}{Whether to use finite sample correction factors.}
  \item{ncsteps}{The maximum number of concentration steps used each iteration of the fasttle algorithm.}
  \item{getalpha}{Argument specifying if the \sQuote{alpha} parameter (roughly the percentage of the sample used for computing the trimmed likelihood) should be estimadted from the data, or if the value of the argument \sQuote{alpha} should be used instead. When set to \dQuote{TwoStep}, \sQuote{alpha} is estimated by a two-step procedure with the value of argument \sQuote{alpha} specifying the size of the samples used in the first step. Otherwise the value of argument \sQuote{alpha} is used directly.}   
    \item{rawMD2Dist}{The assumed reference distribution of the raw MCD squared distances, which is used to find to cutoffs defining the observations kept in one-step reweighted MCD estimates. Alternatives are \sQuote{ChiSq},\sQuote{HardRockeAsF} and \sQuote{HardRockeAdjF}, respectivelly for the usual Chi-squared, and the asymptotic and adjusted scaled F distributions proposed by Hardin and Rocke (2005).}    
    \item{MD2Dist}{The assumed reference distributions used to find cutoffs defining the observations assumed as outliers. Alternatives are \dQuote{ChiSq} and \dQuote{CerioliBetaF} respectivelly for the usual Chi-squared, the Beta and F distributions proposed by Cerioli (2010).}    
    \item{eta}{Nominal size of the null hypothesis that a given observation is not an outlier. Defines the raw MCD Mahalanobis distances cutoff used to choose the observations kept in the reweightening step.}
    \item{multiCmpCor}{Whether a multicomparison correction of the nominal size (eta) for the outliers tests should be performed. Alternatives are: \sQuote{never} -- ignoring the multicomparisons and testing all entities at \sQuote{eta}. \sQuote{always} -- testing all n entitites at 1.- (1.-\sQuote{eta}^(1/n)); and \sQuote{iterstep} -- as sugested by Cerioli (2010), make an initial set of tests using the nominal size 1.- (1-\sQuote{eta}^(1/n)), and if no outliers were detected stop. Otherwise, make a second step testing for outliers at \sQuote{eta}. }   
  \item{getkdblstar}{Argument specifying the size of the initial small (in order to minimize the probability of outliers) subsets. If set to the string \dQuote{Twopplusone} (default) the initial sets have twice the number of interval-value variables plus one which are they are the smaller samples that lead to a non-singular covaraince estimate). Otherwise, an integer with the size of the initial sets.}
  \item{outlin}{The type of outliers to be considered. \dQuote{MidPandLogR} if outliers may be present in both MidPpoints and LogRanges, \dQuote{MidP} if outliers are only present in MidPpoints, or \dQuote{LogR} if outliers are only present in LogRanges.}
  \item{trialmethod}{The method to find a trial subset used to initialize each replication of the fasttle algorithm. The current options are \dQuote{simple} (default) that simply selects \sQuote{kdblstar} observations at random, and \dQuote{Poolm} that divides the original sample into \sQuote{m} non-overlaping subsets, applies the \sQuote{simple trial} and the refinement methods to each one of them, and merges the results into a trial subset.}
  \item{m}{Number of non-overlaping subsets used by the trial method when the argument of \sQuote{trialmethod} is set to 'Poolm'.}
  \item{reweighted}{Should a (Re)weighted estimate of the covariance matrix be used in the computation of the trimmed likelihood or just a \dQuote{raw} covariance estimate; default is (Re)weighting.}
  \item{k2max}{Maximal allowed l2-norm condition number for correlation matrices. Correlation matrices with condition number above k2max are considered to be numerically singular, leading to degenerate results.}
  \item{otpType}{The amount of output returned by fasttle. \cr Current options are \dQuote{SetMD2andEst} (default) which returns an \sQuote{IdtSngNDRE} object with the fasttle estimates,
a vector with the final trimmed subset elements used to compute these estimates and the corresponding robust squared Mahalanobis distances, and \dQuote{SetMD2EstandPrfSt} wich returns an \sQuote{IdtSngNDRE} object with the previous slots plust a list of some performance statistics concerning the algorithm execution.}

}

\value{
A \code{RobEstControl} object
}

\seealso{
  \code{\linkS4class{RobEstControl}}, \code{\link{fasttle}}, \code{\link{RobMxtDEst}}, \code{\link{Roblda}}, \code{\link{Robqda}}
}


\references{
Brito, P., Duarte Silva, A. P.  (2012): "Modelling Interval Data with Normal and Skew-Normal Distributions". Journal of Applied Statistics, Volume 39, Issue 1, 3-20.\cr

Cerioli, A. (2010), Multivariate Outlier Detection with High-Breakdown Estimators.
\emph{Journal of the American Statistical Association} \bold{105} (489), 147--156.\cr 

Duarte Silva, A.P., Filzmoser, P. and Brito, P. (2017), Outlier detection in interval data. \emph{Advances in Data Analysis and Classification}, 1--38.\cr

Hardin, J. and Rocke, A. (2005), The Distribution of Robust Distances.
\emph{Journal of Computational and Graphical Statistics} \bold{14}, 910--927.\cr

Todorov V. and Filzmoser P. (2009), An Object Oriented Framework for Robust Multivariate Analysis. \emph{Journal of Statistical Software} \bold{32}(3), 1--47.
}

