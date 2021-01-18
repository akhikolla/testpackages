\name{bmpt.fia}
\alias{bmpt.fia}

\title{
Compute FIA for MPTs
}
\description{
R-port of the function to compute FIA for MPT models by Wu, Myung, and Batchelder (2010a, 2010b). This function is essentially a copy of the original Matlab code to R (with significant parts moved to C++ and allowing for multicore functionality). Also, the order of input arguments is more R-like.
}
\usage{
bmpt.fia(s, parameters, category, N, ineq0 = NULL, Sample = 2e+05, 
        multicore = FALSE, split = NULL, mConst = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s}{
see Details
}
  \item{parameters}{
see Details
}
  \item{category}{
see Details
}
  \item{N}{
see Details
}
  \item{ineq0}{
see Details
}
  \item{Sample}{
see Details
}
  \item{multicore}{
logical. Should fitting be distributed across several cores? Requires \pkg{snowfall} and initialized cluster. See also below.
}
  \item{split}{
\code{NULL} (the default) or integer specifying in how many separate calls to the C++ workhorse the integrant should be calculated. See below.
}
  \item{mConst}{
A constant which is added in the Monte Carlo integration to avoid numerical underflows and is later subtracted (after appropriate transformation). Should be a power of 2 to avoid unnecessary numerical imprecision. 
}
}
\details{
The following is the original description by Wu, Myung, & Batchelder (2010a) for their Matlab function. All changes to the original document are in squared brackets []:

This function computes the FIA complexity measure, C_FIA, using a Monte Carlo numerical integration algorithm. When inequality is present, sampling from the restricted parameter space is performed by rejection algorithm.

[...] [see References for References]

The following symbols are used in the body of the function:\cr
S denotes number of parameters.\cr
C denotes the number of categories.\cr
M denotes the number of leaves in the tree.

The first input argument \code{s} is related to the string representation of the BMPT model. It can be obtained by replacing all categories in the string by the capital letter C and all branching probabilities by the lower case letter p. 

The second input argument \code{parameters} is a row vector that assigns parameters or constants to the p's in the  string \code{s}. Its length should be the same as the number of p's in \code{s}, and its elements correspond to the p's according to their order in \code{s}. Positive integer elements in \code{parameters} assign parameters to the corresponding p's, with the same integer denoting the same parameter. Constants are assigned to the p's using the negation of their values. 

The [third] input argument \code{category} is a 1 by M vector assigning categories to the C's in the string `s' in the same way \code{parameters} assigns branching probabilities, except that only positive consecutive integers from 1 to \code{J}, the total number of categories, are allowed. 

The [fourth] input argument \code{N} specifies the total sample size.

The [fifth] input argument \code{ineq0} assigns inequality constraints imposed on the parameters. It is a matrix with two columns. Each element denotes a parameter coded in the same way as in \code{parameters}. For each row, the parameter on the left column is constrained to be smaller than that on the right column. The number of rows is determined by the total number of simple inequality constraints of the form \code{theta_1 < theta_2} in the model. [Default is \code{NULL} corresponding to no inequality restrictions.]

The last input argument `Sample' specifies the number of random samples to be drawn in the Monte Carlo algorithm. [Default is 200000.]

[For returned values see Value]

It should be noted that `lnconst' can be computed analytically free of Monte Carlo error on a case by case basis described below. For this reason, the users can calculate \code{C_FIA} [see Wu, Myung & Batchelder, 2010a; Equation 7] by adding \code{(S/2)*ln*(N/(2*pi))}, \code{lnInt} and their hand-calculated \code{lnconst} to minimize the Monte Carlo errors. [In our experience this error is rather low and negligible.]

A sequence of inequalities \code{theta_1 < theta_2 < ... < theta_k} reduces the parameter space to its \code{1/k!}, so in this case \code{lnconst} should be \code{-ln * (k!)}. In general, any combination of inequality constraints specifies a union of subsets of the parameter space, each satisfying some sequence of inequalities. For example, the subspace defined by \code{theta_1 < theta_2} and \code{theta_3 < theta_2} is a union of two subspaces, one satisfying \code{theta_1 < theta_3 < theta_2} and the other \code{theta_3 < theta_1 < theta_2}, so the proportion is given by \code{2 * (1/3!) = 1/3}.

A coding example:\cr
Suppose that for model 1HTM-5c of source monitoring [see Wu et al., 2010a] , the sample sizes of source A, source B and new items are 300, 300 and 400, respectively and the inequality constraint of \code{d_1 < d_2} is imposed. In this case, the six input arguments should be specified as follows:\cr
s = 'ppppCpCCppCCCppCpCCppCCCppCCC';\cr
parameters = c(-.6,-.5,1,2,5,4,5,1,3,5,4,5,4,5); [adapted for R]\cr
ineq0 = matrix(c(2,3), 1,2); [adapted for R]\cr
category = c(1,1,2,1,2,3,5,4,5,4,5,6,7,8,9); [adapted for R]\cr
N = 1000;

Another coding example:\cr
For the pair-clustering model in Batchelder and Riefer (1999, Figure 1), suppose in a pair-clustering experiment there are 300 pairs of words and 100 singletons, the six input arguments should be specified as follows:\cr
s = 'pppCCppCCpCCpCC';
parameters = c(-.75,1,2,3,3,3,3); [adapted for R]\cr
ineq0 = NULL;  [adapted for R]\cr
category = c(1,4,2,3,3,4,5,6); [adapted for R]\cr
N = 400;

[For more examples, see Examples]

Since MPTinR version 1.1.3 the Monte Carlo integration is performed in C++ using \pkg{RcppEigen}. With the default arguments, one instance of the C++ workhorse is called. To call multiple instances of the C++ workhorse, you can use the \code{split} argument (which can be useful to replicate results obtained with \code{multicore = TRUE} as described below). Note, that each time before calling the C++ code, the seed is set (the set of random seeds are generated before calling the function for the first time).

Multicore functionality is achieved via \pkg{snowfall} which needs to be loaded and a cluster initialized via \code{\link[snowfall]{sfInit}} when setting \code{multicore = TRUE}. When \code{split = NULL} (the default), the \code{Samples} will be evenly distributed on the different cores (using \code{\link[snowfall]{sfClusterSplit}}), so that only one instance of the underlying C++ workhorse is called on each core. Setting \code{split} to non-\code{NULL} will produce as many instances (distributed across cores). Note that in order to obtain comparable results (as snowfall uses load balancing), the random seed is set (at each core) before calling each instance of the C++ workhorse. This allows to replicate results obtained via multicore in a non-multicore environment when seting \code{split} appropriately (and \code{set.seed} beforehands).
}
\value{
[A named vector:]

The first output argument \code{CFIA} gives the FIA complexity value of the model.

The second [and third] output argument \code{CI} gives the Monte Carlo confidence interval of \code{CFIA}. [\code{CI.l}, gives the lower, \code{CI.u}, the upper bound of the interval].

The [fourth] output argument \code{lnInt} gives the log integral term in \code{C_FIA} [see Wu, Myung & Batchelder, 2010a; Equation 7] for models without inequality constraints. When inequality constraints are present, \code{lnInt} does not take into account the change in the normalizing constant in the proposal distribution and must be adjusted with the output argument \code{lnconst}.

The [fifth and sixth] output argument [\code{CI.lnint}] gives the Monte Carlo confidence interval of \code{lnInt}. [.l = lower & .u = upper bound of the CI]

When inequality constraints are present, the [seventh] output argument \code{lnconst} serves as an adjustment of `lnInt'. It estimates the logarithm of the proportion of parameter space [0,1]^S that satisfies those inequality constraints, and the log integral term is given by lnInt+lnconst.

The next [two] output argument [\code{CI.lnconst}] give the Monte Carlo confidence interval of `lnconst'. [.l = lower & .u = upper bound of the CI]
}
\references{
Wu, H., Myung, J.I., & Batchelder, W.H. (2010a). Minimum description length model selection of multinomial processing tree models. \emph{Psychonomic Bulletin & Review}, 17, 275-286.

Wu, H., Myung, J.I., & Batchelder, W.H. (2010b). On the minimum description length complexity of multinomial processing trees. \emph{Journal of Mathematical Psychology}, 54, 291-303.
}
\author{
The original Matlab code was written by Hao Wu, Jay I. Myung, and William H. Batchelder.\cr This code was ported to R by Henrik Singmann and David Kellen. RcppEigen was added by Henrik Singmann and Christian Mueller. Multicore functionality was added by Henrik Singmann.
}
\note{
The R version of the code should now (after moving the code to RcppEigen) be considerably faster than the Matlab version of this code.
}
\seealso{
\code{\link{fit.mpt}} for the main function of MPTinR.\cr
\code{\link{get.mpt.fia}} for a convenient wrapper of this function.
}
\examples{
\dontrun{
# The following example is the code for the first example in Wu, Myung & Batchelder (2010a, pp. 280)
# The result should be something like: CFIA = 12.61... or 12.62..., CI = 12.61... - 12.62....
# Executing this command can take a while.

bmpt.fia(s = "ppppCpCCppCCCppCpCCppCCCppCCC", 
	parameters = c(-0.5, -0.5, 3, 2, 5, 1, 5, 4, 2, 5, 1, 5, 1, 5), 
	category = c(1,1,2,1,2,3,5,4,5,4,5,6,7,8,9), 
  N = 1000, ineq0 = matrix(c(4,3),1,2))

bmpt.fia(s = "ppppCpCCppCCCppCpCCppCCCppCCC", 
	parameters = c(-0.5, -0.5, 3, 2, 5, 1, 5, 4, 2, 5, 1, 5, 1, 5), 
	category = c(1,1,2,1,2,3,5,4,5,4,5,6,7,8,9), 
  N = 1000, ineq0 = matrix(c(4,3),1,2), mConst = 2L^8)
}
  }

