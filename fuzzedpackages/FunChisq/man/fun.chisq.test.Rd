\name{fun.chisq.test}
\alias{fun.chisq.test}
\title{
Model-Free Functional Chi-Squared and Exact Tests
}

\description{
Asymptotic chi-squared, normalized chi-squared or exact tests on contingency tables to determine model-free functional dependency of the column variable on the row variable.
}

\usage{
fun.chisq.test(
  x,
  method = c("fchisq", "nfchisq",
             "exact", "exact.qp", "exact.dp", "exact.dqp",
             "default", "normalized", "simulate.p.value"),
  alternative = c("non-constant", "all"), log.p=FALSE,
  index.kind = c("conditional", "unconditional"),
  simulate.nruns = 2000,
  exact.mode.bound=TRUE
)
}

\arguments{
  \item{x}{
  a matrix representing a contingency table. The row variable represents the independent variable or all unique combinations of multiple independent variables. The column variable is the dependent variable.
}

  \item{method}{
  a character string to specify the method to compute the functional chi-squared test statistic and its p-value. The options are \code{"fchisq"} (equivalent to \code{"default"}, the default), \code{"nfchisq"} (equivalent to \code{"normalized"}), \code{"exact"}, \code{"exact.qp"}, \code{"exact.dp"}, \code{"exact.dqp"} or \code{"simulate.p.value"}. See Details.

  Note: \code{"default"} and \code{"normalized"} are deprecated.
}
  \item{alternative}{
  a character string to specify the alternative hypothesis. The options are \code{"non-constant"} (default, non-constant functions) and \code{"all"} (all types of functions including constant ones).
  }
  \item{log.p}{
  logical; if \code{TRUE}, the p-value is given as \code{log(p)}. Taking the log improves the accuracy when p-value is close to zero. The default is \code{FALSE}.
  }
  \item{index.kind}{
  a character string to specify the kind of function index xi.f to be estimated. The options are \code{"conditional"} (default) and \code{"unconditional"}. See Details.
  }
  \item{simulate.nruns}{
   A number to specify the number of tables generated to simulate the null distribution. Default is \code{2000}. Only used when \code{method="simulate.p.value"}.
  }
  \item{exact.mode.bound}{
  logical; if \code{TRUE}, a fast branch-and-bound algorithm is used for the exact functional test (\code{method="exact"}). If \code{FALSE}, a slow brute-force enumeration method is used to provide a  reference for runtime analysis. Both options provide the same exact p-value. The default is \code{TRUE}.
  }
}

\details{

The functional chi-squared test determines whether the column variable is a function of the row variable in contingency table \code{x} \insertCite{zhang2013deciphering,zhang2014nonparametric}{FunChisq}. This function supports three hypothesis testing methods:

When \code{method="fchisq"} (equivalent to \code{"default"}, the default), the test statistic is computed as described in \insertCite{zhang2013deciphering,zhang2014nonparametric}{FunChisq} and the p-value is computed using the chi-squared distribution.

When \code{method="nfchisq"} (equivalent to \code{"normalized"}), the test statistic is obtained by shifting and scaling the original test statistic \insertCite{zhang2013deciphering,zhang2014nonparametric}{FunChisq}; and the p-value is computed using the standard normal distribution \insertCite{Box2005}{FunChisq}. The normalized chi-squared, more conservative on the degrees of freedom, was used by the Best Performer NMSUSongLab in HPN-DREAM (DREAM8) Breast Cancer Network Inference Challenges.

When \code{method="exact"}, \code{"exact.qp"} (quadratic programming) \insertCite{zhong2019eft,zhong2019modelfree}{FunChisq}, \code{"exact.dp"} (dynamic programming) \insertCite{nguyen2018modelfree,Nguyen2020EFT}{FunChisq}, or \code{"exact.dqp"} (dynamic and quadratic programming) \insertCite{nguyen2018modelfree,Nguyen2020EFT}{FunChisq}, an exact functional test is performed. The option of \code{"exact"} uses \code{"exact.dqp"}, the fastest method. All methods compute an exact p-value.

For the \code{"exact.qp"} and \code{"exact.dp"} options, if the sample size is no more than 200 or the average cell count is less than five, and the table size is no more than 10 in either row or column, the exact test will not be called and the asymptotic functional chi-squared test (\code{method="fchisq"}) is used instead.  %When the sample size of an input contingency table is large, the exact functional test and the functional chi-squared test will return similar p-values. When the expected values for all entries in a contingency table are greater than 5, the asymptotic tests will perform similarly with the exact test.
% On a 2-by-2 contingency tables, Fisher's exact test (Fisher, 1922) will be applied.

For \code{"exact.dqp"}, the exact functional test will always be performed.

For 2-by-2 contingency tables, the asymptotic test options (\code{method="fchisq"} or \code{"nfchisq"}) are recommended to test functional dependency, instead of the exact functional test.

When \code{method="simulate.p.value"}, a simulated null distribution is used to calculate \code{p-value}. The null distribution is a multinomial distribution that is the product of two marginal distributions. Like other Monte Carlo based methods, this method is slower but may be more accurate than other methods based on asymptotic distributions.

\code{index.kind} specifies the kind of function index to be computed. If the experimental design controls neither the row nor column marginal sums, \code{index.kind = "unconditional"} is recommended; If the column marginal sums are controlled, \code{index.kind = "conditional"} is recommended. The choice of \code{index.kind} affects only the function index xi.f value, but not the test statistic or p-value.
}

\value{
A list with class "\code{htest}" containing the following components:
\item{statistic}{the functional chi-squared statistic if \code{method = "fchisq"}, \code{"default"}, or \code{"exact"}; or the normalized functional chi-squared statistic if \code{method = "nfchisq"} or \code{"normalized"}.}
\item{parameter}{degrees of freedom for the functional chi-squared statistic.}
\item{p.value}{p-value of the functional test. If \code{method = "fchisq"} (or \code{"default"}), it is computed by an asymptotic chi-squared distribution; if \code{method = "nfchisq"} (or \code{"normalized"}), it is computed by the standard normal distribution; if \code{method} \code{= "exact"}, it is computed by an exact hypergeometric distribution.}
\item{estimate}{an estimate of function index between 0 and 1. The value of 1 indicates a strictly mathematical function. It is asymmetrical with respect to transpose of the input contingency table, different from the symmetrical Cramer's V based on the Pearson's chi-squared test statistic. See \insertCite{Zhong2019FANTOM5,KumarZSLS18}{FunChisq} for the definition of function index.}
}

\references{
\insertAllCited{}
}


\author{
Yang Zhang, Hua Zhong, Hien Nguyen, and Joe Song
}

\seealso{
For data discretization, an option is optimal univariate clustering via package \pkg{Ckmeans.1d.dp}. A second option is joint multivariate discretization via package \pkg{GridOnClusters}.

For symmetrical dependency tests on discrete data, see Pearson's chi-squared test \code{\link[stats]{chisq.test}}, Fisher's exact test \code{\link[stats]{fisher.test}}, and mutual information methods in package \pkg{entropy}.
}

\examples{
\donttest{
# Example 1. Asymptotic functional chi-squared test
x <- matrix(c(20,0,20,0,20,0,5,0,5), 3)
fun.chisq.test(x) # strong functional dependency
fun.chisq.test(t(x)) # weak functional dependency

# Example 2. Normalized functional chi-squared test
x <- matrix(c(8,0,8,0,8,0,2,0,2), 3)
fun.chisq.test(x, method="nfchisq") # strong functional dependency
fun.chisq.test(t(x), method="nfchisq") # weak functional dependency

# Example 3. Exact functional chi-squared test
x <- matrix(c(4,0,4,0,4,0,1,0,1), 3)
fun.chisq.test(x, method="exact") # strong functional dependency
fun.chisq.test(t(x), method="exact") # weak functional dependency

# Example 4. Exact functional chi-squared test on a real data set
#            (Shen et al., 2002)
# x is a contingency table with row variable for p53 mutation and
#   column variable for CIMP
x <- matrix(c(12,26,18,0,8,12), nrow=2, ncol=3, byrow=TRUE)

# Test the functional dependency: p53 mutation -> CIMP
fun.chisq.test(x, method="exact")

# Test the functional dependency CIMP -> p53 mutation
fun.chisq.test(t(x), method="exact")

# Example 5. Asymptotic functional chi-squared test with simulated distribution
x <- matrix(c(20,0,20,0,20,0,5,0,5), 3)
fun.chisq.test(x, method="simulate.p.value")
fun.chisq.test(x, method="simulate.p.value", simulate.n = 1000)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\keyword{htest}
\keyword{nonparametric}

