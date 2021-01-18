\name{FunChisq-package}
\alias{FunChisq-package}
\docType{package}
\title{
Model-Free Functional Chi-Squared and Exact Tests
}
\description{
Statistical hypothesis testing methods for model-free functional dependency using asymptotic chi-squared or exact distributions. Functional chi-squared test statistics \insertCite{zhang2013deciphering,zhang2014nonparametric,nguyen2018modelfree,zhong2019modelfree,zhong2019eft,Nguyen2020EFT}{FunChisq} are asymmetric, functionally optimal, and model-free, unique from other related statistical measures.

Tests in this package reveal evidence for causality based on the causality-by-functionality principle \insertCite{Simon1966}{FunChisq}. The tests require data from two or more variables be formatted as a contingency table. Continuous variables need to be discretized first, for example, using R packages \pkg{Ckmeans.1d.dp} or \pkg{GridOnClusters}.

The package implements asymptotic functional chi-squared tests \insertCite{zhang2013deciphering,zhang2014nonparametric}{FunChisq} and an exact functional test \insertCite{nguyen2018modelfree,zhong2019modelfree,zhong2019eft,Nguyen2020EFT}{FunChisq}. The normalized functional chi-squared test was used by Best Performer NMSUSongLab in HPN-DREAM (DREAM8) Breast Cancer Network Inference Challenges \insertCite{Hill:2016fk}{FunChisq}.

A function index derived from the functional chi-squared offers a new effect size measure for the strength of function dependency. It is asymmetrically functionally optimal, different from the symmetric Cramer's V, also a better alternative to conditional entropy in many aspects.

A simulator is provided to generate functional, dependent non-functional, and independent patterns \insertCite{sharma2017simulating}{FunChisq}.

For continuous data, these tests offer an advantage over regression analysis when a parametric form cannot be reliably assumed for the underlying function. For categorical data, they provide a novel means to assess directional dependency not possible with symmetrical Pearson's chi-squared test, G-test, or Fisher's exact test.
}
\details{
\tabular{ll}{
Package: \tab FunChisq\cr
Type: \tab Package\cr
Current version: \tab 2.5.1\cr
Initial release version: \tab 1.0\cr
Initial release date: \tab 2014-03-08\cr
License: \tab LGPL (>= 3) \cr
}

}

\seealso{
For data discretization, an option is optimal univariate clustering via package \pkg{Ckmeans.1d.dp}. A second option is joint multivariate discretization via package \pkg{GridOnClusters}.

For symmetric dependency tests on discrete data, see Pearson's chi-squared test (\code{\link[stats]{chisq.test}}), Fisher's exact test (\code{\link[stats]{fisher.test}}), mutual information (package \pkg{entropy}), and G-test, implemented in packages \pkg{DescTools} and \pkg{RVAideMemoire}.
}

\author{
Yang Zhang, Hua Zhong, Hien Nguyen, Ruby Sharma, Sajal Kumar, and Joe Song
}

\references{
\insertAllCited{}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\keyword{package}
\keyword{htest}
\keyword{datagen}
\keyword{nonparametric}
