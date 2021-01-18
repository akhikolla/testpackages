\name{test_acetp}
\alias{test_acetp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hypothesis testing of the ACE(t)-p models 
}
\description{
Comparison of different ACE(t)-p models to test a linear or a constant variance component.
}
\usage{
test_acetp(acetp, comp, sim = 100, robust = 0, pe = TRUE, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{acetp}{
An object from the AtCtEtp function. 

}
  \item{comp}{
The component for which linearity or constancy is tested. This component must be specified as splines or linear in the AtCtEtp function. 
}
  \item{sim}{
The number of the bootstrap resampling for approximating the null distribution when testing linearity.
}
  \item{robust}{
An integer indicating the number of different initial values that the function will randomly generate and try in the optimization. The default value is 0.
}

  \item{pe}{
A logical argument indicating whether to use penalized spline model to test linearity. The default value is TRUE.
}

  \item{verbose}{
A logical argument indicating whether to print testing information and results. The default value is TRUE.
}
 
}

\value{
  
  \item{p }{The p-value for the test.}
  \item{llr }{The LRT statistic for testing linearity.}
  \item{llr_sim }{The simulated null distribution of the LRT statistic for testing linearity.}
  \item{chisq }{The chisq statistic for testing a constant or linearity.}
  
%% ...
}

\details{
When pe=TRUE, the linearity is tested under a p-spline framework in which an LRT is performed. Otherwise, a \eqn{\chi^2} test is performed for linearity under a spline framework without penalty on smoothness.
}


\references{
He, L., Sillanpää, M.J., Silventoinen, K., Kaprio, J. and Pitkäniemi, J., 2016. Estimating Modifying Effect of Age on Genetic and Environmental Variance Components in Twin Models. Genetics, 202(4), pp.1313-1328.

He, L., Pitkäniemi, J., Silventoinen, K. and Sillanpää, M.J., 2017. ACEt: An R package for estimating dynamic heritability and comparing twin models. Behavior Genetics, 47(6), pp.620-641.
}

\author{
Liang He
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
# data(data_ace)

# result <- AtCtEtp(data_ace$mz, data_ace$dz, knot_e = 7, knot_c = 5, mod=c('d','d','l'))
# re <- test_acetp(result, comp='e')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

