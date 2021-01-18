\name{acetp_mcmc}
\alias{acetp_mcmc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute CIs for the ACE(t)-p model
}
\description{
Compute the posterior mean and CIs for the ACE(t)-p model using the MCMC methods
}
\usage{
acetp_mcmc(acetp, iter_num = 10000, sd = 0.1, burnin = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{acetp}{
An object from the 'AtCtEtp' function.
}
  \item{iter_num}{
The number of the iterations in the MCMC procedure.
}
  \item{sd}{
The standard error of the normal proposal distribution in the MCMC algorithm. The default value is 0.1.
}
  \item{burnin}{
The number of burn-in, which must be smaller than the number of iteration.
}
}

\value{

  \item{beta_a_mc }{The estimates of the spline coefficients for the A component based on the posterior mean from the MCMC method.}
  \item{beta_c_mc }{The estimates of the spline coefficients for the C component based on the posterior mean from the MCMC method.}
  \item{beta_e_mc }{The estimates of the spline coefficients for the E component based on the posterior mean from the MCMC method.}
  \item{cov_mc }{The posterior covariance matrix of the estimates of the spline coefficients.}
  \item{knots_a }{A vector of the knot positions for the A component.}
  \item{knots_c }{A vector of the knot positions for the C component.}
  \item{knots_e }{A vector of the knot positions for the E component.}

%% ...
}
\references{
He, L., Sillanpää, M.J., Silventoinen, K., Kaprio, J. and Pitkäniemi, J., 2016. Estimating Modifying Effect of Age on Genetic and Environmental Variance Components in Twin Models. Genetics, 202(4), pp.1313-1328.
}
\author{
Liang He
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{

# data(data_ace)

# result <- AtCtEp(data_ace$mz, data_ace$dz, knot_a = 7, knot_c = 7)
# result_mc <- acetp_mcmc(result, iter_num=10000, burnin = 500)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

