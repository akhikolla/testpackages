
\name{FHUBdiscrete}
\alias{FHUBdiscrete}
\title{Compute the Frechet Hoeffding Upper Bound for Given Discrete Marginal Distributions}
\usage{
FHUBdiscrete(marg1, marg2, mu1, mu2, od1 = 0, od2 = 0, binomial.size1 = 1,
  binomial.size2 = 1)
}
\description{
	This function implemented the method of computing the Frechet Hoeffding upper bound for discrete marginals
         described in Nelsen (1987), which can only be applied to discrete marginals. Four commonly used
         marginal distributions were implemented. The distribution "nb" (negative binomial) and "zip"
         (zero-inflated Poisson) are parameterized in terms of the mean and overdispersion,
         see Han and De Oliveira (2016).
}

\arguments{
\item{marg1}{ name of the first discrete marginal distribution. Should be one of the
"poisson", "zip", "nb" or "binomial".}

\item{marg2}{ name of the second discrete marginal distribution. Should be one of the
"poisson", "zip", "nb" or "binomial".}

\item{mu1}{mean of the first marginal distribution. If binomial then it is \eqn{n_1 p_1}.}

\item{mu2}{mean of the second marginal distribution. If binomial then it is \eqn{n_2 p_2}.}

\item{od1}{the overdispersion parameter of the first marginal. Only used when marginal distribution is either
"zip" or "nb".}

\item{od2}{the overdispersion parameter of the second marginal. Only used when marginal distribution is either
"zip" or "nb".}

\item{binomial.size1}{the size parameter (number of trials) when \code{marg1 = "binomial"}.}

\item{binomial.size2}{the size parameter (number of trials) when \code{marg2 = "binomial"}.}
}
\value{
A scalar denoting the Frechet Hoeffding upper bound of the two
        specified marginal.
}


\examples{
\dontrun{

FHUBdiscrete(marg1 = 'nb', marg2 = 'zip',mu1 = 10, mu2 = 2, od1 = 2, od2 = 0.2)
FHUBdiscrete(marg1 = 'binomial', marg2 = 'zip', mu1 = 10, mu2 = 4, binomial.size1 = 25, od2 = 2)
FHUBdiscrete(marg1 = 'binomial', marg2 = 'poisson', mu1 = 0.3, mu2 = 20, binomial.size1 = 1)

NBmu = seq(0.01, 30, by = 0.02)
 fhub <- c()
 for(i in 1:length(NBmu)){
  fhub[i] = FHUBdiscrete(marg1 = 'nb', marg2 = 'nb',mu1 = 10, mu2 = NBmu[i], od1 = 0.2, od2 = 0.2)
}
plot(NBmu, fhub, type='l')
}
}
\author{
Zifei Han \email{hanzifei1@gmail.com}
}
\references{
Nelsen, R. (1987) Discrete bivariate distributions with given marginals and correlation.
             \emph{Communications in Statistics Simulation and Computation, 16:199-208}.

 Han, Z. and De Oliveira, V. (2016) On the correlation structure of Gaussian
copula models for geostatistical count data.
 \emph{Australian and New Zealand Journal of Statistics, 58:47-69}.

Han, Z. and De Oliveira, V. (2018) gcKrig: An R Package for the Analysis of Geostatistical Count Data Using Gaussian Copulas.
\emph{Journal of Statistical Software}, \bold{87}(13), 1--32. \doi{10.18637/jss.v087.i13}.

}

\keyword{Correlation Structure}
