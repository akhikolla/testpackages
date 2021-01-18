\name{bootstrap_corr}
\alias{bootstrap_corr}
\title{
Computation of a confidence interval for the estimated heritability.
}
\description{
This function is based on a non parametric bootstrap technique to compute a confidence interval for the heritability. The strength of this method is that it can deal with correlated observations.
}
\usage{
bootstrap_corr(Y,Z,K,eta_hat,sigma2_hat,level,nb_cores)
}
\arguments{
  \item{Y}{
Vector of observations.
}
  \item{Z}{
Matrix of genetic informations.
}
  \item{K}{
number of subsamples of Y used to apply our bootstrap technique.
}
  \item{eta_hat}{
Estimator of the heritability.
}
  \item{sigma2_hat}{
Estimator of the variance involving the variances of the two random parts of the model.
}
  \item{level}{
Percentage of values which will be removed from the estimated heritabilities to build a confidence interval.
}
  \item{nb_cores}{
Number of cores of the computer. It is used for parallelizing the computations.
}
}
\value{
 \item{CI_up}{Upper bound of the confidence interval for the estimated heritability}
\item{CI_low}{Lower bound of the confidence interval for the estimated heritability}
}

\author{
Anna Bonnet and Celine Levy-Leduc
}

