\name{data_philips}
\alias{data_philips}
\docType{data}
\title{
The philips dataset
%%   ~~ data name/kind ... ~~
}
\description{
A dataset containing measurements of \eqn{d = 9} characteristics
of \eqn{n = 677} diaphragm parts, used in the production of TV sets.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("data_philips")}
\format{
A matrix with \eqn{677} rows and \eqn{9} columns.
}

\source{
The data were provided in 1997 by Gertjan Otten and permission to analyze them was given by Herman Veraa and Frans Van Dommelen at Philips Mecoma in The Netherlands.
}

\examples{
data("data_philips")
}


\references{
Rousseeuw, P.J., and Van Driessen, K. (1999). A fast algorithm  for the Minimum Covariance Determinant estimator. \emph{Technometrics}, \bold{41}, 212--223.
}