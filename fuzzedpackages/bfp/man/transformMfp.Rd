\name{transformMfp}
\alias{transformMfp}
\title{Transform a fitted mfp model into a BayesMfp model}
\description{
Transform a fitted mfp model into a BayesMfp model with the correct powers etc.
to compare with other (true) BayesMfp models fitted to the same (!) data.
}
\usage{
transformMfp(mfpObject, BayesMfpObject)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mfpObject}{the original mfp object}
  \item{BayesMfpObject}{BayesMfp object, from which the first model is
    used for imputation of the powers from mfpObject} 
}
\value{
  A BayesMfp object with the converted model.
}
\keyword{internal}
\keyword{utilities}
