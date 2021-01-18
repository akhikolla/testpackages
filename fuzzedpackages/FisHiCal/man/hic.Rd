\name{hic}
\alias{hic}
\docType{data}
\title{
in-silico 100*100 Hi-C contact matrix with with 5\% noise
}
\description{
in-silico distances were computed from the random configuration  \code{\link{conf}}. A power law model was then used to generate matching in-silico Hi-C frequencies. Random Noise was further added to long range frequencies to mimic a typical situation for Hi-C data. The resulting frequencies were then used to construct this matrix.
}
\usage{data(hic)}

\details{
This data structure is used to illustrate the usage of \code{\link{calibrate}}.
}
\references{
Y. Shavit, F.K. Hamey, P. Lio', FisHiCal: an R package for iterative FISH-based calibration of Hi-C data, 2014 (submitted).
}
\keyword{datasets}
