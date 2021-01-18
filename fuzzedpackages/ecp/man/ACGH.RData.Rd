\name{ACGH}
\alias{ACGH}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
Bladder Tumor Micro-Array Data
}
\description{
Micro-array data for 43 different individuals with a bladder tumor.

}
\usage{data(ACGH)}
\format{
  A list with the following components.
  
  data: The micro-array data for 43 individuals. This information is stored in a 2215 by 43 matrix.
  
  individual: A numeric vector indicating which individuals' mico-array data are present.
  
}
\source{
Bleakley K., Vert J.-P. (2011), The group fused Lasso for multiple change-point detection

N. Stransky, C. Vallot, F. Reyal, I. Bernard-Pierrot, S.G. Diez de Mediana, R. Segraves, Y. de Rycke,
P. Elvin, A. Cassidy, C. Sparaggon, A. Graham, j. Southgate, B. Asselain, Y. Allory, C. C. Addou, D. G.
Albertson, J.-P. Thiery, D. K. Chopin, D. Pinkel, and F. Radvanyi. Regional copy number-independent 
deregulation of transcription in cancer. Nat. Genet., 38(12):1386-1396, Dec 2006
}
\references{
Bleakley K., Vert J.-P. (2011), The group fused Lasso for multiple change-point detection

Nicholas A. James, David S. Matteson (2014). "ecp: An R Package for Nonparametric 
Multiple Change Point Analysis of Multivariate Data.", "Journal of Statistical Software, 
62(7), 1-25", URL "http://www.jstatsoft.org/v62/i07/"
}
\examples{
data(ACGH, package="ecp")
}
\keyword{datasets}
