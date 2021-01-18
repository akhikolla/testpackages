\name{dataCervical}
\alias{dataCervical}
\alias{mirCerv}
\alias{respCerv}
\alias{mirCons}
\docType{data}
\title{
Contains three R-objects, including the data and the binary response

}
\description{The three objects are:
mirCerv: 772 sequenced microRNAs for 56 samples; 
respCerv: binary response coded as healthy and pre-cursor lesions for cervical cancer; 
 and mirCons: conservation status of the microRNAs in three levels coded as only found in humans, found in most mammals, and found in most vertebrates
}
\usage{data(dataCervical)}
\format{
The formats are:
\describe{
\item{mirCerv}{matrix [1:56, 1:772]}
\item{respCerv}{Factor w/ 2 levels "CIN3","Normal"}
\item{mirCons}{Factor w/ 3 levels "NotCons", "Mammals", "Broadly"}
}
}
\details{
This data is used for illustration in the vignette of the \code{\link{gren}} package.
}
\source{
Novianti, P.W., Snoek, B.C., Wilting, S.M., and van de Wiel, M.A. (2017). Better diagnostic signatures from RNAseq data through use of auxiliary co-data. Bioinformatics, 33, 1572–1574.
}
\references{
Münch, M.M., Peeters, C.F.W., van der Vaart, A.W., and van de Wiel, M.A. (2018). Adaptive group-regularized logistic elastic net regression. arXiv:1805.00389v1 [stat.ME].
}
\examples{
data(dataCervical)
str(mirCerv)
str(respCerv)
str(mirCons)
}
\keyword{datasets}
