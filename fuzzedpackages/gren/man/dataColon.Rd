\name{dataColon}
\alias{dataColon}
\alias{mirCol}
\alias{unpenCol}
\alias{respCol}
\alias{mirExpr}
\docType{data}
\title{
Contains four R-objects, including the data and the binary response
}
\description{The four objects are:
mirCol: 2114 sequenced microRNAs for 88 samples; 
unpenCol: 4 clinical covariates for 88 samples;
respCol: binary response coded as progressive disease and benefitted from therapy;
mirExpr: expression levels of microRNAs from previous study coded as not differentially expressed, medium significantly expressed, and highly significantly expressed.
}
\usage{data(dataColon)}
\format{
The formats are:
\describe{
\item{mirCol}{matrix [1:88, 1:2114]}
\item{unpenCol}{data.frame: 88 obs. of 4 variables}
\item{respCol}{Factor w/ 2 levels "Progr","TherBenefit"}
\item{mirExpr}{Factor w/ 3 levels "nonExpr", "medExpr", "highExpr"}
}
}
\details{
This data is used for illustration in the vignette of the \code{\link{gren}} package and the corresponding paper (see references).
}
\source{
Neerincx, M., Poel, D., Sie, D.L.S., van Grieken, N.C.T., Shankaraiah, R.C., van der Wolf - de Lijster, F.S.W., van Waesberghe, J.H.T.M., Burggraaf, J.D., Eijk, P.P., Verhoef, C., Ylstra, B., Meijer, G.A., van de Wiel, M.A., Buffart, T.E., and others. (2018). Combination of a six microRNA expression profile with four clinicopathological factors improves response prediction to systemic treatment in patients with advanced colorectal cancer. Submitted.

Neerincx, M., Sie, D.L.S., van de Wiel, M.A., van Grieken, N.C.T., Burggraaf, J.D., Dekker, H., Eijk, P.P., Ylstra, B., Verhoef, C., Meijer, G.A., Buffart, T.E., and others. (2015). MiR expression profiles of paired primary colorectal cancer and metastases by next-generation sequencing. Oncogenesis, 4, e170.
}
\references{
Münch, M.M., Peeters, C.F.W., van der Vaart, A.W., and van de Wiel, M.A. (2018). Adaptive group-regularized logistic elastic net regression. arXiv:1805.00389v1 [stat.ME].
}
\examples{
data(dataColon)
str(mirCol)
str(unpenCol)
str(respCol)
str(mirExpr)
}
\keyword{datasets}

 
