\name{burn}
\alias{burn}
\docType{data}
\title{
burn dataset
}
\description{
The burn data frame has 154 rows and 17 columns.
}
\usage{data(burn)}
\format{
  A data frame with 154 observations on the following 17 variables.
  \describe{
    \item{\code{Obs}}{Observation number}
    \item{\code{Z1}}{Treatment: 0-routine bathing 1-Body cleansing}
    \item{\code{Z2}}{Gender (0=male 1=female)}
    \item{\code{Z3}}{Race: 0=nonwhite 1=white}
    \item{\code{Z4}}{Percentage of total surface area burned}
    \item{\code{Z5}}{Burn site indicator: head 1=yes, 0=no}
    \item{\code{Z6}}{Burn site indicator: buttock 1=yes, 0=no}
    \item{\code{Z7}}{Burn site indicator: trunk 1=yes, 0=no}
    \item{\code{Z8}}{Burn site indicator: upper leg 1=yes, 0=no}
    \item{\code{Z9}}{Burn site indicator: lower leg 1=yes, 0=no}
    \item{\code{Z10}}{Burn site indicator: respiratory tract 1=yes, 0=no}
    \item{\code{Z11}}{Type of burn: 1=chemical, 2=scald, 3=electric, 4=flame}
    \item{\code{T1}}{Time to excision or on study time}
    \item{\code{D1}}{Excision indicator: 1=yes 0=no}
    \item{\code{T2}}{Time to prophylactic antibiotic treatment or on study time}
    \item{\code{D2}}{Prophylactic antibiotic treatment: 1=yes 0=no}
    \item{\code{T3}}{Time to straphylocous aureaus infection or on study time}
    \item{\code{D3}}{Straphylocous aureaus infection: 1=yes 0=no}     
  }
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
Klein and Moeschberger (1997) Survival Analysis Techniques for Censored and truncated data, \code{Springer}. 

Ichida et al. \code{Stat. Med.} 12 (1993): 301-310.
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(burn)
## maybe str(burn) ;
}
\keyword{datasets}
