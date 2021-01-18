\name{dmi}
\alias{dmi}

\title{Data monotonicity index for missing values}

\description{\sQuote{dmi} calculates a monotonicity index for data with missing values.}

\usage{dmi(Data)}

\arguments{
  \item{Data}{A data frame containing missing values.}
}

\details{The \dfn{data monotonicity index} examines the ratio of missing values with non-monotonicity and complete monotonicity in all variables. To denote full monotonicity with 1 and no monotonicity with 0 this ratio is subtracted from 1.

\deqn{dmi = 1-\frac{\sum_{j=1}^{p}\sum_{i=1}^{n-\sum_{h=1}^{n}I(r_{hj}=0)} \sum_{h=1}^{n} I (r_{hi}=0)}{\sum_{h=1}^{n}\sum_{j=1}^{p} I(r_{hj}=0)}}{dmi = 1 - (sum_{j=1}^{p} sum_{i=1}^{n-sum_h=1^n I(r_hi==0)} sum_{h=1}^{n} I(r_{hj}==0))/(sum_{h=1}^{n} sum_{j=1}^{p} I(r_{hj} == 0))}
}

\value{
  Returns a value between 1 (fully monotone) and 0 (no monotonicity).
}

\author{Florian Meinfelder, Thorsten Schnapp}

\seealso{\code{\link{BBPMM}}, \code{\link[norm]{prelim.norm}}}

\references{
  Harrell, F.E., with contributions from Charles Dupont and many others. (2013) \emph{Hmisc: Harrell Miscellaneous}. R package version 3.13-0. http://CRAN.R-project.org/package=Hmisc
  
  Venables, W. N. & Ripley, B. D. (2002) \emph{Modern Applied Statistics with S}. Fourth Edition. Springer, 
  New York. ISBN 0-387-95457-0
  
  Ported to R by Alvaro A. Novo. Original by Joseph L. Schafer <jls@stat.psu.edu>. (2013). \emph{norm:
  Analysis of multivariate normal datasets with missing values}. R package version 1.0-9.5.
  http://CRAN.R-project.org/package=norm
  
}

\keyword{monotone data pattern}

\examples{

if(!require(MASS)) install.packages("MASS")
library(MASS)  ## see references
data(survey)

## Sorting via 'norm's prelim.norm
if(!require(Hmisc)) install.packages("Hmisc")
library(Hmisc) ## see references
survey.numeric <- asNumericMatrix(survey)

if(!require(norm)) install.packages("norm")
library(norm) ## see references
su.sort    <- prelim.norm(survey.numeric)
new.survey <- survey[order(su.sort$ro),
                     sort(su.sort$nmis,index.return=TRUE)$ix]

## Comparison
dmi(survey)     # original
dmi(new.survey) # sorted


}