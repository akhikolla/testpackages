\docType{data}

\name{teenpov}

\alias{teenpov}

\title{Data from the National Longitudinal Survey of Youth (NLSY).}

\format{The dataset contains the following variables:
\describe{
\item{ID}{a unique subject-identifier.} 
\item{nonpov}{1 if the girl is currenlty not in poverty according to U.S. federal standards.}
\item{hours}{the number of hours that the girl currently works per week.}
\item{inschool}{1 if the girl is currently enrolled in school, 0 otherwise.}
\item{spouse}{1 if the girl is currently living with a spouse, 0 otherwise.}
\item{age}{the girl's current age.}
\item{mother}{1 if the girl currently has at least 1 child, 0 otherwise.}
}
}
\usage{
data(teenpov)
}
\description{
This dataset is a subset of a dataset described by Allison (2009).
The original data by Allison (2009) comes from the National Longitudinal Survey of Youth (NLSY),
and contains information on 1151 teenage girls who were interviewed annually for 
five years, starting in 1979. The teenpov dataset is constructed by first identifying, 
for each girl, the first interview at which the girl is not in poverty,
then excluding all subsequent interviews for that girl. After this exclusion,
the dataset contains 1342 interviews from 401 girls. This procedure is described in greater
detail by Sjolander (2017). 

}
\references{

Allison P. (2009). Fixed effects regression models: quantitative applications in 
the Social sciences. Thousand Oaks, CA: Sage.

Sjolander A. (2017). The case-time-control method for non-binary exposures. 
\emph{Sociological Methodology} \bold{47}(1), 182-211.
}
