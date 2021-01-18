\name{wages}
\alias{wages}
\docType{data}
\title{Wages and Hours}
\description{
The data are from a national sample of 6000 households with a male head earning less than 
USD 15,000 annually in 1966. The data were clasified into 39 demographic groups for 
analysis. The study was undertaken in the context of proposals for a guaranteed 
annual wage (negative income tax). At issue was the response of labor supply 
(average hours) to increasing hourly wages. 
The study was undertaken to estimate this response from available data.
}
\usage{data(wages)}
\format{
    A data frame with 39 observations on the following 10 variables: 
 
    \tabular{ll}{
        \code{HRS} \tab Average hours worked during the year \cr
        \code{RATE} \tab Average hourly wage (USD)\cr
        \code{ERSP} \tab Average yearly earnings of spouse (USD)\cr
        \code{ERNO} \tab Average yearly earnings of other family members (USD)\cr
        \code{NEIN} \tab Average yearly non-earned income\cr
        \code{ASSET} \tab Average family asset holdings (Bank account, etc.) (USD)\cr
        \code{AGE} \tab Average age of respondent\cr
        \code{DEP} \tab Average number of dependents\cr
        \code{RACE} \tab Percent of white respondents\cr
        \code{SCHOOL} \tab Average highest grade of school completed\cr
    }
}
\references{
    D.H. Greenberg and M. Kosters, (1970). Income Guarantees and the Working Poor, The Rand Corporation. 
}
\source{
    DASL library (http://lib.stat.cmu.edu/DASL/Datafiles/wagesdat.html),
        the dataset is not anymore available at this source.
}
\keyword{datasets}
