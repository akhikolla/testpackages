\name{testMod-methods}
\docType{methods}
\alias{testMod-methods}
\alias{testMod}   
\title{ Methods for Function testMod in Package `MAINT.Data'}

\description{Performs statistical likelihood-ratio tests that evaluate the goodness-of-fit of a nested model against a more general one. }


\usage{testMod(ModE,RestMod=ModE@ModelConfig[2]:length(ModE@ModelConfig),FullMod="Next")}

\arguments{
  \item{ModE}{An object of class \code{\linkS4class{IdtE}} representing the estimates of a model fitted to a data set of interval-value variables}
  \item{RestMod}{Indices of the restricted models being evaluated in the NULL hypothesis}
  \item{FullMod}{Either indices of the general models being evaluated in the alternative hypothesis or the strings "Next" (default) or "All". In the former case a Restricted model is always compared against the most parsimonious alternative that encompasses it, and in latter all possible comparisons are performed}

}

\value{An object of class ConfTests with the results of the tests performed}


\keyword{methods}
\keyword{interval data}
\keyword{likelihood ratio test}

\examples{

# Create an Interval-Data object containing the intervals of temperatures by quarter 
# for 899 Chinese meteorological stations.

ChinaT <- IData(ChinaTemp[1:8])

# Estimate by maximum likelihood the parameters of Gaussian models 
# for the Winter (1st and 4th) quarter intervals

ChinaWTE <- mle(ChinaT[,c(1,4)])
cat("China maximum likelhiood estimation results for Winter quarters:\n")
print(ChinaWTE)

# Perform Likelihood-Ratio tests comparing models with consecutive nested Configuration 
testMod(ChinaWTE)

# Perform Likelihood-Ratio tests comparing all possible models 
testMod(ChinaWTE,FullMod="All")

# Compare model with covariance Configuration case 3 (MidPoints independent of LogRanges) 
# against model with covariance Configuration 1 (unrestricted covariance)  
testMod(ChinaWTE,RestMod=3,FullMod=1)
}



