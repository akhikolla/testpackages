\name{docvisits}
\docType{data}
\alias{docvisits}
\title{Demand for Health Care Data}
\description{
  This data set gives the number of doctor visits in the last three
  months for a sample of German male individuals in 1994. The data set
  is taken from Riphahn et al. (2003) and is a subsample of the German
  Socioeconomic Panel (SOEP). In contrast to Riphahn et al. (2003) only
  male individuals from the last wave are considered. See Jochmann
  (2013) for further details.
}
\usage{data(docvisits)}
\format{
  This data frame contains 1812 observations on the following 22 variables:
  \describe{
    \item{docvisits}{number of doctor visits in last 3 months}         
    \item{age}{age}                                              
    \item{agesq}{age squared / 1000}
    \item{age30}{1 if age >= 30}
    \item{age35}{1 if age >= 35}
    \item{age40}{1 if age >= 40}
    \item{age45}{1 if age >= 45}
    \item{age50}{1 if age >= 50}
    \item{age55}{1 if age >= 55}
    \item{age60}{1 if age >= 60}
    \item{health}{health satisfaction, 0 (low) - 10 (high)}         
    \item{handicap}{1 if handicapped, 0 otherwise}                    
    \item{hdegree}{degree of handicap in percentage points}          
    \item{married}{1 if married, 0 otherwise}                        
    \item{schooling}{years of schooling}                               
    \item{hhincome}{household monthly net income, in German marks / 1000}                                 
    \item{children}{1 if children under 16 in the household, 0 otherwise}                                            
    \item{self}{1 if self employed, 0 otherwise}                  
    \item{civil}{1 if civil servant, 0 otherwise}                  
    \item{bluec}{1 if blue collar employee, 0 otherwise}           
    \item{employed}{1 if employed, 0 otherwise}                       
    \item{public}{1 if public health insurance, 0 otherwise}        
    \item{addon}{1 if add-on insurance, 0 otherwise}               
  }
}

\references{
  Jochmann, M. (2013). ``What Belongs Where? Variable Selection for Zero-Inflated Count Models with an Application to the Demand for
  Health Care'', \emph{Computational Statistics}, 28, 1947--1964.
  
  Riphahn, R. T., Wambach, A., Million, A. (2003). 
  ``Incentive Effects in the Demand for Health Care: A Bivariate Panel Count Data Estimation'',
  \emph{Journal of Applied Econometrics}, 18, 387--405.
  
  Wagner, G. G., Frick, J. R., Schupp, J. (2007).
  ``The German Socio-Economic Panel Study (SOEP) -- Scope, Evolution and Enhancements'',
  \emph{Schmollers Jahrbuch}, 127, 139--169.
}
\keyword{datasets}
