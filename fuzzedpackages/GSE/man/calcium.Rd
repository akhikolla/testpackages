\name{calcium}
\alias{calcium}
\docType{data}
\title{Calcium data}
\description{
The Calcium data is from the article by Holcomb and Spalsbury (2005). The dataset used for
class was compiled by Boyd, Delost, and Holcomb (1998) for the use of a study to determine if
significant gender differences existed between subjects 65 years of age and older with regard to
calcium, inorganic phosphorous, and alkaline phosphatase levels. Although the original data from
Boyd, Delost, and Holcomb (1998) had observations needing investigation, Holcomb and Spalsbury
(2005) further massaged the original data to include data problems and issues that have arisen in
other research projects for pedagogical purposes.
}
\usage{data(calcium)}
\format{A data frame with 178 observations on the following 8 variables.
	\tabular{ll}{
		\code{obsno} \tab Patient Observation Number \cr
		\code{age} \tab Age in years \cr
		\code{sex} \tab 1=Male, 2=Female \cr
		\code{alkphos} \tab Alkaline Phosphatase International Units/Liter \cr
		\code{lab} \tab 1=Metpath; 2=Deyor; 3=St. Elizabeth's; 4=CB Rouche; 5=YOH; 6=Horizon \cr
		\code{cammol} \tab Calcium mmol/L \cr
		\code{phosmmol} \tab Inorganic Phosphorus mmol/L \cr
		\code{agegroup} \tab Age group 1=65-69; 2=70-74; 3=75-79; 4=80-84; 5=85-89 Years \cr
	}
}
\references{
Boyd, J., Delost, M., and Holcomb, J., (1998). Calcium, phosphorus, 
and alkaline phosphatase laboratory values of elderly subjects, 
\emph{Clinical Laboratory Science}, 11, 223-227.

Holcomb, J., and Spalsbury, A. (2005), Teaching Students to Use Summary Statistics and Graphics
to Clean and Analyze Data. \emph{Journal of Statistics Education}, 13, Number 3.
}
\source{
  The original data have been taken from the Journal of Statistics Education
  Databases at
  \itemize{
    \item \url{http://www.amstat.org/publications/jse/datasets/calcium.dat.txt},
  }
  the corrected data have been taken from Statlib at
  \itemize{
    \item \url{http://www.amstat.org/publications/jse/datasets/calciumgood.dat.txt}
  }
}
\examples{
\dontrun{
data(calcium)
## remove the categorical variables
calcium.cts <- subset(calcium, select=-c(obsno, sex, lab, agegroup) )
res <- GSE(calcium.cts)
getOutliers(res)
## able to identify majority of the contaminated cases identified 
## in the reference
}
}
\keyword{datasets}
