\name{plot-methods}
\docType{methods}
\alias{plot}
\alias{plot-methods}
\alias{plot,IData-method}
\alias{plot,IData,IData-method}
\alias{plot,IData,missing-method}

\title{Methods for function plot in Package `MAINT.Data'}
\description{S4 methods for function plot. As in the generic plot S3 `graphics' method, these methods plot Interval-valued data contained in \linkS4class{IData} objects.}

\keyword{methods}
\keyword{plot}
\keyword{Interval Data}

\usage{

\S4method{plot}{IData,IData}(x, y, type=c("crosses","rectangles"), append=FALSE, \dots)
\S4method{plot}{IData,missing}(x, casen=NULL, layout=c("vertical","horizontal"), append=FALSE, \dots)

}

\arguments{
  \item{x}{An object of type \linkS4class{IData} representing the values of an Interval-value variable.}

  \item{y}{An object of type \linkS4class{IData} representing the values of a second Interval-value variable, to be displayed along y (vertical) coordinates.}

    \item{type}{What type of plot should de drawn. Alternatives are "crosses" (default) and "rectangles".}

    \item{append}{A boolean flag indicating if the interval-valued variables should be displayed in a new plot, or added to an existing plot.}

    \item{casen}{An optional character string with the case names.}

    \item{layout}{The axes along which the interval-valued variables be displayed. Alternatives are "vertical" (default) and "horizontal".}
    
  \item{\dots}{ Graphical arguments to be passed to methods.}
}

\seealso{
  \linkS4class{IData} 
}

\examples{

\dontrun{

# Create an Interval-Data object containing the Length, Diameter, Height, Whole weight, 
# Shucked weight, Viscera weight (VW), and Shell weight (SeW) of 4177 Abalones, 
# aggregated by sex and age.
# Note: The original micro-data (imported UCI Machine Learning Repository Abalone dataset) 
# is given in the AbaDF data frame, and the corresponding values of the sex by age combinations 
# is represented by the AbUnits factor. 

AbaloneIdt <- AgrMcDt(AbaDF,AbUnits)

# Dispaly a plot of the Length versus the Whole_weight interval variables

plot(AbaloneIdt[,"Length"],AbaloneIdt[,"Whole_weight"])
plot(AbaloneIdt[,"Length"],AbaloneIdt[,"Whole_weight"],type="rectangles")

# Display the Abalone lengths using different colors to distinguish the Abalones age 
# (measured by the number of rings)   

# Create a factor with three levels (Young, Adult and Old) for Abalones with 
# respectively less than 10 rings, between 11 and 18 rings, and more than 18 rings. 

Agestrg <- substring(rownames(AbaloneIdt),first=3)
AbalClass <- factor(ifelse(Agestrg=="1-3"|Agestrg=="4-6"| Agestrg=="7-9","Young",
  ifelse(Agestrg=="10-12"|Agestrg=="13-15"| Agestrg=="16-18","Adult","Old") ) )

plot(AbaloneIdt[AbalClass=="Young","Length"],col="blue",layout="horizontal") 
plot(AbaloneIdt[AbalClass=="Adult","Length"],col="green",layout="horizontal",append=TRUE) 
plot(AbaloneIdt[AbalClass=="Old","Length"],col="red",layout="horizontal",append=TRUE) 
legend("bottomleft",legend=c("Young","Adult","Old"),col=c("blue","green","red"),lty=1)


}

}



