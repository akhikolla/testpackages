\name{plotInfCrt-methods}
\docType{methods}
\alias{plotInfCrt}
\alias{plotInfCrt-methods}
\alias{plotInfCrt,IdtMclust-method}

\title{Information criteria plot.}
\description{Method plotInfCrt displays a plot representing the values of an appropriate information criterion (currently either BIC or AIC) for the models whose results are stored in an IdtMclust-method object. A supplementary short output message prints the values of the chosen criterion for the 'nprin' best models.}

\keyword{methods}
\keyword{BIC}
\keyword{AIC}
\keyword{Interval Data}
\keyword{model based clustering}

\usage{

\S4method{plotInfCrt}{IdtMclust}(object,crt=object@SelCrit,legpos="bottomleft",nprnt=5, \dots)

}

\arguments{
  \item{object}{An object of type \dQuote{IdtMclust} representing the the clusterig results of an Interval-valued data set obtained by the function \dQuote{IdtMclust}.}

  \item{crt}{The information criteria whose values are to be displayed.}

  \item{legpos}{Legend position. Alternatives are \dQuote{bottomleft} (default), \dQuote{bottomright}, \dQuote{topleft} and \dQuote{bottomright}.} 

  \item{nprnt}{Number of solutions for which the value of the information criterio should be printed in an suplmentary short output message.}

  \item{\dots}{ Graphical arguments to be passed to methods.}

}

\seealso{
  \linkS4class{IdtMclust}, \code{\link{Idtmclust}},  \code{\link{pcoordplot}}
}

\examples{

\dontrun{

# Create an Interval-Data object containing the intervals of loan data
# (from the Kaggle Data Science platform) aggregated by loan purpose

LbyPIdt <- IData(LoansbyPurpose_minmaxDt,
                 VarNames=c("ln-inc","ln-revolbal","open-acc","total-acc")) 

#Fit homoscedastic and heteroscedastic mixtures up to Gaussian mixtures with up to seven components

mclustres <- Idtmclust(LbyPIdt,G=1:7,Mxt="HomandHet")

#Compare de model fit according to the BIC

plotInfCrt(mclustres,legpos="bottomleft")

#Display the results of the best three mixtures according to the BIC

summary(mclustres,parameters=TRUE,classification=TRUE)
pcoordplot(mclustres)
summary(mclustres,parameters=TRUE,classification=TRUE,model="HetG2C2")
summary(mclustres,parameters=TRUE,classification=TRUE,model="HomG6C1")
pcoordplot(mclustres,model="HomG6C1")


}

}


