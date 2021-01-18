\name{getIdtOutl}
\alias{getIdtOutl}

\title{Get Interval Data Outliers}

\description{Identifies outliers in a data set of Interval-valued variables}

\usage{getIdtOutl(Idt, IdtE=NULL, muE=NULL, SigE=NULL,
  eta=0.025, Rewind=NULL, m=length(Rewind),
  RefDist=c("ChiSq","HardRockeAdjF","HardRockeAsF","CerioliBetaF"),
  multiCmpCor=c("never","always","iterstep"), 
  outlin=c("MidPandLogR","MidP","LogR"))
}


\arguments{
    \item{Idt}{An IData object representing interval-valued entities.}

    \item{IdtE}{Ao object of class \code{\linkS4class{IdtSngNDRE}} or \code{\linkS4class{IdtSngNDE}} containing mean and covariance estimates.}

    \item{muE}{Vector with the mean estimates used to find Mahalanobis distances. When specified, it overrides the mean estimate supplied in \dQuote{IdtE}.}

    \item{SigE}{Matrix with the covariance estimates used to find Mahalanobis distances. When specified, it overrides the covariance estimate supplied in \dQuote{IdtE}.}

    \item{eta}{Nominal size of the null hypothesis that a given observation is not an outlier.}

    \item{Rewind}{A vector with the subset of entities used to compute trimmed mean and covariance estimates when using a reweighted MCD. Only used when the \sQuote{RefDist} argument is set to \dQuote{CerioliBetaF.}}

    \item{m}{Number of entities used to compute trimmed mean and covariance estimates when using a reweighted MCD. Not used when the \sQuote{RefDist} argument is set to \dQuote{ChiSq.}}

    \item{multiCmpCor}{Whether a multicomparison correction of the nominal size (eta) for the outliers tests should be performed. Alternatives are: \sQuote{never} -- ignoring the multicomparisons and testing all entities at the \sQuote{eta} nominal level. \sQuote{always} -- testing all n entitites at 1.- (1.-\sQuote{eta}^(1/n)); and \sQuote{iterstep} -- use the iterated rule proposed by  Cerioli (2010), \emph{i.e.}, make an initial set of tests using the nominal size 1.- (1-\sQuote{eta}^(1/n)), and if no outliers are detected stop. Otherwise, make a second step testing for outliers at the \sQuote{eta} nominal level.}   

    \item{RefDist}{The assumed reference distributions used to find cutoffs defining the observations assumed as outliers. Alternatives are \dQuote{ChiSq},\dQuote{HardRockeAsF}, \dQuote{HardRockeAdjF} and \dQuote{CerioliBetaF}, respectivelly for the usual Chi-squared, the asymptotic and adjusted scaled F distributions proposed by Hardin and Rocke (2005), and the Beta and F distributions proposed by Cerioli (2010).}    

    \item{outlin}{The type of outliers to be considered. \dQuote{MidPandLogR} if outliers may be present in both MidPpoints and LogRanges, \dQuote{MidP} if outliers are only present in MidPpoints, or \dQuote{LogR} if outliers are only present in LogRanges.}
}

\value{A vector with the indices of the entities identified as outliers.} 

\references{
Cerioli, A. (2010), Multivariate Outlier Detection with High-Breakdown Estimators.
\emph{Journal of the American Statistical Association} \bold{105} (489), 147--156.\cr 

Duarte Silva, A.P., Filzmoser, P. and Brito, P. (2017), Outlier detection in interval data. \emph{Advances in Data Analysis and Classification}, 1--38.\cr

Hardin, J. and Rocke, A. (2005), The Distribution of Robust Distances.
\emph{Journal of Computational and Graphical Statistics} \bold{14}, 910--927.
}

\seealso{ \code{\link{fasttle}}, \code{\link{fulltle}} }

\examples{

\dontrun{

# Create an Interval-Data object containing the intervals for characteristics 
# of 27 cars models.

CarsIdt <- IData(Cars[1:8],VarNames=c("Price","EngineCapacity","TopSpeed","Acceleration"))

# Estimate parameters by the fast trimmed maximum likelihood estimator, 
# using a two-step procedure to select the trimming parameter, a reweighed 
# MCD estimate, and the classical 97.5\% chi-squared quantile cut-offs.
			
Carstle1 <- fulltle(CarsIdt)
			
# Get and display the outliers using the classical 97.5\% chi-squared quantile cut-offs.
		
CarsOtl1 <- getIdtOutl(CarsIdt,Carstle1)
print(CarsOtl1)
plot(CarsOtl1)
			
# Estimate parameters by the fast trimmed maximum likelihood estimator, 
# using a two-step procedure to select the trimming parameter, and a reweighed  
# based on the 97.5\% quantiles of Hardin and Rocke adjusted F distributions.
			
Carstle2 <- fulltle(CarsIdt,rawMD2Dist="HardRockeAdjF")
			
# Get and display the outliers using the 97.5% qunatiles of Cerioli Beta and F distributions.
			
CarsTtl2 <- getIdtOutl(CarsIdt,Carstle2,RefDist="CerioliBetaF")
print(CarsTtl2)
plot(CarsTtl2)


}

}
