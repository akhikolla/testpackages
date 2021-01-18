\name{MANOVA-methods}
\docType{methods}
\alias{MANOVA}
\alias{MANOVA-methods}
\alias{MANOVA,IData-method}
\alias{H0res}
\alias{H0res,IdtMANOVA-method}
\alias{H1res}
\alias{H1res,IdtMANOVA-method}
\alias{show,IdtMANOVA-method}

\title{Methods for Function MANOVA in Package `MAINT.Data'}
\description{Function MANOVA performs MANOVA tests based on likelihood ratios allowing for both Gaussian and Skew-Normal distributions and homoscedastic or heteroscedastic setups. Methods H0res and H1res retrieve the model estimates under the null and alternative hypothesis, and method show displays the  MANOVA results.}

\usage{

MANOVA(Idt, grouping, Model=c("Normal","SKNormal","NrmandSKN"), CovCase=1:4, 
		SelCrit=c("BIC","AIC"), Mxt=c("Hom","Het","Loc","Gen"), 
                CVtol=1.0e-5, k2max=1e6,
		OptCntrl=list(), onerror=c("stop","warning","silentNull"), \dots)

\S4method{H0res}{IdtMANOVA}(object)
\S4method{H1res}{IdtMANOVA}(object)
\S4method{show}{IdtMANOVA}(object)

}

\arguments{
  \item{object}{An object representing a MANOVA analysis on interval-valued entities.}

  \item{Idt}{An IData object representing interval-valued entities.}

  \item{grouping}{Factor indicating the group to which each observation belongs to.}

  \item{Model}{The joint distribution assumed for the MidPoint and LogRanges. Current alternatives are \dQuote{Normal} for Gaussian,
distributions, \dQuote{SKNormal} for Skew-Normal and \dQuote{NrmandSKN} for both Gaussian and Skew-Normal distributions.}

  \item{CovCase}{Configuration of the variance-covariance matrix: a set of integers between 1 and 4.}

   \item{SelCrit}{The model selection criterion.}

	\item{Mxt}{Indicates the type of mixing distributions to be considered. Current alternatives are \dQuote{Hom} (homoscedastic) and \dQuote{Het} (heteroscedastic) for Gaussian models, \dQuote{Loc} (location model -- groups differ only on their location parameters) and \dQuote{Gen} \dQuote{Loc} (general model -- groups differ on all parameters) for Skew-Normal models.}

  \item{CVtol}{Tolerance level for absolute value of the coefficient of variation of non-constant variables. When a MidPoint or LogRange has an absolute value within-groups coefficient of variation below CVtol, it is considered to be a constant.}

  \item{k2max}{Maximal allowed l2-norm condition number for correlation matrices. Correlation matrices with condition number above k2max are considered to be numerically singular, leading to degenerate results.}

  \item{OptCntrl}{List of optional control parameters to passed to the optimization routine. See the documentation of RepLOptim for a description of the available options.}

   \item{onerror}{Indicates whether an error in the optimization algorithm should stop the current call, generate a warning, or return silently a NULL object.}

  \item{\dots}{Other named arguments.}
}

\value{An object of class IdtMANOVA, containing the estimation and test results.}


\keyword{methods}
\keyword{MANOVA}
\keyword{interval data}
\keyword{likelihood ratio tets}


\seealso{ \code{\linkS4class{IdtMANOVA}}, \code{\link{RepLOptim}} }

\examples{

#Create an Interval-Data object containing the intervals of temperatures by quarter 
# for 899 Chinese meteorological stations.
ChinaT <- IData(ChinaTemp[1:8])

#Classical (homoscedastic) MANOVA tests

ManvChina <- MANOVA(ChinaT,ChinaTemp$GeoReg)
cat("China, MANOVA by geografical regions results =\n")
print(ManvChina)

#Heteroscedastic MANOVA tests

HetManvChina <- MANOVA(ChinaT,ChinaTemp$GeoReg,Mxt="Het")
cat("China, heterocedastic MANOVA by geografical regions results =\n")
print(HetManvChina)

#Skew-Normal based MANOVA assuming the the groups differ only according to location parameters
\dontrun{

SKNLocManvChina <- MANOVA(ChinaT,ChinaTemp$GeoReg,Model="SKNormal",Mxt="Loc")
cat("China, Skew-Normal MANOVA (location model) by geografical regions results =\n")
print(SKNLocManvChina)

#Skew-Normal based MANOVA assuming the the groups may differ in all parameters

SKNGenManvChina <- MANOVA(ChinaT,ChinaTemp$GeoReg,Model="SKNormal",Mxt="Gen")
cat("China, Skew-Normal MANOVA (general model)  by geografical regions results =\n")
print(SKNGenManvChina)

}

}

