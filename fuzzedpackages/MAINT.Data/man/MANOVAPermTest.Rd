\name{MANOVAPermTest}
\alias{MANOVAPermTest}

\title{MANOVA permutation test}
\description{Function MANOVAPermTest performs a MANOVA permutation test allowing for both Gaussian and Skew-Normal distributions and homoscedastic or heteroscedastic setups.}

\usage{

MANOVAPermTest(MANOVAres, Idt, grouping, nrep=200,
    Model=c("Normal","SKNormal","NrmandSKN"), CovCase=1:4,
    SelCrit=c("BIC","AIC"), Mxt=c("Hom","Het","Loc","Gen"), CVtol=1.0e-5, k2max=1e6,
    OptCntrl=list(), onerror=c("stop","warning","silentNull"), \dots)

}

\details{Function MANOVAPermTest performs a MANOVA permutation test allowing for both Gaussian and Skew-Normal distributions and homoscedastic or heteroscedastic setups. This test is implemented by simulating the null distribution of the MANOVA likelihood ratio statistic, using many random permutations of the observation group labels. It is intended as an alternative of the classical Chi-squares based MANOVA likelihood ratio tests, when small sample sizes cast doubt on the applicability of the Chi-squared distribution. We note that this test may be computationally intensive, in particular when used for the Skw-Normal model.}


\arguments{
  \item{MANOVAres}{An object representing a MANOVA analysis on interval-valued entities.}

  \item{Idt}{An IData object representing interval-valued entities.}

  \item{grouping}{Factor indicating the group to which each observation belongs to.}

  \item{nrep}{Number of random generated permutations used to approximate the null distribution of the likelihood ratio statistic.}

  \item{Model}{The joint distribution assumed for the MidPoint and LogRanges. Current alternatives are \dQuote{Normal} for Gaussian,
distributions, \dQuote{SKNormal} for Skew-Normal and \dQuote{NrmandSKN} for both Gaussian and Skew-Normal distributions.}

  \item{CovCase}{Configuration of the variance-covariance matrix: a set of integers between 1 and 4.}

   \item{SelCrit}{The model selection criterion.}

   \item{Mxt}{Indicates the type of mixing distributions to be considered. Current alternatives are \dQuote{Hom} (homocedastic) and \dQuote{Het} (heteroscedastic) for Gaussian models, \dQuote{Loc} (location model -- groups differ only on their location parameters) and \dQuote{Gen} \dQuote{Loc} (general model -- groups differ on all parameters) for Skew-Normal models.}

  \item{CVtol}{Tolerance level for absolute value of the coefficient of variation of non-constant variables. When a MidPoint or LogRange has an absolute value within-groups coefficient of variation below CVtol, it is considered to be a constant.}

  \item{k2max}{Maximal allowed l2-norm condition number for correlation matrices. Correlation matrices with condition number above k2max are considered to be numerically singular, leading to degenerate results.}

   \item{OptCntrl}{List of optional control parameters to passed to the optimization routine. See the documentation of RepLOptim for a description of the available options.}

   \item{onerror}{Indicates whether an error in the optimization algorithm should stop the current call, generate a warning, or return silently a NULL object.}

  \item{\dots}{Other named arguments.}
}

\value{the p-value of the MANOVA permutation test.}


\keyword{MANOVA}
\keyword{interval data}
\keyword{likelihood ratio tets}
\keyword{permutation tets}

\examples{

\dontrun{

#Perform a MANOVA of the AbaloneIdt data set, comparing the Abalone variable means 
# according to their age 

# Create an Interval-Data object containing the Length, Diameter, Height, Whole weight, 
# Shucked weight, Viscera weight (VW), and Shell weight (SeW) of 4177 Abalones, 
# aggregated by sex and age.
# Note: The original micro-data (imported UCI Machine Learning Repository Abalone dataset) 
# is given in the AbaDF data frame, and the corresponding values of the sex by age combinations 
# is represented by the AbUnits factor. 

AbaloneIdt <- AgrMcDt(AbaDF,AbUnits)

# Create a factor with three levels (Young, Adult and Old) for Abalones with respectively 
# less than 10 rings, between 11 and 18 rings, and more than 18 rings. 

Agestrg <- substring(rownames(AbaloneIdt),first=3)
AbalClass <- factor(ifelse(Agestrg=="1-3"|Agestrg=="4-6"| Agestrg=="7-9","Young",
  ifelse(Agestrg=="10-12"|Agestrg=="13-15"| Agestrg=="16-18","Adult","Old") ) )

#Perform a classical MANOVA, computing the p-value from the asymptotic Chi-squared distribution 
# of the Wilk's lambda statistic

MANOVAres <- MANOVA(AbaloneIdt,AbalClass)
summary(MANOVAres)

#Find a finite sample p-value of the test statistic, using a permutation test.

MANOVAPermTest(MANOVAres,AbaloneIdt,AbalClass)

}

}


\seealso{ \code{\link{MANOVA}}, \code{\linkS4class{IdtMANOVA}} }

