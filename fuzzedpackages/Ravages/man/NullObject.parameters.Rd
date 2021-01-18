\name{NullObject.parameters}
\alias{NullObject.parameters}
\title{ Null Model for SKAT and burden tests }
\description{ Get the parameters under the null model to peforms burden tests or SKAT}
\usage{
  NullObject.parameters(pheno, RVAT, pheno.type, ref.level, data, formula)
}
\arguments{
  \item{pheno}{ The phenotype of each individual: a factor if \code{pheno.type = "categorial"}, and a numeric vector if \code{pheno.type = "continuous"} }
  \item{RVAT}{ The type of Rare Variant Association Test (RVAT) to perform: should be "burden" or "SKAT"}
  \item{pheno.type}{ The type of phenotype: "categorial" for binary or multinomial traits, or "continuous"}
  \item{ref.level}{ The reference group of individuals for the regression, only needed if \code{RVAT = "burden"} and \code{pheno.type = "categorial"}}
  \item{data}{ Optional, a matrix containing the covariates with one column per covariate and one row per individual }
  \item{formula}{ Optional, an R formula corresponding to the regression model indicating which covariates from \code{data} to include in the model if only some of them are to be included}
}
\details{
  This function gets the parameters under the null model for SKAT or the burden tests.  

  For burden tests, it computes the Log-Likelihood under the null model used to perform the Likelihood Ratio Test.   

  For SKAT, it computes the probabilites for each individual of belonging to each group based on the group sizes and the potential covariates.  

  If \code{formula} is missing, all columns from \code{data} will be included as covariates.
}
\value{
  A list containing different elements depending on the \code{RVAT} performed and the \code{pheno.type}.  

  - if \code{RVAT = "burden"} and \code{pheno.type = "categorial"}:
  \item{group}{A factor containing the group of each individual as given}
  \item{ref.level}{ The reference group of individuals for the regression as given}
  \item{H0.LogLik}{ The Log-Likelihood of the null model}
  \item{covar.toinclude}{ Which covariates to include in the regression, depends on the argument \code{formula} }
  \item{data}{ The \code{data} argument containing covariates, NULL if it was missing}   


  - if \code{RVAT = "burden"} and \code{pheno.type = "continuous"}:
  \item{pheno}{A numeric vector containing the phenotype value for each individual as given}
  \item{covar.toinclude}{ Which covariates to include in the regression, depends on the argument \code{formula} }
  \item{data}{ The \code{data} argument containing covariates, NULL if it was missing}   


  - if \code{RVAT = "SKAT"} and \code{pheno.type = "categorial"}:  
  \item{Pi.data}{ A matrix n.individuals x n.groups containing the probabilities that each individual belong to each group }
  \item{X}{ A matrix containing 1 in the first column for the intercept, and covariates from \code{data} and \code{formula}}
  \item{group}{ A factor containing the group of each individual as given}
  \item{get.moments}{ How to compute moments based on sample size for p-value calculations (only used if \code{get.moments = "size.based"} for a categorial phenotype in \code{SKAT}. }
  \item{P1}{ The vairance-covariance matrix of (Y - Pi_hat)}   

 
  - if \code{RVAT = "SKAT"} and \code{pheno.type = "continuous"}:
  \item{ymp}{ A matrix n.individuals x 1 containing the (y - pi_hat) values, i.e. the residuals from the regression of the phenotype on the potential covariates }
  \item{X}{ A matrix containing 1 in the first column for the intercept, and covariates from \code{data} and \code{formula}}
  \item{pheno}{ The phenotype of each individual as given}
  \item{P1}{ The variance matrix of \code{ymp} }

}
\seealso{ \code{\link{SKAT}}, \code{\link{burden}} }
\examples{
#Random phenotype of 100 individuals 
random.multi.pheno <- sample(1:3, 100, replace = TRUE)
#Random continuous phenotype
random.continuous.pheno <- rnorm(100)
#Random sex covariate
random.covar <- matrix( sample(1:2, prob = c(0.4, 0.6), size = 100, replace = TRUE),
                        ncol = 1 )

#Null Model for burden with a multinomi-category phenotype
#Controls as reference group, no covariates
H0.burden.multi <- NullObject.parameters(pheno = as.factor(random.multi.pheno), 
                                         RVAT = "burden", pheno.type = "categorial", ref.level = 1)
#Null Model for SKAT with a continuous phenotype and a covariate
H0.SKAT.continuous <- NullObject.parameters(pheno = random.continuous.pheno,
                                            RVAT = "SKAT", pheno.type = "continuous",
                                            data = random.covar)
}

