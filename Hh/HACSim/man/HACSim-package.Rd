\name{HACSim-package}
\alias{HACSim-package}
\docType{package}
\title{
\packageTitle{HACSim}
}
\description{
HACSim (\strong{H}aplotype \strong{A}ccumulation \strong{C}urve \strong{Sim}ulator) employs a novel nonparametric stochastic (Monte Carlo) optimization method of iteratively generating species' haplotype accumulation curves through extrapolation to assess sampling completeness based on the approach outlined in Phillips et al. (2015) <doi:10.1515/dna-2015-0008> and Phillips et al. (2019) <doi:10.1002/ece3.4757>.  Any genomic marker can be targeted to assess likely required specimen sample sizes. The method is well-suited to assess sampling sufficiency for DNA barcoding initiatives. The package outputs a number of useful summary statistics of sampling coverage ("Measures of Sampling Closeness"), including an estimate of the likely required sample size (along with desired level confidence intervals) necessary to recover a given number/proportion of observed unique species' haplotypes. 
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{HACSim}

\packageIndices{HACSim}
}

\author{
\packageAuthor{HACSim}

Maintainer: \packageMaintainer{HACSim}
}

\references{
Phillips, J.D., Gwiazdowski, R.A., Ashlock, D. and Hanner, R. (2015). An exploration of sufficient sampling effort to describe intraspecific DNA barcode haplotype diversity: examples from the ray-finned fishes (Chordata: Actinopterygii). DNA Barcodes, 3: 66-73.

Phillips, J.D., Gillis, D.J. and Hanner, R.H. (2019). Incomplete estimates of genetic diversity within species: Implications for DNA barcoding. Ecology and Evolution, 9(5): 2996-3010.

Phillips, J.D., Gillis, D.J. and Hanner, R.H. (2020). HACSim: An R package to estimate intraspecific sample sizes for genetic diversity assessment using haplotype accumulation curves. PeerJ Computer Science
}

\keyword{ package }

\examples{

## Simulate hypothetical species ##

N <- 100 # total number of sampled individuals
Hstar <- 10 # total number of haplotypes
probs <- rep(1/Hstar, Hstar) # equal haplotype frequency distribution

HACSObj <- HACHypothetical(N = N, Hstar = Hstar, 
probs = probs, filename = "output") # outputs a CSV 
# file called "output.csv"

## Simulate hypothetical species - subsampling ##
HACSObj <- HACHypothetical(N = N, Hstar = Hstar, 
probs = probs, perms = 1000, p = 0.95, 
subsample = TRUE, prop = 0.25, conf.level = 0.95, 
filename = "output")

## Simulate hypothetical species and all paramaters changed - subsampling ##
HACSObj <- HACHypothetical(N = N, Hstar = Hstar, probs = probs, 
perms = 10000, p = 0.90, subsample = TRUE, prop = 0.15, 
conf.level = 0.95, filename = "output")

HAC.simrep(HACSObj) # runs a simulation


## Simulate real species ##

\donttest{
## Simulate real species ##
# outputs file called "output.csv"
HACSObj <- HACReal(filename = "output") 

## Simulate real species - subsampling ##
HACSObj <- HACReal(subsample = TRUE, prop = 0.15, 
conf.level = 0.95, filename = "output")

## Simulate real species and all parameters changed - subsampling ##
HACSObj <- HACReal(perms = 10000, p = 0.90, subsample = TRUE, 
prop = 0.15, conf.level = 0.99, filename = "output")

# user prompted to select appropriate FASTA file
HAC.simrep(HACSObj) }

}
