\name{candes}
\alias{candes}

\title{Candidate Description}

\description{
An R-Object is created containing all information describing the individuals, which is usually a sample from the current population and includes the selection candidates. Average kinships and trait values, and the available objective functions and constraints for optimum contribution selection (OCS) are reported.
}

\usage{
candes(phen,  cont=NULL, N=1000, quiet=FALSE, t=NA, bc=NULL, reduce.data=TRUE, ...)
}

\arguments{
 \item{phen}{Data frame with column \code{Indiv} containing animal IDs and possibly \code{Sex} containng sexes, coded as \code{'male'} and \code{'female'}, or \code{NA} if sexes are to be ignored. It also contains column \code{Born} with year of birth if generations are assumed to be overlapping. The other columns may contain traits, e.g. breeding values or native contributions, column \code{Breed} with breed names for multi-breed evaluations, logical column \code{isCandidate} indicating the selection candidates, and columns \code{Sire} and \code{Dam} with IDs of sires and dams.}
 \item{cont}{Data frame frame with column \code{age} (equal to the row number), and columns \code{male}, and \code{female}, containing the contributions of males and females from each age class to the population. It is usually created with function \link{agecont}. The default means that non-overlapping generations are  assumed, so there is only one age class for males and one for females.}
 \item{N}{The population size. A small value accelerates the increase in kinship due to genetic drift. For overlapping generations it can be calculated as \code{N=N0/r0}, where \code{N0} is the number of individuals born each year, and \code{r0<=1} is the  percentage which this age class represents in the population. The default is \code{N=1000}.}
  \item{quiet}{Should the report be suppressed?}
  \item{t}{The time at which the population should be evaluated. The default means that \code{t=max(floor(phen$Born))}.}
\item{bc}{Only needed if multi-breed data is provided. Named vector with breed contributions, with component names being the names of the breeds in \code{phen}.  It contains the proportion of each breed to a hypothetical multi-breed population for which the diversity across breeds should be managed. Alternatively, \code{bc} can be a character string containig the name of a kinship. In this case, optimum contributions of the breeds are determined automatically so that the mean kinship across breeds is minimized.
}
\item{reduce.data}{Logical. Should data from individuals not contributing to the population at time \code{t} be removed from the output?}
 \item{...}{One or more objects of class \code{'matrix'}, \code{'quadFun'}, or \code{'ratioFun'} defining the pairwise kinships and native kinships of individuals.}
}
\value{
List of class candes with the following components:
\item{kinship}{Objects of class \code{'quadFun'}, or \code{'ratioFun'}, one for each additional parameter. These objects define the functions needed to estimate the mean kinships and mean native kinships in the next year or generation.}
\item{phen}{Supplied data frame \code{phen} containing phenotypes, individual IDs, and some appended columns that are needed for OCS. These are

* Column \code{Age} with the ages of the individuals,  

* Column \code{Class} with the \code{Breed x Age x Sex} or \code{Breed x Age} classes to which the individuals belong. 

* Column \code{c0} containing the contribution each individual itself has to the current population. 

* Column \code{c1} containing the contribution each individual itself has to the population in the next year (for overlapping generations) or to the next generation (for non-overlapping generations). In the latter case, \code{c1} contains zeros.

* Column \code{isCandidate} indicating the selection candidates. }
\item{mean}{Data frame containing estimates of the current mean values (at time \code{t}) of the parameters in a population consisting of \code{N} individuals for which the individuals in argument \code{phen} are representative. }
\item{current}{Data frame containing the same values as component \code{mean}, but also some additional information on the parameters.}
\item{bc}{Character vector with optimum breed contributions (see above).}
\item{classes}{Data frame containing the number of individuals in each class (column \code{n}), the contribution of each class to the population in this year/generation (column \code{rcont0}) and in the next year/generation (column \code{rcont1}), and the expected proportion of offspring animals from a given sex have at a particular age.}
\item{breed}{List describing the breeds included in the data set.}
}

\details{
An R-Object is created containing all information describing the individuals, which is usually the current population and includes the selection candidates. Average kinships and trait values are estimated and reported. The weights of Age x Sex classes are in accordance with argument \code{cont}. The available objective functions and constraints for optimum contribution selection are reported.
}

\examples{
data(PedigWithErrors)

Pedig      <- prePed(PedigWithErrors, thisBreed="Hinterwaelder", lastNative=1970, 
                   keep=PedigWithErrors$Born\%in\%1992)
use        <- Pedig$Born \%in\% (1980:1990) & Pedig$Breed=="Hinterwaelder"
Population <- Pedig$Indiv[use]

Pedig$NC   <- pedBreedComp(Pedig, thisBreed="Hinterwaelder")$native
pKin       <- pedIBD(Pedig, keep.only=Population)
pKinatN    <- pedIBDatN(Pedig, thisBreed="Hinterwaelder",  keep.only=Population)
Phen       <- Pedig[Population, ]

### Example 1: Overlapping Generations
### Old individuals contribute only little to the means:

cont <- agecont(Pedig, Population, maxAge=10)
cand <- candes(phen=Phen, pKin=pKin, pKinatN=pKinatN, cont=cont)

cand$current[,c("Name", "Type", "Breed", "Val", "Var")]
#     Name      Type         Breed         Val     Var
#1      BV     trait Hinterwaelder -0.55979308      BV
#2      NC     trait Hinterwaelder  0.56695077      NC
#3    pKin   kinship Hinterwaelder  0.02230896    pKin
#4 pKinatN nat. kin. Hinterwaelder  0.04678453 pKinatN

# BV:      simulated breeding values
# NC:      native genetic contribution computed from pedigree
# pKin:    pedigree-based kinship
# pKinatN: pedigree-based native kinship


### Example 2: Discrete Generations (cont=NULL). 
### Old individuals and young individuals contribute equally to the means:

Phen$Born <- 1
cand <- candes(phen=Phen, pKin=pKin, pKinatN=pKinatN, cont=NULL)

cand$current[,c("Name", "Type", "Breed", "Val", "Var")]

#     Name      Type         Breed         Val     Var
#1      BV     trait Hinterwaelder -0.71910508      BV
#2      NC     trait Hinterwaelder  0.58226604      NC
#3    pKin   kinship Hinterwaelder  0.01979228    pKin
#4 pKinatN nat. kin. Hinterwaelder  0.04053012 pKinatN


### Shorthand:

cand$mean
#          BV       NC       pKin    pKinatN
#1 -0.7191051 0.582266 0.01979228 0.04053012

cand$mean$pKin
#[1] 0.01979228

}

\author{Robin Wellmann}

