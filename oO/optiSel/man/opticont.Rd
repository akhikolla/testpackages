\name{opticont}
\Rdversion{1.1}
\alias{opticont}
\title{Optimum Contributions of Selection Candidates}

\description{
The optimum contributions of selection candidates to the offspring are calculated. The optimization procedure can take into account conflicting breeding goals, which are to achieve genetic gain, to reduce the rate of inbreeding, and to recover the original genetic background of a breed. 

It can be used for overlapping as well as for non-overlapping generations. In the case of overlapping generations, average values of the parameters for the population in the next year will be optimized, whereas for non-overlapping generations, average values of the parameters in the next generation will be optimized. Below, the "next evaluation time" means the next year for populations with overlapping generations, but the next generation for populations with non-overlapping generations. 

Optimization can be done for several breeds or breeding lines simultaneously, which is adviseable if the aim is to increase diversity or genetic distance between them.
}
\usage{
opticont(method, cand, con, bc=NULL,  solver="default", quiet=FALSE, 
         make.definite=FALSE, ...)
}

\arguments{
\item{method}{Character string \code{"min.VAR"}, or \code{"max.VAR"}, whereby \code{VAR} is the name of the variable to be minimized or maximized. Available methods are reported by function \link{candes}.}

\item{cand}{An R-Object containing all information describing the individuals (phenotypes and kinships). These indivdiuals are a sample from the population that includes the selection candidates. It can be created with function \link{candes}. This object also defines whether generations are overlapping or non-overlapping.

* If the aim is to increase genetic distance between breeds, then samples from several breeds are needed.

* If column \code{Sex} of data frame \code{cand$phen} contains \code{NA} for one breed, then the constraint stating that contributions of both sexes must be equal is omitted.}

\item{con}{List defining threshold values for constraints. The components are described in the Details section. If one is missing, then the respective constraint is not applied. Permitted constraint names are reported by function \link{candes}.}

\item{bc}{Named numeric vector with breed contributions, which is only needed if \code{cand$phen} contains individuals from different breeds. It contains the proportion of each breed in a hypothetical multi-breed population for which the diversity across breeds should be managed. The names of the components are the breed names. }

\item{solver}{Name of the solver used for optimization. Available solvers are  \code{"alabama"}, \code{"cccp"}, \code{"cccp2"}, and \code{"slsqp"}. Solver \code{"csdp"} is disabled because the package Rcsdp has been removed from Cran. By default, the solver is chosen automatically. The solvers are the same as for function \link[optiSolve]{solvecop} from package \code{optiSolve}.}

\item{quiet}{If \code{quiet=FALSE} then detailed information is shown.}
\item{make.definite}{Logical variable indicating whether non-positive-semidefinite matrices should be approximated by positive-definite matrices. This is always done for solvers that are known not to convergue otherwise.}

\item{...}{Tuning parameters of the solver. The available parameters depend on the solver and will be printed when function \code{opticont} is used with default values. Definitions of the tuning parameters can be found for \code{alabama}  in \link[alabama]{auglag} and \link[stats]{optim}, for \code{cccp} and \code{cccp2} in \link[cccp]{ctrl}, and for \code{slsqp}  in \link[nloptr]{nl.opts}.}

}
\details{

The optimum contributions of selection candidates to the offspring are calculated. The proportion of offspring that should have a particular selection candidate as parent is twice its optimum contribution.

\bold{Constraints}

Argument \code{con} is a list defining the constraints. Permitted names for the components are displayed by  function \link{candes}. Their meaning is as follows:

\bold{uniform}: Character vector specifying the breeds or sexes for which the contributions are not to be optimized. Within each of these groups it is assumed that all individuals have equal (uniform) contributions. Character string \code{"BREED.female"} means that all females from breed \code{BREED} have equal contributions and thus equal numbers of offspring. Column 'isCandidate' of \code{cand$phen} is ignored for these individuals.

\bold{lb}: Named numeric vector containing lower bounds for the contributions of the selection candidates. The component names are their IDs. By default the lower bound is 0 for all individuals.

\bold{ub}: Named numeric vector containing upper bounds for the contributions of the selection candidates. Their component names are the IDs. By default no upper bound is specified.

\bold{ub.VAR}: Upper bound for the expected mean value of kinship or trait \bold{VAR} in the population at the next evaluation time. Upper bounds for an arbitrary number of different kinships and traits may be provided. If data frame \code{cand$phen} contains individuals from several breeds, the bound refers to the mean value of the kinship or trait in the multi-breed population.

\bold{ub.VAR.BREED}: Upper bound for the expected mean value of kinship or trait \bold{VAR} in the  breed \bold{BREED} at the next evaluation time. Upper bounds for an arbitrary number of different kinships and traits may be provided. 

Note that \bold{VAR} must be replaced by the name of the variable and \bold{BREED} by the name of the breed. For traits, lower bounds can be defined as \bold{lb.VAR} or \bold{lb.VAR.BREED}. Equality constraints can be defined as  \bold{eq.VAR} or \bold{eq.VAR.BREED}.


\bold{Application to multi-breed data}

Optimization can be done for several breeds or breeding lines simultaneously, which is adviseable if the aim is to increase genetic diversity in a multi-breed population, or to increase the genetic distances between breeds or breeding lines. However, for computing the kinship of individuals from different breeds, marker data is needed. 

The multi-breed population referred above is a hypothetical subdivided population consisting of purebred animals from the breeds included in column \code{Breed} of \code{cand$phen}. The proportion of individuals from a given breed in this population is its breed contribution specified in argument \code{bc}. It is not the proportion of individuals of this breed in data frame \code{cand$phen}.

The aim is to minimize or to constrain the average genomic kinship in this multi-breed population. This causes the genetic distance between the breeds to increase, and thus may increase the conservation value of the breeds, or the heterosis effects in crossbred animals.


\bold{Remark}

If the function does not provide a valid result due to numerical problems then try to use another solver, use other optimization parameters, define upper or lower bounds instead of equality constraints, or relax the constraints to ensure that the optimization problem is solvable. 
}

\value{
A list with the following components

\item{parent}{Data frame \code{cand$phen} with some appended columns. Column \code{oc} contains the optimum contributions of the selection candidates, column \code{lb} the lower bounds, and \code{ub} the upper bounds for the contributions.}
\item{info}{Data frame with component \code{valid} indicating if all constraints are fulfilled, component \code{solver} containing the name of the solver used for optimization, and component \code{status} describing the solution as reported by the solver.}
\item{mean}{Data frame containing the expected mean value of each kinship and trait in the population at the next evaluation time.}
\item{bc}{Data frame with breed contributions in the hypothetical multi-breed population used for computing the average kinship across breeds.}
\item{obj.fun}{Named numeric value with value and name of the objective function.}
\item{summary}{Data frame containing one row for each constraint with the value of the constraint in column \code{Val}, and the bound for the constraint in column \code{Bound}. Column \code{OK} states if the constraint is fulfilled, and column \code{Breed} contains the name of the breed to which the constraint applies.  The value of the objective function is shown in the first row.  Additional rows contain the mean values of traits and kinships in the population at the next evaluation time which are not constrained.}
}


\examples{

## For other objective functions and constraints see the vignettes


######################################################
# Example 1: Advanced OCS with overlapping           #
#            generations using pedigree data         #
#   - maximize genetic gain                 (BV)     #
#   - restrict increase of mean kinship     (pKin)   #
#   - restrict increase of native kinship   (pKinatN)#
#   - avoid decrease of native contribution (NC)     #
######################################################

### Define object cand containing all required 
### information on the individuals

data(PedigWithErrors)
Pedig    <- prePed(PedigWithErrors, thisBreed="Hinterwaelder", lastNative=1970, 
                   keep=PedigWithErrors$Born\%in\%1992)
Pedig$NC <- pedBreedComp(Pedig, thisBreed="Hinterwaelder")$native
use      <- Pedig$Born \%in\% (1980:1990) & Pedig$Breed=="Hinterwaelder"
use      <- use & summary(Pedig)$equiGen>=3
cont     <- agecont(Pedig, use, maxAge=10)

Phen     <- Pedig[use, ]
pKin     <- pedIBD(Pedig, keep.only=Phen$Indiv)
pKinatN  <- pedIBDatN(Pedig, thisBreed="Hinterwaelder",  keep.only=Phen$Indiv)
Phen$isCandidate <-  Phen$Born < 1990
cand     <- candes(phen=Phen, pKin=pKin, pKinatN=pKinatN, cont=cont)

### Mean values of the parameters in the population:

cand$mean
#          BV        NC      pKin    pKinatN
#1 -0.5648208 0.5763161 0.02305245 0.0469267


### Define constraints for OCS
### Ne: Effective population size
### L:  Generation interval

Ne   <- 100
L    <- 1/(4*cont$male[1]) + 1/(4*cont$female[1])
con <- list(uniform    = "female",
            ub.pKin    = 1-(1-cand$mean$pKin)*(1-1/(2*Ne))^(1/L),
            ub.pKinatN = 1-(1-cand$mean$pKinatN)*(1-1/(2*Ne))^(1/L),
            lb.NC      = cand$mean$NC)

### Solve the optimization problem

Offspring  <- opticont("max.BV", cand, con, trace=FALSE)

### Expected average values of traits and kinships 
### in the population now and at the next evaluation time

rbind(cand$mean, Offspring$mean)   
#          BV        NC       pKin    pKinatN
#1 -0.5648208 0.5763161 0.02305245 0.04692670
#2 -0.4972679 0.5763177 0.02342014 0.04790944

### Data frame with optimum contributions

Candidate <- Offspring$parent
Candidate[Candidate$oc>0.01, c("Indiv", "Sex", "BV", "NC", "lb", "oc", "ub")] 


######################################################
# Example 2: Advanced OCS with overlapping           #
#            generations using genotype data         #
#   - minimize mean kinship                 (sKin)   #
#   - restrict increase of native kinship   (sKinatN)#
#   - avoid decrease of breeding values     (BV)     #
#   - cause increase of native contribution (NC)     #
######################################################

\dontrun{
### Prepare genotype data

data(map) 
data(Cattle)

### Compute genomic kinship and genomic kinship at native segments
dir     <- system.file("extdata", package = "optiSel")
files   <- file.path(dir, paste("Chr", 1:2, ".phased", sep=""))
sKin    <- segIBD(files, map, minL=1.0)
sKinatN <- segIBDatN(files, Cattle, map, thisBreed="Angler",  minL=1.0)

### Compute migrant contributions of selection candidates 
Haplo   <- haplofreq(files, Cattle, map, thisBreed="Angler", minL=1.0, what="match")
Comp    <- segBreedComp(Haplo$match, map)
Cattle[Comp$Indiv, "NC"] <- Comp$native

Phen  <- Cattle[Cattle$Breed=="Angler",]
cand  <- candes(phen=Phen, sKin=sKin, sKinatN=sKinatN, cont=cont)

### Define constraints for OCS
### Ne: Effective population size
### L:  Generation interval

Ne <- 100 
L  <- 4.7 
con <- list(uniform    = "female",
            ub.sKinatN = 1-(1-cand$mean$sKinatN)*(1-1/(2*Ne))^(1/L),
            lb.NC      = 1.03*cand$mean$NC,
            lb.BV      = cand$mean$BV)

# Compute optimum contributions; the objective is to minimize mean kinship 
Offspring   <- opticont("min.sKin", cand, con=con)

# Check if the optimization problem is solved 
Offspring$info           

# Average values of traits and kinships 
rbind(cand$mean, Offspring$mean)         
#           BV        NC       sKin    sKinatN
#1 -0.07658022 0.4117947 0.05506277 0.07783431
#2 -0.07657951 0.4308061 0.04830328 0.06395410

# Value of the objective function 
Offspring$obj.fun
#      sKin 
#0.04830328 

### Data frame with optimum contributions

Candidate <- Offspring$parent
Candidate[Candidate$oc>0.01, c("Indiv", "Sex", "BV", "NC", "lb", "oc", "ub")] 


#######################################################
# Example 3: Advanced OCS with overlapping            #
#            generations using genotype data          #
#            for multiple breeds or beeding lines     #
#   - Maximize breeding values in all breeds          #
#   - restrict increase of kinships within each breed #
#   - reduce average kinship across breeds            #
#   - restrict increase of native kinship in Angler   #
#   - cause increase of native contribution in Angler #
# by optimizing contributions of males from all breeds#
#######################################################


cand <- candes(phen=Cattle, sKin=sKin, sKinatN.Angler=sKinatN, cont=cont)
L  <- 5
Ne <- 100

con  <- list(uniform          = "female", 
             ub.sKin          = cand$mean$sKin - 0.01/L,
             ub.sKin.Angler   = 1-(1-cand$mean$sKin.Angler)*(1-1/(2*Ne))^(1/L),
             ub.sKin.Holstein = 1-(1-cand$mean$sKin.Holstein)*(1-1/(2*Ne))^(1/L),
             ub.sKin.Rotbunt  = 1-(1-cand$mean$sKin.Rotbunt)*(1-1/(2*Ne))^(1/L),
             ub.sKin.Fleckvieh= 1-(1-cand$mean$sKin.Fleckvieh)*(1-1/(2*Ne))^(1/L),
             ub.sKinatN.Angler= 1-(1-cand$mean$sKinatN.Angler)*(1-1/(2*Ne))^(1/L), 
             lb.NC            = cand$mean$NC + 0.05/L)
            
Offspring <- opticont("max.BV", cand, con, trace=FALSE, solver="slsqp")

Offspring$mean


Candidate <- Offspring$parent[Offspring$parent$Sex=="male", ]
Candidate[Candidate$oc>0.01, c("Indiv", "Sex", "BV", "NC", "lb", "oc", "ub")] 

}

}

\references{
Wellmann, R. (2018). Optimum Contribution Selection and Mate
Allocation for Breeding: The R Package optiSel. submitted
}


\author{Robin Wellmann}
