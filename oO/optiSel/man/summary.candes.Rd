\name{summary.candes}
\Rdversion{1.1}
\alias{summary.candes}
\title{Population Genetic Parameters at Different Times}

\description{
For every time point (age cohort), several population genetic parameters are estimated. These may include the generation interval, the average kinship, the average native kinship, the native effective size (native Ne), and the native genome equivalent (NGE) of the population at that point in time. 
}
\usage{
\method{summary}{candes}(object, tlim=range(object$phen$Born, na.rm=TRUE), 
     histNe=NA, base=tlim[1], df=4, ...)}

\arguments{
\item{object}{R-Object created with function \link{candes} containing phenotypes and kinship information on individuals from the same breed. Data frame \code{cand$phen} has columns \code{Indiv} (with IDs of the individuals),  \code{Born} (with the years-of-birth or generation-numbers), \code{Sex}, \code{I} (average age of the parents at date of birth), and \code{Offspring} (indicating if the individual has offspring). Typically function \code{prePed} is used to create them. For computing the native Ne the individuals must be from different age cohorts.}
\item{tlim}{Numeric vector with 2 components giving the time span for which genetic parameters are to be computed.}
\item{histNe}{The historic effective size of the population assumed for the time between year \code{base} and \code{tlim[1]}, which affects the NGE.}
\item{base}{The base year in which individuals are assumed to be unrelated. The base year affects the NGE. The default is \code{tlim[1]}.} 
\item{df}{Smoothing parameter used for computing the native effective size. The default is \code{df=4}.}
\item{...}{further arguments passed to or from other methods}
}

\details{
For every time point (age cohort), several population genetic parameters are estimated. These may include the generation interval, the average kinship, the average native kinship, the native effective size (native Ne), and the native genome equivalent (NGE) of the population at that point in time. The population at a time t consists of all individuals born between t-I and t, where I is the generation interval. The population genetic parameters are described below. 
}

\value{
A data frame providing for each time point (age cohort) several population genetic parameters. These may include 

\item{\code{t}}{The age cohort, containing e.g. year-of-birth or the generation number. These are the levels of column \code{Born} from data frame \code{cand$phen}.}

\item{\code{I}}{The estimated generation interval at the time when the individuals were born.}

\item{\code{KIN}}{The average kinship \code{KIN} in the population at the time when the individuals were born, where \code{KIN} is the name of a kinship. It is an estimate of the probability that 2 alleles chosen from the population are IBD.
}

\item{\code{NATKIN}}{The average native kinship \code{NATKIN} in the population at the time when the individuals were born, where \code{NATKIN} is the name of a native kinship. It is an estimate of the conditional probability that 2 alleles chosen from the population are IBD, given that both are from native ancestors.}

\item{\code{Ne}}{The native effective size of the population at the time when the individuals were born. The native effective size, quantifies how fast the smoothed native kinship is increasing. The native kinship may decrease for a short time span, in which case the estimate would be NA. Use a smaller value for parameter \code{df} to get a smoother estimate.}

\item{\code{NGE}}{The native genome equivalents of the population at the time when the individuals were born. The NGE estimates the  number of unrelated individuals that would be needed to establish a hypothetical new population that has the same genetic diversity at native alleles as the population under study, whereby the individuals born in the base-year are assumed to be unrelated.}
}


\examples{
data(ExamplePed)
Pedig   <- prePed(ExamplePed, thisBreed="Hinterwaelder", lastNative=1970)
phen    <- Pedig[Pedig$Breed=="Hinterwaelder",]
pKin    <- pedIBD(Pedig)
pKinatN <- pedIBDatN(Pedig, thisBreed="Hinterwaelder")
pop     <- candes(phen=phen, pKin=pKin, pKinatN=pKinatN, quiet=TRUE, reduce.data=FALSE)
Param   <- summary(pop, tlim=c(1970,1995), histNe=150, base=1800, df=4)

plot(Param$t, Param$pKinatN,type="l", ylim=c(0,0.1))
lines(Param$t,Param$pKin, col="red")

plot(Param$t, Param$Ne,  type="l", ylim=c(0,100))
plot(Param$t, Param$NGE, type="l", ylim=c(0,10))
plot(Param$t, Param$I,   type="l", ylim=c(0,10))
}



\author{Robin Wellmann}
