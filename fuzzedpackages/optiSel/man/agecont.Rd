\name{agecont}
\Rdversion{1.1}
\alias{agecont}
\title{Contributions of age cohorts to the population}

\description{
Contributions of age classes to the population are calculated such that the contribution of each age class to the population is proportional to the expected proportion of offspring that is not yet born. 

Note that the contribution of a class to the population is not equal to the proportion of individuals belonging to the class.
}
\usage{
agecont(Pedig, use=Pedig$Born >= quantile(Pedig$Born, 0.75), maxAge=NA)}

\arguments{
\item{Pedig}{Pedigree with colums \code{Indiv}, \code{Sire}, \code{Dam}, and \code{Born}, usually created with function \link{prePed}.}
\item{use}{Logical vector or character vector with IDs indicating 
the individuals from the current population.}
\item{maxAge}{Parents that are more than \code{maxAge} years older than their offspring are ignored. By default, old parents are not ignored.}
}

\details{

Contributions of age classes to the population are calculated such that the contribution of each age class to the population is proportional to the expected proportion of offspring that is not yet born. 


More precisely:

Individuals born in the current year are in age class \code{k=1}. Typically, each age class spans one year. No individual can have offspring in the same age class. Males and females that are not born in the current year are assumed to have equal contributions to the population. Moreover, as stated above, it is assumed that the contribution of each class to the population is proportional to the proportion of offspring from this class that is not yet born when the individuals leaves the class. 

This approach to define contributions has the advantage that it does not need to be known which individuals are still alive and which are removed from the breeding pool. Moreover, it causes old age classes to have a smaller contribution to the population than young age classes.

The contributions are estimated from the ages of the parents when the individuals in vector \code{use} were born. Obviously, the contributions of age classes to the offspring in the next year do not coincide with the contributions of the age classes to the population.

}

\value{
Data frame containing the contributions of all age cohorts to the current population.
}


\examples{

data(PedigWithErrors)
Pedig <- prePed(PedigWithErrors)
use   <- Pedig$Breed=="Hinterwaelder" & !is.na(Pedig$Born)
use   <- use & Pedig$Born>=2000 & Pedig$Born<=2004 

# Calculate the contribution of each age class ##

cont <- agecont(Pedig, use)

# Contribution of each age class to 
# the current population:

head(cont)

# Note: In this case, young males have a higher contribution to the  
# population than young females because they are used for breeding 
# for a shorter time span, i.e. they are culled earlier.

# Males and females (excluding the newborn individuals)
# have equal contributions to the current population:

sum(cont$male[-1])
#[1] 0.3925894

sum(cont$female[-1])
#[1] 0.3925894


# The total contribution of classes to the curent population is equal to 1 

sum(cont$female) + sum(cont$male) 
#[1] 1

# When used for OCS, the contribution of the offspring to the 
# population in the next year is equal to the contribution of the individuals
# born in this year to the current population:

cont$male[1]+cont$female[1]
#[1] 0.2148212

# This is approximately 1/L, where L is the generation interval.

}


