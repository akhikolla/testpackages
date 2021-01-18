\name{pk.estim}
\alias{pk.estim}
\title{
The pharmacokinetic's (PK) measure of exposure.
}

\description{
Estimation of the pharmacokinetic's (PK) measure of exposure. 
}

\usage{
pk.estim(par, t, dose, conc)
}

\arguments{
  \item{par}{
The pharmacokinetic's parameters.
}
  \item{t}{
The time sampling.
}
  \item{dose}{
The doses levels of the drug.
}
  \item{conc}{
The concentration of the drug in blood plasma.
}
}

\references{
Ursino, M., et al, (2017) Dose-finding methods for Phase I clinical trials using pharmacokinetics in small populations, Biometrical Journal, <doi:10.1002/bimj.201600084>.

Toumazi, A., et al, (2018) dfpk: An R-package for Bayesian dose-finding designs using pharmacokinetics (PK) for phase I clinical trials, Computer Methods and Programs in Biomedicine, <doi:10.1016/j.cmpb.2018.01.023>.
}

\author{Artemis Toumazi \email{artemis.toumazi@gmail.com}, 
        Moreno Ursino \email{moreno.ursino@inserm.fr},
        Sarah Zohar \email{sarah.zohar@inserm.fr}
        }

\seealso{
\code{\link{AUC.estim}}, \code{\link{nsim}}
}

\examples{
par <- c(2,10,100) 
conc <- c(20,30,40,50,60,70)
t <- seq(0,24,length.out=48)
dose <- c(12.59972,34.65492,44.69007,60.80685,83.68946,100.37111)
pk.estim(par,t,dose,conc)
}
