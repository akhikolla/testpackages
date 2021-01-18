\name{dfpk-package}
\alias{dfpk-package}
\alias{dfpk}

\docType{package}
\title{Bayesian Dose-Finding Designs using Pharmacokinetics (PK) for Phase I Clinical Trials.
}

\description{
Statistical methods involving PK measures are provided, in the dose allocation process during a Phase I clinical trials. These methods, proposed by Ursino et al, (2017) <doi:10.1002/bimj.201600084>, enter pharmacokinetics (PK) in the dose finding designs in different ways, including covariates models, dependent variable or hierarchical models. This package provides functions to generate data from several scenarios and functions to run simulations which their objective is to determine the maximum tolerated dose (MTD).
}

\details{
The three main functions of the dfpk package are \code{\link{sim.data}}, \code{\link{nsim}} and \code{\link{nextDose}}, for generating PK data based on the input settings, simulating "n" clinical trials and determining the next recommended dose for an ongoing phase I clinical trial based on an enrolled patient's data, respectively.
Subsequently, six dose-finding methods/models can be applied: \code{\link{dtox}}, \code{\link{pktox}}, \code{\link{pkcrm}}, \code{\link{pkcov}}, \code{\link{pkpop}}, and \code{\link{pklogit}} which each one creates a Bayesian model and fits it using \pkg{Stan}.

Since \pkg{dfpk} is based on \pkg{Stan} models, a C++ compiler is required. The program Rtools (available on https://cran.r-project.org/bin/windows/Rtools/) comes with a C++ compiler for Windows while on Mac, you should use Xcode.
For further instructions on how to get the compilers running,
see the prerequisites section on https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started.
}

\author{
Artemis Toumazi \email{artemis.toumazi@gmail.com} , Moreno Ursino \email{moreno.ursino@inserm.fr} and Sarah Zohar \email{sarah.zohar@inserm.fr}

Maintainer: Artemis Toumazi \email{artemis.toumazi@gmail.com}
}

\references{
Ursino, M., et al, (2017) Dose-finding methods for Phase I clinical trials using pharmacokinetics in small populations, Biometrical Journal, <doi:10.1002/bimj.201600084>.

Toumazi, A., et al, (2018) dfpk: An R-package for Bayesian dose-finding designs using pharmacokinetics (PK) for phase I clinical trials, Computer Methods and Programs in Biomedicine, <doi:10.1016/j.cmpb.2018.01.023>.
}

\keyword{package}
