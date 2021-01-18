\name{stosim-package}
\alias{stosim}
\docType{package}
\title{
Stochastic Simulator for Reliability Modeling of Repairable Systems
}
\description{
A toolkit intended for Reliability Availability and Maintainability (RAM) modeling of industrial process systems. 

It is recommended for use with RExcel for data input, holding of the model 
scripts, and ad hoc analysis of reliability parameters.

}

\section{Introduction}{
	\pkg{stosim} provides functions for creating reliability models using observed
	data reduced to probability distributions for failure and repair mechanisms
	on related operations in order to infer expected performance of new systems
	or alteration of existing systems.   
	Models can be assembled from small sub-systems and accumulated to describe
	an entire production plant or refinery.  Stochastic modeling provides an ideal 
	means for study of the performance of product inventory storage and parallel 
	operations as reliability enhancement features.  Time dependent issues such as 
	seasonal variation, and equipment degradation can be accurately assessed for 
	impact on ultimate production capability.  Contractual conditions such as 
	bonus/penalty clauses can be evaluated with realitic statistical projections.
}

\author{
Jacob T. Ormerod

Maintainer: Jacob T. Ormerod <jake@openreliability.org>
}
\references{
  Jones, O.D., R. Maillardet, and A.P. Robinson (2009) An Introduction
  to Scientific Programming and Simulation, Using R.  Chapman And Hall/CRC
  
  Robert, Christian P., G. Casella (2010) Introducing Monte Carlo Methods with R.
  Springer

  Taylor HM,  Karlin S (1998) An Introduction to Stochastic Modeling, 3rd Edition,
  Acadmic Press.
  
  Silkworth, David J. (1998) "Confidence Curves: A Reliability Modelling Technique for
  the Practical Application of Process Unit and Subsystem Failure Data". American
  Institute of Chemical Engineers
  
  Tobias, Paul A., D.C. Trinidade (1986)Applied Reliability.  Van Nostand Reinhold
}

\keyword{ package }



