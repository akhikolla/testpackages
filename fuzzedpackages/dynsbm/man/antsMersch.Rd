\name{antsMersch}
\docType{data}
\alias{antsMersch}
\title{Dynamic contact network of ants}
\description{
Colonies  of the  ant Camponotus fellah were  followed with  a
tracking system that  monitored the individual positions  over days of
observations  and  dynamic  social   interactions  were  deduced  from
physical proximity. The data corresponds to a colony of N=152 ants observed during T=10 days.
Edge weights were binned into K=3 categories corresponding to low, medium and high interaction intensity.
}
\usage{data(antsMersch)}
\format{An array of size 10x152x152.}
\source{http://datadryad.org/resource/doi:10.5061/dryad.8d8h7}
\references{
Mersch, D. P., Crespi, A., Keller, L., \emph{Tracking individuals shows spatial fidelity is a key regulator of ant social organization}, Science, 340(6136), 1090-1093  (2013)
http://dx.doi.org/10.1126/science.1234316

Miele, V and Matias, C, \emph{Revealing the hidden structure of dynamic ecological networks}, Royal Society Open Science (2017)
}
\examples{
\dontrun{
data(antsMersch)

## better to use nstart>>1 starting points
## but estimation can be long;
## better to use nb.cores>1 cores
list.dynsbm <- select.dynsbm(antsMersch, 
				Qmin=1, Qmax=6, edge.type="discrete", K=3, 
				nstart=20, nb.cores=4) 
				

## selection of Q=3 with the 'elbow' method
dynsbm <- list.dynsbm[[3]]

## plotting intra/inter connectivity patterns
connectivity.plot(dynsbm, antsMersch)

## plotting switches between groups
alluvial.plot(dynsbm, timestep.abbrev="D", only.present=F)
}
}
