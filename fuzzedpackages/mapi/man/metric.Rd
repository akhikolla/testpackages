\name{metric}
\alias{metric}
\title{'metric' test dataset}
\usage{
data("metric")
}
\value{
A data.table object with 19900 rows (one per sample pair, non symmetrical) and 3 columns ("ind1", "ind2", "value") 
containing respectively the first sample name, the second sample name and the value of their relation.
}
\description{
The individuals were simulated as described in the \code{\link{samples}} section.
The value is a genetic distance (â estimator in Rousset, 2000) computed between the 
200 simulated samples using SPAGeDi v1.4 (Hardy & Vekemans, 2002).
}
\examples{
data("metric")
}
\references{
Rousset, F. (2000). Genetic differentiation between individuals. Journal of Evolutionary Biology, 13:58--62.

Hardy OJ, Vekemans X (2002) SPAGeDi: a versatile computer program to analyse spatial genetic structure 
at the individual or population levels. Molecular Ecology Notes 2: 618-620.
}
