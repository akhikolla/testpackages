\name{fugeR.save}
\alias{fugeR.save}
\title{Save a fuzzy system into a file.}
\usage{fugeR.save(fuzzySystem, file="")}
\description{Save a fuzzy system into a file.}
\seealso{\code{\link{fugeR.load}}}
\author{Alexandre Bujard, HEIG-VD, Jul'2012}
\arguments{\item{fuzzySystem}{[NULL] The fuzzy system to save.}
\item{file}{[\"\"] A character string naming a file.}}
\examples{##
##
\dontrun{
fis <- fugeR.run (
In,
Out,
generation=100,
population=200,
elitism=40,
verbose=TRUE,
threshold=0.5,
sensiW=1.0,
speciW=1.0,
accuW=0.0,
rmseW=1.0,
maxRules=10,
maxVarPerRule=2,
labelsMf=2
)

fugeR.save( fis, file=\'./myFis.R\' )

savedFis <- fugeR.load( file=\'./myFis.R\' )
}}
