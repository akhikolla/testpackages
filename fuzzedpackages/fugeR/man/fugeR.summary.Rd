\name{fugeR.summary}
\alias{fugeR.summary}
\title{Summarize a fuzzy system.}
\usage{fugeR.summary(fuzzySystem)}
\description{Summarize a fuzzy system.}
\details{Show the text description of a fuzzy system in a human readable form.}
\seealso{\code{\link{fugeR.run}}}
\author{Alexandre Bujard, HEIG-VD, Jul'2012}
\arguments{\item{fuzzySystem}{[NULL] The fuzzy system to show.}}
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

fugeR.summary(fis)
}}
