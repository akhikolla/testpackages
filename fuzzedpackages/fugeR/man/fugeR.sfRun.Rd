\name{fugeR.sfRun}
\alias{fugeR.sfRun}
\title{The parallel version of fugeR...}
\usage{fugeR.sfRun(data, labels, maxRules=4,
    maxVarPerRule=3, labelsMf=2, population=200,
    elitism=NA, mutation=0.01, generation=100,
    sensiW=1, speciW=1, accuW=0, threshold=0.5,
    rmseW=0.2, verbose=FALSE, path, rep=300,
    parallel=FALSE, cpus=1)}
\description{The parallel version of fugeR.run using snowfall.}
\details{The parallel version of fugeR.run. Will launch \code{fugeR.run} a number of times
given as argument. This function use \code{\link{snowfall}} package in order to take benefit
of mutli-core computers.

\code{fugeR.sfRun} sould be used when you want to repeat an experience many times. \cr
This is usefull when you are searching the good parameters (maxRules, macVarPerRule) for a problem.
\code{fugeR.sfRun} will launch \code{fugeR.run} and test the obtained system. It automatically resamples the data
using bootstrapping method. \cr
For example if the argument \code{rep} has the value 1000 and the number of sample in data is 100.
FugeR.sfRun resample the data with replacement with the size of the resample equal to 100
(the size of the original data set) this constitute the training set, the samples that were not picked are taken
to create the validation set. FugeR.run is then called with the training set and the obtained fuzzy systems
is tested on the validation set. If \code{rep} value was 1000, this operation is repeated 1000 times. \cr
FugeR.sfRun saves every systems in the directory specified by \code{path} and return a resume of the performance
obtained by each system on their training and validation set.}
\value{\code{res}, A data.frame of size \code{rep} containing the performance of each fuzzy system on
training and validation set.}
\seealso{\code{\link{fugeR.run}} \code{\link{fugeR.predict}}\code{\link{fugeR.summary}}
\code{\link{fugeR.save}} \code{\link{fugeR.load}}}
\author{Alexandre Bujard, HEIG-VD, Jul'2012}
\arguments{\item{data}{[NULL] Data frame to be used for training (only numeric values are supported).}
\item{labels}{[NULL] Labels of \code{data} (only numeric values are supported).}
\item{maxRules}{[4] Maximum number of rule.}
\item{maxVarPerRule}{[3] Maximum number of input variable per rule.}
\item{labelsMf}{[2] Number of singleton for output variable membership function.}
\item{population}{[200] The population size.}
\item{elitism}{[NA] The number of chromosomes that are kept into the next generation.
By default is about 20\% of the population size.}
\item{mutation}{[0.01] The chance that a gene in the chromosome mutates.}
\item{generation}{[100] The number of generation made by the genetic algorithm.}
\item{sensiW}{[1.0] The weight of the sensitivity in the fitness function.}
\item{speciW}{[1.0] The weight of the specificity in the fitness function.}
\item{accuW}{[0.0] The weight of the accuracy in the fitness function.}
\item{threshold}{[0.5] The threshold to apply in order to calculate sensitivity, specificity and accuracy.}
\item{rmseW}{[0.2] The weight of the "root mean square error" between labels and
values predicted by the fuzzy system.}
\item{verbose}{[FALSE] If true the algorithm will be more verbose. By default False.}
\item{path}{[NULL] THe path where to save the fuzzy systems.}
\item{rep}{[300] Number fuzzy system to find.}
\item{parallel}{[TRUE] Logical value indicating if the function can run in parralel}
\item{cpus}{[1] number of cpus that can be used}}
\examples{##
##
\dontrun{
expResume <- fugeR.sfRun (
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
labelsMf=2,
path=\'./exp\',
rep=100,
parallel=TRUE,
cpus=2
)
}}
