\name{fugeR.run}
\alias{fugeR.run}
\title{R based Fuzzy logic evolutionary algorithm...}
\usage{fugeR.run(data, labels, maxRules=4,
    maxVarPerRule=3, labelsMf=2, population=200,
    elitism=NA, mutation=0.01, generation=100,
    sensiW=1, speciW=1, accuW=0, threshold=0.5,
    rmseW=0.2, verbose=FALSE)}
\description{R based evolutionary algorithm for finding fuzzy systems}
\details{A machine learning algorithm for fuzzy system. \cr
This function use a genetic algorithm in order to construct a fuzzy system able
to fit the values given as \code{labels}. The \code{data} and \code{labels} are used
has learning data.

This is a fuzzy system evolutionnary algorithm. A genetic algorithm is used to find a
fuzzy system able to fit the the data given as labels. \cr
The genetic algorithm generate a random \code{population} of fuzzy system.
At each \code{generation} all the fuzzy system are tested with the data.
Their predictions are then compared with the labels and a "performance" is given at each system.
The top best system (\code{elitsm} are taken without modification for the next generation.
The population is then used to generate the population for the next generation using crossover
and mutation. At the end of the process (at the last generation) the fuzzy system that obtained
the best performance is returned.}
\value{\code{fis}, A list containing the logs of the evolution, the peformances of the best system and
its description.
\item{inputVarIds}{The IDs of the variable used in the fuzzy system}
\item{inputMfIds}{The IDs of the membership function used by each variable in the fuzzy system}
\item{inputMfs}{The values used to caclculate the membership functions}
\item{outputVarIds}{The IDs of each output variable of each rule}
\item{outputMfIds}{The IDs of the membership function used by each output variable}
\item{outputMfs}{the value used to compute the membership functions of the output variables}
\item{fitness}{The fitness value reached by the best fuzzy system}
\item{mse}{The Mean Square Error of the best fuzzy system}
\item{rmse}{The Root Mean Square Error between labels and the prediction made by the best fuzzy system}
\item{accu}{The accuracy of the prediction made by the best fuzzy system
(only if a threshold different of NA was given as argument)}
\item{sensi}{The sensitivity of the prediction made by the best fuzzy system
(only if a threshold different of NA was given as argument)}
\item{speci}{The specificity of the prediction made by the best fuzzy system
(only if a threshold different of NA was given as argument)}
\item{evo}{A list containing the evolution logs}}
\seealso{\code{\link{fugeR.sfRun}} \code{\link{fugeR.predict}}\code{\link{fugeR.summary}}
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
\item{verbose}{[FALSE] If true the algorithm will be more verbose. By default False.}}
\examples{##
## 
\dontrun{
#We use the iris dataset for this example
#We need to convert the output in a numeric format.
data(iris)
OUT <- data.matrix(iris[5])[,1]
fIris <- cbind(iris[1:4], OUT)
In <- fIris[1:4]
Out <- fIris[5]

#Launch the evolution, fugeR.run will return
#the best fuzzy system found during the evolution
fuzzySystem  <- fugeR.run( In,
                   Out,
                   generation=100, # Increase the number of generation for a better accuracy
                   population=100,
                   elitism=20,
                   verbose=TRUE,
                   threshold=NA,
                   sensiW=0.0,
                   speciW=0.0,
                   accuW=0.0,
                   rmseW=1.0,
                   maxRules=5,
                   maxVarPerRule=2,
                   labelsMf=3
)

#Plot the predicton given by the best fuzzy system found during the evolution
prediction <- fugeR.predict(fuzzySystem, In)
plot(prediction[[1]], ylim=c(1,max(unlist(Out))), col='blue', pch=21, axes=FALSE, ann=FALSE)
points(Out[[1]], col="red", pch=21)
axis(1)
axis(2, at=1:3, lab=c('setosa', 'versicolor', 'virginica'))
title(main='Fuzzy system prediction on Iris problem')
title(xlab="Cases")
title(ylab="Specie")
box()
legend(0.0, 3.0, c("Predicted","Actual"), cex=0.8, 
       col=c("blue","red"), pch=c(21,21))

#Display the fuzzy system
fugeR.summary(fuzzySystem)
}}
