\name{fugeR.predict}
\alias{fugeR.predict}
\title{Compute the prediction of a fuzzy system for the given input data.}
\usage{fugeR.predict(fuzzySystem, dataset)}
\description{Compute the prediction of a fuzzy system for the given input data.}
\value{\code{prediction}, A data.frame containing the predictions.}
\note{The dataset must contain the same headers (in the same order) that the \code{data}
used to find the system with fugeR.run.}
\seealso{\code{\link{fugeR.run}}}
\author{Alexandre Bujard, HEIG-VD, Jul'2012}
\arguments{\item{fuzzySystem}{[NULL] The fuzzy system to use for computing the prediction.}
\item{dataset}{[NULL] The data to use.}}
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
