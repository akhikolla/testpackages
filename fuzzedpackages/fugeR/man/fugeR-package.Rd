\name{fugeR-package}
\title{FUzzy GEnetic, a machine learning algorithm to construct prediction model based on fuzzy logic.}
\description{This is an evolutionary algorithm for fuzzy systems, a genetic algorithm is used to construct
 a fuzzy system able to fit the given training data.
 This fuzzy system can then be used as a 
 prediction model, it's composed of fuzzy logic rules that provide a good
 lingustic representation.}
\details{\tabular{ll}{
Package: \tab fugeR\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2012-07-11\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}

This package allow you to use a genetic algorithm
in order to find a fuzzy system that
can be used as a prediction model.

\code{\link{fugeR.run}} find a fuzzy system. \code{\link{fugeR.predict}}
make the prediction for the given input data.}
\alias{fugeR-package}
\docType{package}
\author{Alexandre Bujard \email{alexandre.bujard@gmail.com}}
\references{\url{http://library.epfl.ch/en/theses/?nr=2634}}
\keyword{package}
\keyword{machine}
\keyword{learning}
\keyword{fuzzy}
\keyword{logic}
\keyword{system}
\keyword{genetic}
\keyword{algorithm}
\keyword{model}
\keyword{prediction}
\keyword{expert}
\keyword{decision}
\keyword{data}
\keyword{mining}
\seealso{\code{\link{fugeR.run}}}
\examples{\dontrun{
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
