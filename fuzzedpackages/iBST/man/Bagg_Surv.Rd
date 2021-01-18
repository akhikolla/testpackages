\name{Bagg_Surv}
\alias{Bagg_Surv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bagging improper survival trees
}
\description{
Bagging procedure to aggregate several improper trees using either the pseudo-R2 procedure or the adjusted Logrank procedure. Several scores for variables importance are computed.
}
\usage{
Bagg_Surv(xdata, Y.names, P.names, T.names, method = "R2", args.rpart, 
          args.parallel = list(numWorkers = 1), Bag = 100)
}

\arguments{
  \item{xdata}{
The learning data frame
}
  \item{Y.names}{
A vector of the names of the two variables of interest (the time-to-event is follow by the event indicator)
}
  \item{P.names}{
The names of independant variables acting on the non-susceptible population (the plateau)
}
  \item{T.names}{
The names of independant variables acting on the survival of the susceptible population
}
  \item{method}{
The choosen method (either \code{"LR"} for the Logrank or \code{"R2"} for the proposed pseudo-R2 criterion)
}
  \item{args.rpart}{
The improper survival tree parameters: a list of options that control details of the rpart algorithm.
\code{minbucket}: the minimum number of observations in any terminal <leaf> node; \code{cp}: 
complexity parameter (Any split that does not decrease the overall lack of fit by a factor of
cp is not attempted); \code{maxdepth}:  the maximum depth of any node of the final tree, with
the root node counted as depth 0. ...
See \code{\link[rpart]{rpart.control}} for further details 
}
  \item{args.parallel}{
a list  containing the number of parallel computing arguments: The number of workers, the type of parallelization to achieve, ...  see \code{\link[parallel]{mclapply}} for further details.
}
  \item{Bag}{
The number of Bagging samples to consider 
}
}
\details{
For the Bagging procedure, it is mendatory to set \code{maxcompete = 0} and \code{maxsurrogate = 0} within the \code{args.rpart} arguments. This will ensured the correct calculation of the  importance of variables and also a better computation time.
}
\value{
A list of ten elements
\item{MaxTreeList}{The list of improper survival trees computed during the bagging procedure}
\item{IIS}{The Index Importance Score }
\item{DIIS}{The Depth Index Importance Score }
\item{DEPTH}{The minimum depth importance Score }
\item{IND_OOB}{A list of length \code{Bag} containing the Out Of Bag (OOB) individuals for improper
survival tree model }
\item{IIND_SAMP}{The final list of length \code{Bag} of sample individuals used for each improper survival tree}
\item{IIND_SAMP}{The initial list  of sample individuals used for each improper survival tree at teh beginning}
\item{Bag}{The number of bagging samples retained at the end of the procedure after removing the trees without leaves}
\item{indrpart}{a vector of \code{TRUE} or \code{FALSE} with the value FALSE when the corresponding tree is removed from the final bagged predictor}
\item{Timediff}{The ellapsed time of the Bagging procedure}
}
\references{
Mbogning, C. and Broet, P. (2016). Bagging survival tree procedure for
variable selection and prediction in the presence of nonsusceptible patients.
BMC bioinformatics, 17(1), 1.
}
\author{
Cyprien Mbogning and Philippe Broet
}
\note{
This version of the code allows for the moment only one variable to have an impact on the cured population.The next version will allow more than one variable.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Bagg_pred_Surv}}
}
\examples{
\dontrun{
 data(burn)
 myarg = list(cp = 0, maxcompete = 0, maxsurrogate = 0, maxdepth = 2)
 Y.names = c("T3" ,"D3")
 P.names = 'Z2'
 T.names = c("Z1", paste("Z", 3:11, sep = ''))
 mybag = 40
 set.seed(5000)
 
 burn.BagEssai0 <- Bagg_Surv(burn, Y.names, P.names, T.names, method = "LR", args.rpart = myarg, 
                           args.parallel = list(numWorkers = 1), Bag = mybag)
 
 burn.BagEssai1 <- Bagg_Surv(burn, Y.names, P.names, T.names, method = "R2", args.rpart = myarg, 
                           args.parallel = list(numWorkers = 1), Bag = mybag)

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ documentation }% use one of  RShowDoc("KEYWORDS")
\keyword{ tree }% __ONLY ONE__ keyword per line
\keyword{multivariate}
\keyword{survival}
