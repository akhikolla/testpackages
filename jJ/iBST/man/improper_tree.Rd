\name{improper_tree}
\alias{improper_tree}
\title{
imprper survival tree
}
\description{
Fit an improper survival tree for the mixed population (susceptible and nonsusceptible) using either the 
proposed pseudo R2 criterion or an adjusted Logrank criterion
}
\usage{
improper_tree(xdata, Y.names, P.names, T.names, method = "R2", args.rpart)
}
%- maybe also 'usage' for other objects documented here.
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
}
\value{
An unprunned improper survival tree
}
\references{
Mbogning, C. and Broet, P. (2016). Bagging survival tree procedure for
variable selection and prediction in the presence of nonsusceptible patients.
BMC bioinformatics, 17(1), 1.
}
\author{
Cyprien Mbogning and Philippe Broet
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Bagg_Surv}} \code{\link{Bagg_pred_Surv}} 
}
\examples{
\dontrun{
 data(burn)
 myarg = list(cp = 0, maxcompete = 0, maxsurrogate = 0, maxdepth = 3)
 Y.names = c("T3" ,"D3")
 P.names = 'Z2'
 T.names = c("Z1", paste("Z", 3:11, sep = ''))
 burn.tree <- improper_tree(burn, Y.names, P.names, T.names, method = "R2", args.rpart = myarg)
 plot(burn.tree)
 text(burn.tree, cex = .7)
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
\keyword{ documentation }% use one of  RShowDoc("KEYWORDS")
\keyword{ tree }% __ONLY ONE__ keyword per line
\keyword{multivariate}
\keyword{survival}
