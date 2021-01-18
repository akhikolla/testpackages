\name{permute_select_surv}
\alias{permute_select_surv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
permutation variable selection
}
\description{
Variable selection using the permutation test on several scores of importance: \code{IIS}, \code{DIIS} and \code{DEPTH}. 
}
\usage{
permute_select_surv(xdata, Y.names, P.names, T.names, importance = "IIS", method = "R2",
                    Bag, args.rpart, args.parallel = list(numWorkers = 1), nperm = 50)
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
  \item{importance}{
The importance score to consider: either \code{IIS}, \code{DIIS} or \code{DEPTH}
}
  \item{method}{
The splitting method: either \code{"R2"} for the proposed pseudo-R2 criterion or \code{"LR"} for the adjusted Logrank criterion
}
  \item{Bag}{
The number of Bagging samples to consider
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
  \item{nperm}{
The number of permutation samples to consider for the permutation test
}
}
\details{
Testing weither the importance score is null or not.
}
\value{
A list of five elements:
\item{pvalperm1}{The permutation test P-values ranking in decreasing order}
\item{pvalperm2}{The permutation test P-values ranking in decreasing order considering an approximate gaussian distribution under the null hypothesis}
\item{pvalKS}{The Kolmogorov-Smirnov P-values of the comparisons between the observed importance under the null hypothesis and a theoretical gaussian distribution}
\item{IMPH1}{The observed importance score }
\item{PERMH0}{ A matrix with the importance scores for each permutation sample in each column}
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
 myarg = list(cp = 0, maxcompete = 0, maxsurrogate = 0, maxdepth = 2)
 Y.names = c("T3" ,"D3")
 P.names = 'Z2'
 T.names = c("Z1", paste("Z", 3:11, sep = ''))
 mybag = 40
 set.seed(5000)
 
 data(burn)
 resperm0 <- permute_select_surv(xdata = burn, Y.names, P.names, T.names, method = "LR", 
        Bag = mybag, args.rpart = myarg, args.parallel = list(numWorkers = 1), nperm = 150)
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ documentation }% use one of  RShowDoc("KEYWORDS")
\keyword{ tree }% __ONLY ONE__ keyword per line
\keyword{multivariate}
\keyword{survival}
\keyword{nonparametric}