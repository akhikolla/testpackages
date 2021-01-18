\name{Bagg_pred_Surv}
\alias{Bagg_pred_Surv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bagging survival tree prediction 
}
\description{
Use the Bagging improper survival tree to predict on new features and to evaluate the predictor using 
Out Of Bag Integrated Brier Scores with either the Nelson Aalen estimator or the Breslow estimator. A permutation importance score is also computed using OOB observations.
}
\usage{
Bagg_pred_Surv(xdata, Y.names, P.names, resBag, args.parallel = list(numWorkers = 1), 
               new_data = data.frame(), OOB = FALSE)
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
  \item{resBag}{
The result of the \code{\link{Bagg_Surv}} function
}
  \item{args.parallel}{
a list  containing the number of parallel computing arguments: The number of workers, the type of parallelization to achieve, ...  see \code{\link[parallel]{mclapply}} for further details.
}
  \item{new_data}{
An optional data frame to validate the bagging procedure (the test dataset)
}
  \item{OOB}{
A value of \code{TRUE} or \code{FALSE} with \code{TRUE} indicating  the computation of the  OOB error using the Integrated Brier Score and also the computation of the permutation importance score. 
}
}

\value{
\item{PREDNA}{A matrix with Nelson Aalen predictions on all individuals of the learning sample }
\item{PREDBRE}{A matrix with Breslow predictions on all individuals of the learning sample}
\item{tabhazNAa}{A list of matrix with Nelson Aalen  prediction of each tree of the bagging sequence with the leaf node prediction in each column}
\item{tabhazBRe}{A list of matrix with Breslow  prediction of each tree of the bagging sequence with the leaf node prediction in each column}
\item{OOB }{A value of \code{NULL} if \code{OOB} is \code{FALSE}. A list of twelve elements otherwise:
\code{IBSKM}: The Kaplan-Meier estimation of the Integrated Brier Score; \code{IBSNAOOB}: The OOB error using the Nelson-Aalen estimator; \code{IBSBREOOB}: The OOB error using the Breslow estimator; \code{vimpoobpbpna}: The permutation variable importance using the Nelson-Aalen estimator; \code{vimpoobpbpbre}: The permutation variable importance using the Breslow estimator; \code{oobibspbpna}: The mean OOB error predictor by predictor using the Nelson-Aalen estimator; \code{oobibspbpbre}: The mean OOB error predictor by predictor using the Breslow estimator;  \code{SURVNAOOB}: A matrix with the predicted OOB survival using the Nelson-Aalen estimator; \code{SURVBREOOB}: A matrix with the predicted OOB survival using the Breslow estimator; \code{BSTKM}: The vector of Brier scores using the KM estimator; \code{BSTNAOOB}: The vector of Brier scores using the NA estimator; \code{BSTBREOOB}: The vector of Brier scores using the BRE estimator.
}
\item{Timediff}{The execution time of the prediction  procedure}
\item{TEST}{A value of \code{NULL} if \code{new_data} is not available. A list of seven elements otherwise: \code{IBSNAKMnew}: The IBS using the NA estimator on the new dataset; \code{IBSBRKMnew}: The IBS using the BRE estimator on the new dataset; \code{IBSKMnew}: The IBS using the KM estimator on the new dataset; \code{SURVNAnew}: A matrix of predicted survival on the new dataset using the NA estimator;
\code{SURVBREnew}: A matrix of predicted survival on the new dataset using the BRE estimator; \code{SURV_NAnew}: a vector of survival prediction on the testing dataset using the NA estimator; \code{SURV_BREnew}: a vector of survival prediction on the testing dataset using the BRE estimator.
}
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
\code{\link{Bagg_Surv}}
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

pred0 <- Bagg_pred_Surv(burn, Y.names, P.names, burn.BagEssai0, args.parallel =
                        list(numWorkers = 1), OOB = TRUE) 
 
 
 pred1 <- Bagg_pred_Surv(burn, Y.names, P.names, burn.BagEssai1, args.parallel =
                         list(numWorkers = 1), OOB = TRUE) 
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
\keyword{ documentation }% use one of  RShowDoc("KEYWORDS")
\keyword{ tree }% __ONLY ONE__ keyword per line
\keyword{multivariate}
\keyword{survival}
