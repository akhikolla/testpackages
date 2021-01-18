\name{iBST-package}
\alias{iBST-package}
\alias{iBST}
\docType{package}
\title{
Improper Bagging Survival Tree
}
\description{
             Fit a bagging survival tree on a mixture of population (susceptible and nonsusceptible)
             using either a pseudo R2 criterion or an adjusted Logrank criterion. The predictor is 
             evaluated using the Out Of Bag Integrated Brier Score (IBS) and several scores of importance
             are computed for variable selection. The thressholds values for variable selection are 
             computed using a nonparametric permutation test.
}
\details{
\tabular{ll}{
Package: \tab iBST\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2017-01-30\cr
License: \tab GPL(>=2.0)\cr
}
}
\author{
Cyprien Mbogning and Philippe Broet

Maintainer: \packageMaintainer{iBST}
}
\references{
Mbogning, C. and Broet, P. (2016). Bagging survival tree procedure for
variable selection and prediction in the presence of nonsusceptible patients.
BMC bioinformatics, 17(1), 1.
}

\keyword{ package }
\keyword{ documentation }% use one of  RShowDoc("KEYWORDS")
\keyword{ tree }% __ONLY ONE__ keyword per line
\keyword{multivariate}
\keyword{survival}

\seealso{
\code{\link{Bagg_Surv}}  \code{\link{Bagg_pred_Surv}}  \code{\link{improper_tree}}
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
 
 ## fit an improper survival tree
 burn.tree <- improper_tree(burn, Y.names, P.names, T.names, method = "R2", args.rpart = myarg)
 plot(burn.tree)
 text(burn.tree, cex = .7, xpd = TRUE)
 
 ## fit an improper Bagging survival tree with the adjusted Logrank criterion
 burn.BagEssai0 <- Bagg_Surv(burn, Y.names, P.names, T.names, method = "LR", args.rpart = myarg, 
                           args.parallel = list(numWorkers = 1), Bag = mybag)
 
 ## fit an improper Bagging survival tree with the pseudo R2 criterion
 burn.BagEssai1 <- Bagg_Surv(burn, Y.names, P.names, T.names, method = "R2", args.rpart = myarg, 
                           args.parallel = list(numWorkers = 1), Bag = mybag)

 ## Plot the variable importance scores
 par(mfrow=c(1,3))
barplot(burn.BagEssai1$IIS, main = 'IIS', horiz = TRUE, las = 1,
        cex.names = .8, col = 'lightblue')
barplot(burn.BagEssai1$DIIS, main = 'DIIS', horiz = TRUE, las = 1,
        cex.names = .8, col = 'grey') 
barplot(burn.BagEssai1$DEPTH, main = 'MinDepth', horiz = TRUE, las = 1,
        cex.names = .8, col = 'purple')


 ## evaluation of the Bagging predictors 
pred0 <- Bagg_pred_Surv(burn, Y.names, P.names, burn.BagEssai0, 
                        args.parallel = list(numWorkers = 1), OOB = TRUE) 
 
 
 pred1 <- Bagg_pred_Surv(burn, Y.names, P.names, burn.BagEssai1, 
                         args.parallel = list(numWorkers = 1), OOB = TRUE) 
  
 ## The OOB integrated Brier score using the Breslow estimator
 pred1$OOB$IBSBREOOB
 
 ## The permutation importance score using the Breslow estimator
 pred1$OOB$vimpoobpbpbre
 }
}
