\name{PimpTest}
\alias{PimpTest}
\alias{PimpTest.default}
\alias{print.PimpTest}
\title{
PIMP testing approach
}
\description{
Uses permutations to approximate the null importance distributions for all variables and computes the p-values based on the null importance distribution according to the approach of Altmann et al. (2010).
}
\usage{
\method{PimpTest}{default}(Pimp, para = FALSE, ...)
\method{print}{PimpTest}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Pimp}{ an object of class \code{\link[vita]{PIMP}} }
  \item{para}{ If para is TRUE  the null importance distributions are approximated with Gaussian
                distributions else with empirical cumulative distributions. Default is \code{ para = FALSE}}
 \item{\dots}{ optional parameters, not used }
  \item{x}{for the print method, an \code{PimpTest} object }
}
\details{
The vector \code{perVarImp} of \eqn{S} variable importance measures for every predictor variables from code {\link[vita]{PIMP}} are used to approximate the null importance distributions.
If \code{para} is \code{TRUE} this implementation of the PIMP algorithm fits for each variable a \emph{Gaussian distribution} to the \eqn{S} null importances. If \code{para} is \code{FALSE} the PIMP algorithm uses the empirical distribution of the \eqn{S} null importances.
Given the fitted null importance distribution, the p-value is the probability of observing the \emph{original VarImp} or a larger value.
}
\value{
\item{VarImp }{ the \emph{original permutation variable importance} measures of the random forest. }
\item{PerVarImp }{a matrix, where the l-th row contains the \code{S} permuted VarImp
                  measures for the l-th predictor variable. }
\item{para }{ Was the null distribution approximated by a Gaussian distribution or by the empirical distribution? }
\item{meanPerVarImp }{mean for each row of \code{PerVarImp}. \code{NULL} if \code{ para = FALSE} }
\item{sdPerVarImp }{ standard deviation for each row of \code{PerVarImp}.\code{NULL} if \code{ para = FALSE} }
\item{p.ks.test }{the p-values of the Kolmogorov-Smirnov Tests for each row \code{PerVarImp}. Is the
                null importance distribution significantly different from a normal distribution with the mean(PerVarImp) and
                sd(PerVarImp)? \code{NULL} if \code{ para = FALSE} }
\item{pvalue }{ the p-value is the probability of observing the \code{original VarImp} or a larger value, given the fitted null importance distribution.}
}
\references{
 Breiman L. (2001), \emph{Random Forests}, Machine Learning 45(1),5-32, <doi:10.1023/A:1010933404324>

Altmann A.,Tolosi L., Sander O. and  Lengauer T. (2010),\emph{Permutation importance: a corrected feature importance measure}, Bioinformatics Volume 26 (10), 1340-1347, <doi:10.1093/bioinformatics/btq134>
}

\seealso{
 \code{\link[vita]{PIMP}},  \code{\link[vita]{summary.PimpTest}}
}
\examples{
\donttest{
###############################
#      Regression            #
##############################

## Simulating data
X = replicate(15,rnorm(100))
X = data.frame(X) #"X" can also be a matrix
y = with(X,2*X1 + 1*X2 + 2*X3 + 1*X4 - 2*X5 - 1*X6 - 1*X7 + 2*X8 )

##############################
## Regression with Random Forest:
library("randomForest")
reg.rf = randomForest(X,y,mtry = 3,ntree=500,importance=TRUE)
##############################
## PIMP-Permutation variable importance measure

system.time(pimp.varImp.reg<-PIMP(X,y,reg.rf,S=100, parallel=TRUE, ncores=2))
pimp.t.reg = PimpTest(pimp.varImp.reg)
summary(pimp.t.reg,pless = 0.1)

##############################
#      Classification        #
##############################

## Simulating data
X = replicate(10,rnorm(200))
X= data.frame( X) #"X" can also be a matrix
z  = with(X,2*X1 + 3*X2 + 2*X3 + 1*X4 -
            2*X5 - 2*X6 - 2*X7 + 1*X8 )
pr = 1/(1+exp(-z))         # pass through an inv-logit function
y = as.factor(rbinom(200,1,pr))

##############################
## Classification with Random Forest:
cl.rf = randomForest(X,y,mtry = 3,ntree = 500, importance = TRUE)
##############################
## PIMP-Permutation variable importance measure
system.time(pimp.varImp.cl<-PIMP(X,y,cl.rf,S=100, parallel=TRUE, ncores=2))
pimp.t.cl = PimpTest(pimp.varImp.cl,para = TRUE)
summary(pimp.t.cl,pless = 0.1)
}
}

