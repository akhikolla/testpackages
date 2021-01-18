\name{readData}
\alias{readData}
\alias{EditIn.data}
\alias{EditIn.data-class}
\docType{class}
\title{Read Data and Edit Rules}
\description{
add description about this function
}
\usage{
readData(Y.original, ratio = NULL, range = NULL, balance = NULL, eps.bal = 0.6) 
}
\arguments{
	\item{Y.original}{original dataset of (\eqn{n}, \eqn{p}) dimension with missing and edit-failing values where \eqn{n} is the number of records and \eqn{p} is the number of variables.}
	
	\item{ratio}{ratio edit.}
	
	\item{range}{range restriction.}
		
	\item{balance}{balance edit.}
			
	\item{eps.bal}{threshold for balance edit. Defaults to 0.6.}
}
\details{
\code{Y.original} has \eqn{n} records and \eqn{p} variables. The variable names (column names) of \code{Y.original} are used to specify \code{ratio} edits.

The edit rules are either imported from text files or written by \pkg{editrules} package's syntax. 

A balance edit is considered as two inequality constraints with the threshold, i.e., `A = B' is converted to `-eps.bal < A - B < eps.bal' before computation. 

For accurate computation, nested balances are written as `total variable = sum of component variables'. For example, it is recommended to replace `X1 = X2 + X3' and `X3 = X4 + X5' with `X1 = X2 + X4 + X5' and `X3 = X4 + X5' so that `X3' does not appear both sides of the balance edits.
}
\value{
readData returns an \code{\link{EditIn.data}} object which consists of
\itemize{
 	\item \code{Y.input}: input dataset which replaces NA in \code{Y.original} with -999 and zero values with 0.01.
  \item \code{Edit.editmatrix}: \code{\link{editmatrix}} of edit rules. It can be used for functions of \pkg{editrules} package.
	\item \code{Edit.matrix}:\code{matrix} of edit rules. 
	\item \code{Bound.LU}: range restrictions. For variable X whose range is not specified in \code{range}, the default values are set as max( 0.1min(X), 1e-5 ) for the lower bound and 10max(x) for the upper bound.
	\item \code{ratio}: ratio edits.
	\item \code{n.balance}: number of balance edit, i.e., the row number of \code{balance}.
	\item \code{FaultyRecordID}: record IDs of \code{Y.orig} whose values violate edit rules.
}
}
\references{
Hang J. Kim, Lawrence H. Cox, Alan F. Karr, Jerome P. Reiter and Quanli Wang (2015). "Simultaneous Edit-Imputation for Continuous Microdata", Journal of the American Statistical Association, DOI: 10.1080/01621459.2015.1040881.	\cr\cr
Edwin de Jonge and Mark van der Loo (2013). editrules: R package for parsing, applying, and
  manipulating data cleaning rules. R package version 2.7.2.
}
\seealso{\code{\link{editmatrix}}} 
\examples{
### option 1. import from text files ###

data(NestedEx)
 
D_obs1 = NestedEx$D.obs
Ratio1 = NestedEx$Ratio.edit
Range1 = NestedEx$Range.edit
Balance1 = NestedEx$Balance.edit

data1 = readData(Y.original=D_obs1, ratio=Ratio1, range=Range1, 
balance=Balance1, eps.bal=0.6)

# print(data1$Edit.editmatrix)
# plot(data1$Edit.editmatrix)	  ## function of 'editrules' package

### option 2. Using the syntax of R package 'editrules' ###

data(NestedEx) ; D_obs2 = NestedEx$D.obs

Ratio2 <- editmatrix(c(
 "X1 <= 1096.63*X5", "X1 <= 2980.96*X7", "X1 <= 148.41*X8", "X1 <= 7.39*X9",
 "X5 <= 0.37*X1", "X5 <= 54.60*X7", "X5 <= 2.72*X8", "X5 <= 0.14*X9",
 "X7 <= 0.14*X1", "X7 <= 1.65*X5", "X7 <= 7.39*X8", "X7 <= 0.05*X9",
 "X8 <= 1.65*X1", "X8 <= 54.60*X5", "X8 <= 403.43*X7", "X8 <= 1.65*X9",
 "X9 <= 20.09*X1", "X9 <= 403.43*X5", "X9 <= 13359.73*X7", "X9 <= 148.41*X8"
))
Range2 <- editmatrix(c(
 "X1 >= 2", "X2 <= 1.2e+06", "X11 >= 0.002", "X11 <= 1.2e+04"
))
Balance2 <- editmatrix(c(
 "X1 == X2+X3+X4", "X5 == X6 + 0.4*X10 + 0.6*X11", "X7 == 0.4*X10 + 0.6*X11"
))

data2 = readData(D_obs2, Ratio2, Range2, Balance2)

# print(data2$Edit.editmatrix)
  # Note: data2 is equivalent to data1
# plot(data2$Edit.editmatrix)	  ## function of 'editrules' package
}
% \keyword{internal}