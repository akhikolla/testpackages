\name{flars}
\alias{flars}
\title{
Functional least angle regression.
}
\description{
This is the main function for the functional least angle regression algorithm. Under certain conditions, the function only needs the input of two arguments: \code{x} and \code{y}. This function can do both variable selection and parameter estimation. 
}
\usage{
flars(x,y,method=c('basis','gq','raw'),max_selection,cv=c('gcv'),
      normalize=c('trace','rank','norm','raw'),lasso=TRUE,check=1,
      select=TRUE,VarThreshold=0.1,SignThreshold=0.8,
      control=list())
}
\arguments{
  \item{x}{
The mixed scalar and functional variables. Note that each of the functional variables is expected to be stored in a matrix. Each row of the matrix should represent a sample or a curve. If there is only one functional variable, \code{x} can be a matrix. If there is only scalar variables, \code{x} can be a vector or a matrix. If there are more than one functional variables, or there are mixed functional and scalar variables, \code{x} should be a list. If \code{x} is a list, each item of the list should correspond to one variable. 
}
  \item{y}{
The scalar variable. It can be a matrix or a vector. 
}
  \item{method}{
The representative methods for the functional coefficients. The method could be one of the 'basis', 'gq' and 'raw' for basis function expression, Gaussian quadrature and representative data points, respectively. 
}
\item{max_selection}{
Number of maximum selections when stopping the algorithm. Set a reasonable number for this argument to increase the calculation speed.
}
\item{cv}{
Choise of cross validation. At the moment, the only choice is the generalized cross validation, i.e., \code{cv='gcv'}.
}
\item{lasso}{
Use lasso modification or not. In other words, can variables selected in the former iterations be removed in the later iterations. 
}
\item{check}{
Type of check methods for lasso modification. 1 means variance check, 2 means sign check. \code{check=1} is much better than the other one. 
}
\item{select}{
If \code{TRUE}, the aim is to do selection rather than parameter estimation, and the stopping rule can be used when \code{lasso=TRUE}. If \code{FALSE}, the stopping rule may not work when \code{lasso=TRUE}.
}
\item{VarThreshold}{
Threshold for removing variables based on variation explained. More specifically, one condition to remove a variable is that the variation explained by a variable is less than \code{VarThreshold*Var(y)}. To remove this variable, there is another condition: the variation explained by this variable is less than largest variation it explained in the previous iterations. 
}
\item{SignThreshold}{
This is a similar argument to \code{VarThreshold}. If a functional coefficient has less than \code{SignThreshold} same as that from the previous iteration, the variable is removed. 
}
\item{normalize}{
Choice of normalization methods. This is to remove any effects caused by the different dimensions of functional variables and scalar variables. Currently we have \code{trace}, \code{rank}, \code{norm}, \code{raw}. \code{norm} and \code{raw} are recommended.
}
\item{control}{
list of control elements for the functional coefficients. See \code{\link{fccaGen}} for details.
}
}

\value{
\item{Mu}{Estimated intercept from each of the iterations}
\item{Beta}{Estimated functional coefficients from each of the iterations}
\item{alpha}{Distance along the directions from each of the iterations}
\item{p2_norm}{Normalization constant applied to each of the iterations}
\item{AllIndex}{All the index. If one variable is removed, it will become a negative index.}
\item{index}{All the index at the end of the selection.}
\item{CD}{Stopping rule designed for this algorithm. The algorithm should stop when this value is very small. Normally we can observe an obvious and severe drop of the value.}
\item{resid}{Residual from each of the iteration.}
\item{RowMeans}{Point-wise mean of the functional variables and mean of the scalar variables.}
\item{RowSds}{Point-wise sd of the functional variables and sd of the scalar variables.}
\item{yMean}{Mean of the response variable.}
\item{ySD}{SD of the response variable.}
\item{p0}{The projections obtained from each iteration without normalization.}
\item{cor1}{The maximum correlation obtained from the first iteration.}
\item{lasso}{Weather have lasso step or not.}
\item{df}{The degrees of freedom calculated at the end of each iteration.}
\item{Sigma2Bar}{Estimated $sigma^2$.}
\item{StopStat}{Conventional stopping criteria.}
\item{varSplit}{The variation explained by each of the candidate variables at each iteration.}
\item{SignCheckF}{The proportion of sign changing for each of the candidate variables at each iteration.}
}

\examples{
library(flars)
library(fda)
#### Ex1 ####
## Generate some data.
dataL=data_generation(seed = 1,uncorr = TRUE,nVar = 8,nsamples = 120,
    var_type = 'm',cor_type = 3)

## Do the variable selection
out=flars(dataL$x,dataL$y,method='basis',max_selection=9,
    normalize='norm',lasso=FALSE)

## Check the stopping point with CD
plot(2:length(out$alpha),out$CD) # plot the CD with the iteration number

## In simple problems we can try
(iter=which.max(diff(out$CD))+2)


#### Ex2 ####
## Generate some data.
# dataL=data_generation(seed = 1,uncorr = FALSE,nVar = 8,nsamples = 120,
#      var_type = 'm',cor_type = 3)
## add more variables to the candidate
# for(i in 2:4){
# dataL0=data_generation(seed = i,uncorr = FALSE,nVar = 8,nsamples = 120,
#      var_type = 'm',cor_type = 3) 
# dataL$x=c(dataL$x,dataL0$x)
# }
# names(dataL$x)=paste0('v_',seq(length(dataL$x)))

## Do the variable selection
# out=flars(dataL$x,dataL$y,method='basis',max_selection=9,
#     normalize='norm',lasso=FALSE)

#### Ex3 (small subset of a real data set) ####
data(RealDa, package = 'flars')
out=flars(RealDa$x,RealDa$y,method='basis',max_selection=9,
    normalize='norm',lasso=FALSE)
# out=flars(RealDa$x,RealDa$y,method='basis',max_selection=9,
#     normalize='norm',lasso=TRUE)

## Check the stopping point with CD
plot(2:length(out$alpha),out$CD) # plot the CD with the iteration number
## The value drops to very small compare to others at iteration six and 
###  stays low after that, so the algorithm may stop there. 


}
