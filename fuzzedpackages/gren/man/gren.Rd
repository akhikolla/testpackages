\name{gren}
\alias{gren}
\title{
Group-regularized logistic elastic net regression
}
\description{
Function that estimates a group-regularized elastic net model.
}
\usage{
gren(x, y, m=rep(1, nrow(x)), unpenalized=NULL, partitions=NULL, alpha=0.5, 
     lambda=NULL, intercept=TRUE, monotone=NULL, psel=TRUE, compare=TRUE, 
     posterior=FALSE, nfolds=nrow(x), foldid=NULL, trace=TRUE,
     init=list(lambdag=NULL, mu=NULL, sigma=NULL, chi=NULL, ci=NULL),
     control=list(epsilon=0.001, maxit=500, maxit.opt=1000, maxit.vb=100))
}

\arguments{
\item{x}{
feature data as either \code{numeric} \code{matrix} or \code{data.frame} of \code{numeric} variables.}
\item{y}{
response as either a \code{numeric} with binomial/binary successes of length \code{nrow(x)} or a \code{matrix} of \code{nrow(x)} rows and two columns, where the first column contains the binomial/binary failures and the second column the binomial/binary successes.}
\item{m}{
\code{numeric} of length \code{nrow(x)} that contains the number of Bernoulli trials.}
\item{unpenalized}{
Optional \code{numeric} \code{matrix} or \code{data.frame} of \code{numeric} unpenalized covariates of \code{nrow(x)} rows.}
\item{partitions}{
\code{list} that contains the (possibly multiple) partitions of the data. Every \code{list} object corresponds to one partition, where every partition is a \code{numeric} of length \code{ncol(x)} containing the group ids of the features.}
\item{alpha}{
proportion of L1 penalty as a \code{numeric} of length 1.}
\item{lambda}{
global penalty parameter. The default \code{NULL} will result in estimation by cross-validation.}
\item{intercept}{
\code{logical} to indicate whether an intercept should be included.}
\item{monotone}{
\code{list} of two \code{logical} vectors of length \code{length(partitions)}. The first one \code{monotone} indicates whether the corresponding partition's penalty parameters should be monotonically estimates, the second vector \code{decreasing} indicates whether the monotone penalty parameters are decreasing with group number.}
\item{psel}{
either a \code{numeric} vector that indicates the number of features to select or a \code{logical}. If \code{TRUE} feature selection is done by letting \code{\link{glmnet}} determine the penalty parameter sequence.}
\item{compare}{
\code{logical}, if \code{TRUE}, a regular non-group-regularized model is estimated.}
\item{posterior}{
if \code{TRUE}, the full variational Bayes posterior is returned.}
\item{nfolds}{
\code{numeric} of length 1 with the number of folds used in the cross-validation of the global \code{lambda}. The default is \code{nrow(x)}.}
\item{foldid}{
optional \code{numeric} vector of length \code{nrow(x)} with the fold assignments of the observations.}
\item{trace}{
if \code{TRUE}, progress of the algorithm is printed.}
\item{init}{
optional \code{list} containing the starting values of the iterative algorithm. See Details for more information.}
\item{control}{
a \code{list} of algorithm control parameters. See Details for more information.}
}

\details{
This is the main function of the package that estimates a group-regularized elastic net regression. The elastic net penalty's proportion of L1-norm penalisation is determined by \code{alpha}. \code{alpha} close to 0 implies more ridge-like penalty, while \code{alpha} close to 1 implies lasso-like penalty. The algorithm is a two-step procedure: first, a global lambda penalty is estimates by cross-validation. Next, the groupwise lambda multipliers are estimates by an EM algorithm. The EM algorithm consists of: i) an expectation step in which the expected marginal likelihood of the penalty multipliers is iteratively approximated by a variational Bayes EM algorithm and ii) a maximisation step in which the approximate expected marginal likelihood is maximised with respect to the penalty multipliers. After convergence of the algorithm an (optional) frequentist elastic net model is fit using the estimated penalty multipliers by setting \code{psel=TRUE} or by setting \code{psel} to a \code{numeric} vector. 

The user may speed up the procedure by specifying initial values for the EM algorithm in \code{init}. \code{init} is a \code{list} that contains:
\describe{
\item{\code{lambdag}}{
initial values for \eqn{\lambda_g} in a \code{list} of length \code{length(partitions)}.
}
\item{\code{mu}}{
initial values for the \eqn{\mu_j} in a \code{numeric} vector of length \code{ncol(x) + ncol(unpenalized) + intercept}.
}
\item{\code{chi}}{
initial values for the \eqn{\chi_j} in a \code{numeric} vector of length \code{ncol(x)}.
}
\item{\code{ci}}{
initial values for the \eqn{c_i} in a \code{numeric} vector of length \code{nrow(x)}.
}
\item{\code{sigma}}{
The initial values for the \eqn{\Sigma_{ij}} in a \code{matrix} of \code{numeric}s with \code{ncol(x)} rows and columns.
}
}
\code{control} is a \code{list} with parameters to control the estimation procedure. It consists of the following components:
\describe{
\item{\code{epsilon}}{
\code{numeric} with the relative convergence tolerance. Default is \code{epsilon=0.001}.
}
\item{\code{maxit}}{
\code{numeric} with whole number that gives the maximum number of iterations to update the \code{lambdag}. Default is \code{maxit=500}.
}
\item{\code{maxit.opt}}{
\code{numeric} with whole number that gives the maximum number of iterations to numerically maximise the \code{lambdag}. Maximisation occurs at every iteration. Default is \code{maxit.opt=1000}.
}
\item{\code{maxit.vb}}{
\code{numeric} with whole number that gives the maximum number of iterations to update the variational parameters \code{mu}, \code{sigma}, \code{chi}, and \code{ci}. One full update sequence per iteration. Default is \code{maxit=100}.
}
}
}

\value{
Function returns an S3 \code{list} object of class \code{gren} containing output with the following components:
\item{\code{call}}{
The function call that produced the output.
}
\item{\code{alpha}}{
proportion of L1 penalty as a \code{numeric} of length 1.
}
\item{\code{lambda}}{
global penalty parameter as \code{numeric}. Estimated by cross-validation if \code{lambda=NULL}.
}
\item{\code{lambdag.seq}}{
\code{list} with full sequence of penalty multipliers over iterations.
}
\item{\code{lambdag}}{
\code{list} with final estimates of penalty multipliers.
}
\item{\code{vb.post}}{
\code{list} with variational posterior parameters \eqn{mu_j}, \eqn{sigma_{ij}}, \eqn{c_i}, and \eqn{chi_j}.
}
\item{\code{freq.model}}{
frequentist elastic net model as output of \code{glmnet} call. \code{NULL} if \code{psel=FALSE}.
}
\item{\code{iter}}{
\code{list} with number of iterations of \code{lambdag} estimation, with number of optimisation iterations of \code{lambdag}, and number of variational Bayes iterations.
}
\item{\code{conv}}{
\code{list} of \code{logical}s with convergence of \code{lambdag} sequence, optimisation steps, and variational Bayes iterations.
}
\item{\code{args}}{
\code{list} with input arguments of \code{gren} call.
}
}

\references{
Münch, M.M., Peeters, C.F.W., van der Vaart, A.W., and van de Wiel, M.A. (2018). Adaptive group-regularized logistic elastic net regression. arXiv:1805.00389v1 [stat.ME].
}

\author{
Magnus M. Münch <m.munch@vumc.nl>
}

\seealso{
\code{\link{predict.gren}}, \code{\link{coef.gren}}, \code{\link{cv.gren}}
}

\examples{
\donttest{
## Create data
p <- 1000
n <- 100
set.seed(2018)
x <- matrix(rnorm(n*p), ncol=p, nrow=n)
beta <- c(rnorm(p/2, 0, 0.1), rnorm(p/2, 0, 1))
m <- rep(1, n)
y <- rbinom(n, m, as.numeric(1/(1 + exp(-x \%*\% as.matrix(beta)))))
partitions <- list(groups=rep(c(1, 2), each=p/2))

## estimate model
fit.gren <- gren(x, y, m, partitions=partitions)
}
}