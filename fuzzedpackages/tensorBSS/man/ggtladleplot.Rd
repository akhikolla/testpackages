\name{ggtladleplot}
\alias{ggtladleplot}

\title{
Ladle plot for each mode of an object of class tladle using ggplot2
}
\description{
The ladle plot is a measure for deciding about the number of interesting components. Of interest for the ladle criterion is the minimum.
The function here offers however also to plot other criterion values which are part of the actual ladle criterion. 
}
\usage{
ggtladleplot(x, crit = "gn", type="l", ylab = crit,  
  xlab = "component", main = deparse(substitute(x)), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class ladle.}
  \item{crit}{the criterion to be plotted, options are \code{"gn"}, \code{"fn"}, \code{"phin"} and \code{"lambda"}.}
  \item{type}{plotting type, either lines \code{l} or points \code{p}.}
  \item{ylab}{default ylab value.}
  \item{xlab}{default xlab value.}
  \item{main}{default title.}
  \item{\dots}{other arguments for the plotting functions.}
}
\details{
The main criterion of the ladle is the scaled sum of the eigenvalues and the measure of variation of the eigenvectors up to the component of interest.

The sum is denoted \code{"gn"} and the individual parts are \code{"fn"} for the measure of the eigenvector variation and \code{"phin"} for the scaled eigenvalues.
The last option \code{"lambda"} corresponds to the unscaled eigenvalues yielding then a screeplot.

The plot is drawn separately for each mode of the data.
}

\references{
\cite{Koesner, C, Nordhausen, K. and Virta, J. (2019), Estimating the signal tensor dimension using tensorial PCA. Manuscript.}

\cite{Luo, W. and Li, B. (2016), Combining Eigenvalues and Variation of Eigenvectors for Order Determination, Biometrika, 103. 875--887. <doi:10.1093/biomet/asw051>}
}
\author{
Klaus Nordhausen, Joni Virta
}

\seealso{
\code{\link{tPCAladle}}
}
\examples{
library(ICtest)
n <- 500
sig <- 0.6

Z <- rbind(sqrt(0.7)*rt(n,df=5)*sqrt(3/5),
           sqrt(0.3)*runif(n,-sqrt(3),sqrt(3)),
           sqrt(0.3)*(rchisq(n,df=3)-3)/sqrt(6),
           sqrt(0.9)*(rexp(n)-1),
           sqrt(0.1)*rlogis(n,0,sqrt(3)/pi),
           sqrt(0.5)*(rbeta(n,2,2)-0.5)*sqrt(20)
)

dim(Z) <- c(3, 2, n)

U1 <- rorth(12)[,1:3]
U2 <- rorth(8)[,1:2]
U <- list(U1=U1, U2=U2)
Y <- tensorTransform2(Z,U,1:2)
EPS <- array(rnorm(12*8*n, mean=0, sd=sig), dim=c(12,8,n))
X <- Y + EPS


TEST <- tPCAladle(X, n.boot = 100)
TEST
ggtladleplot(TEST, crit = "gn")
ggtladleplot(TEST, crit = "fn")
ggtladleplot(TEST, crit = "phin")
ggtladleplot(TEST, crit = "lambda")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }% use one of  RShowDoc("KEYWORDS")
