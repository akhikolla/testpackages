\name{s2net-package}
\alias{s2net}
\docType{package}
\title{\packageTitle{s2net}}
\description{ 
\if{html}{\figure{s2net.png}{options: width=200 align='right'}}
  \packageDescription{s2net}
}
\details{
  The DESCRIPTION file: \packageDESCRIPTION{s2net}
  \packageIndices{s2net}

  This package includes a very easy-to-use interface for handling data, with the \code{\link{s2Data}} function. The main function of the package is the \code{\link{s2netR}} function, which is a wrapper for the \code{\link{Rcpp_s2net}} (\code{s2net}) class. 
  }
\author{
  \packageAuthor{s2net}
}
\references{
Laria, J.C., L. Clemmensen (2019). A generalized elastic-net for semi-supervised learning of sparse features.

Sogaard Larsen, J. et. al. (2019).  Semi-supervised covariate shift modelling of spectroscopic data.

Ryan, K. J., & Culp, M. V. (2015). On semi-supervised linear regression in covariate shift problems. \emph{The Journal of Machine Learning Research, 16}(1), 3183-3217.
}
% Optionally other standard keywords, one per line,
% from the file KEYWORDS in the R documentation.
\keyword{package}
\seealso{
  \code{\link{s2Data}},
  \code{\link{s2netR}},
  \code{\link{Rcpp_s2net}}
}
\examples{
data("auto_mpg")
train = s2Data(xL = auto_mpg$P1$xL, yL = auto_mpg$P1$yL,  xU = auto_mpg$P1$xU)

model = s2netR(train, 
                s2Params(lambda1 = 0.1, 
                           lambda2 = 0,
                           gamma1 = 0.1,
                           gamma2 = 100,
                           gamma3 = 0.1))

# here we tell it to transform the valid data as we did with train.
valid = s2Data(auto_mpg$P1$xU, auto_mpg$P1$yU, preprocess = train) 
ypred = predict(model, valid$xL)

\dontrun{
if(require(ggplot2)){
  ggplot() + 
    aes(x = ypred, y = valid$yL) + geom_point() + 
    geom_abline(intercept = 0, slope = 1, linetype = 2)
}
}
}
