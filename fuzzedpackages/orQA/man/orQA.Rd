\name{orQA}
\docType{package}
\title{
  Order restricted quality assessment of microarray titration data
}
\description{
  orQA provides methods for quality assessment microarray titration
  data. These include permutation based trend and shape tests useful for
  the assessment of accuracy and cross platform consistency as well as
  methods for the estimation of variance components under order
  restrictions.
}
\details{
orQA provides methods for quality assessment of microarray titration
data. Exploiting the monotonic nature of such measurements accuracy,
precision and cross-platform agreement can be derived. For the assessment
of accuracy using shape tests see \code{\link{pttest}}. To evaluate
precision using variance component estimates see
\code{\link{est.lme}}. Order restricted inference of monotonic trends
for the purpose of cross-platform comparison is provided by
\code{\link{e2test}}. For an example see below.
}
\author{
Florian Klinglmueller

Maintainer: Florian Klinglmueller <float_at_lefant.net>
}
\references{
    Klinglmueller, F., Tuechler, T., Posch, M. (2010)
  "Cross Platform Comparison Of Microarray Data Using Order
    Restricted Inference" Under Review

  Pinheiro, J., Bates, D., DebRoy, S., Sarkar, D.,R Development Core
    Team (2010)
    "nlme: Linear and Nonlinear Mixed Effects Models"

  Barlow, R. E., Bartholomew, D. J., Bremner, J. M., and Brunk, H. D. (1972)
    "Statistical inference under order restrictions"; Wiley, London.

  Robertson, T., Wright,F. T. and Dykstra, R. L. (1988)
    "Order Restricted Statistical Inference"; Wiley, New York.

  Guo W., Sarkar SK., Peddada SD. (2010)
    "Controlling False Discoveries in Multidimensional Directional
    Decisions, with Applications to Gene Expression Data on Ordered
    Categories" Biometrics
}
\keyword{ package }
\examples{
## Let's assume a titration study with 10 genes 4 titration groups and
## 10 replicates which come from either of 2 biologically different
## samples (i.e. two times 4 technical replicates)

## WARNING: examples are run with a very limited number of permutations 

## Titration groups, random factor
g <- rep(1:4,each=8)
r <- rep(rep(1:2,each=4),4)

## No differences (global null)
nulldata <- matrix(rnorm(320),nc=32)
## Differences between titration levels in each gene
altdata <- t(t(nulldata)+g)

## Accuracy are there any significantly non monotonous trends

res <- pttest(nulldata,g,1000,r) # apply shape test

## no significant trends at all
sigdirPttest(res)

## with alternatives
res <- pttest(altdata,g,1000,r) # apply shape test


## some significant trends 
sigdirPttest(res)

## type of monotonicity (no significant trend, up, down,
##   "anti-monotonous")
table(monotonicity(res))

## Precision - estimate variance components

res <- est.lme(nulldata,g,r)
round(apply(res,2,summary),2)
res <- est.lme(altdata,g,r)
round(apply(res,2,summary),2)


## Trend test and agreement
## some other dataset with some alternatives in both directions
tdir <- sample((1:3)-2,10,rep=TRUE)
altdata2 <- matrix(rnorm(320),nc=32)+ (tdir \%*\% t(g))
res1 <- e2test(altdata,g,1000,r)
res2 <- e2test(altdata2,g,1000,r)

## trinomial coding for directions (-1 down, 0 non sig., 1 up)
sigdir1 <- sigdirE2test(res1)
sigdir2 <- sigdirE2test(res2)

## contingency table of directional decisions
table(sigdir1,sigdir2)

}
