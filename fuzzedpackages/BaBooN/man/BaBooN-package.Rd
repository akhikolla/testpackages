\name{BaBooN-package}
\alias{BaBooN-package}

\docType{package}
\title{
Package for multiple imputation of missing values based on Bayesian Bootstrap with Predictive Mean Matching.
}
\description{
Included are two variants of Bayesian Bootstrap Predictive Mean Matching to multiply impute missing data. The first variant is a variable-by-variable imputation combining sequential regression and Predictive Mean Matching (\acronym{PMM}) that has been extended for unordered categorical data. The Bayesian Bootstrap allows for generating approximately proper multiple imputations. The second variant is also based on PMM, but the focus is on imputing several variables at the same time. The suggestion is to use this variant, if the missing-data pattern resembles a data fusion situation, or any other missing-by-design pattern, where several variables have identical missing-data patterns. Both variants can be run as 'single imputation' versions, in case the analysis objective is of a purely descriptive nature. 
}
\details{
\tabular{ll}{
Package: \tab BaBooN\cr
Type: \tab Package\cr
Version: \tab 0.2-0\cr
Date: \tab 2015-06-15\cr
License: \tab GPL (>= 2)\cr
}

}
\author{
Florian Meinfelder [aut, cre] <florian.meinfelder[AT]uni-bamberg.de>\cr
Thorsten Schnapp [aut] <thorsten.schnapp[AT]uni-bamberg.de>\cr
Maintainer: Florian Meinfelder <florian.meinfelder[AT]uni-bamberg.de>
}
\references{

  Koller-Meinfelder, F. (2009) \emph{Analysis of Incomplete Survey Data
    -- Multiple Imputation Via Bayesian Bootstrap Predictive Mean
    Matching}, doctoral thesis. 
    
  Cowles, M.K. and Carlin, B.P. (1996) Markov Chain Monte Carlo Convergence Diagnostics: A Comparative Review.
  \emph{Journal of the American Statistical Association}, Vol. \bold{91}, pp. 883--904.

  Eddelbuettel, D. and Francois, R. (2011) Rcpp: Seamless R and C++ Integration.
  \emph{Journal of Statistical Software}, Vol. \bold{40}, No. 8, pp. 1--18. URL http://www.jstatsoft.org/v40/i08/.
  
  Eddelbuettel, D. and Sanderson, C. (2014) RcppArmadillo: Accelerating R with high-performance C++ linear algebra.
    \emph{Computational Statistics and Data Analysis}, Vol. \bold{71}, March 2014, pp. 1054--1063.
  
  Harrell, F.E., with contributions from Charles Dupont and many others. (2013) \emph{Hmisc: Harrell
    Miscellaneous}. R package version 3.13-0. http://CRAN.R-project.org/package=Hmisc
  
  Little, R.J.A. (1988) Missing-Data Adjustments in Large
    Surveys, \emph{Journal of Business and Economic Statistics}, Vol. \bold{6},
    No. 3, pp. 287-296.
	
  Plummer, M. and Best, N. and  Cowles, K. and  Vines, K. (2006) CODA: Convergence 
    Diagnosis and Output Analysis for MCMC, \emph{R News}, Vol. \bold{6}, pp. 7--11
	
  Ported to R by Alvaro A. Novo. Original by Joseph L. Schafer <jls@stat.psu.edu>. (2013). \emph{norm:
	Analysis of multivariate normal datasets with missing values}. R package version 1.0-9.5.
	http://CRAN.R-project.org/package=norm

  R Core Team (2015) \emph{R: A language and environment for statistical computing. R
    Foundation for Statistical Computing}, Vienna, Austria. URL
    http://www.R-project.org/.
  
  Raghunathan T.E. and Lepkowski, J.M. and Van Hoewyk, J. and Solenberger, P (2001) A multivariate technique
    for multiply imputing missing values using a sequence of regression models. \emph{Survey Methodology}, Vol. \bold{27}, pp. 85--95.
  
  Rubin DB (1981) The Bayesian Bootstrap. \emph{The Annals of Statistics}, Vol. \bold{9}, pp. 130--134.
  
	 Rubin, D.B. (1987) \emph{Multiple Imputation for Non-Response in Surveys}. New York: John Wiley & Sons, Inc.
  
  Sarkar, Deepayan (2008) \emph{Lattice: Multivariate Data Visualization with R}. New York: Springer.
    
  Van Buuren, S. (2012) \emph{Flexible imputation of missing data}. Boca Raton: CRC Press.
  
  Van Buuren, S. and Brand, J.P.L. and Groothuis-Oudshoorn, C.G.M. and Rubin, D.B. (2006) Fully conditional specification
  in multivariate imputation. \emph{Journal of Statistical Computation and Simulation}, Vol. \bold{76}, No. 12, pp. 1049--1064.
  
  Van Buuren, S. and Groothuis-Oudshoorn, K. (2011) mice: Multivariate Imputation by Chained Equations in R. \emph{Journal of
    Statistical Software}, Vol. \bold{45}, No. 3, pp. 1--67.
    URL http://www.jstatsoft.org/v45/i03/.
   
  Venables, W. N. and Ripley, B. D. (2002) \emph{Modern Applied Statistics with S}. Fourth Edition. New York: Springer.
  
}

\keyword{package}

