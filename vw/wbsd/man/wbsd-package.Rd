\name{wbsd-package}

\alias{wbsd-package}

\alias{wbsd}

\docType{package}

\title{
\packageTitle{wbsd}
}
\encoding{latin1}
\description{
The package \pkg{wbsd} provides functions for testing affine hypotheses on the regression coefficient vector in linear
regression models with heteroskedastic errors based on the wild bootstrap. 
In particular, it provides (i) the function \code{theta} to compute the value of the
diagnostic ``theta'' suggested in the article \enc{Pötscher}{Poetscher} and 
Preinerstorfer (2020); and it provides (ii) the function \code{boot.pval} to obtain
bootstrap p-values. Various test statistics (null-restricted/unrestricted covariance estimator; classical F-test and F-tests based on HC0 - HC4 covariance estimators) are included.
The functions also cover various bootstrap schemes. We refer the user to \enc{Pötscher}{Poetscher} and Preinerstorfer (2020) for:
details concerning the framework, test statistics, bootstrap schemes,
the definition of theta (depending on the test statistic and scheme used), the
underlying theoretical results, and for further references.
}

\references{
\enc{Pötscher}{Poetscher}, B. M. and Preinerstorfer, D. (2020). How Reliable are Bootstrap-based Heteroskedasticity Robust Tests? <arXiv:2005.04089>
}


