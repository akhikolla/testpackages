\name{acrt-package}

\alias{acrt-package}

\alias{acrt}

\docType{package}

\title{
\packageTitle{acrt}
}
\encoding{latin1}
\description{
The package \pkg{acrt} provides functions for testing affine restrictions on 
the regression coefficient vector in linear models with autocorrelated errors.
The methods implemented in \pkg{acrt} are based on the article 
\enc{Pötscher}{Poetscher} and Preinerstorfer (2016).
In particular, the package can be used to compute various autocorrelation robust 
test statistics; to compute critical values that provide size control
when the error process is Gaussian (or, more generally, elliptically symmetric) 
and autocorrelated; and to compute the size of a test that is 
obtained from an autocorrelation robust test statistic and a user-supplied 
critical value.
}

\details{\pkg{acrt} provides three functions: 
\enumerate{
\item The function \code{F.type.test.statistic} can be used to compute
test statistics of the form \eqn{T_w} or \eqn{T_{E, \mathsf{W}}} as defined in 
\enc{Pötscher}{Poetscher} and Preinerstorfer (2016). The class of 
test statistics of the form \eqn{T_w} or \eqn{T_{E, \mathsf{W}}}
includes F-type tests based on covariance estimators with data-independent
bandwidth parameters and without prewhitening as considered in, e.g., 
Newey and West (1987), Andrews (1991), Kiefer and Vogelsang (2002, 2005), cf.
also Preinerstorfer and \enc{Pötscher}{Poetscher} (2016). 
\item The function \code{critical.value} provides an implementation of 
Algorithm 1 in \enc{Pötscher}{Poetscher} and Preinerstorfer (2016), and can be 
used to determine size-controlling critical values for test statistics of the 
form \eqn{T_w} or \eqn{T_{E, \mathsf{W}}}.
\item The function \code{size} provides an 
implementation of Algorithm 2 in \enc{Pötscher}{Poetscher} and Preinerstorfer 
(2016), and can be used to determine the size of tests that are based on 
test statistics of the form \eqn{T_w} or \eqn{T_{E, \mathsf{W}}} and a 
user-supplied critical value (e.g., obtained from asymptotic theory).
}

We refer the user to \enc{Pötscher}{Poetscher} and Preinerstorfer (2016) for 
details concerning the framework, the test statistics, the algorithms, and the
underlying theoretical results.
}

\references{
Andrews, D. W. K. (1991). Heteroskedasticity and autocorrelation consistent 
covariance matrix estimation. \emph{Econometrica}, \strong{59} 817-858.

Kiefer, N. M. and Vogelsang, T. J. (2002). Heteroskedasticity - autocorrelation 
robust standard errors using the Bartlett kernel without truncation. 
\emph{Econometrica}, \strong{70} 2093-2095.

Kiefer, N. M. and Vogelsang, T. J. (2005). A new asymptotic theory for 
heteroskedasticity - autocorrelation robust tests. \emph{Econometric Theory},
\strong{21} 1130-1164.

Newey, W. K. and West, K. D. (1987). A simple, positive semi-definite, 
heteroskedasticity and autocorrelation consistent covariance matrix. 
\emph{Econometrica}, \strong{55} 703-708.

\enc{Pötscher}{Poetscher}, B.M. and Preinerstorfer, D. (2016). Controlling the 
size of autocorrelation robust tests. \url{https://arxiv.org/abs/1612.06127/}

Preinerstorfer, D. and \enc{Pötscher}{Poetscher}, B. M. (2016). On size and 
power of heteroskedasticity and autocorrelation robust tests. 
\emph{Econometric Theory}, \strong{32} 261-358.
}


