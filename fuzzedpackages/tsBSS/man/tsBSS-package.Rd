\name{tsBSS-package}
\alias{tsBSS-package}
\docType{package}
\title{
Blind Source Separation and Supervised Dimension Reduction for Time Series
}

\description{
Different estimators are provided to solve the blind source separation problem for multivariate time series with stochastic volatility (Matilainen, Nordhausen and Oja (2015) <doi:10.1016/j.spl.2015.04.033>; Matilainen, Miettinen, Nordhausen, Oja and Taskinen (2017) <doi:10.17713/ajs.v46i3-4.671>) and supervised dimension reduction problem for multivariate time series (Matilainen, Croux, Nordhausen and Oja (2017) <doi:10.1016/j.ecosta.2017.04.002>). Different functions based on AMUSE and SOBI are also provided for estimating the dimension of the white noise subspace. 
}

\details{
\tabular{ll}{
Package: \tab tsBSS\cr
Type: \tab Package\cr
Version: \tab 0.5.6\cr
Date: \tab 2020-11-17\cr
License: \tab GPL (>= 2)
}
This package contains functions for the blind source separation (BSS) problem for multivariate time series. The methods are designed for time series with stochastic volatility, such as GARCH and SV models.
The main functions of the package for the BSS problem are
\itemize{
\item{\code{\link{FixNA}}}  {Function to solve the BSS problem. Algorithm is an alternative to \code{\link{vSOBI}} algorithm to acommodate stochastic volatility.}
\item{\code{\link{gFOBI}}}  {Function to solve the BSS problem. Algorithm is a generalization of \code{\link[JADE]{FOBI}} designed for time series with stochastic volatility.}
\item{\code{\link{gJADE}}}  {Function to solve the BSS problem. Algorithm is a generalization of \code{\link[JADE]{JADE}} designed for time series with stochastic volatility.}
\item{\code{\link{vSOBI}}}  {Function to solve the BSS problem. Algorithm is a variant of \code{\link[JADE]{SOBI}} algorithm and an alternative to \code{\link{FixNA}} to acommodate stochastic volatility.}
\item{\code{\link{gSOBI}}}  {Function to solve the BSS problem. Algorithm is a combination of \code{\link[JADE]{SOBI}} and \code{\link{vSOBI}} algorithms.}
\item{\code{\link{PVC}}}  {Function to solve the BSS problem. Algorithm is a modified version of Principal Component Volatility Analysis by Hu and Tsay (2011).}
}
The input data can be a numeric matrix or a multivariate time series object of class \code{\link{ts}}, \code{\link[xts]{xts}} or \code{\link[zoo]{zoo}}. For other classes, the \code{\link[tsbox:tsbox-package]{tsbox}} package provides appropriate conversions to and from these classes.

The main function of the package for the supervised dimension reduction is
\itemize{
\item{\code{\link{tssdr}}}  {Function for supervised dimension reduction for multivariate time series. Includes methods TSIR, TSAVE and TSSH.}
}

Methods for ARMA models, such as AMUSE and SOBI, and some non-stationary BSS methods for time series are implemented in the \code{\link[JADE:JADE-package]{JADE}} package. See function \code{\link[dr]{dr}} for methods for supervised dimension reduction for iid observations.

Several functions in this package utilize \code{\link[JADE:rjd]{"rjd"}} (real joint diagonalization) and \code{\link[JADE:rjd]{"frjd"}} (fast rjd) from the \code{\link[JADE:JADE-package]{JADE}} package for joint diagonalization of k real-valued square matrices.

There are several functions for estimating the number of white noise latent series in second-order source separation (SOS) models. The functions are 
\itemize{\item{\code{\link{AMUSEboot}}, \code{\link{AMUSEladle}} and  \code{\link{AMUSEasymp}} which are based on \code{\link[JADE]{AMUSE}}.}

\item{\code{\link{SOBIboot}}, \code{\link{SOBIladle}} and  \code{\link{SOBIasymp}} which are based on \code{\link[JADE]{SOBI}}.}
}
Additionally, there is function \code{\link{lbtest}} for a modified Ljung-Box test and a volatility clustering test for univariate and multivariate time series.

The package also contains a dataset \code{\link{WeeklyReturnsData}}, which has logarithmic returns of exchange rates of 7 currencies against US Dollar.
}

\author{
Markus Matilainen, Christophe Croux, Jari Miettinen, Klaus Nordhausen, Hannu Oja, Sara Taskinen, Joni Virta

Maintainer: Markus Matilainen <markus.matilainen@outlook.com> 
}

\references{
Matilainen, M., Nordhausen, K. and Oja, H. (2015), \emph{New Independent Component Analysis Tools for Time Series}, Statistics & Probability Letters, 105, 80--87.

Matilainen, M., Miettinen, J., Nordhausen, K., Oja, H. and Taskinen, S. (2017), \emph{On Independent Component Analysis with Stochastic Volatility Models}, Austrian Journal of Statistics, 46(3--4), 57--66.

Matilainen M., Croux C., Nordhausen K. and Oja H. (2017), \emph{Supervised Dimension Reduction for Multivariate Time Series}, Econometrics and Statistics, 4, 57--69.

Matilainen M., Croux C., Nordhausen K. and Oja H. (2019), \emph{Sliced Average Variance Estimation for Multivariate Time Series}. Statistics: A Journal of Theoretical and Applied Statistics, 53, 630--655.

Shi, Z., Jiang, Z. and Zhou, F. (2009), \emph{Blind Source Separation with Nonlinear Autocorrelation and Non-Gaussianity}, Journal of Computational and Applied Mathematics, 223(1): 908--915.

Matilainen, M., Nordhausen, K. and Virta, J. (2018), \emph{On the Number of Signals in Multivariate Time Series}. In Deville, Y., Gannot, S., Mason, R., Plumbley, M.D. and  Ward, D. (editors) "International Conference on Latent Variable Analysis and Signal Separation", LNCS 10891, 248--258. Springer, Cham., <doi:10.1007/978-3-319-93764-9_24>.

Nordhausen, K. and Virta, J.(2018), \emph{Ladle Estimator for Time Series Signal Dimension}. In 2018 IEEE Statistical Signal Processing Workshop (SSP), pp. 428--432, <doi:10.1109/SSP.2018.8450695>.

Virta, J. and Nordhausen, K. (2019), \emph{Determining the Signal Dimension in Second Order Source Separation}. To appear in Statistica Sinica, <doi:10.5705/ss.202018.0347>.

Miettinen, M., Matilainen, M., Nordhausen, K. and Taskinen, S. (2020), \emph{Extracting Conditionally Heteroskedastic Components Using Independent Component Analysis}, Journal of Time Series Analysis,41, 293--311.

Hu, Y.-P. and Tsay, R. S. (2014), \emph{Principal Volatility Component Analysis}, Journal of Business & Economic Statistics, 32(2), 153--164.
}

\keyword{package}
\keyword{multivariate}
\keyword{ts}