\name{MTS-package}
\alias{MTS-package}
\alias{MTS}
\docType{package}
\title{
Multivariate Time Series
}
\description{
Multivariate Time Series (MTS) is a general package for analyzing multivariate linear time series and estimating multivariate volatility models. It also handles factor models, constrained factor models, asymptotic principal component analysis commonly used in finance and econometrics, and principal volatility component analysis.  (a) For the multivariate linear time series analysis, the package performs model specification, estimation, model checking, and prediction for many widely used models, including vector AR models, vector MA models, vector ARMA models, seasonal vector ARMA models, VAR models with exogenous variables, multivariate regression models with time series errors, augmented VAR models, and Error-correction VAR models for co-integrated time series. For model specification, the package performs structural specification to overcome the difficulties of identifiability of VARMA models. The methods used for structural specification include Kronecker indices and Scalar Component Models.  (b) For multivariate volatility modeling, the MTS package handles several commonly used models, including multivariate exponentially weighted moving-average volatility, Cholesky decomposition volatility models, dynamic conditional correlation (DCC) models, copula-based volatility models, and low-dimensional BEKK models. The package also considers multiple tests for conditional heteroscedasticity, including rank-based statistics.  (c) Finally, the MTS package also performs forecasting using diffusion index, transfer function analysis, Bayesian estimation of VAR models, and multivariate time series analysis with missing values.Users can also use the package to simulate VARMA models, to compute impulse response functions of a fitted VARMA model, and to calculate theoretical cross-covariance matrices of a given VARMA model. 
}
\details{
\tabular{ll}{
Package: \tab MTS\cr
Type: \tab Package\cr
License: \tab Artistic License 2.0\cr
}
}
\author{
Ruey S. Tsay and David Wood
}