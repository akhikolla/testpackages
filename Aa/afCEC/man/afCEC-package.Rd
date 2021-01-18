\name{afCEC-package}
\alias{afCEC-package}
\docType{package}
\title{
Active Function Cross-Entropy Clustering
}
\description{
Active Function Cross-Entropy Clustering (afCEC) partitions the n-dimensional data into the clusters by finding the parameters of the mixed generalized multivariate normal distribution, that optimally approximates the scattering of the data in the n-dimensional space, whose density function is of the form:
\cr
\cr
p_1*N(mi_1,^sigma_1,sigma_1,f_1)+...+p_k*N(mi_k,^sigma_k,sigma_k,f_k)
\cr
\cr
The above-mentioned generalization is performed by introducing so called "f-adapted gaussian densities" (i.e. the ordinary gaussian densities adapted by the "active function"). Additionally, the Active Function Cross-Entropy Clustering performs the automatic reduction of the unnecessary clusters. It's implemented in the form of a customizable function \code{\link{afCEC}}.
}
\details{
\tabular{ll}{
Package: \tab afCEC\cr
Type: \tab Package\cr
Version: \tab 1.0.2\cr
Date: \tab 2018-01-04\cr
License: \tab GPL (>= 2) \cr
}
}
\author{
Krzysztof Byrski
}

\seealso{
\code{\link{afCEC}.}
}
