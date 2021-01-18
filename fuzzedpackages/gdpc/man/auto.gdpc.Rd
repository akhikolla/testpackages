\name{auto.gdpc}
\alias{auto.gdpc}
\title{Automatic Fitting of Generalized Dynamic Principal Components}
\description{
Computes Generalized Dynamic Principal Components. The number of components can be supplied by the user or chosen automatically so that a given proportion of variance is explained. The number of lags is chosen automatically using one of the following criteria: Leave-one-out cross-validation, an AIC type criterion, a BIC type criterion or a criterion based on a proposal of Bai and Ng (2002). See Peña, Smucler and Yohai (2020) for more details. 
}


\usage{auto.gdpc(Z, crit = 'LOO', normalize = 1, auto_comp = TRUE, expl_var = 0.9,
          num_comp = 5, tol = 1e-4, k_max = 10,
          niter_max = 500, ncores = 1, verbose = FALSE)}
\arguments{
  \item{Z}{Data matrix. Each column is a different time series.}
  \item{crit}{A string specifying the criterion to be used. Options are 'LOO', 'AIC', 'BIC' and 'BNG'. Default is 'LOO'. See Details below.}
  \item{normalize}{Integer. Either 1, 2 or 3. Indicates whether the data should be standardized. Default is 1. See Details below.}
  \item{auto_comp}{Logical. If TRUE compute components until the proportion of explained variance is equal to expl_var, otherwise use num_comp components. Default is TRUE.}
  \item{expl_var}{A number between 0 and 1. Desired proportion of explained variance (only used if auto_comp==TRUE). Default is 0.9.}
  \item{num_comp}{Integer. Number of components to be computed (only used if auto_comp==FALSE). Default is 5.}
  \item{tol}{Relative precision. Default is 1e-4.}
  \item{k_max}{Integer. Maximum possible number of lags. Default is 10.}
  \item{niter_max}{Integer. Maximum number of iterations. Default is 500.}
  \item{ncores}{Integer. Number of cores to be used for parallel computations. Default is 1.}
  \item{verbose}{Logical. Should progress be reported? Default is FALSE.}
}

\details{

Suppose the data matrix consists of \eqn{m} series of length \eqn{T}.
Let \eqn{\bold{f}} be the dynamic principal component defined using \eqn{k} lags, let \eqn{R} be the corresponding matrix of residuals and let \eqn{\Sigma = (R^{\prime} R) / T}.

If crit = 'LOO' the number of lags is chosen among \eqn{0,\dots, k_{max}} as the value \eqn{k} that minimizes the leave-one-out (LOO) cross-validation mean squared error, given by
\deqn{ LOO = \frac{1}{T m}\sum\limits_{i=1}^{m}\sum\limits_{t=1}^{T}\frac{R_{t,i}^{2}}{(1-h_{t,t})^{2}},}
where \eqn{h_{t,t}} are the diagonal elements of the hat matrix \eqn{H = F(F^{\prime} F)^{-1} F^{\prime} }, with \eqn{F} being the \eqn{T \times (k+2)} matrix with rows \eqn{(f_{t-k}, f_{t-k+1}, \dots, f_{t}, 1)}.

If crit = 'AIC' the number of lags is chosen among \eqn{0,\dots, k_{max}} as the value \eqn{k} that minimizes the following AIC type criterion
\deqn{ AIC = T \log(trace(\Sigma)) + 2 m (k+2) .}

If crit = 'BIC' the number of lags is chosen among \eqn{0,\dots, k_{max}} as the value \eqn{k} that minimizes the following BIC type criterion
\deqn{ BIC = T \log(trace(\Sigma)) + m (k+2) \log(T) .}

If crit = 'BNG' the number of lags is chosen among \eqn{0,\dots, k_{max}} as the value \eqn{k} that minimizes the following criterion
\deqn{ BNG = \min(T, m) \log(trace(\Sigma)) + (k+1) \log(\min(T, m)).}
This is an adaptation of a criterion proposed by Bai and Ng (2002).

For problems of relatively small dimension, say \eqn{T \geq m 10}, 'AIC' can can give better results than the
default 'LOO'.

If normalize = 1, the data is analyzed in the original units, without mean and variance standarization. If normalize = 2, the data is standardized to zero mean and unit variance before computing the principal components, but the intercepts and loadings are those needed to reconstruct the original series. If normalize = 3 the data are standardized as in normalize = 2, but the intercepts and the loadings are those needed to reconstruct the standardized series. Default is normalize = 1.

}

\value{
  An object of class \code{gdpcs}, that is, a list of length equal to the number of computed components. The i-th entry of this list is an object of class \code{gdpc}, that is, a list with entries
   \item{expart}{Proportion of the variance explained by the first i components.}
  \item{mse}{Mean squared error of the reconstruction using the first i components.}
  \item{crit}{The value of the criterion of the reconstruction, according to what the user specified.}
  \item{k}{Number of lags chosen.}
  \item{alpha}{Vector of intercepts corresponding to f.}
  \item{beta}{Matrix of loadings corresponding to f. Column number \eqn{k} is the vector of \eqn{k-1} lag loadings.}
  \item{f}{Coordinates of the i-th dynamic principal component corresponding to the periods \eqn{1,\dots,T}.}
  \item{initial_f}{Coordinates of the i-th dynamic principal component corresponding to the periods \eqn{-k+1,\dots,0}. Only for the case \eqn{k>0}, otherwise 0.}
      \item{call}{The matched call.}
      \item{conv}{Logical. Did the iterations converge?}
      \item{niter}{Integer. Number of iterations.}
  \code{components}, \code{fitted}, \code{plot} and \code{print} methods are available for this class.
}
\author{
  Daniel Peña, Ezequiel Smucler, Victor Yohai
}

\references{
Bai J. and Ng S. (2002). “Determining the Number of Factors in Approximate Factor Models.”
Econometrica, 70(1), 191–221.

Peña D., Smucler E. and Yohai V.J. (2020). “{gdpc}: An {R} Package for Generalized Dynamic Principal Components.” Journal of Statistical Software, to appear.
}

\seealso{
\code{\link{gdpc}}, \code{\link{plot.gdpc}}, \code{\link{plot.gdpcs}}, \code{\link{fitted.gdpcs}}, \code{\link{components.gdpcs}}
}
\examples{
T <- 200 #length of series
m <- 200 #number of series
set.seed(1234)
f <- rnorm(T + 1)
x <- matrix(0, T, m)
u <- matrix(rnorm(T * m), T, m)
for (i in 1:m) {
    x[, i] <- 10 * sin(2 * pi * (i/m)) * f[1:T] + 10 * cos(2 * pi * (i/m)) * f[2:(T + 1)] + u[, i]
}
#Choose number of lags using the LOO criterion.
#k_max=3 to keep computation time low
autofit <- auto.gdpc(x, k_max = 3)
autofit
fit_val <- fitted(autofit, 1) #Get fitted values
resid <- x - fit_val #Residuals
plot(autofit, which_comp = 1) #Plot component
}

\keyword{ts}