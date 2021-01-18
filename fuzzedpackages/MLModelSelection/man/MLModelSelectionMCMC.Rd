\name{MLModelSelectionMCMC}
\alias{MLModelSelectionMCMC}
\title{Model estimation for multivariate longitudinal models.}
\description{Using MCMC procedure to generate posterior samples and provide AIC, BIC, DIC, MPL, MSPE, and predicted values.}
\usage{MLModelSelectionMCMC(Num.of.iterations, list.Data, list.InitialValues, list.HyperPara, 
	list.UpdatePara, list.TuningPara)}
\arguments{
\item{Num.of.iterations}{Number of iterations.}
\item{list.Data}{\emph{List} of data set containing response \eqn{Y}, design matrix \eqn{X}, avialable time points for each subject, GARP model, and ISD model.}
\item{list.InitialValues}{\emph{List} of initial values for parameters.}
\item{list.HyperPara}{\emph{List} of given hyperparameters in priors.}
\item{list.UpdatePara}{Determine which parameter will be updated.}
\item{list.TuningPara}{Provide turning parameters in proposal distributions.}
}

\details{
We set the subject \eqn{i} (\eqn{i=1, \ldots, N}) has \eqn{K} continuous responses at each time point \eqn{t} (\eqn{t=1, \ldots, n_i}). Assume that the measurement times are common across subjects, but not necessarily equally-spaced. Let \eqn{{y}_{it} = (y_{it1}, \ldots, y_{itK})} denote the response vector containing \eqn{K} continuous responses for \eqn{i}th subject at time \eqn{t} along with a \eqn{p\times 1} vectof of covariates, \eqn{{x}_{it} = (x_{it1}, \ldots, x_{itp})}. An efficient Gibbs sampling algorithm is developed for model estimation in the multivariate longitudinal model given by
\deqn{
y_{i1k}  = {x}'_{it}{\beta}_k + e_{i1k}, t=1;
}
\deqn{
y_{itk}  = {x}'_{it}{\beta}_k + \sum_{g=1}^K\sum_{j=1}^{t-1} \phi_{itj, kg} (y_{ijg}-x'_{ij}{\beta}_g)+ e_{itk},  t\geq 2,
}
where \eqn{{\beta}_k = (\beta_{k1}, \ldots, \beta_{kp})'} is a vector of regression coefficients of length \eqn{p}, \eqn{\phi_{itj, kg}} is a generalized autoregressive parameter (GARP) to explain the serial dependence of responses across time. Moreover, 
\deqn{
\phi_{itj, kg} = \alpha_{kg} \mathbf{1}\{|t-j|=1\} ,\;  \log(\sigma_{itk}) = \lambda_{k0} + \lambda_{k1} h_{it}, \; \log\left(\frac{\omega_{ilm}}{\pi-\omega_{ilm}}\right) = \nu_l + \nu_m.
}
The priors for the parameters in the model given by
\deqn{
{\beta} \sim \mathcal{N}(0, \sigma_\beta^2 I);
}
\deqn{
{\lambda}_k \sim \mathcal{N}(0, \sigma_\lambda^2 I);
}
\deqn{
{\nu}_k \sim \mathcal{N}(0, \sigma_\nu^2 I),  \quad k=1, \ldots, K,
}
where \eqn{\sigma_\beta^2}, \eqn{\sigma_\lambda^2}, and \eqn{\sigma_\nu^2} are prespecified values. For \eqn{k, g = 1, \ldots, K} and \eqn{m=1, \ldots, a}, we further assume
\deqn{
\alpha_{kgm} \sim \delta_{kgm} \mathcal{N}(0, \sigma^2_\delta) + (1-\delta_{kgm})\eta_0,
}
where \eqn{\sigma^2_\delta} is prespecified value and \eqn{\eta_0} is the point mass at 0.

}


\value{Lists of posterior samples, parameters estimates, AIC, BIC, DIC, MPL, MSPE, and predicted values are returned}

\references{
Keunbaik Lee et al. (2015) Estimation of covariance matrix of multivariate longitudinal data using modified Choleksky and hypersphere decompositions. \emph{Biometrics}. \bold{75-86, 2020}. \doi{doi:10.1111/biom.13113}.
}

\author{
Kuo-Jung Lee
}
\note{
We'll provide the reference for details of the model and the algorithm for performing model estimation whenever the manuscript is accepted. 
}

\examples{
\donttest{
library(MASS)
library(MLModelSelection)


AR.Order = 6 #denote \phi_{itj, kg} = \alpha_{kg} \mathbf{1}\{|t-j|=1\} 
ISD.Model = 1 #denote \log(\sigma_{itk}) = \lambda_{k0} + \lambda_{k1} h_{it}

data(SimulatedData)

N = dim(SimulatedData$Y)[1] # the number of subjects
T = dim(SimulatedData$Y)[2] # time points
K = dim(SimulatedData$Y)[3] # the number of attributes
P = dim(SimulatedData$X)[3] # the number of covariates
M = AR.Order  # the demension of alpha
nlamb = ISD.Model + 1 # the dimension of lambda

Data = list(Y = SimulatedData$Y, X = SimulatedData$X, 
	TimePointsAvailable = SimulatedData$TimePointsAvailable, 
	AR.Order = AR.Order, ISD.Model = ISD.Model)

beta.ini = matrix(rnorm(P*K), P, K)
delta.ini = array(rbinom(K*K*M, 1, 0.1), c(K, K, M)) 
alpha.ini = array(runif(K*K*M, -1, 1), c(K, K, M))
lambda.ini = matrix(rnorm(nlamb*K), K, nlamb, byrow=T)
nu.ini = rnorm(K)


InitialValues = list(beta = beta.ini, delta = delta.ini, alpha = alpha.ini, 
	lambda = lambda.ini, nu = nu.ini)

# Hyperparameters in priors
sigma2.beta = 1
sigma2.alpha = 10
sigma2.lambda = 0.01
sigma2.nu = 0.01

# Whehter the parameter will be updated
UpdateBeta = TRUE
UpdateDelta = TRUE
UpdateAlpha = TRUE
UpdateLambda = TRUE	
UpdateNu = TRUE


HyperPara = list(sigma2.beta = sigma2.beta, sigma2.alpha=sigma2.alpha, 
	sigma2.lambda=sigma2.lambda, sigma2.nu=sigma2.nu)


UpdatePara = list(UpdateBeta = UpdateBeta, UpdateAlpha = UpdateAlpha, UpdateDelta = UpdateDelta, 
	              UpdateLambda = UpdateLambda, UpdateNu = UpdateNu)

# Tuning parameters in proposal distribution within MCMC
TuningPara = list(TuningAlpha = 0.01, TuningLambda = 0.005, TuningNu = 0.005)

num.of.iter = 100

start.time <- Sys.time()

PosteriorSamplesEstimation = MLModelSelectionMCMC(num.of.iter, Data, InitialValues, 
	HyperPara, UpdatePara, TuningPara)

end.time <- Sys.time()

cat("Estimate of beta\n")
print(PosteriorSamplesEstimation$PosteriorEstimates$beta.mean)

}
}