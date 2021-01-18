library(testthat)
library(MLModelSelection)
context("A simulated study for MLModelSelectionMCMC()")



test_that("A simulation study", {
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

})

