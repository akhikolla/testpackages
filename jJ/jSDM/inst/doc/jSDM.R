## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(
	fig.align = "center",
	fig.width = 6, fig.height = 6,
	cache = FALSE,
	collapse = TRUE,
	comment = "#>",
	highlight = TRUE
)

## ----frog-picture, echo=FALSE, out.width=500, out.height=400, fig.cap="(ref:cap-frog)"----
knitr::include_graphics("figures/Litoria_ewingii.jpg")

## ----libraries-----------------------------------------------------------
# Load libraries
library(jSDM)

## ----frogs-data----------------------------------------------------------
# frogs data
data(frogs, package="jSDM")
head(frogs)

## ----arranging-data------------------------------------------------------
# data.obs
PA_frogs <- frogs[,4:12]

# Normalized continuous variables
Env_frogs <- cbind(scale(frogs[,1]),frogs[,2],scale(frogs[,3]))
colnames(Env_frogs) <- colnames(frogs[,1:3])

## ----jSDM, cache=TRUE----------------------------------------------------
mod_jSDM_block_frogs <- jSDM_probit_block (
  # Response variable 
  presence_site_sp = as.matrix(PA_frogs), 
  # Explanatory variables 
  site_suitability = ~.,   
  site_data = as.data.frame(Env_frogs), n_latent=2,
  # Chains
  burnin=20000, mcmc=5000, thin=5,
  # Starting values
  alpha_start=0, beta_start=0,
  lambda_start=0, W_start=0,
  V_alpha_start=1, 
  # Priors
  shape=0.5, rate=0.0005,
  mu_beta=0, V_beta=1.0E6,
  mu_lambda=0, V_lambda=10,
  # Various 
  seed=1234, verbose=1)

## ----plot-results--------------------------------------------------------
## alpha_i of the first two sites
plot(coda::as.mcmc(mod_jSDM_block_frogs$mcmc.alpha[,1:2]))

## Valpha
par(mfrow=c(1,2))
coda::traceplot(mod_jSDM_block_frogs$mcmc.Valpha, main="V_alpha")
coda::densplot(mod_jSDM_block_frogs$mcmc.Valpha, main="V_alpha")

np <- nrow(mod_jSDM_block_frogs$model_spec$beta_start)

## beta_j of the first two species
par(mfrow=c(np,2))
for (j in 1:2) {
  for (p in 1:np) {
      coda::traceplot(coda::as.mcmc(mod_jSDM_block_frogs$mcmc.sp[[paste0("sp_",j)]][,p]))
      coda::densplot(coda::as.mcmc(mod_jSDM_block_frogs$mcmc.sp[[paste0("sp_",j)]][,p]), 
main = paste(colnames(mod_jSDM_block_frogs$mcmc.sp[[paste0("sp_",j)]])[p],
", species : ",j))
  }
}

## lambda_j of the first two species
n_latent <- mod_jSDM_block_frogs$model_spec$n_latent
par(mfrow=c(n_latent*2,2))
for (j in 1:2) {
  for (l in 1:n_latent) {
      coda::traceplot(coda::as.mcmc(mod_jSDM_block_frogs$mcmc.sp[[paste0("sp_",j)]][,np+l]))
      coda::densplot(coda::as.mcmc(mod_jSDM_block_frogs$mcmc.sp[[paste0("sp_",j)]][,np+l]), 
      main = paste(colnames(mod_jSDM_block_frogs$mcmc.sp[[paste0("sp_",j)]])[np+l],
      ", species : ",j))
    }
  }

## Latent variables W_i for the first two sites
par(mfrow=c(1,2))
for (l in 1:n_latent) {
  plot(mod_jSDM_block_frogs$mcmc.latent[[paste0("lv_",l)]][,1:2], main = paste0("Latent variable W_", l))
}

## probit_theta
par (mfrow=c(1,1))
hist(mod_jSDM_block_frogs$probit_theta_pred, main = "Predicted probit theta", xlab ="predicted probit theta")

## Deviance
plot(mod_jSDM_block_frogs$mcmc.Deviance,main = "Deviance")


## ----correlation-matrix--------------------------------------------------
plot_residual_cor(mod_jSDM_block_frogs)

## ----predictions---------------------------------------------------------
# Sites and species concerned by predictions :
## 50 sites among the 104
Id_sites <- sample.int(nrow(PA_frogs), 50)
## All species 
Id_species <- colnames(PA_frogs)
# Simulate new observations of covariates on those sites 
simdata <- matrix(nrow=50, ncol = ncol(mod_jSDM_block_frogs$model_spec$site_data))
colnames(simdata) <- colnames(mod_jSDM_block_frogs$model_spec$site_data)
rownames(simdata) <- Id_sites
simdata <- as.data.frame(simdata)
simdata$Covariate_1 <- rnorm(50)
simdata$Covariate_3 <- rnorm(50)
simdata$Covariate_2 <- rbinom(50,1,0.5)

# Predictions 
theta_pred <- predict.jSDM(mod_jSDM_block_frogs, newdata=simdata, Id_species=Id_species,
						   Id_sites=Id_sites, type="mean")
hist(theta_pred, main="Predicted theta with simulated data", xlab="predicted theta")

