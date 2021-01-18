## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
## loading and looking at the data
welfare = read.csv2('welfare.csv')
head(welfare)

## ----setup--------------------------------------------------------------------
library(hdpGLM)

## ---- results='hide'----------------------------------------------------------
## estimating the model
mcmc = list(burn.in=10, ## MCMC burn-in period
            n.iter =500) ## number of MCMC iterations to keep
mod = hdpGLM(support ~ inequality + income + ideology, data=welfare,
             mcmc=mcmc)


## -----------------------------------------------------------------------------
## printing the outcome
summary(mod)

## -----------------------------------------------------------------------------
welfare_clustered = hdpGLM_classify(welfare, mod)
head(welfare_clustered)
tail(welfare_clustered)

## ---- fig.width=7.2, fig.height=5---------------------------------------------
plot(mod, separate=T, ncols=4)

## -----------------------------------------------------------------------------
## loading and looking at the data
welfare = read.csv2('welfare2.csv')
head(welfare)
tail(welfare)

## ---- results='hide'----------------------------------------------------------
## estimating the model
mcmc = list(burn.in=1, ## MCMC burn-in period
            n.iter =50) ## number of MCMC iterations to keep
mod = hdpGLM(support ~ inequality + income + ideology, 
             support ~ gap,
	     data=welfare, mcmc=mcmc)

## -----------------------------------------------------------------------------
summary(mod)

## ---- fig.width=7.2, fig.height=7---------------------------------------------
plot_tau(mod)

## ---- fig.width=7.2, fig.height=5---------------------------------------------
plot_pexp_beta(mod, smooth.line=TRUE, ncol.beta=2)

