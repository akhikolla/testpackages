## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  screenshot.force = FALSE, 
  echo = TRUE,
  rows.print = 5,
  message = FALSE, 
  warning = FALSE)

## ----requirement--------------------------------------------------------------
library(PLNmodels)
library(ggplot2)
library(corrplot)

## ----data_load----------------------------------------------------------------
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)

## ----simple PLNPCA------------------------------------------------------------
PCA_models <- PLNPCA(
  Abundance ~ 1 + offset(log(Offset)),
  data  = trichoptera, 
  ranks = 1:5
)

## ----show nocov---------------------------------------------------------------
PCA_models

## ----collection criteria------------------------------------------------------
PCA_models$criteria %>% knitr::kable()

## ----convergence criteria-----------------------------------------------------
PCA_models$convergence  %>% knitr::kable()

## ----plot nocov, fig.width=7, fig.height=5------------------------------------
plot(PCA_models)

## ----model extraction---------------------------------------------------------
myPCA_ICL <- getBestModel(PCA_models, "ICL") 
myPCA_BIC <- getModel(PCA_models, 3) # getBestModel(PCA_models, "BIC")  is equivalent here 

## ----map, fig.width=8, fig.height=8-------------------------------------------
plot(myPCA_ICL, ind_cols = trichoptera$Group)

## ----regression---------------------------------------------------------------
coef(myPCA_ICL) %>% head() %>% knitr::kable()

## ----sigma, fig.width=7-------------------------------------------------------
sigma(myPCA_ICL) %>% corrplot(is.corr = FALSE)

## ----rotation-----------------------------------------------------------------
myPCA_ICL$rotation %>% head() %>% knitr::kable()

## ----scores-------------------------------------------------------------------
myPCA_ICL$scores %>% head() %>% knitr::kable()

## ----show PLNPCAfit-----------------------------------------------------------
myPCA_ICL

## ----cov----------------------------------------------------------------------
PCA_models_cov <- 
  PLNPCA(
    Abundance ~ 1 + offset(log(Offset)) + Temperature + Wind + Cloudiness,
    data  = trichoptera,
    ranks = 1:4
  )

## ----extraction cov, fig.width=7, fig.height=7--------------------------------
plot(PCA_models_cov)
myPCA_cov <- getBestModel(PCA_models_cov, "ICL")

## ----maps, fig.height=4, fig.width=7------------------------------------------
gridExtra::grid.arrange(
  plot(myPCA_cov, map = "individual", ind_cols = trichoptera$Group, plot = FALSE),
  plot(myPCA_cov, map = "variable", plot = FALSE),
  ncol = 2
)

## ----fitted, fig.cap = "fitted value vs. observation", fig.dim=c(7,5)---------
data.frame(
  fitted   = as.vector(fitted(myPCA_cov)),
  observed = as.vector(trichoptera$Abundance)
) %>% 
  ggplot(aes(x = observed, y = fitted)) + 
    geom_point(size = .5, alpha =.25 ) + 
    scale_x_log10(limits = c(1,1000)) + 
    scale_y_log10(limits = c(1,1000)) + 
    theme_bw() + annotation_logticks()

