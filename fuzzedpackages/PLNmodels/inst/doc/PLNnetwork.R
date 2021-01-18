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

## ----data_load----------------------------------------------------------------
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)

## ----simple PLNnetwork--------------------------------------------------------
network_models <- PLNnetwork(Abundance ~ 1 + offset(log(Offset)), data = trichoptera)

## ----show---------------------------------------------------------------------
network_models

## ----collection criteria------------------------------------------------------
network_models$criteria %>% head() %>% knitr::kable()

## ----convergence criteria-----------------------------------------------------
network_models$convergence %>% head() %>% knitr::kable()

## ----diagnostic, fig.width=7, fig.height=5------------------------------------
plot(network_models, "diagnostic")

## ----plot, fig.width=7, fig.height=5------------------------------------------
plot(network_models)

## ----path_coeff, fig.width=7, fig.height=7------------------------------------
coefficient_path(network_models, corr = TRUE) %>% 
  ggplot(aes(x = Penalty, y = Coeff, group = Edge, colour = Edge)) + 
    geom_line(show.legend = FALSE) +  coord_trans(x="log10") + theme_bw()

## ----extract models-----------------------------------------------------------
model_pen <- getModel(network_models, network_models$penalties[20]) # give some sparsity
model_BIC <- getBestModel(network_models, "BIC")   # if no criteria is specified, the best BIC is used

## ----extract models stars-----------------------------------------------------
model_StARS <- getBestModel(network_models, "StARS") # if StARS is requested, stability selection is performed if needed 

## ----plot stability, fig.width=7, fig.height=5--------------------------------
plot(network_models, "stability")

## ----show/print---------------------------------------------------------------
model_StARS

## ----extract------------------------------------------------------------------
my_graph <- plot(model_StARS, plot = FALSE)
my_graph

## ----stars_network, fig.width=7, fig.height=7---------------------------------
plot(model_StARS)
plot(model_StARS, type = "support", output = "corrplot")

## ----fitted, fig.cap = "fitted value vs. observation", fig.dim=c(7,5)---------
data.frame(
  fitted   = as.vector(fitted(model_StARS)),
  observed = as.vector(trichoptera$Abundance)
) %>% 
  ggplot(aes(x = observed, y = fitted)) + 
    geom_point(size = .5, alpha =.25 ) + 
    scale_x_log10(limits = c(1,1000)) + 
    scale_y_log10(limits = c(1,1000)) + 
    theme_bw() + annotation_logticks()

