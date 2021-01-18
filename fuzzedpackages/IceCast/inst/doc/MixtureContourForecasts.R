## ----load package,  message = FALSE, warning = FALSE---------------------
library("IceCast")

## ---- fit weights--------------------------------------------------------
weight <- fit_weights(mod1 = clim_9_2005_2007, mod2 = ppe_9_2005_2007,
                 obs = obs_9_2005_2007, prop_area = prop_area) 

## ---- apply weights, fig.height = 5, fig.width  = 5----------------------
MCF_9_2008 <- wght_mod(w = weight, mod1 = clim_9_2008, mod2 = ppe_9_2008)

## ---- load packages,  message = FALSE, warning = FALSE-------------------
library("fields")
library("viridis")

## ---- fig.height = 10, fig.width  = 6------------------------------------
par(mfrow = c(2, 1), oma = c(0, 0, 0, 4))
image.plot(MCF_9_2008, main = "Mixture Contour Forecast, September 2008", 
           col = viridis(8), zlim = c(0, 1), xaxt = "n", yaxt = "n")
image.plot(obs_9_2008, main = "Observed Sea Ice, September 2008", 
           col = viridis(8), zlim = c(0, 1), xaxt = "n", yaxt = "n")

