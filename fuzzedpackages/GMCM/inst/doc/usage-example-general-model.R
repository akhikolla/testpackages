## ----knit-int, echo=FALSE, include=FALSE---------------------------------
set.seed(7869670)
knitr::opts_knit$set(self.contained = TRUE)

## ----load-packages, include=TRUE-----------------------------------------
#install.packages("GMCM")  # Uncomment to install the GMCM package
library("GMCM")

## ----sim-data, include=TRUE, echo=TRUE-----------------------------------
true_theta <- rtheta(m = 3, d = 2)
plot(true_theta)
ds <- SimulateGMCMData(n = 1000, theta = true_theta)
str(ds)

## ----data-plot-----------------------------------------------------------
plot(ds$u, col = ds$K)

## ----select-data, include=TRUE, echo=TRUE--------------------------------
uhat <- Uhat(ds$u)
plot(uhat)

## ----show-initial-params, include=TRUE, echo=TRUE------------------------
init_theta <- choose.theta(uhat, m = 3)
print(init_theta)

## ----fit_model, error=TRUE-----------------------------------------------
est_theta <- fit.full.GMCM(u = uhat,  # Ranking function is applied automatically
                           theta = init_theta,
                           method = "NM",
                           max.ite = 5000,
                           verbose = FALSE)
print(est_theta)

## ----compute_probs-------------------------------------------------------
membership <- classify(uhat, est_theta)
str(membership)

## ----post_prob-----------------------------------------------------------
post_prob <- get.prob(uhat, est_theta)  # Compute component probabilities

## ----classes_table-------------------------------------------------------
table(membership)

## ----plot_results--------------------------------------------------------
par(mfrow = c(1,2))
plot(uhat, col = membership, asp = 1) # Plot of estimated copula values
z <- GMCM:::qgmm.marginal(uhat, theta = est_theta) # Estimate latent process
plot(z, col = membership, asp = 1) # Plot of estimated latent process

## ----plot_theta----------------------------------------------------------
plot(est_theta)

## ----session-info--------------------------------------------------------
sessionInfo()

## ----citation, echo=FALSE, results='asis'--------------------------------
cites <- lapply(c("GMCM", "knitr", "rmarkdown"), citation)
fmt_cites <- unlist(lapply(cites, format, style = "text"))
cat(paste0("  ", seq_along(fmt_cites), ". ", fmt_cites, "\n"))

