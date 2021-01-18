## ----knit-int, echo=FALSE, include=FALSE---------------------------------
library("knitr")
set.seed(5828993)
opts_knit$set(self.contained = TRUE)

## ----load-packages, include=TRUE-----------------------------------------
#install.packages("GMCM")  # Uncomment to install the GMCM package
library("GMCM")

## ----load-data, include=TRUE, echo=TRUE----------------------------------
x <- get(data("u133VsExon"))
head(x, n = 4)

## ------------------------------------------------------------------------
x <- x[1:5000, ]

## ----preprocess-data-----------------------------------------------------
u <- Uhat(1 - x)

## ------------------------------------------------------------------------
par(mfrow = c(1,2))  # Visualizing P-values and the ranked P-values
plot(x, cex = 0.5, pch = 4, col = "tomato", main = "P-values",
     xlab = "P-value (U133)", ylab = "P-value (Exon)")
plot(u, cex = 0.5, pch = 4, col = "tomato", main = "Ranked P-values",
     xlab = "rank(1-P) (U133)", ylab = "rank(1-P) (Exon)")

## ----show-initial-params, include=TRUE, echo=TRUE------------------------
init_par <- c(pie1 = 0.6, mu = 1, sigma = 1, rho = 0.2)

## ----fit_model, error=TRUE-----------------------------------------------
par <- fit.meta.GMCM(u = u,
                     init.par = init_par,
                     method = "NM",
                     max.ite = 1000,
                     verbose = FALSE)
print(par)

## ----compute_probs-------------------------------------------------------
meta_IDR_thres <- 0.05
out <- get.IDR(x, par = par,
               threshold = meta_IDR_thres) # Compute IDR
str(out)

out <- get.IDR(u, par = par, threshold = meta_IDR_thres)
below <- out[["IDR"]] < meta_IDR_thres
out$l <- sum(below)
out$Khat <- ifelse(below, 2, 1)

## ----classes_table-------------------------------------------------------
table(out$Khat)

## ----plot_results--------------------------------------------------------
plot(x, col = out$Khat, asp = 1) # Plot of raw values
plot(u, col = out$Khat, asp = 1) # Plot of copula values
z <- GMCM:::qgmm.marginal(u, theta = meta2full(par, d = ncol(u))) # Estimate latent process
plot(z, col = out$Khat, asp = 1) # Plot of estimated latent process

## ----plot_theta----------------------------------------------------------
plot(meta2full(par, d = ncol(u)))

## ----session-info--------------------------------------------------------
sessionInfo()

## ----citation, echo=FALSE, results='asis'--------------------------------
cites <- lapply(c("GMCM", "knitr", "rmarkdown"), citation)
fmt_cites <- unlist(lapply(cites, format, style = "text"))
cat(paste0("  ", seq_along(fmt_cites), ". ", fmt_cites, "\n"))

