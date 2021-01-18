## ----setup, include = FALSE---------------------------------------------------
library(knitr)
opts_chunk$set(
  collapse = TRUE,
  comment = ""
)
options(knitr.kable.NA = '')

## ---- warning=FALSE, message=FALSE--------------------------------------------
library(DNAtools)

## -----------------------------------------------------------------------------
data(dbExample, package = "DNAtools")
kable(head(dbExample)[,2:7])

## -----------------------------------------------------------------------------
allele_freqs <- lapply(1:10, function(x){
  al_freq <- table(c(dbExample[[x*2]], dbExample[[1+x*2]]))/(2*nrow(dbExample))
  al_freq[sort.list(as.numeric(names(al_freq)))]
})
names(allele_freqs) <- sub("\\.1", "", names(dbExample)[(1:10)*2])

## ---- results = 'asis', echo=FALSE, fig.align="left"--------------------------
for(l in head(names(allele_freqs))){
  # cat(paste0("**",l,"**\n\n"))
  cat("\n")
  cat(kable(t(c(setNames(NA, paste0(l,":")), allele_freqs[[l]]))), sep = "\n")
}

## -----------------------------------------------------------------------------
m <- 3
P3_D16 <- Pnm_locus(m = m, theta = 0, alleleProbs = allele_freqs$D16S539)

## ---- results='asis'----------------------------------------------------------
cat(kable(t(setNames(P3_D16, paste0("P(n=",1:(2*m), "|m=3)")))), sep = "\n")

## -----------------------------------------------------------------------------
ms <- 1:6
Ps_D16 <- matrix(NA, ncol = length(ms), nrow = max(ms)*2, 
                 dimnames = list(alleles = 1:(max(ms)*2), noContrib = ms))
for(m in seq_along(ms)){
  Ps_D16[1:(ms[m]*2), m] <- Pnm_locus(m = ms[m], theta = 0, alleleProbs = allele_freqs$D16S539)
}

## ---- results='asis'----------------------------------------------------------
dimnames(Ps_D16) <- list(paste0("n=", rownames(Ps_D16),""), 
                         paste0("P(n|m=", colnames(Ps_D16),")"))
kable(Ps_D16, row.names = TRUE)

## ---- fig.width=7, fig.asp=0.62-----------------------------------------------
par(mar = c(4,5,0,0))
plot(c(1,max(ms)*2), range(Ps_D16, na.rm = TRUE), type = "n",
     xlab = "no. of allele", ylab = "Probability")
for (m in seq_along(ms)) {
  lines(1:(2*ms[m]), Ps_D16[1:(ms[m]*2),m], col = m)
}
legend("topright", bty = "n", title = "m", legend = ms, col = seq_along(ms), lty = 1)

## -----------------------------------------------------------------------------
Ps_D16_cumsum <- apply(Ps_D16, 2, cumsum)

## ---- results='asis', echo=FALSE----------------------------------------------
nm <- dim(Ps_D16_cumsum)
dimnames(Ps_D16_cumsum) <- list(paste0("n'=",1:nm[1]), 
                                paste0("P(n&le;n'|m=",1:nm[2],")"))
kable(Ps_D16_cumsum, row.names = TRUE, escape = TRUE)

## ---- results='asis',echo=FALSE-----------------------------------------------
Pm_D16 <- setNames(rep(NA, max(ms)-1), paste0("P(n&le;",(1:(max(ms)-1))*2,"|m=",2:max(ms),")"))
for(m in 2:max(ms)) Pm_D16[m-1] <- Ps_D16_cumsum[2*(m-1),m]
kable(t(Pm_D16))

## -----------------------------------------------------------------------------
no_contrib <- 3
Pnm_all_tab <- Pnm_all(m = no_contrib, theta = 0, probs = allele_freqs, locuswise = TRUE)

## ---- echo = FALSE------------------------------------------------------------
kable(Pnm_all_tab, row.names = TRUE, col.names = paste0("P(n=",1:(2*no_contrib),"|m=",no_contrib,")"))

## -----------------------------------------------------------------------------
locus_Ps <- lapply(ms, Pnm_all, theta = 0, probs = allele_freqs, locuswise = TRUE)
locus_Ps_cumsum <- lapply(locus_Ps, apply, 1, cumsum)
all_Ps_cumsum <- lapply(locus_Ps_cumsum, apply, 1, prod)
all_Ps_table <- matrix(NA, ncol = length(ms), nrow = max(ms)*2, 
                       dimnames = list(n0 = 1:(max(ms)*2), m = ms))
for(m in seq_along(ms)) all_Ps_table[1:(ms[m]*2), m] <- all_Ps_cumsum[[m]]

## ---- echo=FALSE--------------------------------------------------------------
nm <- dim(all_Ps_table)
dimnames(all_Ps_table) <- list(paste0("n~0~=",1:nm[1]), paste0("P(n&le;n~0~|m=",1:nm[2],")"))
kable(all_Ps_table, row.names = TRUE, escape = FALSE)

## -----------------------------------------------------------------------------
Pm_all <- setNames(rep(NA, max(ms)-1), paste0("P(n&le;",(1:(max(ms)-1))*2,"|m=",2:max(ms),")"))
for(m in 2:max(ms)) Pm_all[m-1] <- all_Ps_table[2*(m-1),m]

## ----results='asis',echo=FALSE------------------------------------------------
kable(t(round(Pm_all, 4)))

## ---- fig.width=7, fig.asp=0.62-----------------------------------------------
par(mar = c(4,5,0,0))
P3_D16_t0.03 <- Pnm_locus(m = 3, theta = 0.03, alleleProbs = allele_freqs$D16S539)
P3_D16_t0.1 <- Pnm_locus(m = 3, theta = 0.1, alleleProbs = allele_freqs$D16S539)
plot(P3_D16_t0.03, type = "o", pch = 16, ylab = "Probability", xlab = "no. of alleles", col = 2)
points(P3_D16, type = "o", lty = 2, col = 1)
points(P3_D16_t0.1, type = "o", lty = 3, pch = 2, col = 3)
legend("topleft", bty = "n", title = expression(theta), col = 1:3,
       legend = c("0.00", "0.03", "0.10"), lty = c(2, 1, 3), pch = c(1, 16, 2))

## ---- fig.width=7, fig.asp=0.62-----------------------------------------------
par(mar = c(4,5,0,0))
P3_all <- Pnm_all(m = 3, theta = 0, probs = allele_freqs, locuswise = FALSE)
plot(P3_all, xlab = "no. of alleles", ylab = "Probability", type = "o")

## ---- fig.width=7, fig.asp=0.62, echo=FALSE-----------------------------------
par(mar = c(4,5,0,0))
m_prior <- c(0.2,0.35,0.3,0.2,0.1,0.025)
m_prior <- m_prior/sum(m_prior)
barplot(m_prior, names.arg = seq_along(m_prior), ylab = "Prior probability,"~italic(P(m)), xlab = expression(italic(m)))

## -----------------------------------------------------------------------------
n_vWA <- 4
P_vWA <- lapply(ms, Pnm_locus, theta = 0, alleleProbs = allele_freqs$vWA)
P_vWA_n <- sapply(P_vWA, function(p) ifelse(length(p) < n_vWA, 0, p[n_vWA]))
Pn_vWA <- (m_prior * P_vWA_n)/sum(m_prior * P_vWA_n)

## ---- echo=FALSE--------------------------------------------------------------
kable(rbind("P(m)" = m_prior, "P(n~vWA~=4|m)" = P_vWA_n, "P(m|n~vWA~=4)" = Pn_vWA), col.names = paste0("m=",paste(ms)))

## -----------------------------------------------------------------------------
get_pn <- function(p, n){
  pn <- split(as.data.frame(p), n)
  p_n <- unlist(lapply(names(pn), function(n_l) if(ncol(pn[[n_l]])<n_l) rep(0,nrow(pn[[n_l]])) else pn[[n_l]][,n_l]))
  prod(p_n)
}

n_loci <- rep(4, length(allele_freqs))
P_all <- lapply(ms, Pnm_all, theta = 0, probs = allele_freqs, locuswise = TRUE)
P_all_n <- sapply(P_all, get_pn, n_loci)
Pn_all <- (m_prior * P_all_n)/sum(m_prior * P_all_n)

## ---- echo=FALSE--------------------------------------------------------------
kable(t(setNames(Pn_all, paste0("P(m=", seq_along(Pn_all),"|n=4)"))))

