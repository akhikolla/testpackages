## ----knitr_options, include=FALSE---------------------------------------------
knitr::opts_chunk$set(fig.width=7, fig.height=4.5, dev="svg")
old_digits <- options("digits")$digits
options(digits=3)

## ----load_package-------------------------------------------------------------
library(lineup2)
gastroc <- lineup2ex$gastroc
islet <- lineup2ex$islet

## ----align_rows---------------------------------------------------------------
aligned <- align_matrix_rows(gastroc, islet)

## ----align_cols---------------------------------------------------------------
aligned <- align_matrix_cols(aligned[[1]], aligned[[2]])

## ----calc_paired_corr---------------------------------------------------------
paired_corr <- corr_betw_matrices(aligned[[1]], aligned[[2]])

## ----to_reset_par, include=FALSE----------------------------------------------
old_mar <- par("mar")
old_mfrow <- par("mfrow")
old_las <- par("las")

## ----plot_paired_corr---------------------------------------------------------
par(mar=c(4.1, 5.1, 1.1, 1.1), las=1)
plot(sort(paired_corr), ylab="Correlations between column pairs, sorted", pch=16)
abline(h=0, lty=2, col="gray60")

## ----select_top_genes---------------------------------------------------------
selected_genes <- which(paired_corr > 0.75)

## ----sample_correlations------------------------------------------------------
corr_betw_samples <- corr_betw_matrices(t(gastroc[,selected_genes]),
                                        t(islet[,selected_genes]), what="all")

## ----hist_self_nonself, fig.height=6, dev="png"-------------------------------
hist_self_nonself(corr_betw_samples, xlabel="correlation")

## ----self_values--------------------------------------------------------------
self <- get_self(corr_betw_samples)

## ----missing_self-------------------------------------------------------------
self[is.na(self)]

## ----self_sorted--------------------------------------------------------------
sort(self)[1:8]

## ----best_by_row_and_col------------------------------------------------------
best_byrow <- get_best(corr_betw_samples, get_min=FALSE, dimension="row")
best_bycol <- get_best(corr_betw_samples, get_min=FALSE, dimension="column")

## ----plot_self_v_best---------------------------------------------------------
par(mfrow=c(1,2), mar=c(4.1, 4.1, 1.6, 0.6), las=1)
plot(self, best_byrow, xlab="Self-self correlation", ylab="Best islet correlation",
     main="Gastroc samples", pch=16, xlim=c(-0.3, 1.0), ylim=c(0.6, 1.0))
label <- best_byrow-self > 0.2
text(self[label]+0.03, best_byrow[label], names(self)[label],
     adj=c(0,0.5), cex=0.8)

plot(self, best_bycol, xlab="Self-self correlation", ylab="Best gastroc correlation",
     main="Islet samples", pch=16, xlim=c(-0.3, 1.0), ylim=c(0.6, 1.0))
label <- best_bycol-self > 0.3
text(self[label]+0.03, best_bycol[label], names(self)[label],
     adj=c(0,0.5), cex=0.8)
label <- best_bycol-self > 0.2 & best_bycol-self < 0.3
text(self[label]-0.03, best_bycol[label], names(self)[label],
     adj=c(1,0.5), cex=0.8)

## ----problems_byrow-----------------------------------------------------------
get_problems(corr_betw_samples, threshold=0.2, get_min=FALSE)

## ----problems_bycol-----------------------------------------------------------
get_problems(corr_betw_samples, threshold=0.2, get_min=FALSE, dimension="column")

## ----which_best_islet3296-----------------------------------------------------
which_best(corr_betw_samples, get_min=FALSE, dimension="column")["Mouse3296"]

## ----plot_six_samples, fig.height=6.5-----------------------------------------
samples <- sort(c(names(self)[best_bycol-self > 0.2], "Mouse3295"))
par(mfrow=c(3,2), las=1, mar=c(4.1, 4.1, 2.1, 0.6))
rn <- rownames(corr_betw_samples)
green <- "#2ecc40"
red <- "#ff4136"
for(sample in samples) {
    plot(corr_betw_samples[,sample], xlab="Index", ylab="Correlation",
         main=sample, pch=16, ylim=range(corr_betw_samples))
    mx <- which.max(corr_betw_samples[,sample])
    points(mx, corr_betw_samples[mx,sample], pch=16, col=red)
    text(mx-8, corr_betw_samples[mx, sample], names(mx), adj=c(1, 0.5))
    points(which(rn==sample), corr_betw_samples[sample,sample], pch=16, col=green)
}

## ----second_best_by_row_and_col-----------------------------------------------
secbest_byrow <- get_2ndbest(corr_betw_samples, get_min=FALSE, dimension="row")
secbest_bycol <- get_2ndbest(corr_betw_samples, get_min=FALSE, dimension="column")

## ----sec_best_vs_best---------------------------------------------------------
par(mfrow=c(1,2), mar=c(4.1, 4.1, 1.6, 0.6), las=1)
plot(best_byrow, secbest_byrow, xlab="Best islet correlation",
     ylab="Second best islet correlation",
     main="Gastroc samples", pch=16, xlim=c(0.5, 1), ylim=c(0.5, 1))
label <- best_byrow-self > 0.2
points(best_byrow[label], secbest_byrow[label], pch=16, col=red)
abline(0,1, lty=2)

plot(best_bycol, secbest_bycol, xlab="Best gastroc correlation",
     ylab="Second best gastroc correlation",
     main="Islet samples", pch=16, xlim=c(0.5, 1), ylim=c(0.5, 1))
label <- best_bycol-self > 0.2
points(best_bycol[label], secbest_bycol[label], pch=16, col=red)
abline(0,1, lty=2)

## ----reset_par_and_options, include=FALSE-------------------------------------
par(mar=old_mar, mfrow=old_mfrow, las=old_las)
options(digits=old_digits)

