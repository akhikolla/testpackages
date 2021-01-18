## ----setup, include = FALSE----------------------------------------------
library(gasper)
library(rwavelet)
if (as.numeric(R.version$minor)>6){
  RNGkind(sample.kind = "Rounding")
}
set.seed(0)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=3,
  fig.height=3, 
  fig.align="center"
)

## ------------------------------------------------------------------------
graphname <- "netz4504"
groupname <- "AG-Monien"
download_graph(graphname,groupname)
attributes(netz4504)

## ----fig.show='hold'-----------------------------------------------------
f <- rnorm(nrow(netz4504$xy))
plot_graph(netz4504)
plot_signal(netz4504, f,size = f)

## ------------------------------------------------------------------------
A <- full(rlogo$sA)
L <- laplacian_mat(A)
val1 <- eigendec(L)
evalues <- val1$evalues
evectors <- val1$evectors
lmax <- max(evalues)
b <- 2
kmax <- floor(log(lmax)/log(b)) + 2
tf <- tight_frame(evalues, evectors, b=b)

## ------------------------------------------------------------------------
x1 <- rlogo$xy[,1]
x2 <- rlogo$xy[,2]
n <- length(x1)
f <- randsignal(0.01, 3, A)
sigma <- 0.01
noise <- rnorm(n, sd = sigma)
tilde_f <- f + noise

## ----fig.show='hold'-----------------------------------------------------
plot_signal(rlogo, f, size = 3)
plot_signal(rlogo, tilde_f, size = 3)

## ------------------------------------------------------------------------
wcn <- analysis(tilde_f,tf)
wcf <- analysis(f,tf)

## ------------------------------------------------------------------------
diagWWt <- colSums(t(tf)^2)
thresh <- sort(abs(wcn))
opt_thresh <- SURE_MSEthresh(wcn, 
                           wcf, 
                           thresh, 
                           diagWWt, 
                           b=2, 
                           sigma, 
                           NA,
                           policy = "dependent")


## ------------------------------------------------------------------------
plot(thresh,opt_thresh$res$MSE,type="l",xlab = "t",ylab = "risk")
lines(thresh,opt_thresh$res$SURE-n*sigma^2,col="red")
legend("topleft", legend=c("MSE", "SURE"),
       col=c("black", "red"),lty=1)

## ------------------------------------------------------------------------
hatf_oracle <- synthesis(opt_thresh$wc[,opt_thresh$min[1]], tf)
hatf_SURE  <- synthesis(opt_thresh$wc[,opt_thresh$min[2]], tf)

SNR(f,tilde_f)
SNR(f,hatf_oracle)
SNR(f,hatf_SURE)

