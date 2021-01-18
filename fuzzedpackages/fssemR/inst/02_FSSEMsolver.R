##' FSSEMsolver
seed = as.numeric(Sys.time())
N  = 100                                                              # sample size. 500 sample better than 200 sample, very very
Ng = 100                                                              # gene number
Nk = 100 * 3                                                          # eQTL number
Ns = 100 / Ng                                                          # sparse ratio
sigma2 = 0.01                                                         # sigma2
set.seed(seed)
library(fssemR)

data = randomFSSEMdata(n = N, p = Ng, k = Nk, sparse = Ns, df = 0.3, sigma2 = sigma2, u = 5, type = "DG", nhub = 1, dag = T)
## data = randomFSSEMdata(n = N, p = Ng, k = Nk, sparse = Ns, df = 0.3, sigma2 = sigma2, u = 5, type = "ER", nhub = 1)
gamma = cv.multiRegression(data$Data$X, data$Data$Y, data$Data$Sk, ngamma = 50, nfold = 5, N, Ng, Nk)
fit   = multiRegression(data$Data$X, data$Data$Y, data$Data$Sk, gamma, N, Ng, Nk, trans = FALSE)
Xs    = data$Data$X
Ys    = data$Data$Y
Sk    = data$Data$Sk

## cross-validation
cvfitc <- cv.multiFSSEMiPALM(Xs = Xs, Ys = Ys, Bs = fit$Bs, Fs = fit$Fs, Sk = Sk,
                             sigma2 = fit$sigma2, nlambda = 10, nrho = 10,
                             nfold = 5, p = Ng, q = Nk, wt = T)

system.time(fitc <<- multiFSSEMiPALM(Xs = Xs, Ys = Ys, Bs = fit$Bs, Fs = fit$Fs, Sk = Sk,
                            sigma2 = fit$sigma2, lambda = cvfitc$lambda, rho = cvfitc$rho,
                            Wl = inverseB(fit$Bs), Wf = flinvB(fit$Bs),
                            p = Ng, maxit = 100, threshold = 1e-5, sparse = T, verbose = T, trans = T, strict = F))


(TPR(fitc$Bs[[1]], data$Vars$B[[1]]) + TPR(fitc$Bs[[2]], data$Vars$B[[2]])) / 2
(FDR(fitc$Bs[[1]], data$Vars$B[[1]]) + FDR(fitc$Bs[[2]], data$Vars$B[[2]])) / 2
TPR(fitc$Bs[[1]] - fitc$Bs[[2]], data$Vars$B[[1]] - data$Vars$B[[2]])
FDR(fitc$Bs[[1]] - fitc$Bs[[2]], data$Vars$B[[1]] - data$Vars$B[[2]])


fitm <- opt.multiFSSEMiPALM(Xs = Xs, Ys = Ys, Bs = fit$Bs, Fs = fit$Fs, Sk = Sk,
                            sigma2 = fit$sigma2, nlambda = 10, nrho = 10,
                            p = Ng, q = Nk, wt = T)

fitc0 <- fitm$fit

(TPR(fitc0$Bs[[1]], data$Vars$B[[1]]) + TPR(fitc0$Bs[[2]], data$Vars$B[[2]])) / 2
(FDR(fitc0$Bs[[1]], data$Vars$B[[1]]) + FDR(fitc0$Bs[[2]], data$Vars$B[[2]])) / 2
TPR(fitc0$Bs[[1]] - fitc0$Bs[[2]], data$Vars$B[[1]] - data$Vars$B[[2]])
FDR(fitc0$Bs[[1]] - fitc0$Bs[[2]], data$Vars$B[[1]] - data$Vars$B[[2]])




