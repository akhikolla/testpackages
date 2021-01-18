## ----setting, echo = FALSE----------------------------------------------------
knitr::opts_chunk$set(fig.width = 6, fig.height = 4) 

## ----install, eval = FALSE----------------------------------------------------
#  ## on CRAN
#  ## install.packages("Compack")
#  library(Compack)

## ----Fcomp_Model--------------------------------------------------------------
library(Compack)
df_beta = 5
p = 30
beta_C_true = matrix(0, nrow = p, ncol = df_beta)
beta_C_true[1, ] <- c(-0.5, -0.5, -0.5 , -1, -1)
beta_C_true[2, ] <- c(0.8, 0.8,  0.7,  0.6,  0.6)
beta_C_true[3, ] <- c(-0.8, -0.8 , 0.4 , 1 , 1)
beta_C_true[4, ] <- c(0.5, 0.5, -0.6  ,-0.6, -0.6)

n_train = 50
n_test = 30
Data <- Fcomp_Model(n = n_train, p = p, m = 0, intercept = TRUE,
                    SNR = 4, sigma = 3, rho_X = 0.6, rho_T = 0.2,
                    df_beta = df_beta, n_T = 20, obs_spar = 1, theta.add = FALSE,
                    beta_C = as.vector(t(beta_C_true)))
arg_list <- as.list(Data$call)[-1]
arg_list$n <- n_test
Test <- do.call(Fcomp_Model, arg_list)

## ----FuncompCL----------------------------------------------------------------
m1 <- FuncompCGL(y = Data$data$y, X = Data$data$Comp, 
                 Zc = Data$data$Zc, intercept = Data$data$intercept, 
                 k = df_beta)
plot(m1, xlab = "-log")
betamat <- coef(m1)
predmat <- predict(m1, Data$data$Comp, Data$data$Zc)

## ----cv_cgl-------------------------------------------------------------------
k_list = c(4,5,6)
nfolds = 10
foldid <- sample(rep(seq(nfolds), length = n_train))
## cv_cgl: Constrained group lasso
cv_cgl <-  cv.FuncompCGL(y = Data$data$y, X = Data$data$Comp,
                         Zc = Data$data$Zc, intercept = Data$data$intercept,
                         k = k_list, foldid = foldid,
                         keep = TRUE)
plot(cv_cgl,k = k_list)
cv_cgl$Ftrim[c("lam.min", "lam.1se")]
beta <-  coef(cv_cgl, trim = FALSE, s = "lam.min")
k_opt <- cv_cgl$Ftrim$lam.min['df']
plot(cv_cgl$FuncompCGL.fit[[as.character(k_opt)]])
m1 <- ifelse(is.null(ncol(Data$data$Zc)), 0, ncol(Data$data$Zc))
m1 <- m1 + Data$data$intercept
if(k_opt == df_beta) {
  plot(Data$beta, col = "red", pch = 19,
       ylim = range(c(range(Data$beta), range(beta))))
  abline(v= seq(from = 0, to = (p*df_beta), by = df_beta ))
  abline(h = 0)
  points(beta)
  if(m1 > 0) points(p*df_beta + 1:m1, tail(Data$beta, m1),
                    col = "blue", pch = 19)
} else {
  plot(beta, ylim = range(c(range(Data$beta), range(beta))) )
  abline(v= seq(from = 0, to = (p*k_opt), by = k_opt ))
  abline(h = 0, col = "red")
  if(m1 > 0) points(p*k_opt + 1:m1, tail(Data$beta, m1),
                    col = "blue", pch = 19)
}
beta_C <- matrix(beta[1:(p*k_opt)], byrow = TRUE, nrow = p)
## satisfies zero-sum constraints
cat("colSums:", colSums(beta_C))
Nonzero <- (1:p)[apply(beta_C, 1, function(x) max(abs(x)) >0)]
cat("selected groups:", Nonzero)

sseq <- Data$basis.info[, 1]
beta_curve_true <- Data$basis.info[, -1] %*%  t(beta_C_true)
Nonzero_true <- (1:p)[apply(beta_C_true, 1, function(x) max(abs(x)) >0)]
matplot(sseq, beta_curve_true, type = "l", ylim = range(beta_curve_true),
        ylab = "True coeffcients curves", xlab = "TIME")
abline(a = 0, b = 0, col = "grey", lwd = 2)
text(0, beta_curve_true[1, Nonzero_true], labels = Nonzero_true)

beta_curve <- splines::bs(sseq, df = k_opt, intercept = TRUE) %*% t(beta_C)
matplot(sseq, beta_curve, type = "l", ylim = range(beta_curve_true),
        ylab = "Estimated coefficient curves", xlab = "TIME")
abline(a = 0, b = 0, col = "grey", lwd = 2)
text(0, beta_curve[1, Nonzero], labels = Nonzero)

## set a thresholding for variable selection via cross-validation model
## example: cut by average L2-norm for estiamted coefficient curves
Curve_L2 <- colSums(beta_curve^2)
Curve_L2 <- Curve_L2 - colSums(beta_curve[c(1, nrow(beta_curve)), ]^2) / 2
Curve_L2 <- Curve_L2 * (Data$basis.info[2,1] - Data$basis.info[1,1])
Curve_L2 <- sqrt(Curve_L2)
plot(Curve_L2, xlab = "Component index", ylab = "L2-norm for coefficient curves")
cutoff <- sum(Curve_L2) / p
Nonzero_cut <- (1:p)[which(Curve_L2 >= cutoff)]
cat("selected groups after thresholding cut-off:", Nonzero_cut)
y_hat <- predict(cv_cgl, Data$data$Comp, Data$data$Zc, s = "lam.min")
MSE <- sum((drop(Data$data$y) - y_hat)^2) / n_train
y_hat <- predict(cv_cgl, Test$data$Comp, Test$data$Zc, s = "lam.min")
PRE <- sum((drop(Test$data$y) - y_hat)^2) / n_test
cgl_result <- list(cv.result = cv_cgl, beta = beta,
                   Nonzero = c("Original" = Nonzero, "Cut" = Nonzero_cut),
                   MSE = MSE, PRE = PRE)

## ----cv_naive-----------------------------------------------------------------
## set mu_raio = 0 to identifying without linear constraints,
## no outer_loop for Lagrange augmented multiplier
cv_naive <-  cv.FuncompCGL(y = Data$data$y, X = Data$data$Comp,
                           Zc = Data$data$Zc, intercept = Data$data$intercept,
                           k = k_list, foldid = foldid, keep = TRUE,
                           mu_ratio = 0)
plot(cv_naive, k = k_list)
beta <-  coef(cv_naive, trim = FALSE, s = "lam.min")
k_opt <- cv_naive$Ftrim$lam.min['df']
beta_C <- matrix(beta[1:(p*k_opt)], byrow = TRUE, nrow = p)
## does NOT satisfy zero-sum constraints
cat("colSums:", colSums(beta_C))
Nonzero <- (1:p)[apply(beta_C, 1, function(x) max(abs(x)) >0)]

beta_curve <- splines::bs(sseq, df = k_opt, intercept = TRUE) %*% t(beta_C)
Curve_L2 <- colSums(beta_curve^2) - colSums(beta_curve[c(1, nrow(beta_curve)), ]^2) / 2
Curve_L2 <- sqrt(Curve_L2 * (Data$basis.info[2,1] - Data$basis.info[1,1]))
cutoff <- sum(Curve_L2) / p
Nonzero_cut <- (1:p)[which(Curve_L2 >= cutoff)]
MSE <- sum((drop(Data$data$y) - predict(cv_naive, Data$data$Comp, Data$data$Zc, s = "lam.min"))^2) / n_train
PRE <- sum((drop(Test$data$y) - predict(cv_naive, Test$data$Comp, Test$data$Zc, s = "lam.min"))^2) / n_test
naive_result <- list(cv.result = cv_naive, beta = beta,
                     Nonzero = c("Original" = Nonzero, "Cut" = Nonzero_cut),
                     MSE = MSE, PRE = PRE)


## ----cv_base------------------------------------------------------------------
## mu_ratio is set to 0 automatically once ref is set to a integer
ref = sample(1:p, 1)
cv_base <- cv.FuncompCGL(y = Data$data$y, X = Data$data$Comp,
                         Zc = Data$data$Zc, intercept = Data$data$intercept,
                         k = k_list, foldid = foldid, keep = TRUE,
                         ref = ref)
plot(cv_base, k = k_list)
beta <-  coef(cv_base, trim = FALSE, s = "lam.min")
k_opt <- cv_base$Ftrim$lam.min['df']
beta_C <- matrix(beta[1:(p*k_opt)], byrow = TRUE, nrow = p)
## satisfies zero-sum constraints
cat("colSums:", colSums(beta_C))
Nonzero <- (1:p)[apply(beta_C, 1, function(x) max(abs(x)) >0)]
beta_curve <- splines::bs(sseq, df = k_opt, intercept = TRUE) %*% t(beta_C)
Curve_L2 <- colSums(beta_curve^2) - colSums(beta_curve[c(1, nrow(beta_curve)), ]^2) / 2
Curve_L2 <- sqrt(Curve_L2 * (Data$basis.info[2,1] - Data$basis.info[1,1]))
cutoff <- sum(Curve_L2) / p
Nonzero_cut <- (1:p)[which(Curve_L2 >= cutoff)]
MSE <- sum((drop(Data$data$y) - predict(cv_base, Data$data$Comp, Data$data$Zc, s = "lam.min"))^2) / n_train
PRE <- sum((drop(Test$data$y) - predict(cv_base, Test$data$Comp, Test$data$Zc, s = "lam.min"))^2) / n_test
base_result <- list(cv.result = cv_base, beta = beta,
                    Nonzero = c("Original" = Nonzero, "Cut" = Nonzero_cut),
                    MSE = MSE, PRE = PRE)

## ----GIC_cgl------------------------------------------------------------------
GIC_cgl <- GIC.FuncompCGL(y = Data$data$y, X = Data$data$Comp,
                          Zc = Data$data$Zc, intercept = Data$data$intercept,
                          k = k_list)
beta <- coef(GIC_cgl)
plot(GIC_cgl)
y_hat <- predict(GIC_cgl, Znew = Test$data$Comp, Zcnew = Test$data$Zc)
plot(Test$data$y, y_hat, xlab = "Observed response", ylab = "Predicted response")

## ----GIC_naive----------------------------------------------------------------
GIC_naive <- GIC.FuncompCGL(y = Data$data$y, X = Data$data$Comp,
                            Zc = Data$data$Zc, intercept = Data$data$intercept,
                            k = k_list, mu_ratio = 0)
beta <- coef(GIC_naive)
plot(GIC_naive)
y_hat <- predict(GIC_naive, Znew = Test$data$Comp, Zcnew = Test$data$Zc)
plot(Test$data$y, y_hat, xlab = "Observed response", ylab = "Predicted response")

## ----GIC_base-----------------------------------------------------------------
GIC_base <- GIC.FuncompCGL(y = Data$data$y, X = Data$data$Comp,
                            Zc = Data$data$Zc, intercept = Data$data$intercept,
                            k = k_list, ref = ref)
beta <- coef(GIC_base)
plot(GIC_base)
y_hat <- predict(GIC_base, Znew = Test$data$Comp, Zcnew = Test$data$Zc)
plot(Test$data$y, y_hat, xlab = "Observed response", ylab = "Predicted response")

## ----comp_Model---------------------------------------------------------------
library(Compack)
p = 30
n = 50
beta = c(1, -0.8, 0.6, 0, 0, -1.5, -0.5, 1.2)
beta = c(beta, rep(0, times = p - length(beta)))
Comp_data = comp_Model(n = n, p = p, beta = beta, intercept = FALSE)
Comp_data2 = comp_Model(n = n, p = p, beta = Comp_data$beta, intercept = FALSE)

## ----compCL-------------------------------------------------------------------
m1 <- compCL(y = Comp_data$y, Z = Comp_data$X.comp,
             Zc = Comp_data$Zc, intercept = Comp_data$intercept)
plot(m1, label = TRUE)
coef(m1)[1:10, 90:100]

## ----cv.compCL----------------------------------------------------------------
cvm1 <- cv.compCL(y = Comp_data$y, Z = Comp_data$X.comp,
                  Zc = Comp_data$Zc, intercept = Comp_data$intercept)
plot(cvm1, xlab = "-log")
beta_est <- coef(cvm1, s = "lam.min")
head(beta_est, 10)
sum(beta_est[1:p]) # satisfies zero-sum constraint
y_hat <- predict(cvm1, Comp_data2$X.comp, Comp_data2$Zc, s = "lam.min")

## ----GIC.compCL---------------------------------------------------------------
GICm1 <- GIC.compCL(y = Comp_data$y, Z = Comp_data$X.comp,
                    Zc = Comp_data$Zc, intercept = Comp_data$intercept)
plot(GICm1, xlab = "-log")
beta_est <- coef(GICm1, s = "lam.min")
head(beta_est, 10)
sum(beta_est[1:p]) # satisfies zero-sum constraint
y_hat <- predict(GICm1, Comp_data2$X.comp, Comp_data2$Zc, s = "lam.min")

