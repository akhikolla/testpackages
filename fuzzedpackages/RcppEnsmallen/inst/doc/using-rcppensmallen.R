## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----data-generation----------------------------------------------------------
n <- 1e6
beta <- c(-2, 1.5, 3, 8.2, 6.6)
p <- length(beta)
X <- cbind(1, matrix(rnorm(n), ncol = p - 1))
y <- X %*% beta + rnorm(n / (p - 1))

## ----data-estimation, eval = FALSE--------------------------------------------
#  coefs_lbfgs <- lin_reg_lbfgs(X, y)
#  coefs_lm <- lm.fit(X, y)$coefficients

## ----echo = FALSE-------------------------------------------------------------
coefs_lbfgs <- RcppEnsmallen::lin_reg_lbfgs(X, y)
coefs_lm <- lm.fit(X, y)$coefficients
compare_coefs = cbind(coefs_lbfgs, coefs_lm)
colnames(compare_coefs) = c("LBFGS", "LM")
rownames(compare_coefs) = paste0("Beta", seq_len(nrow(compare_coefs)))
knitr::kable(compare_coefs, longtable = FALSE, caption = "Comparison of Estimated Coefficients")

