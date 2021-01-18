# If resid is a matrix of B residuals, outputs a B vector
# Expects x_alt to be scaled and a vector, and resid to be scaled and centred
# gives the correlation between resid and x_alt
# Used in pvals_all with the normal option
cov_test_twoside <- function(x_alt, resid, mc.cores=NULL, x_null=NULL) {
  if (is.matrix(resid) && (ncol(resid) > 1)) {
    return(colMeans(x_alt * resid))
  } else {
    resid <- as.numeric(resid)
    return(mean(x_alt*resid))
  }
}

rf_test <- function(x_alt, resid, mc.cores=NULL, x_null=NULL) {
  mean((resid - predict(randomForest::randomForest(x_alt, resid), x_alt))^2)
  #mean((resid-randomForest(x_alt, resid)$predicted)^2)
}

rf_test_S <- function(x_alt, resid, mc.cores=NULL, x_null=NULL) {
  tol <- 1e-3
  n <- nrow(x_alt)
  # Find S_hat
  abs_cor <- abs(colMeans(resid*x_alt))
  S_hat <- which(abs_cor > max(abs_cor) - tol)
  x_alt <- cbind(1, x_alt[, S_hat])
  # out of bag error
  mean((resid - randomForest::randomForest(x_alt, resid)$predicted)^2)
  #mean((resid - predict(randomForest(x_alt, resid), x_alt))^2)
}

# test_aux will be a sequence of lambda values
lasso_test <- function(x_alt, resid, test_aux, mc.cores=1L, x_null=NULL) {
  if (missing(test_aux)) {
    glmnet_out <- glmnet(x_alt, resid, standardize = FALSE)
    return(list("test" = colMeans((resid - glmnet::predict.glmnet(glmnet_out, newx=x_alt))^2),
                "test_aux" = glmnet_out$lambda))
  } else {
    # resid will now be resid_sim, an n by B matrix of residuals
    B <- ncol(resid)
    rep_fun <- function(b) {
      glmnet_out <- glmnet(x_alt, resid[, b], lambda=test_aux, standardize = FALSE)
      return(colMeans((resid[, b] - predict.glmnet(glmnet_out, newx=x_alt))^2))
    }
    test_sim <- parallel::mclapply(1:B, rep_fun, mc.cores=mc.cores)
    return(test_sim)
  }
}

# Assumes resid is from a fit to a design x with scaled cols
# x_null contains the original x (x_null)
lasso_test_S_inp <- function(x_alt, resid, test_aux, mc.cores=1L, x_null=NULL, S_hat_inp) {
  tol <- 1e-3
  if (missing(test_aux)) {
    S_hat <- S_hat_inp
    x_alt <- cbind(rep(1, nrow(x_alt)), x_alt[, S_hat, drop=FALSE])
    glmnet_out <- glmnet(x_alt, resid, standardize = FALSE)
    return(list("test" = colMeans((resid - predict.glmnet(glmnet_out, newx=x_alt))^2),
                "test_aux" = glmnet_out$lambda))
  } else {
    # resid will now be resid_sim, an n by B matrix of residuals
    B <- ncol(resid)
    rep_fun <- function(b) {
      S_hat <- S_hat_inp[[b]]
      x_alt <- cbind(rep(1, nrow(x_alt)), x_alt[, S_hat, drop=FALSE])
      glmnet_out <- glmnet(x_alt, resid[, b], lambda=test_aux, standardize=FALSE)
      return(colMeans((resid[, b] - predict.glmnet(glmnet_out, newx=x_alt))^2))
    }
    test_sim <- parallel::mclapply(1:B, rep_fun, mc.cores=mc.cores)
    return(test_sim)
  }
}

# To test for heteroscedasticity
# Similar to Breusch-Pagan test
sq_lasso_test <- function(x_alt, resid, test_aux, mc.cores=1L, x_null=NULL) {
  if (missing(test_aux)) {
    resid <- resid^2
    resid <- Scale(resid)
    return(lasso_test(x_alt, resid))
  } else {
    n <- nrow(resid)
    resid <- resid^2
    resid <- scale(resid) * sqrt(n/(n-1))
    return(lasso_test(x_alt, resid, test_aux))
  }
}

sq_lasso_test_S <- function(x_alt, resid, test_aux, mc.cores=1L, x_null=NULL) {
  tol <- 2e-3
  if (missing(test_aux)) {
    abs_cor <- abs(colMeans(resid*x_alt))
    S_hat <- which(abs_cor > max(abs_cor) - tol)

    resid <- resid^2
    resid <- Scale(resid)
    return(lasso_test_S_inp(x_alt, resid, S_hat_inp=S_hat))
  } else {
    S_hat_inp <- lapply(1:ncol(resid), function(b) {
      abs_cor <- abs(colMeans(resid[, b]*x_alt))
      return(which(abs_cor > max(abs_cor) - tol))
    })
    n <- nrow(resid)
    resid <- resid^2
    resid <- scale(resid) * sqrt(n/(n-1))
    return(lasso_test_S_inp(x_alt, resid, test_aux, mc.cores=mc.cores, S_hat_inp=S_hat_inp))
  }
}
