NBCEE <- function (X, Y, U, omega, niter = 5000, nburn = 500, nthin = 10, 
    maxmodelX = NA, maxmodelY = NA, family.X = "gaussian") 
{
    n = length(Y)
    n_cov = ncol(as.matrix(U))
    if (is.na(maxmodelX)) 
        maxmodelX = min(niter + nburn, 2^n_cov)
    if (is.na(maxmodelY)) 
        maxmodelY = min(niter + nburn, 2^n_cov)
    alpha_X = as.numeric(bic.glm(y = X, x = U, glm.family = family.X, 
        OR = 1.0000001)$which[1, ])
    alpha_Y = as.numeric(bic.glm(y = Y, x = cbind(X, U), glm.family = "gaussian", 
        OR = 1.0000001, prior.param = c(1, rep(0.5, n_cov)))$which[1, 
        2:(n_cov + 1)])
    if (sum(alpha_X) == 0) {
        model_X = glm(X ~ 1, family = family.X)
    } else {
        model_X = glm(X ~ U[, alpha_X == 1], family = family.X)
    }
    models_X = matrix(0, nrow = 1 + floor((niter-1)/nthin), ncol = n_cov)
    tested_models_X = matrix(-1, nrow = maxmodelX, ncol = 1)
    pX_tested = matrix(NA, nrow = maxmodelX, ncol = 1)
    bic_x = BIC(model_X)
    if (sum(alpha_Y) == 0) {
        model_Y0 = lm(Y ~ X)
    } else {
        model_Y0 = lm(Y ~ X + U[, alpha_Y == 1])
    }
    bic_y = BIC(model_Y0)
    betas = numeric(1 + floor((niter-1)/nthin))
    models_Y = matrix(0, nrow = 1 + floor((niter-1)/nthin), ncol = n_cov)
    tested_models_Y = matrix(-1, nrow = maxmodelY, ncol = 1)
    pY_tested = matrix(NA, nrow = maxmodelY, ncol = 1)
    pY_tested1 = matrix(NA, nrow = maxmodelY, ncol = n_cov)
    j = 1
    k_X = 1
    k_Y = 1
    sy = sd(Y)
    su = apply(U, 2, sd)
    tested_models_Y[k_Y] = paste0(alpha_Y, collapse = "")
    a = rep(NA, n_cov)
    if (sum(alpha_X != 0)) {
        missing_causes = which(alpha_X == 1)
        for (d in missing_causes) {
            if (alpha_Y[d] == 1) {
                indice = sum(alpha_Y[(1:n_cov) <= d]) + 2
                a[d] = coef(model_Y0)[indice]
            }
            else {
                machin = numeric(n_cov)
                machin[d] = 1
                alpha_Ym = alpha_Y + machin
                model_Ym = lm(Y ~ X + U[, alpha_Ym == 1])
                indice = sum(alpha_Ym[(1:n_cov) <= d]) + 2
                a[d] = coef(model_Ym)[indice]
            }
        }
        pY_tested[k_Y] = bic_y
        pY_tested1[k_Y, ] = a
        k_Y = k_Y + 1
    }
    for (i in 1:(nburn + niter)) {
        alpha_X_0 = alpha_X
        alpha_X_1 = alpha_X
        temp = sample(n_cov, 1)
        alpha_X_1[temp] = (alpha_X_1[temp] + 1)%%2
        bic_x_0 = bic_x
        char_alpha_X_1 = paste0(alpha_X_1, collapse = "")
        tested = which(tested_models_X == char_alpha_X_1)
        if (length(tested)) {
            bic_x_1 = pX_tested[tested]
        }
        else {
            if (sum(alpha_X_1) == 0) {
                bic_x_1 = BIC(glm(X ~ 1, family = family.X))
            }
            else {
                model_X = glm(X ~ U[, alpha_X_1 == 1], family = family.X)
                bic_x_1 = BIC(model_X)
            }
            tested_models_X[k_X] = char_alpha_X_1
            pX_tested[k_X] = bic_x_1
            k_X = k_X + 1
        }
        ratio = exp(-bic_x_1/2 + bic_x_0/2)
        if (sample(c(1, 0), size = 1, prob = c(min(1, ratio), 
            1 - min(1, ratio)))) {
            alpha_X = alpha_X_1
            bic_x = bic_x_1
        }
        alpha_Y_0 = alpha_Y
        alpha_Y_1 = alpha_Y
        temp = sample(n_cov, 1)
        alpha_Y_1[temp] = (alpha_Y_1[temp] + 1)%%2
        bic_y_0 = bic_y
        tested = which(tested_models_Y == paste0(alpha_Y_0, collapse = ""))
        prior0 = rep(1, n_cov)
        if (sum(alpha_X != 0)) {
            a_0 = pY_tested1[tested, ]
            missing_causes = which(alpha_X == 1 & is.na(a_0))
            if (length(missing_causes != 0)) {
                if (sum(alpha_Y_0) == 0) {
                  model_Y0 = lm(Y ~ X)
                }
                else {
                  model_Y0 = lm(Y ~ X + U[, alpha_Y_0 == 1])
                }
                dunno = missing_causes[which(alpha_Y_0[missing_causes] == 
                  1)]
                if (length(dunno)) {
                  presents = sapply(dunno, function(x) {
                    sum(alpha_Y_0[(1:n_cov) <= x]) + 2
                  })
                  a_0[dunno] = coef(model_Y0)[presents]
                }
                loop = missing_causes[!(missing_causes %in% dunno)]
                for (d in loop) {
                  model_Ym = update(model_Y0, . ~ . + U[, d])
                  a_0[d] = coef(model_Ym)[sum(alpha_Y_0) + 3]
                }
                pY_tested1[tested, ] = a_0
            }
            a_0[!is.na(a_0)] = omega * (a_0[!is.na(a_0)] * su[!is.na(a_0)]/sy)^2
            prior0[alpha_X == 1 & alpha_Y_0 == 1] = a_0[alpha_X == 
                1 & alpha_Y_0 == 1]/(a_0[alpha_X == 1 & alpha_Y_0 == 
                1] + 1)
            prior0[alpha_X == 1 & alpha_Y_0 == 0] = 1/(a_0[alpha_X == 
                1 & alpha_Y_0 == 0] + 1)
        }
        char_alpha_Y_1 = paste0(alpha_Y_1, collapse = "")
        tested = which(tested_models_Y == char_alpha_Y_1)
        if (length(tested)) {
            prior1 = rep(1, n_cov)
            if (sum(alpha_X != 0)) {
                a_1 = pY_tested1[tested, ]
                missing_causes = which(alpha_X == 1 & is.na(a_1))
                if (length(missing_causes) != 0) {
                  if (sum(alpha_Y_1) == 0) {
                    model_Y1 = lm(Y ~ X)
                  }
                  else {
                    model_Y1 = lm(Y ~ X + U[, alpha_Y_1 == 1])
                  }
                  dunno = missing_causes[which(alpha_Y_1[missing_causes] == 
                    1)]
                  if (length(dunno)) {
                    presents = sapply(dunno, function(x) {
                      sum(alpha_Y_1[(1:n_cov) <= x]) + 2
                    })
                    a_1[dunno] = coef(model_Y1)[presents]
                  }
                  loop = missing_causes[!(missing_causes %in% 
                    dunno)]
                  for (d in loop) {
                    model_Ym = update(model_Y1, . ~ . + U[, d])
                    a_1[d] = coef(model_Ym)[sum(alpha_Y_1) + 
                      3]
                  }
                  pY_tested1[tested, ] = a_1
                }
                a_1[!is.na(a_1)] = omega * (a_1[!is.na(a_1)] * 
                  su[!is.na(a_1)]/sy)^2
                prior1[alpha_X == 1 & alpha_Y_1 == 1] = a_1[alpha_X == 
                  1 & alpha_Y_1 == 1]/(a_1[alpha_X == 1 & alpha_Y_1 == 
                  1] + 1)
                prior1[alpha_X == 1 & alpha_Y_1 == 0] = 1/(a_1[alpha_X == 
                  1 & alpha_Y_1 == 0] + 1)
            }
            bic_y_1 = pY_tested[tested]
        }
        else {
            if (sum(alpha_Y_1) == 0) {
                model_Y1 = lm(Y ~ X)
            }
            else {
                model_Y1 = lm(Y ~ X + U[, alpha_Y_1 == 1])
            }
            bic_y_1 = BIC(model_Y1)
            prior1 = rep(1, n_cov)
            if (sum(alpha_X != 0)) {
                a_1 = rep(NA, n_cov)
                missing_causes = which(alpha_X == 1)
                dunno = missing_causes[which(alpha_Y_1[missing_causes] == 
                  1)]
                if (length(dunno)) {
                  presents = sapply(dunno, function(x) {
                    sum(alpha_Y_1[(1:n_cov) <= x]) + 2
                  })
                  a_1[dunno] = coef(model_Y1)[presents]
                }
                loop = missing_causes[!(missing_causes %in% dunno)]
                for (d in loop) {
                  model_Ym = update(model_Y1, . ~ . + U[, d])
                  a_1[d] = coef(model_Ym)[sum(alpha_Y_1) + 3]
                }
                pY_tested1[k_Y, ] = a_1
            }
            a_1[!is.na(a_1)] = omega * (a_1[!is.na(a_1)] * su[!is.na(a_1)]/sy)^2
            prior1[alpha_X == 1 & alpha_Y_1 == 1] = a_1[alpha_X == 
                1 & alpha_Y_1 == 1]/(a_1[alpha_X == 1 & alpha_Y_1 == 
                1] + 1)
            prior1[alpha_X == 1 & alpha_Y_1 == 0] = 1/(a_1[alpha_X == 
                1 & alpha_Y_1 == 0] + 1)
            tested_models_Y[k_Y] = char_alpha_Y_1
            pY_tested[k_Y] = bic_y_1
            k_Y = k_Y + 1
        }
        ratio_P_alpha = prod(prior1/prior0)
        if (ratio_P_alpha == 0) {
            ratio = 0
        }
        else if (ratio_P_alpha == Inf) {
            ratio = 1
        }
        else {
            ratio = exp(-bic_y_1/2 + bic_y_0/2) * ratio_P_alpha
        }
        if (sample(c(1, 0), size = 1, prob = c(min(1, ratio), 
            1 - min(1, ratio)))) {
            alpha_Y = alpha_Y_1
            bic_y = bic_y_1
        }
        if (i > nburn & (i - nburn)%%nthin == 1) {
            if (sum(alpha_Y) != 0) {
                model_b = lm(Y ~ X + U[, alpha_Y == 1])
            }
            else {
                model_b = lm(Y ~ X)
            }
            beta_alphaY = rnorm(1, mean = model_b$coefficients[2], 
                sd = sqrt(vcov(model_b)[2, 2]))
            betas[j] = beta_alphaY
            models_X[j, ] = alpha_X
            models_Y[j, ] = alpha_Y
            j = j + 1
        }
    }
    return(list(betas = betas, models.X = models_X, models.Y = models_Y))
}