

# this function calls on other functions in order to return the frequentist estimates
# and parametric bootstrapped confidence intervals, sampling from a multivariate normal distribution

freqFun_para <- function(data, n.boot, estimates, interval, omega.freq.method,
                         item.dropped, alpha.int.analytic, omega.int.analytic, pairwise, callback = function(){},
                         parametric = T){
  p <- ncol(data)
  n <- nrow(data)
  if (pairwise) {
    cc <- cov(data, use = "pairwise.complete.obs")
  } else{
    cc <- cov(data)
  }
  res <- list()
  res$covsamp <- NULL
  boot_cov <- NULL
  if (("alpha" %in% estimates & !alpha.int.analytic) |
      "lambda2" %in% estimates | "lambda4" %in% estimates | "lambda6" %in% estimates |
      "glb" %in% estimates | ("omega" %in% estimates & omega.freq.method == "pfa")){

    boot_cov <- array(0, c(n.boot, p, p))
    for (i in 1:n.boot){
      boot_data <- MASS::mvrnorm(n, colMeans(data, na.rm = T), cc)
      boot_cov[i, , ] <- cov(boot_data)
      callback()
    }

    res$covsamp <- boot_cov
  }
  if (item.dropped){
    Ctmp <- array(0, c(p, p - 1, p - 1))
    for (i in 1:p){
      Ctmp[i, , ] <- cc[-i, -i]
    }
  }

  if ("alpha" %in% estimates){
    res$est$freq_alpha <- applyalpha(cc)
    if (alpha.int.analytic){
      int <- ciAlpha(1 - interval, n, cc)
      res$conf$low$freq_alpha <- int[1]
      res$conf$up$freq_alpha <- int[2]
    } else{
      alpha_obj <- apply(boot_cov, 1, applyalpha, callback)
      if (length(unique(round(alpha_obj, 4))) == 1){
        res$conf$low$freq_alpha <- 1
        res$conf$up$freq_alpha <- 1
      } else{
        res$conf$low$freq_alpha <- quantile(alpha_obj, probs = (1 - interval)/2, na.rm = T)
        res$conf$up$freq_alpha <- quantile(alpha_obj, probs = interval + (1 - interval)/2, na.rm = T)
      }
      res$boot$alpha <- alpha_obj
    }
    if (item.dropped){
      res$ifitem$alpha <- apply(Ctmp, 1, applyalpha)
    }
  }

  if ("lambda2" %in% estimates){
    res$est$freq_lambda2 <- applylambda2(cc)
    lambda2_obj <- apply(boot_cov, 1, applylambda2, callback)
    if (length(unique(round(lambda2_obj, 4))) == 1){
      res$conf$low$freq_lambda2 <- NA
      res$conf$up$freq_lambda2 <- NA
    } else{
      res$conf$low$freq_lambda2 <- quantile(lambda2_obj, probs = (1 - interval)/2, na.rm = T)
      res$conf$up$freq_lambda2 <- quantile(lambda2_obj, probs = interval + (1 - interval)/2, na.rm = T)
    }
    res$boot$lambda2 <- lambda2_obj
    if (item.dropped){
      res$ifitem$lambda2 <- apply(Ctmp, 1, applylambda2)
    }
  }

  if ("lambda4" %in% estimates){
    res$est$freq_lambda4 <- applylambda4_nocpp(cc)
    lambda4_obj <- apply(boot_cov, 1, applylambda4_nocpp, callback)
    if (length(unique(round(lambda4_obj, 4))) == 1){
      res$conf$low$freq_lambda4 <- NA
      res$conf$up$freq_lambda4 <- NA
    } else{
      res$conf$low$freq_lambda4 <- quantile(lambda4_obj, probs = (1 - interval)/2, na.rm = T)
      res$conf$up$freq_lambda4 <- quantile(lambda4_obj, probs = interval + (1 - interval)/2, na.rm = T)
    }
    res$boot$lambda4 <- lambda4_obj
    if (item.dropped){
      res$ifitem$lambda4 <- apply(Ctmp, 1, applylambda4_nocpp)
    }
  }

  if ("lambda6" %in% estimates){
    res$est$freq_lambda6 <- applylambda6(cc)
    lambda6_obj <- apply(boot_cov, 1, applylambda6, callback)
    if (length(unique(round(lambda6_obj, 4))) == 1){
      res$conf$low$freq_lambda6 <- NA
      res$conf$up$freq_lambda6 <- NA
    } else{
      res$conf$low$freq_lambda6 <- quantile(lambda6_obj, probs = (1 - interval)/2, na.rm = T)
      res$conf$up$freq_lambda6 <- quantile(lambda6_obj, probs = interval + (1 - interval)/2, na.rm = T)
    }
    res$boot$lambda6 <- lambda6_obj
    if (item.dropped){
      res$ifitem$lambda6 <- apply(Ctmp, 1, applylambda6)
    }
  }

  if ("glb" %in% estimates){
    res$est$freq_glb <- glbOnArray_custom(cc)
    glb_obj <- glbOnArray_custom(boot_cov, callback)
    if (length(unique(round(glb_obj, 4))) == 1){
      res$conf$low$freq_glb <- NA
      res$conf$up$freq_glb <- NA
    } else{
      res$conf$low$freq_glb <- quantile(glb_obj, probs = (1 - interval)/2, na.rm = T)
      res$conf$up$freq_glb <- quantile(glb_obj, probs = interval + (1 - interval)/2, na.rm = T)
    }
    res$boot$glb <- glb_obj
    if (item.dropped){
      res$ifitem$glb <- glbOnArray_custom(Ctmp)
    }
  }

  #omega --------------------------------------------------------------------------
  if ("omega" %in% estimates){
    if (omega.freq.method == "cfa"){
      out <- omegaFreqData(data, interval, omega.int.analytic, pairwise, n.boot, callback, parametric)
      res$fit.object <- out$fit.object
      if (is.null(res$fit.object)) {
        if (is.null(boot_cov)) {
          boot_cov <- array(0, c(n.boot, p, p))
          for (i in 1:n.boot){
            boot_data<- MASS::mvrnorm(n, colMeans(data, na.rm = T), cc)
            boot_cov[i, , ] <- cov(boot_data)
            callback()
          }
        }
        res$est$freq_omega <- applyomega_pfa(cc)
        omega_obj <- apply(boot_cov, 1, applyomega_pfa)
        if (length(unique(round(omega_obj, 4))) == 1){
          res$conf$low$freq_omega <- NA
          res$conf$up$freq_omega <- NA
        }
        else{
          res$conf$low$freq_omega <- quantile(omega_obj, probs = (1 - interval)/2, na.rm = T)
          res$conf$up$freq_omega <- quantile(omega_obj, probs = interval + (1 - interval)/2, na.rm = T)
        }
        res$boot$omega <- omega_obj
        res$omega.error <- TRUE
        res$omega.pfa <- TRUE

        if (item.dropped){
          res$ifitem$omega <- apply(Ctmp, 1, applyomega_pfa)
        }
      } else {
        res$est$freq_omega <- out$omega
        res$loadings <- out$loadings
        res$resid_var <- out$errors
        res$conf$low$freq_omega <- out$omega_low
        res$conf$up$freq_omega <- out$omega_up
        res$omega_fit <- out$indices
        res$boot$omega <- out$omega_boot


        if (item.dropped){
          res$ifitem$omega <- numeric(p)
          for (i in 1:p) {
            dtmp <- data[, -i]
            res$ifitem$omega[i] <- applyomega_cfa_data(dtmp, interval = interval, pairwise = pairwise)
          }
          if (any(is.na(res$ifitem$omega))) {
            res$ifitem$omega <- apply(Ctmp, 1, applyomega_pfa)
            res$omega.item.error <- TRUE
          }
        }
      }
    } else if (omega.freq.method == "pfa"){
      res$est$freq_omega <- applyomega_pfa(cc)
      omega_obj <- apply(boot_cov, 1, applyomega_pfa, callback)
      res$omega.pfa <- TRUE
      if (length(unique(round(omega_obj, 4))) == 1){
        res$conf$low$freq_omega <- NA
        res$conf$up$freq_omega <- NA
      }
      else{
        res$conf$low$freq_omega <- quantile(omega_obj, probs = (1 - interval)/2, na.rm = T)
        res$conf$up$freq_omega <- quantile(omega_obj, probs = interval + (1 - interval)/2, na.rm = T)
      }
      res$boot$omega <- omega_obj
      if (item.dropped){
        res$ifitem$omega <- apply(Ctmp, 1, applyomega_pfa)
      }
    }
  }

  return(res)
}
