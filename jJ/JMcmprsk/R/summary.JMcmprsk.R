##' Joint modelling for longitutal and censored data with competing risks
##' @title Extraction of standard error and 95\% confidence interval of longitudinal/survival sub-model fixed effects
##' @param object  The JMcmprsk object returned by either jmo or jmc function.
##' @param coeff The coefficients returned by the JMcmprsk object. Results of longitudinal/survival sub-model will be printed out when typing "longitudinal" or "survival". 
##' @param digits The number of digits to print out.
##' @param ... further arguments passed to or from other methods.
##' @seealso \code{\link{jmc}}
##' @return Return standard errors of parameters with variable names
##' @references
##' \itemize{
##' \item Elashoff, Robert M., Gang Li, and Ning Li. "A joint model for longitudinal measurements and survival data in the presence of multiple failure types." Biometrics 64.3 (2008): 762-771.
##' }
##' @export
summary.JMcmprsk <- 
  function (object, coeff=c("longitudinal", "survival"), digits=4, ...) {
    if (!inherits(object, "JMcmprsk"))
      stop("Use only with 'JMcmprsk' objects.\n")
    if (object$type == "jmo") {
      if (coeff == "longitudinal") {
        ##Estimates of betas
        Estimate <- object$betas
        SE <- object$se_betas
        LowerLimit <- Estimate - 1.96 * SE
        UpperLimit <- Estimate + 1.96 * SE
        zval = (Estimate/SE)
        pval = 2 * pnorm(-abs(zval))
        out <- data.frame(Estimate, SE, LowerLimit, UpperLimit, pval)
        out <- cbind(rownames(out), out)
        rownames(out) <- NULL
        names(out) <- c("Longitudinal", "coef", "SE", "95%Lower", "95%Upper", "p-values")
        
        
        ##Estimates of alphas
        Estimate <- t(object$alphamatrix)
        Estimate <- reshape2::melt(Estimate)
        SE <- object$se_alphas
        LowerLimit <- Estimate[, 3] - 1.96 * SE
        UpperLimit <- Estimate[, 3] + 1.96 * SE
        zval = (Estimate[, 3]/SE)
        pval = 2 * pnorm(-abs(zval))
        outalphas <- data.frame(Estimate, SE, LowerLimit, UpperLimit, pval)
        outalphas[, 1] <- paste(outalphas[, 1], (outalphas[, 2] + 1), sep = "_")
        outalphas <- outalphas[, -2]
        names(outalphas) <- c("Longitudinal", "coef", "SE", "95%Lower", "95%Upper", "p-values")
        
        ##Estimates of thetas
        Estimate <- object$thetas
        SE <- object$se_thetas
        LowerLimit <- Estimate - 1.96 * SE
        UpperLimit <- Estimate + 1.96 * SE
        zval = (Estimate/SE)
        pval = 2 * pnorm(-abs(zval))
        thetaname <- NULL
        for (i in 1:length(Estimate)) {
          thetaname <- c(thetaname, paste0("theta", i))
        }
        outthetas <- data.frame(thetaname, Estimate, SE, LowerLimit, UpperLimit, pval)
        names(outthetas) <- c("Longitudinal", "coef", "SE", "95%Lower", "95%Upper", "p-values")
        out <- rbind(out, outalphas, outthetas)
        
        ##round p-val (TBD Hong Wang)
        out[, 2:ncol(out)] <- round(out[, 2:ncol(out)], digits = digits)
        
        return(out)
      } else if (coeff == "survival") {
        Estimate <- t(object$gamma_matrix)
        Estimate <- reshape2::melt(Estimate)
        SE <- t(object$se_gamma_matrix)
        SE <- reshape2::melt(SE)
        LowerLimit <- Estimate[, 3] - 1.96 * SE[, 3]
        UpperLimit <- Estimate[, 3] + 1.96 * SE[, 3]
        zval = (Estimate[, 3]/SE[, 3])
        pval = 2 * pnorm(-abs(zval))
        out <- data.frame(Estimate, exp(Estimate[, 3]), SE[, 3], LowerLimit, UpperLimit, pval)
        out[, 1] <- paste(out[, 1], out[, 2], sep = "_")
        out <- out[, -2]
        names(out) <- c("Survival", "coef", "exp(coef)", "SE(coef)", "95%Lower", "95%Upper", "p-values")
        ##round p-val (TBD Hong Wang)
        out[, 2:ncol(out)] <- round(out[, 2:ncol(out)], digits = digits)
        return(out)
      }  else {
        stop("Unexpected arguments! Must choose one of the following options: longitudinal, survival")
      }
    } else if (object$type == "jmc") {
      if (coeff == "longitudinal") {
        ##Estimates of betas
        Estimate <- object$betas
        SE <- object$se_betas
        LowerLimit <- Estimate - 1.96 * SE
        UpperLimit <- Estimate + 1.96 * SE
        zval = (Estimate/SE)
        pval = 2 * pnorm(-abs(zval))
        out <- data.frame(Estimate, SE, LowerLimit, UpperLimit, pval)
        out <- cbind(rownames(out), out)
        names(out)[1] <- "Longitudinal"
        rownames(out) <- NULL
        names(out)[4:6] <- c("95%Lower", "95%Upper", "p-values")
        ##round p-val (TBD Hong Wang)
        out[, 2:ncol(out)] <- round(out[, 2:ncol(out)], digits = digits)
        return(out)
      } else if (coeff == "survival") {
        Estimate <- t(object$gamma_matrix)
        Estimate <- reshape2::melt(Estimate)
        SE <- t(object$se_gamma_matrix)
        SE <- reshape2::melt(SE)
        LowerLimit <- Estimate[, 3] - 1.96 * SE[, 3]
        UpperLimit <- Estimate[, 3] + 1.96 * SE[, 3]
        zval = (Estimate[, 3]/SE[, 3])
        pval = 2 * pnorm(-abs(zval))
        out <- data.frame(Estimate, exp(Estimate[, 3]), SE[, 3], LowerLimit, UpperLimit, pval)
        out[, 1] <- paste(out[, 1], out[, 2], sep = "_")
        out <- out[, -2]
        names(out) <- c("Survival", "coef", "exp(coef)", "SE(coef)", "95%Lower", "95%Upper", "p-values")
        ##round p-val (TBD Hong Wang)
        out[, 2:ncol(out)] <- round(out[, 2:ncol(out)], digits = digits)
        return(out)
      } else {
        stop("Unexpected arguments! Must choose one of the following options: longitudinal, survival")
      }
    }
  }