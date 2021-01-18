##' Print contents of JMcmprsk object.
##'
##'
##' @title Print JMcmprsk
##' @param x Object of class 'JMcmprsk'.
##' @param digits  The desired number of digits after the decimal point. Default is 4.
##' @param ... Further arguments passed to or from other methods.
##' @seealso \code{\link{jmc}}
##' @author Hong Wang
##' @export
print.JMcmprsk <- function(x,digits=4, ...) {
  if (!inherits(x, "JMcmprsk")) 
    stop("Not a legitimate \"JMcmprsk\" object")
  cat("\nCall:\n", sprintf(format(paste(deparse(x$call, width.cutoff = 500), collapse=""))), "\n\n")
  if (x$type == "jmc") {
    
    ##need to add function call (Hong)
    cat("Data Summary:\n")
    cat("Number of observations:", x$SummaryInfo$Numobs, "\n")
    cat("Number of groups:", x$k, "\n\n")
    cat("Proportion of competing risks: \n")
    for (i in 1:2) {
      cat("Risk", i, ":", x$SummaryInfo$PropComp[i+1], "%\n")
    }
    cat("\nNumerical intergration:\n")
    cat("Method: standard Guass-Hermite quadrature\n")
    cat("Number of quadrature points: ", x$point, "\n")
    cat("\nModel Type: joint modeling of longitudinal continuous and competing risks data", "\n\n")
    cat("Model summary:\n")
    cat("Longitudinal process: linear mixed effects model\n")
    cat("Event process: cause-specific Cox proportional hazard model with unspecified baseline hazard\n\n")
    cat("Loglikelihood: ", x$loglike, "\n\n")
    cat("Longitudinal sub-model fixed effects: ", 
        sprintf(format(paste(deparse(x$SummaryInfo$LongitudinalSubmodel, width.cutoff = 500), collapse=""))), "\n")
    cat("\n                  Estimate   Std. Error       95% CI            Pr(>|Z|)    \n")
    cat("Longitudinal:                \n")
    cat(" Fixed effects:                 \n")
    
    for (i in 1:length(x$betas)) {
      #beta = paste0("beta", i)
      beta=names(x$betas)[i]
      uppsd = x$betas[i] + 1.96 * x$se_betas[i]
      lowersd = x$betas[i] - 1.96 * x$se_betas[i]
      zval = (x$betas[i]/x$se_betas[i])
      pval = 2 * pnorm(-abs(zval))
	 
	  
	  cat(" ",formatC(beta,width=14,flag="-"),  formatC(x$betas[i], digits = digits, width=10,format = "f",flag="-"))	  
	 
      cat("  ",  formatC(x$se_betas[i], digits = digits, width=10, format = "f",flag="-"))
	  
      cat("  (", formatC(lowersd, digits = digits, format = "f",flag=" ") ,",",formatC(uppsd, digits = digits, format = "f",flag=" "),")", sep = "")
      cat("     ", formatC(pval, digits = digits, width=10,format = "f",flag="-"))
      cat("\n")
    }
    cat("Random effects:                 \n")
    cat(" ",formatC("sigma^2",width=14,flag="-"), sprintf("% 1.4f", x$sigma2_val))
    cat("     ", sprintf("% 1.4f", x$se_sigma2_val))
    
    uppsd = x$sigma2_val + 1.96 * x$se_sigma2_val
    lowersd = x$sigma2_val - 1.96 * x$se_sigma2_val
    zval = (x$sigma2_val/x$se_sigma2_val)
    pval = 2 * pnorm(-abs(zval))
    cat("  (", formatC(lowersd, digits = digits, format = "f",flag=" ") ,",",formatC(uppsd, digits = digits, format = "f",flag=" "),")", sep = "")
    cat("     ", formatC(pval, digits = digits, width=10,format = "f",flag="-"))
    cat("\n")  
        
    
    cat("\nSurvival sub-model fixed effects: ", 
        sprintf(format(paste(deparse(x$SummaryInfo$SurvivalSubmodel, width.cutoff = 500), collapse=""))), "\n")
    cat("\n                  Estimate   Std. Error       95% CI           Pr(>|Z|)    \n")
    cat("Survival:                \n")
    cat(" Fixed effects:                 \n")
    for (i in 1:dim(x$gamma_matrix)[1]) for (j in 1:dim(x$gamma_matrix)[2]) {
      gamma = paste0("gamma", i, j)
      gamma=colnames(x$gamma_matrix)[j]
      gamma=paste0(gamma,'_', i)
      gammaval = x$gamma_matrix[i, j]
      stdgammaval = x$se_gamma_matrix[i, j]
      uppsd = gammaval + 1.96 * stdgammaval
      lowersd = gammaval - 1.96 * stdgammaval
      zval = (gammaval/stdgammaval)
      pval = 2 * pnorm(-abs(zval))
      cat(" ",formatC(gamma,width=14,flag="-"),formatC(gammaval, digits = digits, width=10, format = "f",flag="-"))
	  
      cat("  ",  formatC(stdgammaval, digits = digits, width=10, format = "f",flag="-"))	  
      cat("  (", formatC(lowersd, digits = digits, format = "f",flag=" ") ,",",formatC(uppsd, digits = digits, format = "f",flag=" "),")", sep = "")
      cat("     ", formatC(pval, digits = digits, width=10,format = "f",flag="-"))
      cat("\n")
    }
    cat("\nAssociation parameter:                 \n")
	cat(" ",formatC("v2",width=14,flag="-"), formatC(x$v_estimate, digits = digits, width=10, format = "f",flag="-"))
	cat("  ",  formatC(x$se_v_estimate, digits = digits, width=10, format = "f",flag="-"))
	      
    uppsd = x$v_estimate + 1.96 * x$se_v_estimate
    lowersd = x$v_estimate - 1.96 * x$se_v_estimate
    zval = (x$v_estimate/x$se_v_estimate)
    pval = 2 * pnorm(-abs(zval))
    cat("  (", formatC(lowersd, digits = digits, format = "f",flag=" ") ,",",formatC(uppsd, digits = digits, format = "f",flag=" "),")", sep = "")
    cat("     ", formatC(pval, digits = digits, width=10,format = "f",flag="-"))
    cat("\n")
    
    # cat(' Random effects: \n')

    # restore a matrix from its uppertri
    sd_sigmamatrix = matrix(0, dim(x$sigma_matrix)[1], dim(x$sigma_matrix)[1])
    sd_sigmamatrix[lower.tri(sd_sigmamatrix, diag = TRUE)] <- x$se_sigma
    sd_sigmamatrix <- t(sd_sigmamatrix)
    
    # print sigmabii
    cat("\nRandom effects:                 \n")
    for (i in 1:(dim(x$sigma_matrix)[1] - 1)) # for (j in 1:dim(x$sigma_matrix)[2])
    {
      sigma = paste0("sigma_b", i, i)
      sigmaval = x$sigma_matrix[i, i]
      stdsigmaval = sd_sigmamatrix[i, i]
      uppsd = sigmaval + 1.96 * stdsigmaval
      lowersd = sigmaval - 1.96 * stdsigmaval
      zval = (sigmaval/stdsigmaval)
      pval = 2 * pnorm(-abs(zval))
	  cat(" ",formatC(sigma,width=14,flag="-"), formatC(sigmaval, digits = digits, format = "f",flag="-"))
      cat("      ", 	  formatC(stdsigmaval, digits = digits,width=10, format = "f",flag="-"))
      cat("  (", formatC(lowersd, digits = digits, format = "f",flag=" ") ,",",formatC(uppsd, digits = digits, format = "f",flag=" "),")", sep = "")
      cat("     ", formatC(pval, digits = digits, width=10,format = "f",flag="-"))
      cat("\n")
    }
    # print sigmau
    si = dim(x$sigma_matrix)[1]
    sigmaval = x$sigma_matrix[si, si]
    stdsigmaval = x$se_sigma[(si * (si + 1))/2]
    uppsd = sigmaval + 1.96 * stdsigmaval
    lowersd = sigmaval - 1.96 * stdsigmaval
    zval = (sigmaval/stdsigmaval)
    pval = 2 * pnorm(-abs(zval))
    cat(" ",formatC("sigma_u",width=14,flag="-"), formatC(sigmaval, digits = digits, format = "f",flag="-"))
	cat("      ", formatC(stdsigmaval, digits = digits, width=10,format = "f",flag="-"))
	cat("  (", formatC(lowersd, digits = digits, format = "f",flag=" ") ,",",formatC(uppsd, digits = digits, format = "f",flag=" "),")", sep = "")
    cat("     ", formatC(pval, digits = digits, width=10,format = "f",flag="-"))
    cat("\n")
    
    cat(" Covariance:                 \n")
    
    
    # print sigmabii
    
    for (i in 1:(dim(x$sigma_matrix)[1] - 1)) for (j in (i + 
      1):(dim(x$sigma_matrix)[2])) {
      if (j < dim(x$sigma_matrix)[2]) {
        sigma = paste0("sigma_b", i, j)
      } else {
        sigma = paste0("sigma_b", i, "u")
      }
      
      sigmaval = x$sigma_matrix[i, j]
      stdsigmaval = sd_sigmamatrix[i, j]
      uppsd = sigmaval + 1.96 * stdsigmaval
      lowersd = sigmaval - 1.96 * stdsigmaval
      zval = (sigmaval/stdsigmaval)
      pval = 2 * pnorm(-abs(zval))
      cat(" ",formatC(sigma,width=14,flag="-"),formatC(sigmaval, digits = digits, format = "f",flag="-"))
      cat("     ", formatC(stdsigmaval, digits = digits, format = "f",flag="-"))
	  
      cat("      (", formatC(lowersd, digits = digits, format = "f",flag=" ") ,",",formatC(uppsd, digits = digits, format = "f",flag=" "),")", sep = "")
       cat("     ", formatC(pval, digits = digits, width=10,format = "f",flag="-"))
	   
      cat("\n")
      
    }
    
  }
  if (x$type == "jmo") {
   
    cat("Data Summary:\n")
    cat("Number of observations:", x$SummaryInfo$Numobs, "\n")
    cat("Number of groups:", x$k, "\n\n")
    cat("Proportion of competing risks: \n")
    for (i in 1:2) {
      cat("Risk", i, ":", x$SummaryInfo$PropComp[i+1], "%\n")
    }
    cat("\nNumerical intergration:\n")
    cat("Method: Standard Guass-Hermite quadrature\n")
    cat("Number of quadrature points: ", x$point, "\n")
    cat("\nModel Type: joint modeling of longitudinal ordinal and competing risks data", "\n\n")
    cat("Model summary:\n")
    cat("Longitudinal process: partial proportional odds model\n")
    cat("Event process: cause-specific Cox proportional hazard model with unspecified baseline hazard\n\n")
    cat("Loglikelihood: ", x$loglike, "\n\n")
    cat("Longitudinal sub-model proportional odds: ", 
        sprintf(format(paste(deparse(x$SummaryInfo$LongitudinalSubmodel, width.cutoff = 500), collapse=""))), "\n")
    cat("Longitudinal sub-model non-proportional odds:", sprintf(x$SummaryInfo$LongNP), "\n")
    cat("\n                  Estimate   Std. Error       95% CI            Pr(>|Z|)    \n")
    cat("Longitudinal:                \n")
    cat(" Fixed effects:                 \n")
    cat("  proportional odds:\n")
    for (i in 1:length(x$betas)) {
      #beta = paste0("beta", i)
      beta=names(x$betas)[i]
      uppsd = x$betas[i] + 1.96 * x$se_betas[i]
      lowersd = x$betas[i] - 1.96 * x$se_betas[i]
      zval = (x$betas[i]/x$se_betas[i])
      pval = 2 * pnorm(-abs(zval))
     
	  cat(" ",formatC(beta,width=14,flag="-"),  formatC(x$betas[i], digits = digits, width=10,format = "f",flag="-"))	  
	 
      cat("  ",  formatC(x$se_betas[i], digits = digits, width=10, format = "f",flag="-"))
	  
      cat("  (", formatC(lowersd, digits = digits, format = "f",flag=" ") ,",",formatC(uppsd, digits = digits, format = "f",flag=" "),")", sep = "")
      cat("     ", formatC(pval, digits = digits, width=10,format = "f",flag="-"))
      cat("\n")  
	  
	  
    }
    cat("  Non-proportional odds:\n")
      for (i in 1:dim(x$alphamatrix)[1]) for (j in 1:dim(x$alphamatrix)[2]) {
      #alpha = paste0("alpha", i+1, j)
      alpha=colnames(x$alphamatrix)[j]
      alpha=paste0(alpha,'_', i+1)
      
      alphaval = x$alphamatrix[i, j]
      stdalphaval = x$se_alphas[(i-1)*dim(x$alphamatrix)[2]+j]
      uppsd = alphaval + 1.96 * stdalphaval
      lowersd = alphaval - 1.96 * stdalphaval
      zval = (alphaval/stdalphaval)
      pval = 2 * pnorm(-abs(zval))

	  cat(" ",formatC(alpha,width=14,flag="-"), formatC(alphaval, digits = digits, width=10,format = "f",flag="-"))	  
	 
      cat("  ",  formatC(stdalphaval, digits = digits, width=10, format = "f",flag="-"))
	  
      cat("  (", formatC(lowersd, digits = digits, format = "f",flag=" ") ,",",formatC(uppsd, digits = digits, format = "f",flag=" "),")", sep = "")
      cat("     ", formatC(pval, digits = digits, width=10,format = "f",flag="-"))
      cat("\n")
	  
	  
	  
	  
    }
    cat("  Logit-specific intercepts: \n")
	for (i in 1:length(x$thetas)) {
      theta = paste0("theta", i)
      uppsd = x$thetas[i] + 1.96 * x$se_thetas[i]
      lowersd = x$thetas[i] - 1.96 * x$se_thetas[i]
      zval = (x$thetas[i]/x$se_thetas[i])
      pval = 2 * pnorm(-abs(zval))
     	  
	  cat(" ",formatC(theta,width=14,flag="-"), formatC(x$thetas[i], digits = digits, width=10,format = "f",flag="-"))	  
	 
      cat("  ",  formatC(x$se_thetas[i], digits = digits, width=10, format = "f",flag="-"))
	  
      cat("  (", formatC(lowersd, digits = digits, format = "f",flag=" ") ,",",formatC(uppsd, digits = digits, format = "f",flag=" "),")", sep = "")
      cat("     ", formatC(pval, digits = digits, width=10,format = "f",flag="-"))
      cat("\n")
	  
	  
	  
    }
    
    cat("\nSurvival sub-model fixed effects: ", 
        sprintf(format(paste(deparse(x$SummaryInfo$SurvivalSubmodel, width.cutoff = 500), collapse=""))), "\n")
    cat("\n                  Estimate   Std. Error       95% CI            Pr(>|Z|)    \n")
    cat("Survival:                \n")
    cat(" Fixed effects:                 \n")
    for (i in 1:dim(x$gamma_matrix)[1]) for (j in 1:dim(x$gamma_matrix)[2]) {
      gamma = paste0("gamma", i, j)
      gamma=colnames(x$gamma_matrix)[j]
      gamma=paste0(gamma,'_', i)
      gammaval = x$gamma_matrix[i, j]
      stdgammaval = x$se_gamma_matrix[i, j]
      uppsd = gammaval + 1.96 * stdgammaval
      lowersd = gammaval - 1.96 * stdgammaval
      zval = (gammaval/stdgammaval)
      pval = 2 * pnorm(-abs(zval))
      
      cat(" ",formatC(gamma,width=14,flag="-"),formatC(gammaval, digits = digits, width=10, format = "f",flag="-"))
	  
      cat("  ",  formatC(stdgammaval, digits = digits, width=10, format = "f",flag="-"))	  
      cat("  (", formatC(lowersd, digits = digits, format = "f",flag=" ") ,",",formatC(uppsd, digits = digits, format = "f",flag=" "),")", sep = "")
      cat("     ", formatC(pval, digits = digits, width=10,format = "f",flag="-"))
      cat("\n")
	  
	  
    }
    cat("\nAssociation prameter:                 \n")
    
	cat(" ",formatC("v2",width=14,flag="-"), formatC(x$v_estimate, digits = digits, width=10, format = "f",flag="-"))
	cat("  ",  formatC(x$se_v_estimate, digits = digits, width=10, format = "f",flag="-"))
	
    uppsd = x$v_estimate + 1.96 * x$se_v_estimate
    lowersd = x$v_estimate - 1.96 * x$se_v_estimate
    zval = (x$v_estimate/x$se_v_estimate)
    pval = 2 * pnorm(-abs(zval))
    cat("  (", formatC(lowersd, digits = digits, format = "f",flag=" ") ,",",formatC(uppsd, digits = digits, format = "f",flag=" "),")", sep = "")
    cat("     ", formatC(pval, digits = digits, width=10,format = "f",flag="-"))
    cat("\n")
    
    # cat(' Random effects: \n')
    
    # restore a matrix from its uppertri
    sd_sigmamatrix = matrix(0, dim(x$sigma_matrix)[1], dim(x$sigma_matrix)[1])
    sd_sigmamatrix[lower.tri(sd_sigmamatrix, diag = TRUE)] <- x$se_sigma
    sd_sigmamatrix <- t(sd_sigmamatrix)
    cat("\nRandom effects:                 \n")
    # print sigmabii
    for (i in 1:(dim(x$sigma_matrix)[1] - 1)) # for (j in 1:dim(x$sigma_matrix)[2])
    {
      sigma = paste0("sigma_b", i, i)
      sigmaval = x$sigma_matrix[i, i]
      stdsigmaval = sd_sigmamatrix[i, i]
      uppsd = sigmaval + 1.96 * stdsigmaval
      lowersd = sigmaval - 1.96 * stdsigmaval
      zval = (sigmaval/stdsigmaval)
      pval = 2 * pnorm(-abs(zval))
		  
	  cat(" ",formatC(sigma,width=14,flag="-"), formatC(sigmaval, digits = digits, format = "f",flag="-"))
      cat("      ", 	  formatC(stdsigmaval, digits = digits,width=10, format = "f",flag="-"))
      cat("  (", formatC(lowersd, digits = digits, format = "f",flag=" ") ,",",formatC(uppsd, digits = digits, format = "f",flag=" "),")", sep = "")
      cat("     ", formatC(pval, digits = digits, width=10,format = "f",flag="-"))
      cat("\n")
	  
    }
    # print sigmau
    si = dim(x$sigma_matrix)[1]
    sigmaval = x$sigma_matrix[si, si]
    stdsigmaval = x$se_sigma[(si * (si + 1))/2]
    uppsd = sigmaval + 1.96 * stdsigmaval
    lowersd = sigmaval - 1.96 * stdsigmaval
    zval = (sigmaval/stdsigmaval)
    pval = 2 * pnorm(-abs(zval))
     
	cat(" ",formatC("sigma_u",width=14,flag="-"), formatC(sigmaval, digits = digits, format = "f",flag="-"))
	cat("      ", formatC(stdsigmaval, digits = digits, width=10,format = "f",flag="-"))
	cat("  (", formatC(lowersd, digits = digits, format = "f",flag=" ") ,",",formatC(uppsd, digits = digits, format = "f",flag=" "),")", sep = "")
    cat("     ", formatC(pval, digits = digits, width=10,format = "f",flag="-"))
    cat("\n")
	
	
	
    cat(" Covariance:                 \n")
    
    
    # print sigmabii
    
    for (i in 1:(dim(x$sigma_matrix)[1] - 1)) for (j in (i + 
                                                         1):(dim(x$sigma_matrix)[2])) {
      if (j < dim(x$sigma_matrix)[2]) {
        sigma = paste0("sigma_b", i, j)
      } else {
        sigma = paste0("sigma_b", i, "u")
      }
      
      sigmaval = x$sigma_matrix[i, j]
      stdsigmaval = sd_sigmamatrix[i, j]
      uppsd = sigmaval + 1.96 * stdsigmaval
      lowersd = sigmaval - 1.96 * stdsigmaval
      zval = (sigmaval/stdsigmaval)
      pval = 2 * pnorm(-abs(zval))
      
	  cat(" ",formatC(sigma,width=14,flag="-"),formatC(sigmaval, digits = digits, format = "f",flag="-"))
      cat("     ", formatC(stdsigmaval, digits = digits, format = "f",flag="-"))
	  
      cat("      (", formatC(lowersd, digits = digits, format = "f",flag=" ") ,",",formatC(uppsd, digits = digits, format = "f",flag=" "),")", sep = "")
       cat("     ", formatC(pval, digits = digits, width=10,format = "f",flag="-"))
	   
      cat("\n")
    }
    
    
  }
}

