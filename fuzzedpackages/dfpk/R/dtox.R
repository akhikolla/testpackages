#' @import ggplot2
#' @import rstan
#' @import Rcpp
#' @import methods
#' @import stats
#' @useDynLib dfpk, .registration = TRUE
#' @export
dtox <-
function(y, doses, x, theta, prob = 0.9, options=list(nchains = 4, niter = 4000, nadapt = 0.8), 
         betapriors = c(0, 16.71, 0, 6.43), thetaL=NULL, auc = NULL, deltaAUC = NULL, p0 = NULL, L = NULL, CI = TRUE){
        
        checking1 <- function(x,target,error){
            sum(x>(target+error))/length(x)              
        }
        
        num <- length(x)  	# how many patients
        dose1 <- log(doses[x])
        # For STAN model
        
        data_s <- list(N=num, y=y, dose=dose1, beta0mean=c(betapriors[1], betapriors[2]), beta1mean=c(betapriors[3], betapriors[4]))
        sm_lrDtox <- stanmodels$cdf_reg_dtox
        reg1 <- sampling(sm_lrDtox, data=data_s, iter=options$niter, chains=options$nchains, control = list(adapt_delta = options$nadapt))
        a1 = get_posterior_mean(reg1)
        sampl1 <- extract(reg1)
        
        beta <- a1[1:2, options$nchains + 1]
        
        beta0 <- -beta[1]
        beta1 <- beta[2]
        pnew <- pnorm(beta0 + beta1*log(doses))
        
        Beta0 <- -sampl1$beta0
        Beta1 <- sampl1$beta1
        pstim_sum <- matrix(0, ncol = options$nchains*options$niter/2, nrow = length(doses))
        p_sum <- NULL 
        

        for(i in 1:ncol(pstim_sum)){
            pstim_sum[1,i] <- pnorm(Beta0[i] + Beta1[i]*log(doses[1])) 
        }

        #######################
        #### Stopping Rule ####
        #######################

        pstop <- checking1(pstim_sum[1,], target=theta, error=0)
        stoptox <- (pstop >= prob)
        stoptrial <- stoptox


        if(CI == "TRUE"){
            p_sum <- summary(pstim_sum[1,])
            for(o in 2:length(doses)){
                for(i in 1:ncol(pstim_sum)){
                    pstim_sum[o,i] <- pnorm(Beta0[i] + Beta1[i]*log(doses[o])) 
                }
                p_sum <- rbind(p_sum, summary(pstim_sum[o,]))
            }
        }else{
            p_sum <- NULL
        }
        
        
        # check if we will stop the trial or not
        
        if (stoptrial){
            newDose = NA
            message("The trial stopped based on the stopping rule \n \n")
        }else{                                      # if we don't stopped
            newDose <- order(abs(pnew-theta))[1]
        }
        
        parameters <- beta 
        names(parameters) <- c("beta0", "beta1")
        list(newDose = newDose, pstim = pnew, p_sum = p_sum, parameters = parameters)
    }
