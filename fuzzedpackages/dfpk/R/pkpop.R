#' @import ggplot2
#' @import rstan
#' @import Rcpp
#' @import methods
#' @import stats
#' @useDynLib dfpk, .registration = TRUE
#' @export
pkpop <-
function(y, auc, doses, x, theta, prob = 0.9, options = list(nchains = 4, niter = 4000, nadapt = 0.8), 
         betapriors = c(10, 10000, 10, 5), thetaL=NULL, p0 = NULL, L = NULL, deltaAUC = NULL, CI = TRUE){
        
        checking1 <- function(x,target,error){
            sum(x>(target+error))/length(x)              
        }
        
        num <- length(x)                           # how many patients
        dose1 <- cbind(rep(1,num), log(doses[x]))
        mu1 <- -log(betapriors[1])
        
        # For STAN model
        data_s <- list(N=num, auc=log(auc),dose=dose1, mu = mu1, beta0=betapriors[2])
        sm_lrauc <- stanmodels$reg_auc
        reg1 <- sampling(sm_lrauc, data=data_s, iter=options$niter, chains=options$nchains,
                         control = list(adapt_delta = options$nadapt))
        
        a1 = get_posterior_mean(reg1)
        sampl1 <- extract(reg1)
        
        beta1 <- a1[1:2,options$nchains+1]
        mu <- beta1[1] + beta1[2]*log(doses)
        nu <- a1[3,options$nchains+1]
        
        mu1 <- mu[x]
        
        # For STAN model
        data_s <- list(N=num,y=y,dose=mu1, beta3mean = betapriors[3], beta4mean = betapriors[4])
        sm_lrPkpop <- stanmodels$logit_reg_pkpop
        reg2 <- sampling(sm_lrPkpop, data=data_s, iter=options$niter, chains=options$nchains,
                         control = list(adapt_delta = options$nadapt))
        a2 = get_posterior_mean(reg2)
        sampl2 <- extract(reg2)
        
        Beta <- a2[1:2,options$nchains+1]
        
        # Computation probability
        pstim = invlogit(-Beta[1] + Beta[2]*mu)
        
        b1 <- sampl1$b[,1]
        b2 <- sampl1$b[,2]
        n <- sampl1$sigma 
        Beta1 <- -sampl2$beta3
        Beta2 <- sampl2$beta4
        
        pstim_sum <- matrix(0, ncol = options$nchains*options$niter/2, nrow = length(doses))
        p_sum <- NULL
        m <- NULL

        m <- sampl1$b[,1] + sampl1$b[,2]*log(doses[1])
        for(i in 1:ncol(pstim_sum)){
            pstim_sum[1,i] <- invlogit(Beta1[i] + Beta2[i]*m[i])
        }

        #######################
        #### Stopping Rule ####
        #######################

        pstop <-  checking1(pstim_sum[1,], target=theta, error=0)
        stoptox <- (pstop >= prob)
        stoptrial <- stoptox



        if(CI == "TRUE"){
            p_sum <- summary(pstim_sum[1,])
            for(o in 2:length(doses)){
                m <- sampl1$b[,1] + sampl1$b[,2]*log(doses[o])
                for(i in 1:ncol(pstim_sum)){
                    pstim_sum[o,i] <- invlogit(Beta1[i] + Beta2[i]*m[i])
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
        }else{                                          # if we not stop
            newDose <- order((abs(pstim-theta)))[1]
        }
        
        
        # newDose = order((abs(pstim-theta)))[1]
        parameters <- c(beta1,nu,Beta)
        names(parameters) <- c("beta0", "beta1", "nu", "beta3", "beta4")
        list(newDose = newDose, pstim = pstim, p_sum=p_sum, parameters = parameters)
    }
