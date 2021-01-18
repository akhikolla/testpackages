#' @import ggplot2
#' @import rstan
#' @import Rcpp
#' @import methods
#' @import stats
#' @useDynLib dfpk, .registration = TRUE
#' @export
pktox <-
function(y, auc, doses, x, theta, prob = 0.9, options = list(nchains = 4, niter = 4000, nadapt = 0.8), 
         betapriors = c(10, 10000, 20, 10), thetaL=NULL, p0=NULL, L=NULL, deltaAUC=NULL, CI = TRUE){
        
        checking1 <- function(x,target,error){
            sum(x>(target+error))/length(x)             
        }
        
        f <- function(v,lambda,parmt){
            pnorm(lambda[1]+lambda[2]*v)*dnorm(v,parmt[1],parmt[2])
        }
        
        f2 <- function(v,lambda1, lambda2, parmt1, parmt2){
            pnorm(lambda1+lambda2*v)*dnorm(v,parmt1,parmt2)
        }
        
        num <- length(x)           # how many patients
        dose1 <- cbind(rep(1,num), log(doses[x]))
        mu1 <- -log(betapriors[1])
        
        # For STAN model
        data_s <- list(N=num, auc=log(auc), dose=dose1, mu = mu1, beta0=betapriors[2])
        sm_lrauc <- stanmodels$reg_auc
        reg1 <- sampling(sm_lrauc, data=data_s, iter=options$niter, chains=options$nchains,
                         control = list(adapt_delta = options$nadapt))
        sampl1 <- extract(reg1)
        a1=get_posterior_mean(reg1)
        beta1 <- a1[1:2,options$nchains+1]
        nu <- a1[3,options$nchains+1]
        auc1 <- log(auc)
        
        # For STAN model
        data_s <- list(N=num, y=y, dose=auc1, beta2mean = betapriors[3], beta3mean = betapriors[4])
        ## auc1 is log(AUCs)
        sm_lr <- stanmodels$cdf_reg_pktox
        reg2 <- sampling(sm_lr, data=data_s, iter=options$niter, chains=options$nchains, control = list(adapt_delta = options$nadapt))
        sampl2 <- extract(reg2)
        a2 = get_posterior_mean(reg2)
        
        Beta <- a2[1:2,options$nchains+1]
        Beta2 <- Beta[1]
        Beta3 <- Beta[2]

        # Computation probability 
        pstim <- NULL
        pstim_sum <- matrix(0, ncol = options$nchains*options$niter/2, nrow = length(doses))
        p_sum <- NULL
        for (o in 1:length(doses)){
            parmt = c(a1[1,options$nchains+1] + a1[2,options$nchains+1]*log(doses[o]), a1[3,options$nchains+1])
            pstim <- c(pstim, integrate(f,-Inf, Inf, lambda=c(-Beta2, Beta3), parmt=parmt)$value)
        }

        ## Calculating the credible interval for one sample to check the stopping rule ##

        parmt1 = sampl1$b[,1] + sampl1$b[,2]*log(doses[1])
        parmt2 = sampl1$sigma
        for(i in 1:ncol(pstim_sum)){
            pstim_sum[1,i] <- integrate(f2,-Inf, Inf, lambda1 = -sampl2$beta2[i], lambda2 = sampl2$beta3[i], 
                                        parmt1 = parmt1[i], parmt2 = parmt2[i])$value
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
                parmt1 = sampl1$b[,1] + sampl1$b[,2]*log(doses[o])
                parmt2 = sampl1$sigma
                for(i in 1:ncol(pstim_sum)){
                    pstim_sum[o,i] <- integrate(f2,-Inf, Inf, lambda1 = -sampl2$beta2[i], lambda2 = sampl2$beta3[i], 
                                                parmt1 = parmt1[i], parmt2 = parmt2[i])$value
                }
                p_sum <- rbind(p_sum, summary(pstim_sum[o,]))
            }
        }else{
            p_sum <- NULL
        }
        
        
        # check if we will stop the trial or not
        
        if (stoptrial){
            newDose= NA 
            message("The trial stopped based on the stopping rule \n \n")
        }else{                # if we not stop
            newDose = order(abs(pstim-theta))[1]
        }
        
        # newDose = order(abs(pstim-theta))[1]
        parameters <- c(beta1, nu, Beta)
        names(parameters) <- c("beta0", "beta1", "nu", "beta2", "beta3")
        list(newDose = newDose, pstim = pstim, p_sum = p_sum, parameters = parameters)
    }
