#' @import ggplot2
#' @import rstan
#' @import Rcpp
#' @import dfcrm
#' @import methods
#' @import stats
#' @useDynLib dfpk, .registration = TRUE
#' @export
pkcrm <-
function(y, auc, doses, x, theta, p0, L, prob = 0.9, options = list(nchains = 4, niter = 4000, nadapt = 0.8), 
         betapriors = c(10, 10000), thetaL=NULL, deltaAUC = NULL, CI = TRUE){
        
        checking1 <- function(x,target,error){
            sum(x>(target+error))/length(x)              
        }
        
        num <- length(x)    		     # how many patients
        dose1 <- cbind(rep(1,num), log(doses[x]))
        mu1 <- -log(betapriors[1])
        
        # For STAN
        data_s <- list(N=num, auc=log(auc), dose=dose1, mu = mu1, beta0=betapriors[2])
        sm_lrauc <- stanmodels$reg_auc
        reg1 <- sampling(sm_lrauc, data=data_s, iter=options$niter, chains=options$nchains,
                         control = list(adapt_delta = options$nadapt))
        a1=get_posterior_mean(reg1)
        sampl1 <- extract(reg1)
        
        beta1 <- c(a1[1,options$nchains+1], a1[2,options$nchains+1])
        nu <- a1[3,options$nchains+1]
        mu <- beta1[1] + beta1[2]*log(doses)
        
        results_crm <- crm(p0,theta,y,x)$mtd
        
        # Computation probability
        p_new <- round(1-pnorm((L-mu)/sqrt(nu)),options$nchains+1)
        
        ## Posterior probabilities
        b1 <- sampl1$b[,1]
        b2 <- sampl1$b[,2]
        n <- sampl1$sigma
        m <- NULL
        
        pstim_sum <- matrix(0, ncol = options$nchains*options$niter/2, nrow = length(doses))
        p_sum <- NULL
        
        m <- b1 + b2*log(doses[1])
        for(i in 1:ncol(pstim_sum)){
            pstim_sum[1,i] <- round(1-pnorm((L-m[i])/sqrt(n[i])), options$nchains+1)
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
                m <- b1 + b2*log(doses[o])
                for(i in 1:ncol(pstim_sum)){
                    pstim_sum[o,i] <- round(1-pnorm((L-m[i])/sqrt(n[i])), options$nchains+1)
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
        }else{             # if we not stop
            if(is.null(thetaL) == FALSE){
                result_safety <- order(abs(p_new - thetaL))[1]
                newDose = min(results_crm, result_safety)
            }else{
                result_safety <- order(abs(p_new - theta))[1]
                newDose = min(results_crm, result_safety)
            }
        }
        
        parameters <- c(beta1, nu)
        names(parameters) <- c("beta0", "beta1", "nu")
        list(newDose = newDose, pstim=p_new, p_sum = p_sum, parameters = parameters)
    }
