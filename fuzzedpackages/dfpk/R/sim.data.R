#' @useDynLib dfpk, .registration = TRUE
#' @export
sim.data <-
function(PKparameters,omegaIIV,omegaAlpha,sigma,doses,limitTox,timeSampling, N, TR, seed=190591){
    
    pk.model <- function(t,dose,ka,CL,V){
        dose*ka/V/(ka-CL/V)* (exp(-CL/V*t)-exp(-ka*t))  
    }
    
    resScen <- list()

    for (tr in 1:TR){

        set.seed(seed + tr)
        nPK <- length(timeSampling)      		
        # doses <- exp(qnorm(preal)*sqrt(omegaIIV^2+omega_a^2) + log(limit_tox) + log(PKparameters[2]))  
        preal <- pnorm((log(doses)-log(limitTox)-log(PKparameters[2]))/sqrt(omegaIIV^2+omegaAlpha^2)) 
        parameter <- NULL               			
        sens_AUC <- NULL                			
        alphatot <- NULL			 			
        tox <- NULL

        tab <- c(0,timeSampling)	 
        npar <- length(PKparameters) - 1             

        for(i in 1:N){
            ipar <- (PKparameters[2:3])*exp(rnorm(npar, sd=omegaIIV))
            ipar <- c(PKparameters[1],ipar)
            parameter <- rbind(parameter,ipar)
            alfa <- exp(rnorm(1, sd=omegaAlpha))
            alphatot <- c(alphatot, alfa)
            for(j in doses){
                concen <- pk.model(timeSampling, dose=j, ka = ipar[1], CL = ipar[2], V = ipar[3]) 
                concPred <- concen*(1+rnorm(nPK, sd=sigma))   	# real concentrations + predictions of them
                tab <- rbind(tab,c(i, concPred))     					
                CL <- ipar[2]		   
                sens <- alfa*j/CL          
                sens_AUC <- c(sens_AUC, sens)  
            }
        }
        
        for(i in 1:N){
            row.names(parameter)[i] <- paste("pid", i)
            colnames(parameter) <- c("ka", "CL", "V")
        }
       
        #####################################
        ########## Concentration  ###########
        #####################################
        
        concentration <- NULL
        for(k in 1:length(doses)){
            conc <- pk.model(timeSampling, dose=doses[k], ka = PKparameters[1], CL = PKparameters[2], V = PKparameters[3])
            concentration <- cbind(concentration,conc)
        }

        for(i in 1:length(sens_AUC)){
            if(sens_AUC[i] >= limitTox){
                tox[i] = 1
            }else{
                tox[i] = 0
            }
        }
        tox <- matrix(tox, ncol = length(doses), byrow = T)
        rownames(tox) <- paste("pid", 1:N)
        colnames(tox) <- paste("dose", 1:length(doses))

        res <- new("scen", PKparameters=PKparameters, nPK= nPK, time=timeSampling, idtr=tr, N=N, doses=doses, preal = preal,
                    limitTox=limitTox, omegaIIV=omegaIIV, omegaAlpha=omegaAlpha, conc=concentration, 
                    concPred=concPred, tox=tox, tab=tab, parameters=parameter, alphaAUC=sens_AUC)

        resScen[[tr]] <- res
    }
    return(resScen)
}


