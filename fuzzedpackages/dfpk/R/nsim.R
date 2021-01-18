#' @import rstan
#' @useDynLib dfpk, .registration = TRUE
#' @export
nsim <- 
function(doses, N, cohort, icon, theta, model, simulatedData, TR=length(simulatedData), prob = 0.9, AUCmethod = 2, options = list(nchains = 4, niter = 4000, 
         nadapt = 0.8), betapriors = NULL, thetaL=NULL, p0 = 0, L = 0, CI = FALSE, seed = 190591){

    if(TR > length(simulatedData)) stop("Subscript out of bounds. TR must be less or equal to the length of simulatedData")
        
        model1 = NULL
        eval(parse(text = paste("model1 =", model, sep="")))
        MTD = NULL
        doseLevels = NULL
        toxicity = NULL
        AUC_s = NULL
        AUCd = NULL
        pstim1 = list()
        pstim3 = list()
        pstim_mean = list()
        seed_trial <- matrix(NA, nrow = 1, ncol = TR + 1)
        rownames(seed_trial) <- "Seed"
        colnames(seed_trial) <- c("Initial", paste("Trial ", 1:TR, sep = ""))
        seed_trial[1, 1] <- seed
        sel = rep(0,length(doses)+1)

        if (model == "pktox" & is.null(betapriors)){
        	betapriors = c(10, 10000, 20, 10)
        }else if(model == "pkcrm" & is.null(betapriors)){
        	betapriors = c(10, 10000)
        }else if (model == "pkpop" & is.null(betapriors)){
        	betapriors = c(10, 10000, 10, 5)
        }else if (model == "dtox" & is.null(betapriors)){
        	betapriors = c(0, 16.71, 0, 6.43)
        }else if(model == "pkcov" & is.null(betapriors)){
        	betapriors = c(-14.76, 0, 3.23 + 5)
        }else if (model == "pklogit" & is.null(betapriors)){
        	betapriors = c(10, 10000, 20, 10)
        }

        for (tr in 1:TR){
            set.seed(seed + tr)
            ndos <- length(doses)
            tox <- simulatedData[[tr]]@tox         
            stab <- simulatedData[[tr]]@tab 
            n_pk <- simulatedData[[tr]]@nPK        
            doses <- simulatedData[[tr]]@doses
            preal <- simulatedData[[tr]]@preal
            x <- rep(1,cohort)
            y <- tox[cbind(1:length(x),x)]  
            M = N/cohort
            nd = rep(0,length(doses))
            
            for (i in 1:length(x)){
                eval(parse(text = paste("conc",i," <- as.vector(stab[((i-1)*ndos +x[i] +1), 2:(n_pk +1)])", sep= "")))
                eval(parse(text = paste("conc",i," <- conc",i,"[icon]", sep = ""))) 
                nd[x[i]] <- nd[x[i]] + 1
            }
            
            time1 <- as.vector(stab[1, 2:(n_pk +1)])
            time1 <- time1[icon]
            
            AUCs <- NULL
            for (i in 1:length(x)){
                eval(parse(text = paste("AUCs <- c(AUCs, AUC.estim(conc=conc",i,", t=time1, dose=doses[x[",i,"]], method = AUCmethod))", sep="")))
            }
            
            pstim_auctox = matrix(0, length(doses)*cohort)
            pstim_Q1 = matrix(0, length(doses)*cohort)
            pstim_Q3 = matrix(0, length(doses)*cohort)

            AUCpop <- rep(0, length(doses))
            for(s in which(nd!=0)){
                AUCpop[s] = mean(AUCs[which(x==s)])
            }
            
            deltaAUC <- (log(AUCs) - log(AUCpop[x]))
            
            stage1 = TRUE
            for (i in 2:M) {
                j= (cohort*(i-1) + 1) : (cohort*i)   # position
                ### starting dose until toxicity
                if (stage1){
                    x <- c(x,rep(min((max(x)+1),length(doses)), cohort))             
                    y <- c(y, tox[cbind(j,x[j])])
                    for (k in j) {
                        conci <- as.vector(stab[((k-1)*ndos + x[k] +1), 2:(n_pk +1)])
                        conci <- conci[icon]
                        AUCs <- c(AUCs, AUC.estim(conc=conci, t=time1, dose=doses[x[k]], method = AUCmethod))
                        nd[x[k]] <- nd[x[k]] + 1 
                    }
                    pstim_auctox = cbind(pstim_auctox, rep(0,length(doses)))
                    if(CI == TRUE){
                        pstim_Q1 = cbind(pstim_Q1, rep(0,length(doses)))
                        pstim_Q3 = cbind(pstim_Q3, rep(0,length(doses)))
                    }else{
                        pstim_Q1 = NULL
                        pstim_Q3 = NULL
                    }
                    
                    for(s in which(nd!=0)){
                        AUCpop[s] = mean(AUCs[which(x==s)])
                    }
                    deltaAUC <- (log(AUCs) - log(AUCpop[x]))
                    
                    if (any(y == "1")) {stage1 <- FALSE}
                } else {
                    
                    results <- model1(y=y, auc = AUCs, doses = doses, x=x, theta=theta, prob = prob, betapriors = betapriors, 
                                     thetaL=thetaL, options = options, p0 = p0, L = L, deltaAUC = deltaAUC, CI = CI)
                    
                    if (is.na(results$newDose) == "TRUE") break

                    newdose <- min(results$newDose, max(x) + 1)     
                    # Check on the skipping dose
                    x <- c(x,rep(newdose,cohort))
                    y <- c(y, tox[cbind(j,x[j])])    
                    for (k in j) {
                        conci <- as.vector(stab[((k-1)*ndos +x[k] +1), 2:(n_pk +1)])
                        conci <- conci[icon]
                        AUCs <- c(AUCs, AUC.estim(conc=conci, t=time1, dose=doses[x[k]], method = AUCmethod))
                        nd[x[k]] <- nd[x[k]] + 1
                    }
                    pstim_auctox = cbind(pstim_auctox, results$pstim)
                    pstim_Q1 = cbind(pstim_Q1, results$p_sum[,2])
                    pstim_Q3 = cbind(pstim_Q3, results$p_sum[,5])
                    for(s in which(nd!=0)){
                        AUCpop[s] = mean(AUCs[which(x==s)]) 
                    }
                    deltaAUC <- (log(AUCs) - log(AUCpop[x]))
                }
            }
            trial <- paste('trial:',tr,sep='')
            # check if we stopped before
            if (length(x) < N){
                nstop <- N-length(x)
                MtD = 0
                mtd = results$newDose 
                if (is.na(mtd) == "TRUE"){
                    mtd = 0
                }
                MTD = c(MTD, mtd)
                sel[MtD+1] <- sel[MtD+1] + 1
                doseLevels = rbind(doseLevels,c(x,rep(0,nstop)))
                toxicity = rbind(toxicity,c(y,rep(NA,nstop)))
                AUC_s = rbind(AUC_s, c(AUCs,rep(NA,nstop)))
                AUCd = rbind(AUCd, c(deltaAUC,rep(NA,nstop)))
                pstim1[[trial]] = pstim_Q1
                pstim3[[trial]] = pstim_Q3
                pstim_mean[[trial]] = pstim_auctox
            }else{
                MtD = model1(y = y, auc = AUCs, doses = doses, x = x, theta = theta, prob = prob, betapriors = betapriors, 
                            thetaL = thetaL, options = options, p0 = p0, L = L, deltaAUC = deltaAUC, CI = CI)$newDose
                sel[MtD+1] <- sel[MtD+1] + 1
                MTD = c(MTD, MtD)
                doseLevels = rbind(doseLevels, x)
                toxicity = rbind(toxicity, y)
                AUC_s = rbind(AUC_s, AUCs)
                AUCd = rbind(AUCd, deltaAUC)
                pstim1[[trial]] <- pstim_Q1
                pstim3[[trial]] <- pstim_Q3
                pstim_mean[[trial]] <- pstim_auctox
            }
                nchains = options$nchains
                niter = options$niter
                nadapt = options$nadapt
                pid = c(1:N)
                seed_trial[1, tr+1] = seed+tr
        }

        if(TR == 1){
            newDose = MtD
        }else{
            newDose = sel/TR
        }
        
        new("dosefinding", pid = pid, N = N, time = time1, doses = doses, conc = conci, p0 = p0,
             L = L,  nchains = options$nchains, niter = options$niter, nadapt = options$nadapt, newDose = newDose, MTD = MTD, MtD=MtD,
             theta = theta, doseLevels = doseLevels, toxicity = toxicity, AUCs = AUC_s, TR = TR, preal = preal, 
             pstim  = pstim_mean, pstimQ1 = pstim1, pstimQ3 = pstim3, model = model, seed = seed_trial)

}
