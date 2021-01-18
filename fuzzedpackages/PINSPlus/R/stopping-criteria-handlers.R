DiffDevianceStoppingHandler <- function(kMin, kMax, origS, iterMin = 20, iterMax = 200, madMin = 1e-03, msdMin = 1e-06, onExcute){
    listAUC <- list()
    clistAUC <- list()
    isStop <- rep(F, kMax)
    isStop[1] <- T
    
    for (k in kMin:kMax) {
        listAUC[k] <- list(c())
        clistAUC[k] <- list(c())
    }

    aucMeanSquareDeviation <- function(arr) sum((arr - mean(arr))^2)/length(arr)
    aucMeanAbsoluteDeviation <- function(arr) sum(abs(arr - mean(arr)))/length(arr)

    function(iter = NULL, k, auc = NULL, type = 0,...){
        if (isStop[k]) return(1)
        stop = 0
        
        if (type == 0){
            
            listAUC[[k]] <<- c(listAUC[[k]], auc)
            m <- as.matrix(table(listAUC[[k]]))
            meanAUC <- sum(m[,1]^2*as.numeric(row.names(m)))/sum(m[,1]^2)
            clistAUC[[k]] <<- c(clistAUC[[k]], meanAUC)
            
            onExcute(k = k, AUCs = listAUC[[k]])
            
            lengthAUC = length(listAUC[[k]])
            
            start = iter - 20
            if (iter >= iterMin){
                if (aucMeanAbsoluteDeviation(clistAUC[[k]][start:iter]) <= madMin & aucMeanSquareDeviation(clistAUC[[k]][start:iter]) <= msdMin){
                    if (clistAUC[[k]][iter] == 1){
                        stop <- 2 # if stop = 2, stop all clustering 
                    }
                    else {
                        stop <- 1
                        isStop[k] <<- T
                        
                        currentAUCs = lapply(clistAUC, function(l) l[length(l)])
                        ma <- max(unlist(currentAUCs))
                        mk <- min(which(unlist(currentAUCs) == ma))
                        
                        if (k == mk){
                            #if (sum(isStop[-k]) > 0){
                                stop <- 0
                                isStop[k] <<- F
                            #}
                        }
                    }
                }
            }
        }
        else {
            currentAUCs = lapply(clistAUC, function(l) l[length(l)])
            ma <- max(unlist(currentAUCs))
            mk <- min(which(unlist(currentAUCs) == ma))
            
            if (currentAUCs[[k]] == ma){
                if (sum(isStop[-k]) == kMax - 1) stop <- 1
            }
            else if (k > mk){
                listAUC.tmp <- c(listAUC[[k]], rep(1,10))
                m.tmp <- as.matrix(table(listAUC.tmp))
                meanAUC.tmp <- sum(m.tmp[,1]^2*as.numeric(row.names(m.tmp)))/sum(m.tmp[,1]^2)
                if (meanAUC.tmp < ma){
                    stop <- 1
                }
            }
        }
        
        stop
    }
}