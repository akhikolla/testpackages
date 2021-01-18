#' An internal function to compute sigma square analytically
#' @param t1 t1
#' @param t2 t2
#' @param b b (see Yashin et. al, 2007)
#' @return sigma_square (see Akushevich et. al, 2005)
sigma_sq <- function(t1, t2, b) {
    # t2 = t_{j}, t1 = t_{j-1}
    ans <- as.numeric(b) %*% as.numeric(t2-t1)
    ans
}

#' An internal function to compute m from 
#' @param y Current value of Y
#' @param t1 t1
#' @param t2 t2
#' @param a a (see Yashin et. al, 2007)
#' @param f1 f1 (see Yashin et. al, 2007)
#' @return m m (see Yashin et. al, 2007)
m <- function(y, t1, t2, a, f1) {
    # y = y_{j-1}, t1 = t_{j-1}, t2 = t_{j}
    ans <- y + as.numeric(a) %*% t(as.numeric(y) - as.numeric(f1)) %*% as.numeric((t2 - t1))
    ans
}

#' An internal function to compute mu
#' @param y Current value of y
#' @param mu0 mu0 (see Yashin et. al, 2007)
#' @param b b (see Yashin et. al, 2007)
#' @param Q Q (see Yashin et. al, 2007)
#' @param theta theta (see Yashin et. al, 2007)
#' @param tt t (time)
#' @return mu Next value of mu
mu <- function(y, mu0, b, Q, theta, tt) {
    ans <- (mu0 + t(y) %*% b + t(y) %*% Q %*% y)*exp(theta*tt)
    ans
}

#mu <- function(tt, y1, gamma1, f, f1, mu0, theta, Q) {
#  hf <- f - y1
#  hf1 <- f1 - y1
#  #if(gomp) {
#  #  mu0Ht = mu0H*exp(thetaH*t);
#  #} else {
#  #  mu0Ht = mu0H;
#  #}
#  #print(t)
#  #print(theta)
#  mu0Ht <- mu0*exp(theta*tt)
#  QH_gamma1 <- Q %*% gamma1
#  mu <- mu0Ht + (t(hf) %*% Q) %*% hf + sum(diag(QH_gamma1))
#  mu
#}

#' An internal function to compute next value of physiological variable Y
#' @param y1 y1
#' @param t1 t1
#' @param t2 t2
#' @param b b (see Yashin et. al, 2007)
#' @param a a (see Yashin et. al, 2007)
#' @param f1 f1 (see Yashin et. al, 2007)
#' @return y.next Next value of y
getNextY.cont2 <- function(y1, t1, t2, b, a, f1) {
    ssd <- sigma_sq(t1, t2, b)
    mean <- m(as.numeric(y1), t1, t2, a, f1)
    dt <- as.numeric(t2-t1)
    #y2 <- sapply(1:length(y1), function(i) {rnorm(1,mean=mean[i], sd=sqrt(b[i])/dt)})
    y2 <- sapply(1:length(y1), function(i) {rnorm(1,mean=mean[i], sd=sqrt(ssd[i]))})
    #y2 <- sapply(1:length(y1), function(i) {rnorm(1,mean=mean[i], sd=b[i])})
    y2
}

#' An internal function to compute the next value of physiological variable Y
#' based on discrete-time model (Akushevich et. al., 2005)
#' @param y1 y1
#' @param u u (see Akushevich et. al, 2005)
#' @param R R (see Akushevich et. al, 2005)
#' @param Sigma Sigma (see Akushevich et. al, 2005)
#' @return y.next Next value of y
getNextY.discr <- function(y1, u, R, Sigma) {
    eps<-matrix(nrow=dim(R)[1], ncol=1)
    #eps[,1] <- sapply(1:length(Sigma), function(i) {rnorm(1, mean=0.0, sd=sqrt(Sigma[i]))})
    eps[,1] <- sapply(1:length(Sigma), function(i) {rnorm(1, mean=0.0, sd=Sigma[i])})
    y2 <- u + R %*% y1 + eps
    return(y2)
}

getNextY.discr2 <- function(y1, u, R, Sigma) {
    y2 <- sapply(1:length(Sigma), function(i) {rnorm(1, mean=getNextY.discr.m(y1, u, R)[i], sd=sqrt(Sigma[i])/2)})
    y2
}


#' An internal function to compute next m based on dicrete-time model 
#' @param y1 y1
#' @param u u
#' @param R R
#' @return m Next value of m (see Yashin et. al, 2007)
getNextY.discr.m <- function(y1, u, R) {
    m <- u + R %*% y1
    return(m)
}

#' An internal function to compute previous m based on discrete-time model
#' @param y2 y2
#' @param u u
#' @param R R
#' @return m Next value of m (see Yashin et. al, 2007)
getPrevY.discr.m <- function(y2, u, R) {
  m <- solve(R) %*% (y2 - u)
  m
}

#' An internal function to compute previous value of
#' physiological variable Y based on 
#' discrete-time model
#' @param y2 y2
#' @param u u
#' @param R R
#' @param Sigma Sigma
#' @return y1 Previous value of y
getPrevY.discr <- function(y2, u, R, Sigma) {
    eps<-matrix(nrow=dim(R)[1], ncol=1)
    eps[,1] <- sapply(1:dim(eps)[1], function(i) {rnorm(1, mean=0.0, sd=Sigma[i])})
    
    y1 <- solve(R) %*% (y2 - u - eps)
    y1
}


#' An internal function to compute m and gamma based on 
#' continuous-time model (Yashin et. al., 2007)
#' @param tt tt - time
#' @param y y 
#' @param a a (see Yashin et. al, 2007)
#' @param f1 f1 (see Yashin et. al, 2007)
#' @param Q Q (see Yashin et. al, 2007)
#' @param f f (see Yashin et. al, 2007)
#' @param b b (see Yashin et. al, 2007)
#' @param theta theta
#' @return list(m, gamma) Next values of m and gamma (see Yashin et. al, 2007)
func1 <- function(tt, y, a, f1, Q, f, b, theta) {
    hf <- f - as.numeric(y[[1]])
    hf1 <- f1 - as.numeric(y[[1]])
    
    dm <- -1.00 * (a %*% as.matrix(hf1)) + (2.00 * as.numeric(y[[2]])) %*% Q %*% hf
    dgamma <- a %*% as.numeric(y[[2]]) + as.numeric(y[[2]]) %*% t(a) + b %*% t(b) - 2.00 * ((as.numeric(y[[2]]) %*% Q) %*% as.numeric(y[[2]]))
    
    res <- list(m=dm, gamma=dgamma)
    return(res)
}

#' An internal function to compute next Y based on 
#' continous-time model (Yashin et. al., 2007)
#' @param y1 y1
#' @param t1 t1
#' @param t2 t2
#' @param a a (see Yashin et. al, 2007)
#' @param f1 f1 (see Yashin et. al, 2007)
#' @param Q Q (see Yashin et. al, 2007)
#' @param f f (see Yashin et. al, 2007)
#' @param b b (see Yashin et. al, 2007)
#' @param theta theta (see Yashin et. al, 2007)
#' @return y.next Next value of Y
getNextY.cont <- function(y1, t1, t2, a, f1, Q, f, b, theta) {
    nsteps <- 4
    tdiff <- t2-t1
    h <- tdiff/nsteps
    gamma1 <- matrix(nrow=dim(a)[1], ncol=dim(a)[1], 0) # set gamma1 to zero matrix
    #gamma1 <- matrix(nrow=1, ncol=1, 0) # set gamma1 to zero matrix
    tt <- t1
  
    out <- list()
    out[[1]] <- y1
    out[[2]] <- gamma1
  
  
    for(j in 1:nsteps) {
        #Runge-Kutta method:
        yn <- out
        k1 <- func1(tt, yn, a, f1, Q, f, b, theta)
    
        yn[[1]] <- out[[1]] + h/2 * k1[[1]]
        yn[[2]] <- out[[2]] + h/2 * k1[[2]]
        k2 <- func1(tt+h/2, yn, a, f1, Q, f, b, theta)
    
        yn[[1]] <- out[[1]] + h/2 * k2[[1]]
        yn[[2]] <- out[[2]] + h/2 * k2[[2]]
        k3 <- func1(tt+h/2, yn, a, f1, Q, f, b, theta)
    
        yn[[1]] <- out[[1]] + h * k3[[1]]
        yn[[2]] <- out[[2]] + h * k3[[2]]
        k4 <- func1(tt+h, yn, a, f1, Q, f, b, theta)
    
        out[[1]] <- out[[1]] + h/6 * (k1[[1]] + 2*k2[[1]] + 2*k3[[1]] + k4[[1]])
        out[[2]] <- out[[2]] + h/6 * (k1[[2]] + 2*k2[[2]] + 2*k3[[2]] + k4[[2]])
    
        tt <- tt + h
    }
  
    m2 <- out[[1]]
    gamma2 <- out[[2]]
  
    # New y2:
    y2 <- matrix(nrow=dim(m2)[1], ncol=1, 0)
    for(ii in 1:dim(m2)[1]) {
        y2[ii,1] <- rnorm(1, mean=m2[ii,1], sd=sqrt(gamma2[ii,ii])) 
        #y2[ii,1] <- rnorm(1, mean=m2[ii,1], sd=gamma2[ii,ii]) 
    }
    y2
}

#'Multiple Data Imputation with SPM
#'@param x A longitudinal dataset with missing observations
#'@param id A name (text) or index (numeric) of ID column. Default: 1
#'@param case A case status column name (text) or index (numeric). Default: 2
#'@param t1 A t1 (or t if short format is used) column name (text) or index (numeric). Default: 3
#'@param t2 A t2 column name (if long format is used) (text) or index (numeric). Default: 4
#'@param covariates A list of covariate column names or indices. Default: 5
#'@param minp Number of imputations. Default: 5
#'@param theta_range A range of parameter theta used for optimization, default: seq(0.01, 0.15, by=0.001).
#'@return A list(imputed, imputations)
#'@return imputed An imputed dataset.
#'@return imputations Temporary imputed datasets used in multiple imputaitons.
#'@export
#'@examples \dontrun{
#'library(stpm) 
#'##Data preparation ##
#'data <- simdata_discr(N=1000, dt = 2)
#'miss.id <- sample(x=dim(data)[1], size=round(dim(data)[1]/4)) # ~25% missing data
#'incomplete.data <- data
#'incomplete.data[miss.id,5] <- NA
#'incomplete.data[miss.id-1,6] <- NA
#'## End of data preparation ##
#'
#'# Estimate parameters from the complete dataset #
#'p <- spm_discrete(data, theta_range = seq(0.075, 0.09, by=0.001))
#'p
#'
#'##### Multiple imputation with SPM #####
#'imp.data <- spm.impute(x=incomplete.data, 
#'                       minp=5, 
#'                       theta_range=seq(0.075, 0.09, by=0.001))$imputed
#'head(imp.data)
#'## Estimate SPM parameters from imputed data and compare them to the p ##
#'pp.test <- spm_discrete(imp.data, theta_range = seq(0.075, 0.09, by=0.001))
#'pp.test
#'}
spm.impute <- function(x, 
                       id=1, 
                       case=2,
                       t1=3, 
                       t2=3, 
                       covariates=4,
                       minp=5, 
                       theta_range=seq(0.01, 0.2, by=0.001)) 
{
    
    # Check input parameters for correctness
    if(class(x) != "data.frame") {
        stop("Class of dataset must be a 'data.frame'.")
    }
  
    datasets <- list() # To hold imputed datasets
    
    col.id.ind <- ifelse(class(id)=="character", get.column.index(x, id), id)
    if(col.id.ind == 0) stop(paste("Column",id, "not found in data table!"))
    
    col.status.ind <- ifelse(class(case)=="character", get.column.index(x, case), case)
    if(col.status.ind == 0) stop(paste("Column",case, "not found in data table!"))
    
    col.age.ind <- ifelse(class(t1)=="character", get.column.index(x, t1), t1)
    if(col.age.ind == 0) stop(paste("Column",t1, "not found in data table!"))
    
    col.age.event.ind <- ifelse(class(t2)=="character", get.column.index(x, t2), t2)
    if(col.age.event.ind == 0) stop(paste("Column", t2, "not found in data table!"))
    
    col.covar.ind <- c()
    for(c in covariates) {
        c.ind <- ifelse(class(c)=="character", get.column.index(x, c), c)
        if(c.ind == 0) stop(paste("Column", c, "not found in data table!"))
        col.covar.ind <- c(col.covar.ind, c.ind)
    } 
    
    if(length(which(is.na(x[, col.id.ind]))) != 0) {
        stop("Dataset contains NAs in ID column.")
    }
    # Remove rows with id = NA/NULL
    #x <- x[which(!is.na(x[, col.id.ind])), ]
    
    # Constructing working dataset
    dataset.work <- data.frame(x[, c(col.id.ind, col.status.ind, col.age.ind, col.covar.ind)])
    
    Ncol <- dim(dataset.work)[2]
    
    # Estimate parameters from raw data
    dataset <- prepare_data_cont(x, 
                                col.id.ind=col.id.ind, 
                                col.status.ind=col.status.ind,
                                col.age.ind=col.age.ind, 
                                col.age.event.ind=col.age.event.ind, 
                                col.covar.ind=covariates, 
                                dt=1, 
                                verbose=FALSE)
    
    pp <- spm_discrete(dataset, theta_range = theta_range)
    
    # Individual IDs
    ids <- unique(dataset[,1])
    #print(pp)
    for(m in 1:minp) {
      
        ####
        # for u:
        pp.dmodel.u <- pp$dmodel$u
        for(i in 1:length(pp$dmodel$u)) {
            mu <- rnorm(10000, mean=pp$dmodel$u[i], sd=sqrt(pp$dmodel$u.std.err[i]))
            z <- rnorm(10000, mean=mu, sd=pp$dmodel$u.std.err[i])
            pp.dmodel.u[i] <- mean(z)
        }
        # For R:
        pp.dmodel.R <- pp$dmodel$R
        for(i in dim(pp$dmodel$R)[1]) {
            for(j in dim(pp$dmodel$R)[2]) {
                mu <- rnorm(10000, mean=pp$dmodel$R[i,j], sd=sqrt(pp$dmodel$R.std.err[i,j]))
                z <- rnorm(10000, mean=mu, sd=pp$dmodel$R.std.err[i,j])
                pp.dmodel.R[i,j] <- mean(z)
            }
        }
        ####
        ## Temporary working dataset:
        dat <- dataset.work
        
        for(k in ids) {
            ########## Forward #########
            # Retrieve records for single subject
            df <- dat[which(dat[,1] == k), ]
            
            Nrec <- dim(df)[1]
      
            if(Nrec == 1) {
                row.cur <- df
                for(j in 4:Ncol) {
                    l <- j-3
                    if(is.na(row.cur[j])) {
                        y.start <- rnorm(1, mean = pp$cmodel$f[l], sd=pp$dmodel$Sigma[l])
                        row.cur[j] <- y.start
                    }
                }
                
                df <- row.cur
                dat[which(dat[,1] == k), ] <- df
                
                next
            }
      
            # Check the first row #
            row.cur <- df[1, ]
            row.next <- df[2, ]
  
            for(j in 4:Ncol) {
                l <- j-3
                if(is.na(row.cur[j]) & !is.na(row.next[j])) {
                    y2 <- row.next[j]
                    # Only one-dimension
                    y1 <- getPrevY.discr(as.matrix(y2), pp.dmodel.u[l], as.matrix(pp.dmodel.R[l,l]), pp$dmodel$Sigma[l])
                    row.cur[j] <- y1
                } else if(is.na(row.cur[j]) & is.na(row.next[j])) {
                    y.start <- rnorm(1, mean = pp$cmodel$f1[l], sd=pp$cmodel$b[l])
                    row.cur[j] <- y.start
                }
            }
            
            df[1, ] <- row.cur
            
      
            #### Main imputation loop ####
        
            for(i in 2:Nrec) {
                row.prev <- df[i-1,]
                row.cur <- df[i,]
                
                if(i != Nrec & Nrec >2) {row.next <- df[i+1,]}
                
                y1 <- row.prev[4:Ncol]
                
                if(any(is.na(row.cur))) {
                    y.next <- getNextY.discr.m(as.matrix(as.numeric(y1)), as.matrix(pp.dmodel.u), pp.dmodel.R)
                    row.cur.na <- which(is.na(row.cur))
                    row.cur[row.cur.na] <- y.next[row.cur.na - 3]
                }
            
                df[i, ] <- row.cur
            }
      
            ### Last record in a dataset ###
            row.cur <- df[Nrec, ]
            row.prev <- df[Nrec-1,]
            if(any(is.na(row.cur[4:Ncol])) & row.cur[2] == 0) {
                y1 <- row.prev
                y.next <- getNextY.discr.m(as.matrix(as.numeric(y1)), as.matrix(pp.dmodel.u), pp.dmodel.R)
                row.cur.na <- which(is.na(row.cur))
                row.cur[row.cur.na] <- y.next[row.cur.na - 3]
            }
      
            df[Nrec, ] <- row.cur
            
            
            tryCatch({
                dat[which(dat[,1] == k), ] <- df
            },error=function(e) {
                print(e)
                print("dat:")
                print(dat[which(dat[,1] == k), ])
                print("df:")
                print(df)
                stop()
            }, warning=function(w){
                print(w)
                print("dat:")
                print(dat[which(dat[,1] == k), ])
                print("df:")
                print(df)
                dat[which(dat[,1] == k), ] <- df
            })
        }
        #################################################################################
        datasets[[m]] <- dat
    }
  
    ### Summarizing imputed datasets together by averaging missing values, i.e. 'completion' ###
    final.dataset <- x
    for(j in 4:Ncol) {
        data.tmp <- matrix(nrow = dim(dataset.work)[1], ncol=0)
        for(i in 1:minp) {
            data.tmp <- cbind(data.tmp, datasets[[i]][,j])
        }
        data.tmp.2.mean <- apply(X = data.tmp, FUN = mean, MARGIN = 1, na.rm=T)
        dataset.work[,j] <- data.tmp.2.mean
    }
    
    colnames.x <- colnames(x)
    final.dataset[, col.covar.ind] <- dataset.work[,4:Ncol]
    colnames(final.dataset) <- colnames(x)
    
    res <- list(imputed=final.dataset, imputations=datasets, data.long.format=dataset)
    res
}