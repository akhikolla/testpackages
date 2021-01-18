### basad.R
### Main fucntion of basad
###
###
### Author: Qingyan Xiang

basad <- function(x = NULL,
                  y = NULL,
                  K = -1,
                  df = 5,
                  nburn = 1000,
                  niter = 1000,
                  alternative  = FALSE,
                  verbose = FALSE,
                  nsplit = 20,
                  tau0 = NULL,
                  tau1 = NULL,
                  prior.dist = "Gauss",
                  select.cri = "median",
                  BIC.maxsize = 20){

###--------------------------------
###Preproccessing
###--------------------------------

    #######Check the data
	if (is.null(y) | is.null(x) )
		stop("x and or y is missing")
    if( nrow(x) != length(y) )
        stop("number of observations does not math")
    if( is.factor(y) )
        stop("y should be continuous for basad regression")
    if( any( is.na(x)) )
        stop("NA not permitted in x")
    if( any( is.na(y)) )
        stop("NA not permitted in y")
    
    
    ########Prepare the data
	Y <- y
	X <- x
    X <- scale(X)     ## scale the X input matrix
    X <- cbind( rep(1, nrow(X)), X)   ##add the intercept
    beta.names <- colnames(X)
    if (length(unique(beta.names)) != ncol(X))
		colnames(X) <- beta.names <- paste("x.", 0:(ncol(X)-1), sep = "")

    
	colnames(X)[1] <- beta.names[1] <- "intercept"			  
	
    X <- data.matrix(X)
    Y <- data.matrix(Y)
    Y <- Y - mean(Y)  ## scale the Y vector
    
  
###----------------------------------
###Set parameters and Initial values
###----------------------------------
  
    ###Set the dimension and sample sizes
    p = dim(X)[2]-1;  n = dim(X)[1]

    ########calculate the prior probability of q
    if( K > 2 ){
        choicep <-function(x){
            return(x - K + qnorm(0.9)*sqrt(x*(1- x/p)))
        }
        cp = uniroot(choicep, c(1,K))$root
        pr = cp/p
    }
    else
        pr = -1
    
    pr0 = 0.1
    
    ##########calculate the ols estimate of sig
    Bols = sapply(seq(1:(p+1)), function(j) summary(lm(Y~ 0+X[,j]))$coef[,1])
    B0 = Bols
    if(p > n) {m = round(p - (n/2))}
    if(p <=n) m = round(p/2)
    ind = which(abs(B0)> sort(abs(B0))[m+1])
    sighat = sum((summary(lm(Y~0+X[,union(1,ind)]))$residuals)^2)/(n-length(union(1,ind)))

    B0 = rep(0,(p+1))
    Z0 = array(0,(p+1))
    sig = sighat
    
    ##########make a proper nsplit in low dimensional cases
    if( p > 30 )
        nsplit = nsplit
    else
        nsplit = 1
        
    if( nsplit > p )
        stop("Number of splits can not exceed dimensions")
        
	cat("Algorithms running:",  "\n" )
###---------------------------------
###RUN THE CORE
###---------------------------------
###The tau0 coefficient

    if( !is.null(tau0) )
        s0 = tau0
    else
        s0 = 1/n

    if( prior.dist == "t"){
        
        nu = df
        fvalue = dt(sqrt(2.1*log(p+1)), df = nu)
        
        if( !is.null(tau1) )
            s1 = tau1
        else
            s1  = max(100*s0, pr0*s0/((1 - pr0)*fvalue))
        
		res <- .Call( 'basadFuncScale', X, Y, Z0, B0, sig, pr, n, p, nu, s0, s1, nburn, niter, nsplit, alternative, priorType = 0, PACKAGE = 'basad' )
    }
    else if( prior.dist == "Laplace"){
    
        fvalue = dlaplace(sqrt(2.1*log(p+1)))
        if( !is.null(tau1)  )
            s1 = tau1
        else
            s1  = max(100*s0, pr0*s0/((1 - pr0)*fvalue));
        
	    res <- .Call( 'basadFuncScale', X, Y, Z0, B0, sig, pr, n, p, lambda = 1, s0, s1, nburn, niter, nsplit, alternative, priorType = 1, PACKAGE = 'basad' )
    }
    else if( prior.dist == "Gauss"){
    
        fvalue = dnorm(sqrt(2.1*log(p+1)))
        if( !is.null(tau1)   )
            s1 = tau1
        else
            s1  = max(100*s0, pr0*s0/((1 - pr0)*fvalue));
        
        res <- .Call( 'basadFunctionG', X, Y, Z0, B0, sig, pr, n, p, s0, s1, nburn, niter, nsplit, alternative,  PACKAGE = 'basad')
    }
    else{
       stop("No such prior type, please choose from 'Gauss', 't', 'Laplace'")
    }

###--------------------------------
###SUMMARY
###--------------------------------
    allZ <- res$Z
    allB <- res$B

    ZZ <- res$Z[(nburn+1):(nburn+niter), ]
    BB <- res$B[(nburn+1):(nburn+niter), ]
    
    ####posterior Z and B vector
    Z <- apply(ZZ, 2, mean)
    B <- apply(BB, 2, mean)
    
    Zsort <- sort.int(Z, decreasing = TRUE, index.return = TRUE )
    
    modelIdx <- c()
    ####return the results if selection criteria is BIC
    
    if( p > BIC.maxsize )
        BICsize = BIC.maxsize
    else
        BICsize = p
    
    if( select.cri == "BIC" ){
        print( BICsize )
        bic <- numeric(BICsize)
        for( i in 1:BICsize ){
            idx <- Zsort$ix[ 1:i ]
            Bbic <- numeric(p + 1)
            Bbic[idx] <- B[idx]
            bic[i] <- n * log( ( t( Y - X %*% Bbic ) %*% ( Y - X %*% Bbic ) ) / n  ) + log(n) * i
        }
        minIdx <- which.min(bic)
        modelIdx <-  Zsort$ix[1:minIdx]
    }
    else{
    ####median probability model vector
        modelIdx <- which( Z > 0.5 )
        select.cri = "Median Probability Model"
    }
    
    ##rerturn all variables
    basad.all <- round( B[ Zsort$ix ], 4)
    posterior.all <- round( Z[ Zsort$ix ], 4)
    basad.all <- data.frame(basad.all, posterior.all)
    
    rownames( basad.all ) <- beta.names[  Zsort$ix ]
    colnames( basad.all ) <- c( "estimated values", "posterior prob" )

    ##return slected variables
	basad.select <- round(B[ modelIdx[-1] ] , 4)
    posterior.select <- round( Z[ modelIdx[-1] ], 4)
	basad.select <- data.frame( basad.select, posterior.select )
	rownames( basad.select ) <- beta.names[ modelIdx[-1] ]
	colnames( basad.select ) <- c("estimated values", "posterior prob")
    
    
    ##return selected model
    modelZ <- numeric(p+1)
    modelZ[modelIdx] <- 1
    names( modelZ ) <- beta.names
    
    ##return the original vector
    names( B ) <- beta.names
	names( Z ) <- beta.names
    
###--------------------------------
###PRINT VERBOSE DETAIL
###--------------------------------

verboseList <- list(
    c(n),
    c(p+1),
    c(nburn),
    c(niter),
    c(nsplit),
	c(alternative),
    c(select.cri)
)

if( verbose ){
cat("-----------------------------------", "\n")
cat("Sample size                      :", verboseList[[1]], "\n" )
cat("Dimension                        :", verboseList[[2]], "\n" )
cat("Burn-in length                   :", verboseList[[3]], "\n" )
cat("Iteration length                 :", verboseList[[4]], "\n" )
cat("Block updating split sizes       :", verboseList[[5]], "\n" )
cat("Alternative fast sampling        :", verboseList[[6]], "\n" )
cat("Model selection criteria         :", verboseList[[7]], "\n" )
cat("\n\n")
cat("-----> Selected variables:\n")
print(round(basad.select, 4))
cat("-----------------------------------", "\n")
}

###----------------------------------
###RETURN RESULTS
###----------------------------------
    out <- list(
        all.var = basad.all,
        select.var = basad.select,
        beta.names = beta.names,
		verbose = verboseList,
		posteriorZ = Z,
        modelSize = length(modelIdx) - 1 ,
        model.index = modelIdx,
        modelZ = modelZ,
		est.B = B,
        allB = allB,
        allZ = allZ,
		x = X,
		y = Y
       )
    
    class(out) <- "basad"
    return(out)
}

