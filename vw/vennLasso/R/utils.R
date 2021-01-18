
#' @importFrom stats sd
checkWhichVarsMissingForWhichConditions <- function(x, group.indicators) {
    # checks if any variables are missing for an entire disease condition
    # returns logical matrix of dimension (num.vars)x(num.groups) with TRUE
    # in position i,j indicating that variable i is missing in all of group j
    t(apply(is.na(x), 2,
            function(x.bool) apply(group.indicators, 2, function(gr) all(x.bool[which(gr == 1)]) ) ) )
}

checkWhichVarsMissingForAllButOneCondition <- function(x, group.indicators) {
    # checks if any variables are missing for all entire disease condition
    # returns logical matrix of dimension (num.vars)x(num.groups) with TRUE
    # in position i,j indicating that variable i is missing in all of group j
    t(apply(is.na(x), 2,
            function(x.bool) apply(group.indicators, 2, function(gr) all(x.bool[which(gr != 1)]) ) ) )
}

checkWhichVarsZeroSDForWhichConditions <- function(x, group.indicators, combin.mat) {
    # checks if any variables are missing for an entire disease condition
    # returns logical matrix of dimension (num.vars)x(num.groups) with TRUE
    # in position i,j indicating that variable i is missing in all of group j
    gr.symbols <- apply(group.indicators, 1, function(xx) paste(xx, collapse = ","))
    combins.unique <- apply(combin.mat, 1, function(xx) paste(xx, collapse = ","))
    ret <- t(apply(x, 2,
                   function(xj) sapply(combins.unique, function(gr) sd(xj[which(gr == gr.symbols)], na.rm = TRUE) == 0 ) ) )
    ret[is.na(ret)] <- TRUE
    ret
}


# taken from glmnet
cvcompute <- function(mat,weights,foldid,nlams)
{
    ###Computes the weighted mean and SD within folds, and hence the se of the mean
    wisum  <- tapply(weights,foldid,sum)
    nfolds <- max(foldid)
    outmat <- matrix(NA,nfolds,ncol(mat))
    good   <- matrix(0,nfolds,ncol(mat))
    mat[is.infinite(mat)] <- NA #just in case some infinities crept in
    for(i in seq(nfolds))
    {
        mati <- mat[foldid==i,,drop=FALSE]
        wi   <- weights[foldid==i]
        outmat[i,] <- apply(mati, 2, weighted.mean, w = wi, na.rm = TRUE)
        good[i,seq(nlams[i])] <- 1
    }
    N <- apply(good,2,sum)
    list(cvraw = outmat, weights = wisum, N = N)
}

# taken from glmnet
getmin <- function(lambda,cvm,cvsd)
{
    cvmin <- min(cvm,na.rm=TRUE)
    idmin <- cvm <= cvmin
    lambda.min <- max(lambda[idmin], na.rm = TRUE)
    idmin <- match(lambda.min,lambda)
    semin <- (cvm + cvsd)[idmin]
    idmin <- cvm <= semin
    lambda.1se <- max(lambda[idmin], na.rm = TRUE)
    list(lambda.min = lambda.min,
         lambda.1se = lambda.1se)
}

# taken from glmnet
lambdaInterp <- function(lambda,s)
{
    ### lambda is the index sequence that is produced by the model
    ### s is the new vector at which evaluations are required.
    ### the value is a vector of left and right indices, and a vector of fractions.
    ### the new values are interpolated bewteen the two using the fraction
    ### Note: lambda decreases. you take:
    ### sfrac*left+(1-sfrac*right)

    if(length(lambda)==1){# degenerate case of only one lambda
        nums=length(s)
        left=rep(1,nums)
        right=left
        sfrac=rep(1,nums)
    }
    else{
        s[s > max(lambda)] = max(lambda)
        s[s < min(lambda)] = min(lambda)
        k=length(lambda)
        if (lambda[1] > lambda[2]) {
            sfrac <- (lambda[1]-s)/(lambda[1] - lambda[k])
            lambda <- (lambda[1] - lambda)/(lambda[1] - lambda[k])
        } else {
            sfrac <- (lambda[k]-s)/(lambda[k] - lambda[1])
            lambda <- (lambda[k] - lambda)/(lambda[k] - lambda[1])
        }
        coord <- approx(lambda, seq(lambda), sfrac)$y
        left <- floor(coord)
        right <- ceiling(coord)
        sfrac=(sfrac-lambda[right])/(lambda[left] - lambda[right])
        sfrac[left==right]=1
    }
    list(left  = left,
         right = right,
         frac  = sfrac)
}

nonzeroCoeff <- function (beta, bystep = FALSE)
{
    ### bystep = FALSE means which variables were ever nonzero
    ### bystep = TRUE means which variables are nonzero for each step
    nr=nrow(beta)
    if (nr == 1) {#degenerate case
        if (bystep)
            apply(beta, 2, function(x) if (abs(x) > 0)
                1
                else NULL)
        else {
            if (any(abs(beta) > 0))
                1
            else NULL
        }
    } else 
    {
        beta  <- abs(beta)>0 # this is sparse
        which <- seq(nr)
        ones  <- rep(1, ncol(beta))
        nz    <- as.vector((beta %*% ones)>0)
        which <- which[nz]
        if (bystep) 
        {
            if(length(which)>0)
            {
                beta <- as.matrix(beta[which,,drop=FALSE])
                nzel <- function(x, which) if (any(x))
                    which[x]
                else NULL
                which <- apply(beta, 2, nzel, which)
                if(!is.list(which)) which <- data.frame(which)# apply can return a matrix!!
                which
            }
            else{
                dn    <- dimnames(beta)[[2]]
                which <- vector("list",length(dn))
                names(which) <- dn
                which
            }

        }
        else which
    }
}


# modified from glmnet
auc <- function(y,prob,w)
{
    if(missing(w))
    {
        rprob <- rank(prob)
        n1    <- sum(y)
        n0    <- length(y) - n1
        u     <- sum(rprob[y==1]) - n1 * (n1 + 1) / 2
        
        exp(log(u) - log(n1) - log(n0))
    }
    else
    {
        rprob <- runif(length(prob))
        op    <- order(prob,rprob)#randomize ties
        y     <- y[op]
        w     <- w[op]
        cw    <- cumsum(w)
        w1    <- w[y==1]
        cw1   <- cumsum(w1)
        wauc  <- log(sum(w1*(cw[y==1] - cw1)))
        sumw1 <- cw1[length(cw1)]
        sumw2 <- cw[length(cw)] - sumw1
        
        exp(wauc - log(sumw1) - log(sumw2))
    }
}
# taken from glmnet
auc.mat <- function(y,prob,weights=rep(1,nrow(y)))
{
    Weights <- as.vector(weights*y)
    ny      <- nrow(y)
    Y       <- rep(c(0,1),c(ny,ny))
    Prob    <- c(prob,prob)
    auc(Y, Prob, Weights)
}


# taken from glmnet
error.bars <- function(x, upper, lower, width = 0.02, ...)
{
    xlim <- range(x)
    barw <- diff(xlim) * width
    segments(x, upper, x, lower, col = 8, lty = 5, lwd = 0.5, ...)
    segments(x - barw, upper, x + barw, upper, col = "grey50", lwd = 1, ...)
    segments(x - barw, lower, x + barw, lower, col = "grey50", lwd = 1, ...)
    range(upper, lower)
}

