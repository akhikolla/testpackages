
fDMA <- function (y,x,alpha,lambda,initvar,W=NULL,initial.period=NULL,V.meth=NULL,kappa=NULL,
                  gprob=NULL,omega=NULL,model=NULL,parallel=NULL,m.prior=NULL,mods.incl=NULL,
                  DOW=NULL,DOW.nmods=NULL,DOW.type=NULL,DOW.limit.nmods=NULL,progress.info=NULL,
                  forced.models=NULL,forbidden.models=NULL,forced.variables=NULL,
                  bm=NULL,small.c=NULL,fcores=NULL,mods.check=NULL,red.size=NULL,
                  av=NULL)
{

### requires "forecast", "parallel", "stats" and "xts" packages

### y - a numeric or a column matrix of a dependent variable,
###     if y is a xts object, then plots will have time index on the x axis

### x - a matrix of independent variables (drivers), different columns correspond to different variables

### alpha - a forgetting factor between 0 and 1 used in probabilities estimations

### lambda - a forgetting factor between 0 and 1 used in variance approximations,
###          or numeric vector of forgetting factors between 0 and 1

### initvar - initial variance in the state space equation

### W - a method for setting the initial values of variance for the models equations,
###     W = "reg" corresponds to the method based on the linear regression as in the paper by Raftery et al. (2010),
###     alternatively an arbitrary positive number can be specified,
###     by default the method of Raftery et al. (2010) is used

### initial.period - a number of observation since which MSE and MAE are computed,
###                  by default the whole sample is used, i.e., initial.period = 1

### V.meth - a method for the state space equation variance updating,
###          V.meth = "rec" corresponds to the recursive moment estimator, as in the paper by Raftery et al. (2010),
###          V.meth = "ewma" corresponds to the exponentially weighted moving average,
###          by default V.meth = "rec" is used

### kappa - a parameter in the exponentially weighted moving average, between 0 and 1,
###         used if V.meth = "ewma"

### gprob - a matrix of Google probabilities, columns should correspond to columns of x

### omega - a parameter between 0 and 1 used in probabilities estimations,
###         used if gprob is specified

### model - model = "dma" for Dynamic Model Averaging, model = "dms" for Dynamic Model Selection,
###         or model = "med" for Median Probability Model,
###         by default model = "dma" is used

### parallel - indicate whether parallel computations should be used,
###            by default parallel = FALSE

### m.prior - a parameter for general model prior (Mitchell and Beauchamp, 1988),
###           by default m.prior = 0.5, which corresponds to the uniform distribution

### mods.incl - a matrix indicating which models should be used for estimation,
###             the first column indicates inclusion of a constant,
###             by default all possible models with a constant are used

### DOW - a threshold for Dynamic Occam's Window (Onorante and Raftery, 2016),
###       should be a number between 0 and 1,
###       if DOW = 0, then no Dynamic Occam's Window is applied,
###       by default DOW = 0

### DOW.nmods - initial number of models for Dynamic Occam's Window,
###             should be less than the number of all possible models and larger than or equal to 2,
###             they are randomly chosen,
###             if DOW.nmods = 0, then initially models with exactly one variable are taken,
###             by default DOW.nmods = 0

### DOW.type - DOW.type = "r" for DMA-R (Onorante and Raftery, 2016),
###            DOW.type = "e" for DMA-E,
###            by default DOW.type = "r"

### bm - indicate whether benchmark forecast should be computed,
###      by default bm = FALSE

### DOW.limit.nmods - the limit of number of models used in Dynamic Occam's Window

### progress.info - applicable only if Dynamic Occam's Window is used,
###                 if progress.info = TRUE number of round and number of models are printed,
###                 by default progress.info = FALSE

### small.c - a small constant added to posterior model probabilities,
###           to prevent possible reduction them to 0 due to computational issues,
###           by default small.c is taken as in small constant as in Eq. (17)
###           in Raftery et al. (2010)

### forced.models - matrix of models which if Dynamic Occam's Window is used
###                 have to be always included, similar as mods.incl,
###                 by default forced.models = NULL

### forbidden.models - matrix of models which should not be used
###                    in Dynamic Occam's Window, similar as mods.incl,
###                    by default forbidden.models = NULL

### forced.variables - vector of variables which must be present in each model
###                    in Dynamic Occam's Window, first slot indicates constant,
###                    by default forced.variables = NULL

### fcores - used only if parallel = TRUE, otherwise ignored,
###          indicates the number of cores that should not be used,
###          by default fcores = 1

### mods.check - indicates if models specified by the user should be checked,
###              by default mods.check = FALSE

### red.size - indicates if outcomes should be reduced to save memory,
###            by default red.size = FALSE 

### av - av = "dma" corresponds to the original DMA averaging scheme,
###      av = "mse" corresponds to averaging based on Mean Squared Error,
###      av = "hr1" corresponds to averaging based on Hit Ratio, assuming time-series are in levels,
###      av = "hr2" corresponds to averaging based on Hit Ratio, assuming time-series represent changes,
###      by default av = "dma"

###################################### checking initial parameters

if (missing(y)) { stop("please, specify y") }
if (missing(x)) { stop("please, specify x") }
if (is.xts(x)) { x <- as.matrix(x) }
if (is.xts(y)) { y <- as.matrix(y) }
if (! (is.numeric(y) || is.matrix(y))) { stop("y must be numeric or matrix") }
if (is.matrix(y) && ! (ncol(y) == 1)) { stop("y must be a one column matrix") }
if (! is.matrix(x)) { stop("x must be a matrix") }
if (is.null(colnames(x)))
  {
    colnames(x) <- colnames(x, do.NULL = FALSE, prefix = "X")
    warning('column names of x were automatically created')
  }
if (anyNA(colnames(x))) { stop("x must have column names") }
if (is.matrix(y) && is.null(colnames(y)))
  {
    warning('column name of y was automatically created')
    colnames(y) <- colnames(y, do.NULL = FALSE, prefix = "Y")
  }
if (is.matrix(y) && anyNA(colnames(y)))
  {
    warning('column name of y was automatically created')
    colnames(y) <- "Y1"
  }
if (! length(y) == nrow(x)) { stop("y and x must have the same number of observations") }
if (anyNA(y)) { stop("missing values in y") }
if (anyNA(x)) { stop("missing values in x") }
if (missing(alpha)) { stop("please, specify alpha") }
if (! missing(alpha) && ! is.numeric(alpha)) { stop("alpha must be numeric") }
if ((alpha <= 0) || (alpha > 1)) { stop("alpha must be greater than 0, and less than or equal to 1") }
if (missing(lambda)) { stop("please, specify lambda") }
if (! missing(lambda) && ! is.numeric(lambda)) { stop("lambda must be numeric or numeric vector") }
if ((any(lambda <= 0)) || (any(lambda > 1))) { stop("lambda must be greater than 0, and less than or equal to 1") }
if (missing(initvar)) { stop("please, specify initvar (i.e., initial variance)") }
if (! missing(initvar) && ! is.numeric(initvar)) { stop("initvar must be numeric") }
if (initvar <= 0) { stop("variance (initvar) must be positive") }
if (is.null(W)) { W <- "reg" }
if (! ((W == "reg") || is.numeric(W)) ) { stop("please, specify correct W (i.e., initial variance)") }
if (is.numeric(W) && (W <= 0)) { stop("variance (W) must be positive") }
if (W == "reg") { W <- NULL }
if (is.null(initial.period)) { initial.period <- 1 }
if (! is.numeric(initial.period)) { stop("initial.period must be numeric") }
if ((initial.period <= 0) || (initial.period > length(y))) { stop("initial.period must be greater than or equal to 1, and less than the number of observations") }
if (is.null(V.meth)) { V.meth <- "rec" }
if (V.meth == "rec" && !is.null(kappa)) { stop("kappa is used only if V.meth is set to ''ewma''") }
if (! V.meth %in% c("rec","ewma")) { stop("please, specify correct V.meth") }
if (V.meth == "ewma" && is.null(kappa)) { stop("please, specify kappa") }
if (V.meth == "ewma" && ! is.numeric(kappa)) { stop("kappa must be numeric") }
if ( (V.meth == "ewma") && ( (kappa < 0) || (kappa > 1) ) ) { stop("kappa must be between 0 and 1") }
if ( ! is.null(gprob) && ! is.matrix(gprob) ) { stop("gprob must be a matrix") }
if ( ! is.null(gprob) && !(length(y) == nrow(gprob)) ) { warning("time-series of gprob and x differ in length") }
if (! is.null(gprob))
  {
    gprob.ind <- nrow(x) - nrow(gprob) + 1
  }
else
  {
    gprob.ind <- nrow(x) + 1
  }
if ( ! is.null(gprob) && !(ncol(x) == ncol(gprob)) ) { stop("gprob and x must have the same number of columns") }
if ( ! is.null(gprob) && anyNA(gprob) ) { stop("missing values in gprob") }
if ( (! is.null(gprob)) && ( (length(gprob[gprob<0]) > 0) || (length(gprob[gprob>1]) > 0) ) )
  {
    stop("values of gprob must be greater than or equal to 0, and less than or equal to 1")
  }
if ( (! is.null(gprob)) &&  (is.null(omega) ) ) { stop("please, specify omega") }
if (! is.null(omega) && ! is.numeric(omega)) { stop("omega must be numeric") }
if (! is.null(omega) && ( (omega < 0) || (omega > 1) ) ) { stop("omega must be greater than or equal to 0, and less than or equal to 1") }
if (is.null(model)) { model <- "dma" }
if (! model %in% c("dma","dms","med")) { stop("please, specify correct model: ''dma'', ''dms'', or ''med''") }
if (is.null(parallel)) { parallel <- FALSE }
if (! is.logical(parallel)) { stop("parallel must be logical, i.e., TRUE or FALSE") }
if (is.null(m.prior)) { m.prior <- 0.5 }
if (! is.null(m.prior) && ! is.numeric(m.prior)) { stop("m.prior must be numeric") }
if ((m.prior <= 0) || (m.prior >= 1)) { stop("m.prior must be greater than 0, and less than 1") }
if (is.null(mods.incl))
  {
    all.mods <- TRUE
    mods.incl <- expand.grid(rep.int(list(0:1), ncol(x)))
    mods.incl <- as.matrix(cbind(rep.int(1,nrow(mods.incl)),mods.incl))
  }
else
  {
    all.mods <- FALSE
  }
  
if (is.null(mods.check)) { mods.check <- FALSE }
if (mods.check == TRUE)
  {
    if ((! is.null (mods.incl)) && (! is.matrix(mods.incl))) { stop("mods.incl must be a matrix") }
    if (is.matrix(mods.incl) && (! (ncol(mods.incl) == (ncol(x)+1)))) { stop("columns of mods.incl do not correspond to variables specified by columns of x") }
    if (is.matrix(mods.incl) && (length(mods.incl[!(mods.incl %in% c(0,1))]) > 0)) { stop("mods.incl should contain only 0 and 1") }
    if (is.matrix(mods.incl) && any(duplicated(mods.incl))) { stop("mods.incl contain duplicated models") }
    if (is.matrix(mods.incl))
      {
        test <- FALSE
        test.row <- rep.int(0,ncol(mods.incl))
        for (i in 1:nrow(mods.incl))
          {
            if (identical(test.row,mods.incl[i,])) { test <- TRUE }
          }
        if (test == TRUE) { stop("mods.incl contain a model with no variables") }
      }
  }
  
if ( (all.mods == TRUE ) && (ncol(x)>29) ) { stop("max number of variables, in case of using all possible models, is 29") }
if ( (all.mods == FALSE ) && (nrow(mods.incl)>2^29) ) { stop("max number of models is 2^29") }

if (all.mods == FALSE && nrow(mods.incl) == 2^ncol(x)) { all.mods <- TRUE }
if (is.null(DOW))
  {
    threshold <- 0
  }
if (all.mods == FALSE && model == "med") { stop("Median Probability Model can be applied only if all possible models with a constant are used, i.e, mods.incl is not specified") }
if (!is.null(DOW) && !is.numeric(DOW)) { stop("DOW must be numeric") }
if (!is.null(DOW) && ((DOW < 0) || (DOW > 1))) { stop("DOW must be between 0 and 1") }
if (!is.null(DOW))
  {
    threshold <- DOW
  }
if (!is.null(DOW) && is.null(DOW.nmods)) { DOW.nmods <- 0 }
if (!is.null(DOW.nmods) && !is.numeric(DOW.nmods)) { stop("DOW.nmods must be numeric") }
if (!is.null(DOW.nmods) && ((DOW.nmods < 2 && !DOW.nmods == 0) || DOW.nmods > nrow(mods.incl))) { stop("DOW.nmods must be greater than or equal to 2, and less than or equal to the number of all possible models") }
if (is.null(DOW.type)) { DOW.type <- "r" }
if (! DOW.type %in% c("r","e")) { stop("please, specify correct DOW.type: ''r'', or ''e''") }
if (!is.null(DOW) && DOW.type %in% c("r","e") && !model == "dma") { stop("Dynamic Occam's Window can be applied only to Dynamic Model Averaging, i.e., model must be ''dma''") }
if (is.null(bm)) { bm <- FALSE }
if (! is.logical(bm)) { stop("bm must be logical, i.e., TRUE or FALSE") }
if (!is.null(DOW.limit.nmods) && !is.numeric(DOW.limit.nmods)) { stop("DOW.limit.nmods must be numeric") }
if (!is.null(DOW.limit.nmods) && DOW.limit.nmods < 2) { stop("DOW.limit.nmods must be greater than or equal to 2") }
rm(all.mods)
if (length(lambda)>1 && model=="med") stop("multiple lambdas cannot be used for Median Probability Model")
if (length(lambda)>1 && !threshold==0) stop("multiple lambdas cannot be used with Dynamic Occam's Window")
if (is.null(progress.info)) { progress.info <- FALSE }
if (! is.logical(progress.info)) { stop("progress.info must be logical, i.e., TRUE or FALSE") }
if (! is.null(small.c) && ! is.numeric(small.c)) { stop("small.c must be a (small) number") }
if (! is.null(small.c) && (small.c<0)) { stop("small.c must be positive") }
if (! is.null(forced.models)) { forced.models <- forced.models[!duplicated(forced.models),,drop=FALSE] }
if (! is.null(forbidden.models)) { forbidden.models <- forbidden.models[!duplicated(forbidden.models),,drop=FALSE] }
if (is.null(fcores)) { fcores <- 1 }
if (is.null(red.size)) { red.size <- FALSE }

lambda <- unique(lambda)
lambda <- sort(lambda,decreasing=TRUE)

mods.incl <- matrix(as.logical(mods.incl), dim(mods.incl)) 
if (! is.null(forced.models)) { forced.models <- matrix(as.logical(forced.models), dim(forced.models)) }
if (! is.null(forbidden.models)) { forbidden.models <- matrix(as.logical(forbidden.models), dim(forbidden.models)) }
if (! is.null(forced.variables)) { forced.variables <- as.logical(forced.variables) }
if (is.null(av)) { av <- "dma" } 

###################################################################
###################################################################
################################################ basic DMA function

f.basicDMA <- function (y,x,alpha,lambda,initvar,W=NULL,kappa=NULL,
                        gprob=NULL,omega=NULL,model=NULL,parallel=NULL,
                        m.prior=NULL,mods.incl=NULL,small.c=NULL)

  {
    ################################ setting initial values for DMA

    len <- length(lambda)

    exp.lambda <- vector()

      ### small constant as in Eq. (17) in Raftery et al. (2010)
    if (is.null(small.c))
      {
        c <- 0.001 * (1/(len*(2^ncol(x))))
      }
    else
      {
        c <- small.c
      }

      ### initialization of \pi_{0|0} as in Raftery et al. (2010)
    if (m.prior == 0.5)
      {
        pi1 <- as.vector(rep.int(1/(nrow(mods.incl)*len),len*nrow(mods.incl)))
      }
    else
      {
        pi1 <- vector()
        for (i in 1:nrow(mods.incl))
          {
             pi1[i] <- m.prior^(sum(mods.incl[i,])) * (1-m.prior)^(ncol(mods.incl)-sum(mods.incl[i,]))
          }
        pi1 <- rep(pi1,len)
        pi1 <- pi1 / sum(pi1)
      }

      ### yhat.all - predictions from all regression models used in averaging
      ### post     - posterior predictive model probabilities
      ###
      ### in DMA all models are averaged, in DMS the model with highest
      ### posterior predictive model probability is chosen,
      ### for MED see Barbieri and Berger (2004)

    if (model == "dma")
      {
       if (red.size==FALSE)
        {
          post <- as.vector(rep.int(0,len*nrow(mods.incl)))
          yhat.all <- matrix(0,ncol=len*nrow(mods.incl),nrow=1)
        }
       else
        {
          post <- as.vector(rep.int(0,len*ncol(mods.incl)))
          yhat.all <- NA
          vm <- as.matrix(apply(mods.incl,1,sum))
          exp.var <- vector()
        }
      }
    if (model == "dms" || model == "med")
      {
        post <- vector()
        p.incl <- matrix(0,ncol=ncol(mods.incl),nrow=1)
      }

    ydma <- vector()
    thetas.exp <- as.vector(rep.int(0, ncol(mods.incl)))

    ###############################################################
    ############################ recursive estimation of sub-models

    f.param <- function(i)
      {
        i.old <- i
        i <- i.old %% nrow(mods.incl)
        if (i==0) { i <- nrow(mods.incl) }
        i.l <- 1 + ( (i.old - 1) %/% nrow(mods.incl) )
        if (length(which(mods.incl[i,-1,drop=FALSE]==1))>0)
          {
            xx <- x[,which(mods.incl[i,-1,drop=FALSE]==1),drop=FALSE]
            if (mods.incl[i,1]==1)
              {
                c.incl <- TRUE
              }
            else
              {
                c.incl <- FALSE
              }
          }
         else
          {
            xx <- matrix(,ncol=0,nrow=length(y))
            c.incl <- TRUE
          }
        out <- tvp(y=y,x=xx,V=initvar,lambda=lambda[i.l],W=W,kappa=kappa,c=c.incl)
        out[[4]] <- NULL
        return(out)
      }

    if (parallel == TRUE)
      {
        ncores <- detectCores() - fcores
        cl <- makeCluster(ncores)
        registerDoParallel(cl)
        est.models <- foreach(i=seq(len*nrow(mods.incl)),.packages=c("xts","fDMA")) %dopar% 
                        {
                          f.param(i)
                        }
        stopCluster(cl)
        rm(cl,ncores)
      }
    else
      {
        est.models <- lapply(seq(len*nrow(mods.incl)),f.param)
      }

    ############################################## models averaging

      ### modification as in Koop and Onorante (2014)
    f.pi2.g <- function(i)
      {
        if (length(which(mods.incl[i,-1,drop=FALSE]==1))>0)
          {
            p1 <- gprob[t-gprob.ind+1,which(mods.incl[i,-1,drop=FALSE]==1)]
          }
        else
          {
            p1 <- 1
          }
        if (length(which(mods.incl[i,-1,drop=FALSE]==0))>0)
          {
            p2 <- 1-gprob[t-gprob.ind+1,which(mods.incl[i,-1,drop=FALSE]==0)]
          }
        else
          {
            p2 <- 1
          }
        p1 <- prod(p1)
        p2 <- prod(p2)
        if (p1==0) { p1 <- 0.001 * (1/(2^ncol(mods.incl))) }
        if (p2==0) { p2 <- 0.001 * (1/(2^ncol(mods.incl))) }
        return( p1*p2 )
      }

    f.pdens <- function(i)
      {
        return(.subset2(.subset2(est.models,i),3)[t])
      }
    f.mse <- function(i)
      {
        out.mse <- mean((y[1:t] - .subset2(.subset2(est.models,i),1)[1:t])^2)
        if (out.mse==0) { out.mse <- 0.001 * (1/(2^ncol(mods.incl))) }
        return(1/out.mse)
      }
    f.hr1 <- function(i)
      {
        out.hr <- hit.ratio(y=y[1:t],y.hat=.subset2(.subset2(est.models,i),1)[1:t],d=FALSE)
        if (out.hr==0) { out.hr <- 0.001 * (1/(2^ncol(mods.incl))) }
        return(out.hr)
      }
    f.hr2 <- function(i)
      {
        out.hr <- hit.ratio(y=y[1:t],y.hat=.subset2(.subset2(est.models,i),1)[1:t],d=TRUE)
        if (out.hr==0) { out.hr <- 0.001 * (1/(2^ncol(mods.incl))) }
        return(out.hr)
      }

    f.yhat <- function(i)
      {
        return(.subset2(.subset2(est.models,i),1)[t])
      }

    f.thetas <- function(i)
      {
        i.old <- i
        i <- i.old %% nrow(mods.incl)
        if (i==0) { i <- nrow(mods.incl) }

        theta.i.tmp <- as.vector(mods.incl[i,])
        theta.i.tmp[theta.i.tmp==1] <- .subset2(.subset2(est.models,i.old),2)[t,]
        return(theta.i.tmp)
      }

    for (t in 1:nrow(x))
      {
        if (t<gprob.ind)
          {
            pi2 <- (pi1^alpha + c) / (sum((pi1)^alpha  + c))
          }
        else
          {
            pi2.g <- unlist(lapply(seq(nrow(mods.incl)),f.pi2.g))
            pi2.g <- pi2.g / sum(pi2.g)
            pi2.g <- pi2.g / len
            pi2.g <- rep(pi2.g,len)

            pi1sum <- sum((pi1)^alpha  + c)

            pi2 <- omega * ( (pi1^alpha + c) / pi1sum ) + (1-omega) * pi2.g

            rm(pi1sum,pi2.g)
          }

        if (model == "dma")
          {
            if (red.size==FALSE) 
              { 
                post <- rbind(post,pi2) 
              }
            else
              {
                f.split.pi2 <- function(i.col)
                  {
                    col.s <- ( ( i.col - 1 ) * nrow(mods.incl) ) + 1
                    col.e <- col.s + nrow(mods.incl) - 1
                    return(pi2[col.s:col.e])
                  }
                pi2.temp <- lapply(1:len,f.split.pi2)
                pi2.temp <- Reduce('+', pi2.temp)
                pi2.temp <- as.vector(pi2.temp)
                post <- rbind(post,pi2.temp %*% mods.incl)
                
                exp.var[t] <- pi2.temp %*% vm
                rm(pi2.temp)
              }
          }
        if (model == "dms")
          {
            j.m <- which.max(pi2)
            post[t] <- pi2[j.m]
            j.m.new <- j.m %% nrow(mods.incl)
            if (j.m.new==0) { j.m.new <- nrow(mods.incl) }
            p.incl <- rbind(p.incl,mods.incl[j.m.new,])
          }
        if (model == "med")
          {
            j.m <- as.vector(pi2 %*% mods.incl)
            j.m1 <- which(j.m >= 0.5)
            j.m <- as.vector(rep.int(0, ncol(mods.incl)))
            j.m[j.m1] <- 1
            j.m <- which(apply(mods.incl, 1, function(x) all(x == j.m)))
            rm(j.m1)
            post[t] <- pi2[j.m]
            p.incl <- rbind(p.incl,mods.incl[j.m,])
          }

        yhat <- unlist(lapply(seq(len*nrow(mods.incl)),f.yhat))

        if (model == "dma")
          {
            if (red.size==FALSE) { yhat.all <- rbind(yhat.all,yhat) }
          }

        if (model == "dma")
          {
            ydma[t] <- crossprod(pi2,yhat)
            
            thetas <- t(sapply(seq(len*nrow(mods.incl)),f.thetas))
            thetas.exp <- rbind(thetas.exp,pi2 %*% thetas)

            exp.lambda[t] <- as.numeric(crossprod(pi2,rep(lambda,each=nrow(mods.incl))))
          }
        if (model == "dms" || model == "med")
          {
            ydma[t] <- yhat[j.m]

            thetas <- f.thetas(j.m)
            thetas.exp <- rbind(thetas.exp,thetas)

            exp.lambda[t] <- lambda[1 + (j.m - 1) %/% nrow(mods.incl)]

            rm(j.m)
          }

        if (av=="dma") { pdens <- unlist(lapply(seq(len*nrow(mods.incl)),f.pdens)) }
        if (av=="mse") { pdens <- unlist(lapply(seq(len*nrow(mods.incl)),f.mse)) } 
        if (av=="hr1") 
          { 
            if (t==1)
              {
                pdens <- rep(1,len*nrow(mods.incl))
              }
            else
              {
                pdens <- unlist(lapply(seq(len*nrow(mods.incl)),f.hr1)) 
              }
          } 
        if (av=="hr2") { pdens <- unlist(lapply(seq(len*nrow(mods.incl)),f.hr2)) } 
        
        pi1 <- (pi2 * pdens) / as.numeric(crossprod(pi2,pdens))
      }

    rm(est.models)
    
    ######################################################## output

    thetas.exp <- thetas.exp[-1,,drop=FALSE]

    if (model == "dma")
      {
        post <- as.matrix(post[-1,,drop=FALSE])

        if (red.size==FALSE)
          {
            yhat.all <- yhat.all[-1,,drop=FALSE]

            f.split <- function(i.col)
              {
                col.s <- ( ( i.col - 1 ) * nrow(mods.incl) ) + 1
                col.e <- col.s + nrow(mods.incl) - 1
                return(post[,col.s:col.e])
              }

            post <- lapply(1:len,f.split)
            post <- Reduce('+', post)
            post <- as.matrix(post)
          }
          
        if (red.size==FALSE) { exp.var <- NA }
        
        return(list(post,yhat.all,thetas.exp,ydma,pdens,NA    ,pi1,thetas,yhat,pi2,exp.lambda, NA, exp.var))
      }
    if (model == "dms" || model == "med")
      {
        p.incl <- p.incl[-1,,drop=FALSE]
        return(list(post,NA      ,thetas.exp,ydma,pdens,p.incl,pi1,thetas,yhat,pi2,exp.lambda))
      }

  }

################################################## end of basic DMA
###################################################################
###################################################################

if (threshold==0)
  {
    comput <- f.basicDMA(y=y,x=x,alpha=alpha,lambda=lambda,initvar=initvar,W=W,kappa=kappa,
                         gprob=gprob,omega=omega,model=model,parallel=parallel,
                         m.prior=m.prior,mods.incl=mods.incl,small.c=small.c)
  }
else
  {
    ##### Dynamic Occam's Window as in Onorante and Raftery (2016)
    ##############################################################

    ydma.dow <- vector()
    post.dow <- as.vector(rep.int(0, ncol(mods.incl)))
    thetas.exp.dow <- as.vector(rep.int(0,ncol(mods.incl)))
    exp.var.dow <- vector()

      ### selection of initial models

    if (!(DOW.nmods == 0))
      {
        mods.incl <- mods.incl[sample.int(nrow(mods.incl),size=DOW.nmods),]
      }
    else
      {
        mods.incl <- as.matrix(rbind(rep.int(0,ncol(x)),diag(1,nrow=ncol(x),ncol=ncol(x))))
        mods.incl <- as.matrix(cbind(rep.int(1,ncol(x)+1),mods.incl))
        mods.incl <- matrix(as.logical(mods.incl), dim(mods.incl)) 
      }
    init.mods <- mods.incl
    nms <- vector()

    for (T in 1:nrow(x))
      {
        nms[T] <- nrow(mods.incl)

        if (progress.info==TRUE)
          {
            print(paste('round:',T))
            print(paste('models:',nms[T]))
          }

        if (parallel==TRUE && nms[T]>2^10)
          {
            dopar <- TRUE
          }
        else
          {
            dopar <- FALSE
          }
        T.dma <- f.basicDMA(y=y[1:T,,drop=FALSE],x=x[1:T,,drop=FALSE],alpha=alpha,lambda=lambda,initvar=initvar,W=W,kappa=kappa,
                            gprob=gprob,omega=omega,model=model,parallel=dopar,
                            m.prior=m.prior,mods.incl=mods.incl,small.c=small.c)
        yhat <- T.dma[[9]]
        thetas <- T.dma[[8]]
        T.dma[[8]] <- NA
        pi1 <- T.dma[[7]]

        if (DOW.type == "r")
          {
            pi2.temp <- T.dma[[10]]
            pi2.temp[which(pi1<(threshold*max(pi1)))] <- 0
            pi2.temp <- pi2.temp / sum(pi2.temp)
            ydma.dow[T] <- crossprod(pi2.temp,yhat)

            post.dow <- rbind(post.dow,pi2.temp %*% mods.incl)
            exp.var.dow[T] <- pi2.temp %*% (as.matrix(apply(mods.incl,1,sum)))

            thetas <- pi2.temp %*% thetas
            thetas.exp.dow <- rbind(thetas.exp.dow,thetas)

            rm(pi2.temp)
          }
        if (DOW.type == "e")
          {
            ydma.dow[T] <- crossprod(T.dma[[10]],yhat)

            post.dow <- rbind(post.dow,T.dma[[10]] %*% mods.incl)
            exp.var.dow[T] <- T.dma[[10]] %*% (as.matrix(apply(mods.incl,1,sum)))

            thetas <- T.dma[[10]] %*% thetas
            thetas.exp.dow <- rbind(thetas.exp.dow,thetas)
          }

          ### models' space update

        mod.red <- mods.incl[which(pi1>=(threshold*max(pi1))),,drop=FALSE]

        if (!is.matrix(mod.red)) { mod.red <- t(as.matrix(mod.red)) }

        if (!is.null(DOW.limit.nmods))
          {
            if (nrow(mod.red)>DOW.limit.nmods)
              {
                ind.red <- sort(pi1,decreasing=TRUE,index.return=TRUE)$ix
                ind.red <- ind.red[1:DOW.limit.nmods]
                mod.red <- mods.incl[ind.red,,drop=FALSE]
                rm(ind.red)
              }
          }

        if (!is.matrix(mod.red)) { mod.red <- t(as.matrix(mod.red)) }

        mod.exp <- mod.red
        unm <- matrix(TRUE,nrow(mod.red),1)

        for (i.col in 2:ncol(mod.red))
          {
            mod.exp.temp <- mod.red
            mod.exp.temp[,i.col] <- xor(mod.red[,i.col],unm)
            mod.exp <- rbind(mod.exp,mod.exp.temp)
          }

        if (! is.null(forced.models))
          {
            mod.exp <- rbind(mod.exp,forced.models)
          }

        if (! is.null(forced.variables))
          {
            for (i.var in which(forced.variables==1))
              {
                mod.exp[,i.var] <- TRUE
              }
          }

        mod.exp <- unique(mod.exp)

        if (! is.null(forbidden.models))
          {
            mod.exp <- rbind(forbidden.models,mod.exp)
            mod.exp <- mod.exp[!duplicated(mod.exp),,drop=FALSE]
            mod.exp <- mod.exp[-(1:nrow(forbidden.models)),,drop=FALSE]
          }

        if (T<nrow(x)) { mods.incl <- mod.exp }

        rm(mod.red,mod.exp,mod.exp.temp)
      }

    comput <- list(post.dow[-1,,drop=FALSE],NA,thetas.exp.dow[-1,,drop=FALSE],ydma.dow,T.dma[[5]])
    comput[[11]] <- T.dma[[11]]
    comput[[12]] <- exp.var.dow

    rm(post.dow,ydma.dow,thetas.exp.dow,T.dma)
  }

############################################################ output
###################################################################

ydma  <- comput[[4]]
comput[[4]] <- NA
pdens <- comput[[5]]
comput[[5]] <- NA

if (length(lambda)>1)
  {
    exp.lambda <- comput[[11]]
  }
else
  {
    exp.lambda <- rep.int(lambda,nrow(x))
  }

if (model == "dma")
  {
    if (threshold==0) 
      { 
        post <- comput[[1]] 
      }

    yhat.all <- comput[[2]]

    if (threshold==0)
      {
        if (red.size==FALSE)
          {
            colnames(yhat.all) <- seq(1,length(lambda)*nrow(mods.incl))
            rownames(yhat.all) <- rownames(x)
          }
      }
  }

thetas.exp <- as.matrix(comput[[3]])
colnames(thetas.exp) <- c("const",colnames(x))

if (model == "dma")
  {
   if (threshold==0)
    {
      if (red.size==FALSE) 
        { 
          post.inc <- post %*% mods.incl 
        }
      else
        {
          post.inc <- post
          post <- NA
        }
    }
  else
    {
      post.inc <- comput[[1]]
    }
  }
if (model == "dms" || model == "med")
  {
    post.inc <- comput[[6]]
  }
colnames(post.inc) <- c("const",colnames(x))

if (model == "dma")
  {
    if (threshold==0)
      {
        if (red.size==FALSE) 
          { 
            exp.var <- post %*% (as.matrix(apply(mods.incl,1,sum))) 
          } 
        else 
          { 
            exp.var <- as.matrix(comput[[13]])
          }
      }
    else
      {
        exp.var <- as.matrix(comput[[12]])
      }
  }
if (model == "dms" || model == "med")
  {
    exp.var <- as.matrix(rowSums(post.inc))
  }

if (model == "dms" || model == "med")
  {
    post <- as.matrix(comput[[1]])
    yhat.all <- NA
  }

mse <- (mean((y[initial.period:length(y)] - ydma[initial.period:length(y)])^2))^(1/2)
mae <- mean(abs(y[initial.period:length(y)] - ydma[initial.period:length(y)]))

naive.mse <- (mean((diff(as.numeric(y[initial.period:length(as.numeric(y))])))^2))^(1/2)
naive.mae <- mean(abs(diff(as.numeric(y[initial.period:length(as.numeric(y))]))))

if (bm == TRUE)
  {
    arima <- auto.arima(as.numeric(y))$residuals[initial.period:length(as.numeric(y))]
    arima.mse <- (mean((arima)^2))^(1/2)
    arima.mae <- mean(abs(arima))
  }
else
  {
    arima.mse <- NA
    arima.mae <- NA
  }

benchmarks <- rbind(cbind(naive.mse,arima.mse),cbind(naive.mae,arima.mae))
rownames(benchmarks) <- c("RMSE", "MAE")
colnames(benchmarks) <- c("naive", "auto ARIMA")

if (model == "dma")
  {
    mod.type <- "DMA"
  }
if (model == "dms")
  {
    mod.type <- "DMS"
  }
if (model == "med")
  {
    mod.type <- "MED"
  }

if (is.null(kappa)) { kappa <- NA }
if (is.null(omega)) { omega <- NA }
if (is.null(W)) { W <- "reg" }
if (threshold==0)
  {
    DOW.nmods <- NA
    DOW.type <- NA
  }

if (length(lambda)>1) { lambda <- paste(lambda,collapse=" ") }

if (threshold==0)
  {
    temp <- list(ydma, post.inc, mse, mae, mods.incl+0, post, exp.var, thetas.exp,
                 cbind(alpha,lambda,initvar,mod.type,W,initial.period,V.meth,kappa,omega,m.prior,threshold,DOW.nmods,DOW.type),
                 yhat.all, y, benchmarks, NA,         NA, pdens,exp.lambda)
  }
else
  {
    colnames(init.mods) <- colnames(post.inc)
    rownames(init.mods) <- seq(1,nrow(init.mods))
    temp <- list(ydma, post.inc, mse, mae, mods.incl+0,   NA, exp.var, thetas.exp,
                 cbind(alpha,lambda,initvar,mod.type,W,initial.period,V.meth,kappa,omega,m.prior,threshold,DOW.nmods,DOW.type),
                 yhat.all, y, benchmarks, init.mods+0, nms, pdens,exp.lambda)
  }

rownames(temp[[2]]) <- rownames(x)
colnames(temp[[5]]) <- c("const", colnames(x))
if (threshold==0 & (!any(is.na(temp[[6]])))) { rownames(temp[[6]]) <- rownames(x) }
if (model == "dma")
  {
    if (threshold==0 & (!any(is.na(temp[[6]])))) { colnames(temp[[6]]) <- seq(1,nrow(mods.incl)) }
  }
if (model == "dms" || model == "med")
  {
    colnames(temp[[6]]) <- "mod. prob."
  }
rownames(temp[[7]]) <- rownames(x)
rownames(temp[[8]]) <- rownames(x)
colnames(temp[[9]]) <- c("alpha","lambda","initvar","model type","W","initial period","V.meth","kappa","omega","m.prior","DOW threshold","DOW.nmods","DOW.type")
rownames(temp[[9]]) <- "parameters"
names(temp) <- c("y.hat","post.incl","RMSE","MAE","models","post.mod","exp.var","exp.coef.","parameters",
                 "yhat.all.mods","y","benchmarks","DOW.init.mods","DOW.n.mods.t","p.dens.","exp.lambda")

class(temp) <- "dma"

return(temp)

}
