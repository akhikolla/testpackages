#' Cross Validation for the vennLasso
#'
#' @param x input matrix or SparseMatrix of dimension nobs x nvars. Each row is an observation,
#' each column corresponds to a covariate
#' @param y numeric response vector of length nobs
#' @param groups A list of length equal to the number of groups containing vectors of integers
#' indicating the variable IDs for each group. For example, groups=list(c(1,2), c(2,3), c(3,4,5)) specifies
#' that Group 1 contains variables 1 and 2, Group 2 contains variables 2 and 3, and Group 3 contains
#' variables 3, 4, and 5. Can also be a matrix of 0s and 1s with the number of columns equal to the
#' number of groups and the number of rows equal to the number of variables. A value of 1 in row i and
#' column j indicates that variable i is in group j and 0 indicates that variable i is not in group j.
#' @param lambda A user-specified sequence of lambda values. Left unspecified, the a sequence of lambda values is
#' automatically computed, ranging uniformly on the log scale over the relevant range of lambda values.
#' @param compute.se logical flag. If \code{TRUE}, standard errors will be computed, otherwise if \code{FALSE} they will not
#' @param conf.int value between 0 and 1 indicating the level of the confidence intervals to be computed. For example
#' if \code{conf.int = 0.95}, 95 percent confidence intervals will be computed.
#' @param type.measure One of \code{c("mse","deviance","class","auc","mae","brier")} indicating measure to evaluate for cross-validation. The default is \code{type.measure = "deviance"}, 
#' which uses squared-error for gaussian models (a.k.a \code{type.measure = "mse"} there), deviance for logistic
#' regression. \code{type.measure = "class"} applies to binomial only. \code{type.measure = "auc"} is for two-class logistic 
#' regression only. \code{type.measure = "mse"} or \code{type.measure = "mae"} (mean absolute error) can be used by all models;
#' they measure the deviation from the fitted mean to the response. \code{type.measure = "brier"} is for models with 
#' \code{family = "coxph"} and will compute the Brier score.
#' @param nfolds number of folds for cross-validation. default is 10. 3 is smallest value allowed. 
#' @param foldid an optional vector of values between 1 and nfold specifying which fold each observation belongs to.
#' @param grouped Like in \pkg{glmnet}, this is an experimental argument, with default \code{TRUE}, and can be ignored by most users. 
#' For all models, this refers to computing nfolds separate statistics, and then using their mean and estimated standard 
#' error to describe the CV curve. If \code{grouped = FALSE}, an error matrix is built up at the observation level from the 
#' predictions from the \code{nfold} fits, and then summarized (does not apply to \code{type.measure = "auc"}). 
#' @param keep If \code{keep = TRUE}, a prevalidated list of arrasy is returned containing fitted values for each observation 
#' and each value of lambda for each model. This means these fits are computed with this observation and the rest of its
#' fold omitted. The folid vector is also returned. Default is \code{keep = FALSE}
#' @param parallel If TRUE, use parallel foreach to fit each fold. Must register parallel before hand, such as \pkg{doMC}.
#' @param ... parameters to be passed to vennLasso
#' @importFrom stats predict
#' @importFrom stats stepfun
#' @importFrom stats weighted.mean
#' @import foreach
#' @return An object with S3 class "cv.vennLasso"
#'
#'
#' @import Rcpp
#'
#' @export
#' @examples
#' 
#' library(Matrix)
#' 
#' set.seed(123)
#' n.obs <- 150
#' n.vars <- 25
#'
#' true.beta.mat <- array(NA, dim = c(3, n.vars))
#' true.beta.mat[1,] <- c(-0.5, -1, 0, 0, 2, rep(0, n.vars - 5))
#' true.beta.mat[2,] <- c(0.5, 0.5, -0.5, -0.5, 1, -1, rep(0, n.vars - 6))
#' true.beta.mat[3,] <- c(0, 0, 1, 1, -1, rep(0, n.vars - 5))
#' rownames(true.beta.mat) <- c("1,0", "1,1", "0,1")
#' true.beta <- as.vector(t(true.beta.mat))
#'
#' x.sub1 <- matrix(rnorm(n.obs * n.vars), n.obs, n.vars)
#' x.sub2 <- matrix(rnorm(n.obs * n.vars), n.obs, n.vars)
#' x.sub3 <- matrix(rnorm(n.obs * n.vars), n.obs, n.vars)
#'
#' x <- as.matrix(rbind(x.sub1, x.sub2, x.sub3))
#'
#' conditions <- as.matrix(cbind(c(rep(1, 2 * n.obs), rep(0, n.obs)),
#'                               c(rep(0, n.obs), rep(1, 2 * n.obs))))
#'
#' y <- rnorm(n.obs * 3, sd = 3) + drop(as.matrix(bdiag(x.sub1, x.sub2, x.sub3)) %*% true.beta)
#'
#' fit <- cv.vennLasso(x = x, y = y, groups = conditions, nfolds = 3)
#'
#' fitted.coef <- predict(fit$vennLasso.fit, type = "coefficients", s = fit$lambda.min)
#' (true.coef <- true.beta.mat[match(dimnames(fit$vennLasso.fit$beta)[[1]], 
#'                                   rownames(true.beta.mat)),])
#' round(fitted.coef, 2)
#'
#' ## effects smaller for logistic regression
#' \dontrun{
#' true.beta.mat <- true.beta.mat / 2
#' true.beta <- true.beta / 2
#' # logistic regression example#'
#' y <- rbinom(n.obs * 3, 1, 
#'        prob = 1 / (1 + exp(-drop(as.matrix(bdiag(x.sub1, x.sub2, x.sub3)) %*% true.beta))))
#'
#' bfit <- cv.vennLasso(x = x, y = y, groups = conditions, family = "binomial",
#'                      nfolds = 3)
#'
#' fitted.coef <- predict(bfit$vennLasso.fit, type = "coefficients", s = bfit$lambda.min)
#' (true.coef <- true.beta.mat[match(dimnames(bfit$vennLasso.fit$beta)[[1]], 
#'                                   rownames(true.beta.mat)),])
#' round(fitted.coef, 2)
#' }
#'
cv.vennLasso <- function(x, y,
                         groups,
                         lambda       = NULL,
                         compute.se   = FALSE,
                         conf.int     = NULL,
                         type.measure = c("mse","deviance","class","auc","mae","brier"),
                         nfolds       = 10,
                         foldid,
                         grouped      = TRUE,
                         keep         = FALSE,
                         parallel     = FALSE,
                         ...)
{
  if(missing(type.measure)) type.measure <- "default"
  else type.measure <- match.arg(type.measure)
  if(!is.null(lambda)&&length(lambda)<2) stop("Need more than one value of lambda for cv.vennLasso")
  N <- nrow(x)
  #if(missing(weights))weights=rep(1.0,N)else weights=as.double(weights)
  ###Fit the model once to get dimensions etc of output
  y <- drop(y) # we dont like matrix responses unless we need them
  ###Next we construct a call, that could recreate a vennLasso object - tricky
  ### This if for predict, exact=TRUE
  vennLasso.call <- match.call(expand.dots = TRUE)
  which <- match(c("type.measure","nfolds","foldid","grouped","keep"), names(vennLasso.call), FALSE)

  if(any(which)) vennLasso.call <- vennLasso.call[-which]

  vennLasso.call[[1]] <- as.name("vennLasso")
  vennLasso.object <- vennLasso(x, y, groups, 
                                lambda = lambda,
                                compute.se = compute.se, 
                                conf.int = conf.int, ...)
  vennLasso.object$call=vennLasso.call
  #is.offset=vennLasso.object$offset
  #lambda=vennLasso.object$lambda
  #if(inherits(vennLasso.object,"multnet")){
  #  nz=predict(vennLasso.object,type="nonzero")
  #  nz=sapply(nz,function(x)sapply(x,length))
  #  nz=ceiling(apply(nz,1,median))
  #}
  #nz=sapply(predict(vennLasso.object,type="nonzero"),length)
  #if(missing(foldid)) foldid=sample(rep(seq(nfolds),length=N)) else nfolds=max(foldid)
  if(missing(foldid)) 
  {
    foldid <- sample(rep(seq(nfolds), length = N))
    foldid <- vector(length = N)
    for (c in 1:length(vennLasso.object$data.indices)) 
    {
      foldid[vennLasso.object$data.indices[[c]]] <-
        sample(rep(seq(nfolds), length = length(vennLasso.object$data.indices[[c]])))
    }
  } else 
  {
    nfolds <- max(foldid)
  }

  if(nfolds < 3) stop("nfolds must be bigger than 3; nfolds=10 recommended")
  outlist <- as.list(seq(nfolds))
  ###Now fit the nfold models and store them
  ###First try and do it using foreach if parallel is TRUE
  if (parallel) 
  {
    outlist = foreach (i=seq(nfolds), .packages=c("vennLasso")) %dopar% {
      which <- foldid==i
      if(is.matrix(y)) y_sub <- y[!which,] else y_sub <- y[!which]
      #if(is.offset)offset_sub=as.matrix(offset)[!which,]
      #else offset_sub=NULL
      vennLasso(x[!which,,drop=FALSE], 
                y_sub, 
                groups[!which,,drop=FALSE], 
                lambda = lambda,
                compute.se = FALSE, 
                conf.int = NULL, ...)
    }
  } else
  {
    for(i in seq(nfolds))
    {
      which <- foldid==i
      if(is.matrix(y))y_sub=y[!which,]else y_sub=y[!which]
      #if(is.offset)offset_sub=as.matrix(offset)[!which,]
      #else offset_sub=NULL
      #cat("fitting fold", i, "\n")
      outlist[[i]] <- vennLasso(x[!which,,drop=FALSE], 
                                y_sub, 
                                groups[!which,,drop=FALSE], 
                                lambda = lambda,
                                compute.se = FALSE, 
                                conf.int = NULL, ...)
    }
  }
  ###What to do depends on the type.measure and the model fit
  fun     <- paste("cv",class(vennLasso.object)[[2]],sep=".")
  cvstuff <- do.call(fun,
                     list(outlist,
                          vennLasso.object$lambda,
                          x,
                          y,
                          groups,
                          foldid,
                          type.measure,
                          grouped,
                          keep))
  cvm     <- cvstuff$cvm
  cvsd    <- cvstuff$cvsd
  cvname  <- cvstuff$name

  out=list(lambda        = vennLasso.object$lambda,
           cvm           = cvm,
           cvsd          = cvsd,
           cvup          = cvm + cvsd,
           cvlo          = cvm - cvsd,
           name          = cvname,
           vennLasso.fit = vennLasso.object)
  if(keep) out <- c(out,list(fit.preval=cvstuff$fit.preval,foldid=foldid))
  lamin = if(type.measure=="auc") getmin(vennLasso.object$lambda,-cvm,cvsd)
  else getmin(vennLasso.object$lambda,cvm,cvsd)
  obj <- c(out,as.list(lamin))
  class(obj) <- "cv.vennLasso"
  obj
}

cv.venncoxph <- function(outlist,lambda,x,y,groups,foldid,type.measure,grouped,keep=FALSE)
{
    typenames <- c(deviance="Coxph Deviance")
    if(type.measure == "default") type.measure <- "deviance"
    
    if(!match(type.measure,c("deviance","brier"),FALSE))
    {
        warning("Only 'deviance' or 'brier' available for coxph model; as default 'deviance' will be used")
        type.measure="deviance"
    }

    N      <- nrow(y)
    nfolds <- max(foldid)
    if( (N/nfolds < 10) && type.measure == "auc")
    {
        warning("Too few (< 10) observations per fold for type.measure='auc' in cv.lognet; changed to type.measure='deviance'. Alternatively, use smaller value for nfolds",call.=FALSE)
        type.measure="deviance"
    }
    if( (N/nfolds <3) && grouped)
    {
        warning("Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold",call.=FALSE)
        grouped <- FALSE
    }


    #if(!is.null(offset)){
    #  is.offset=TRUE
    #  offset=drop(offset)
    #}else is.offset=FALSE
    devmat <- matrix(NA, nfolds, length(lambda))
    nlams  <- double(nfolds)

    logpl <- function(link, yyy)
    {
        t <- yyy[,1]
        c <- yyy[,2]

        index.t <- order(t)

        link[index.t,] <- apply(link[index.t,], 2, function (tt)
        {
            if (max(tt) > 100)
            {
                kkk <- floor(max(tt) / 100)
                tt  <- tt - 100*kkk
            } else 
            {
                tt <- tt
            }
        })

        link.order.exp <- exp(link[index.t,])
        c.order        <- c[index.t]
        #       r <- apply(exp(link),2,function(tt){rev(cumsum(rev(tt[index.t])))})
        logL <- apply(link.order.exp,2,function(tt){r <- rev(cumsum(rev(tt)))
        sum(log(tt[c.order==1])-log(r[c.order==1]))})
    }

    brier <- function(survive,yytest,yytrain,cenFun,cenTime)
    {
        index.test  <- order(yytest[,1])
        ytest       <- yytest[index.test]
        survive     <- survive[index.test,]

        index.train <- order(yytrain[,1])
        ytrain      <- yytrain[index.train]

        nobs        <- dim(survive)[1]
        brierScore  <- array(0,c(nobs,1))
        sub         <- 100 * nobs

        for (i in 1:nobs)
        {
            if (length(survive[i,])<length(ytrain[,1]))
            {
                survm <- c(rep(1,times = length(ytrain[,1])-length(survive[i,])), survive[i,])
            } else 
            {
                survm <- survive[i,]
            }
            surFun  <- stepfun(ytrain[,1],c(1,survm), f=0)
            proFun1 <- stepfun(ytest[i,1],c(1,0), f=0)
            proFun2 <- stepfun(ytest[i,1],c(0,1), f=0)
            knots   <- sort(unique(c(ytrain[,1],ytest[,1],cenTime)))

            Fun1 <- function(t) { (surFun(t))^2*proFun2(t)}
            Fun2 <- function(t) { (1-surFun(t))^2*proFun1(t)/cenFun(t)}
            if (ytest[i,2] == 1)
            {
                brierScore[i] <- brierScore[i] + integrate.step(Fun1,ytest[i,1],max(ytest[,1]),knots) / cenFun(ytest[i,1]) + integrate.step(Fun2,0,ytest[i,1], knots)
            } else 
            {
                brierScore[i] <- brierScore[i] + integrate.step(Fun2,0,ytest[i,1], knots)
            }
    #         for (j in 1:(dim(ytrain)[1]-1)){
    #             if (ytest[i,1]>ytrain[j,1] && ytest[]){
    #                 brierScore[i] <- brierScore[i] + (1-survive[i,j])^2*(min(ytrain[j+1,1],ytest[i,1])-ytrain[j,1]) / cenFun(ytrain[j,1])
    #             }
    #             if (j == dim(ytrain)[1]-1 & ytest[i,1]>ytrain[j+1,1]){
    #                 brierScore[i] <- brierScore[i] + (1-survive[i,j+1])^2*(ytest[i,1]-ytrain[j+1,1]) / cenFun(ytrain[j,1]-10^(-6))
    #             }
    #         }
    #         print(6)
    #
    #         for (j in 2:(dim(ytrain)[1])){
    #             if (ytest[i,1]<ytrain[j,1] && ytest[i,2]==1){
    #                 brierScore[i] <- brierScore[i] + (survive[i,j-1])^2*(ytrain[j,1]-max(ytrain[j-1,1],ytest[i,1])) / cenFun(ytest[i,1])
    #             }
    #         }
    #         print(7)
    #
        }
        brierScore
    }

    integrate.step <- function(Fun,lower,upper,knots)
    {
        index.lower <- which(knots==lower)
        index.upper <- which(knots==upper)-1
        res <- 0
        for(i in index.lower:index.upper)
        {
            res <- res + Fun(knots[i+1] - 0.01) * (knots[i+1]-knots[i])
        }
        res
    }


    for(i in seq(nfolds))
    {
        which  <- foldid==i
        fitobj <- outlist[[i]]
        #if(is.offset)off_sub=offset[which]
        nlami  <- length(outlist[[i]]$lam)

        if(type.measure == "deviance")
        {
            linkout = predict(fitobj, newx = x[!which,,drop=FALSE], group.mat = groups[!which,,drop=FALSE], type="link")
            print(anyNA(linkout))
            link = predict(fitobj, newx = x, group.mat = groups, type="link")
            print(anyNA(link))
            devmat[i,seq(nlami)] = -logpl(link,y) + logpl(linkout,y[!which,,drop=FALSE])
        } else 
        {
            surv = predict(fitobj, newx = x[which,], group.mat = groups[which,,drop=FALSE], type="survival")
            cenSurv = Surv(time = y[which,1], event = 1-y[which,2])
            cenSurvFun = c(1,1, survfit(cenSurv~1)$surv)
            cenTime = c(0, survfit(cenSurv~1)$time)
            cenFun = stepfun(cenTime, cenSurvFun, f=0)

            for (j in 1:nlami)
            {
                devmat[i,j] <- mean(brier(surv[[j]],y[which,],y[!which,],cenFun,cenTime))
                #print(devmat[i,j])
                #print(surv[[j]][6,100])
            }
        }
        nlams[i] <- nlami
    }

    print(nfolds)
    ###If auc we behave differently
        ##extract weights and normalize to sum to 1
        #ywt=apply(y,1,sum)
        #y=y/ywt
        #weights=weights*ywt
        cvraw=devmat
    #   if(grouped){
    #        cvob=cvcompute(cvraw,rep(1, nrow(y)),foldid,nlams)
    #        cvraw=cvob$cvraw;weights=cvob$weights;N=cvob$N
    #    }
    #cvm=apply(cvraw,2,weighted.mean,w=weights,na.rm=TRUE)
    cvm  <- apply(cvraw, 2, mean)
    #cvsd=sqrt(apply(scale(cvraw,cvm,FALSE)^2,2,weighted.mean,w=rep(1, nrow(y)),na.rm=TRUE)/(N-1))
    cvsd <- apply(cvraw, 2, sd)
    out <- list(cvm  = cvm,
                cvsd = cvsd,
                name = typenames[type.measure])
    if(keep)out$fit.preval <- devmat
    out

}

cv.vennbinomial <- function(outlist,
                            lambda,
                            x,
                            y,
                            groups,
                            foldid,
                            type.measure,
                            grouped,
                            keep = FALSE)
{
  typenames  <- c(mse = "Mean-Squared Error",
                  mae = "Mean Absolute Error",
                  deviance = "Binomial Deviance",
                  auc = "AUC",
                  class = "Misclassification Error")
  if(type.measure == "default") type.measure="deviance"
  if(!match(type.measure, c("mse", "mae", "deviance", "auc", "class"), FALSE))
  {
    warning("Only 'deviance', 'class', 'auc', 'mse' or 'mae'  available for binomial models; 'deviance' used")
    type.measure <- "deviance"
  }

  ###These are hard coded in the Fortran, so we do that here too
  prob_min <- 1e-5
  prob_max <- 1-prob_min
  ###Turn y into a matrix
  nc <- dim(y)
  if (is.null(nc)) 
  {
    y    <- as.factor(y)
    ntab <- table(y)
    nc   <- as.integer(length(ntab))
    y    <- diag(nc)[as.numeric(y), ]
  }
  N      <- nrow(y)
  nfolds <- max(foldid)
  if( (N/nfolds <10) && type.measure == "auc")
  {
    warning("Too few (< 10) observations per fold for type.measure='auc' in cv.lognet; changed to type.measure='deviance'. Alternatively, use smaller value for nfolds",call.=FALSE)
    type.measure <- "deviance"
  }
  if( (N/nfolds <3)&&grouped)
  {
    warning("Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold",call.=FALSE)
    grouped <- FALSE
  }


  #if(!is.null(offset)){
  #  is.offset=TRUE
  #  offset=drop(offset)
  #}else is.offset=FALSE
  predmat <- matrix(NA, nrow(y), length(lambda))
  nlams   <- double(nfolds)
  for(i in seq(nfolds))
  {
    which    <- foldid == i
    fitobj   <- outlist[[i]]
    #if(is.offset)off_sub=offset[which]
    preds    <- predict(fitobj, newx = x[which,,drop=FALSE], group.mat = groups[which,,drop=FALSE], type="response")
    nlami    <- length(outlist[[i]]$lam)
    predmat[which,seq(nlami)] <- preds
    nlams[i] <- nlami
  }
  ###If auc we behave differently
  if(type.measure=="auc")
  {
    cvraw <- matrix(NA,nfolds, length(lambda))
    good  <- matrix(0,nfolds,  length(lambda))
    for(i in seq(nfolds))
    {
      good[i, seq(nlams[i])] <- 1
      which <- foldid == i
      for(j in seq(nlams[i]))
      {
        #cvraw[i,j]=auc.mat(y[which,],predmat[which,j],weights[which])
        cvraw[i,j] <- auc.mat(y[which,],predmat[which,j], rep(1, sum(which)))
      }
    }
    N <- apply(good, 2, sum)
    #weights=tapply(weights,foldid,sum)
  } else
  {
    ##extract weights and normalize to sum to 1
    #ywt=apply(y,1,sum)
    #y=y/ywt
    #weights=weights*ywt

    N     <- nrow(y) - apply(is.na(predmat), 2, sum)
    cvraw <- switch(type.measure,
                    "mse"     = (y[,1] - (1 - predmat)) ^ 2 +(y[,2] - predmat) ^ 2,
                    "mae"     = abs(y[,1] - (1 - predmat)) + abs(y[,2] - predmat),
                    "deviance"= {
                        predmat=pmin(pmax(predmat,prob_min),prob_max)
                        lp=y[,1]*log(1-predmat)+y[,2]*log(predmat)
                        ly=log(y)
                        ly[y==0]=0
                        ly=drop((y*ly)%*%c(1,1))
                        2*(ly-lp)
                    },
                    "class"=y[,1]*(predmat>.5) +y[,2]*(predmat<=.5)
    )
    if(grouped)
    {
      cvob    <- cvcompute(cvraw,rep(1, nrow(y)), foldid, nlams)
      cvraw   <- cvob$cvraw
      weights <- cvob$weights
      N       <- cvob$N
    }
  }
  #cvm=apply(cvraw,2,weighted.mean,w=weights,na.rm=TRUE)
  cvm  <- apply(cvraw, 2, mean, w = weights, na.rm = TRUE)
  #cvsd=sqrt(apply(scale(cvraw,cvm,FALSE)^2,2,weighted.mean,w=rep(1, nrow(y)),na.rm=TRUE)/(N-1))
  cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE) ^ 2, 2, mean, na.rm = TRUE) / (N - 1))
  out  <- list(cvm  = cvm,
               cvsd = cvsd,
               name = typenames[type.measure])
  if(keep)out$fit.preval <- predmat
  out

}

cv.venngaussian <- function(outlist,lambda,x,y,groups,foldid,type.measure,grouped,keep=FALSE){
  typenames <- c(deviance = "Mean-Squared Error",
                 mse      = "Mean-Squared Error",
                 mae      = "Mean Absolute Error")
  if(type.measure == "default") type.measure <- "mse"
  if(!match(type.measure,c("mse","mae","deviance"),FALSE))
  {
    warning("Only 'mse', 'deviance' or 'mae'  available for Gaussian models; 'mse' used")
    type.measure <- "mse"
  }
  #if(!is.null(offset))y=y-drop(offset)
  predmat <- matrix(NA,length(y),length(lambda))
  nfolds  <- max(foldid)
  nlams   <- double(nfolds)
  for(i in seq(nfolds))
  {
    which    <- foldid == i
    fitobj   <- outlist[[i]]
    #fitobj$offset=FALSE
    preds    <- predict(fitobj, x[which,,drop=FALSE], group.mat = groups[which,,drop=FALSE], type = "response")
    nlami    <- length(outlist[[i]]$lambda)
    predmat[which,seq(nlami)] <- preds
    nlams[i] <- nlami
  }

  N <- length(y) - apply(is.na(predmat), 2, sum)
  cvraw <- switch(type.measure,
                  "mse"      =(y - predmat) ^ 2,
                  "deviance" =(y - predmat) ^ 2,
                  "mae"      =abs(y - predmat)
  )
  if( (length(y)/nfolds <3) && grouped)
  {
    warning("Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold",call.=FALSE)
    grouped <- FALSE
  }
  if(grouped)
  {
    cvob    <- cvcompute(cvraw, rep(1, length(y)), foldid, nlams)
    cvraw   <- cvob$cvraw
    weights <- cvob$weights
    N       <- cvob$N
  }

  #cvm=apply(cvraw,2,weighted.mean,w=weights,na.rm=TRUE)
  cvm  <- apply(cvraw, 2, mean, na.rm = TRUE)
  #cvsd=sqrt(apply(scale(cvraw,cvm,FALSE)^2,2,weighted.mean,w=weights,na.rm=TRUE)/(N-1))
  cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE) ^ 2, 2, mean, na.rm = TRUE)/(N-1))
  out  <- list(cvm  = cvm,
               cvsd = cvsd,
               name = typenames[type.measure])
  if(keep)out$fit.preval <- predmat
  out
}



