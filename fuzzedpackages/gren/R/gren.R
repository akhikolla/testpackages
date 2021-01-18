# group-regularized elastic net function
  gren <- function(x, y, m=rep(1, nrow(x)), unpenalized=NULL, partitions=NULL, 
                   alpha=0.5, lambda=NULL, intercept=TRUE, monotone=NULL, 
                   psel=TRUE, compare=TRUE, posterior=FALSE, nfolds=nrow(x), 
                   foldid=NULL, trace=TRUE,
                   init=list(lambdag=NULL, mu=NULL, sigma=NULL, 
                             chi=NULL, ci=NULL),
                   control=list(epsilon=0.001, maxit=500, maxit.opt=1000, 
                                maxit.vb=100)) {
  
  # save argument list and call
  argum <- formals(gren)
  fit.call <- match.call()
  
  # change some input for convenience
  if(is.data.frame(x)) {
    names.xr <- colnames(x)
    x <- as.matrix(x)
  } else {
    names.xr <- NULL
  }
  if(is.vector(partitions) & is.atomic(partitions)) {
    partitions <- list(partition1=partitions)
  }
  
  # check input
  if(!is.numeric(x) | !(is.numeric(y) | is.factor(y))) {
    stop("only numerical input data is supported at the moment")
  } else if(!is.null(ncol(y))) {
    if(ncol(y)!=2) {
      stop("y is either a vector or matrix with two columns")
    }
  } else if(ifelse(is.null(ncol(y)), length(y), nrow(y))!=length(m) | 
            nrow(x)!=length(m)) {
    stop("number of observations in y, m, and x not equal")
  } else if(!is.null(unpenalized) & !is.numeric(unpenalized) & 
            !is.data.frame(unpenalized)) {
    stop("only unpenalized data in a matrix or data.frame is supported")
  } else if(ifelse(is.null(ncol(unpenalized)), length(unpenalized), 
                   nrow(unpenalized))!=nrow(x) & !is.null(unpenalized)) {
    stop("number of observations in unpenalized not equal to number of 
         observations in y, x, and m")
  } else if(!is.list(partitions)) {
    stop("partitions should be either a list of partitions or one partition as 
         a numeric vector")
  } else if(any(lapply(partitions, length)!=ncol(x))) {
    stop("all partitions should be vectors of length ncol(x), containing the
         group identifiers of the features")
  } else if(!is.numeric(alpha) | length(alpha)!=1 | (alpha < 0) | (alpha > 1)) {
    stop("alpha should be of length one and a value between 0 and 1")
  } else if(!is.null(lambda) & !(is.vector(lambda) & is.atomic(lambda))) {
    stop("lambda should be either NULL or a numeric vector")
  } else if(!is.logical(intercept)) {
    stop("logical should be either TRUE or FALSE")
  } else if(!is.null(monotone) & 
            !(is.logical(monotone) & length(monotone)==length(partitions)) &
            !(is.list(monotone) & length(monotone)==2 & 
              is.logical(unlist(monotone)) & 
              length(unlist(monotone))==2*length(partitions))) {
    stop("monotone should be a vector of TRUEs and FALSEs of the same length
         as partitions or NULL")
  } else if(!is.null(psel) & !(is.vector(psel) & is.atomic(psel))) {
    stop("psel is either NULL or a vector of non-negative whole numbers")
  } else if(!is.logical(compare)) {
    stop("compare is either TRUE or FALSE")
  } else if(!is.logical(posterior)) {
    stop("posterior is either TRUE or FALSE")
  } else if(!is.numeric(nfolds) | length(nfolds)!=1) {
    stop("nfolds should be a whole number")
  } else if(!is.null(foldid) & !(is.numeric(foldid) & length(foldid==nrow(x)))) {
    stop("foldid is either NULL or a vector of length nrow(x)")
  } else if(!is.null(foldid) & !all(foldid %in% c(1:nrow(x)))) {
    stop("foldid must be a vector of length n, containing only whole numbers
         from 1 to nrow(x)")
  } else if(!is.logical(trace)) {
    stop("trace is either TRUE or FALSE")
  } else if(!is.list(init)) {
    stop("init must be a list")
  } else if(!is.null(init$lambdag) & !is.list(init$lambdag)) {
    stop("lambdag must be a list")
  } else if(!is.null(init$mu) & !is.numeric(init$mu) & 
            !is.numeric(init$chi) & !is.null(init$ci) & !is.numeric(init$ci)) {
    stop("mu, chi, and ci must be numeric vectors")
  } else if(!is.null(init$sigma) & !is.numeric(init$sigma)) {
    stop("sigma must be a numeric matrix")
  } else if(!is.list(control)) {
    stop("control must be a list")
  } else if(!is.numeric(control$epsilon) | length(control$epsilon)!=1) {
    stop("epsilon must be a non-negative number")
  } else if(!is.numeric(control$maxit) | length(control$maxit)!=1) {
    stop("maxit must be a non-negative whole number")
  } else if(!is.numeric(control$maxit.opt) | length(control$maxit.opt)!=1) {
    stop("maxit.opt must be a non-negative whole number")
  } else if(!is.numeric(control$maxit.vb) | length(control$maxit.vb)!=1) {
    stop("maxit.vb must be a non-negative whole number")
  }
  
  # set some auxiliary variables
  if(is.null(unpenalized)) {
    xu <- matrix(1, nrow=2)
  } else if(is.data.frame(unpenalized)) {
    xu <- model.matrix(as.formula(paste("~", paste(colnames(unpenalized), 
                                                   collapse="+"))), 
                       data=unpenalized)
    if(!intercept) {
      xu <- xu[, -1]
    } 
  } else if(intercept) {
    xu <- as.matrix(cbind(1, unpenalized))
  }
  xr <- x
  xr <- scale(xr)
  if(is.null(unpenalized)) {
    x <- xr
  } else {
    x <- cbind(xu, xr)  
    if(intercept) {
      x <- x[, -1]
    }
  }
  n <- nrow(x)
  r <- ncol(xr)
  u <- ifelse(is.null(unpenalized), 0, ncol(xu) - intercept)
  p <- r + u
  unpenalized <- !is.null(unpenalized)
  if(is.factor(y)) {
    y <- as.numeric(y) - 1
  }
  if(is.null(ncol(y))) {
    ymat <- cbind(m - y, y)
  } else {
    ymat <- y
    y <- y[, 2]
  }
  
  # if no penalty parameter lambda is given we estimate it by cross-validation
  if(is.null(lambda)) {
    if(is.null(foldid)) {
      rest <- n %% nfolds
      foldsize <- c(rep(n %/% nfolds + as.numeric(rest!=0), times=rest),
                    rep(n %/% nfolds, times=nfolds - rest))
      foldid <- sample(rep(1:nfolds, times=foldsize))
    }
    
    if(trace) {cat("\r", "Estimating global lambda by cross-validation", 
                   sep="")}
    srt <- proc.time()[3]
    cv.fit <- cv.glmnet(x, y, family="binomial", alpha=alpha, standardize=FALSE,
                        intercept=intercept, foldid=foldid, grouped=FALSE,
                        penalty.factor=c(rep(0, u), rep(1, r)))
    lambda <- cv.fit$lambda.min
    cv.time <- proc.time()[3] - srt
    
    if(trace) {cat("\n", "Estimated global lambda is ", round(lambda, 2), 
                   " in ", round(cv.time, 2), " seconds", sep="")}
  }
  lambda1 <- lambda*alpha*n*2
  lambda2 <- lambda*(1 - alpha)*n
  
  # auxiliary variables for the partitions
  if(is.null(names(partitions))) {
    names(partitions) <- paste("partition", 1:length(partitions), sep="")
  }
  partnames <- names(partitions)
  nparts <- length(partitions)
  if(is.null(monotone)) {
    monotone <- list(monotone=rep(FALSE, nparts), decreasing=rep(FALSE, nparts))
  } else if(is.logical(monotone)) {
    monotone <- list(monotone=monotone, decreasing=rep(FALSE, nparts))
  } else if(is.list(monotone) & length(monotone)==1) {
    monotone$decreasing <- rep(FALSE, nparts)
  }
  names(monotone) <- c("monotone", "decreasing")
  numparts <- lapply(partitions, is.numeric)
  groupnames <- lapply(partitions, function(part) {
    if(is.factor(part)) {
      levels(part)
    } else {
      as.character(sort(unique(part)))
    }})
  partitions <- lapply(partitions, as.numeric)
  sizes <- lapply(partitions, function(part) {rle(sort(part))$lengths})
  G <- lapply(partitions, function(part) {length(unique(part))})
  
  # objects used in optimisation of penalty parameters
  intsec <- do.call("paste", c(partitions, sep=" "))
  uintsec <- unique(intsec)
  intsizes <- sapply(uintsec, function(int) {sum(intsec==int)})
  partsmat <- unique(do.call("cbind", partitions))
  partsind <- rep(1:nparts, times=unlist(G))
  
  # starting values for lambdag, lagrange multiplier s
  if(is.null(control$lambdag.start)) {
    lambdag.start <- lapply(G, function(gpart) {rep(1, gpart)})
  }
  lambdagnew <- lambdag <- lambdagold <- lambdagseq <- lambdag.start
  lambdamultvec <- lambdamultvecold <- apply(sapply(1:nparts, function(part) {
    lambdagnew[[part]][partitions[[part]]]}), 1, prod)
  s <- 0
  
  # fixed parameters (dont change with VB iterations)
  kappa <- y - m/2
  phi <- 0.25*lambda1^2/lambda2
  
  # starting values for the model parameters
  if(is.null(init$mu) | is.null(init$sigma) | is.null(init$chi) | 
     is.null(init$ci)) {
    fit.start <- glmnet(x=x, y=ymat, family="binomial", alpha=0, lambda=lambda, 
                        standardize=FALSE, intercept=intercept, 
                        penalty.factor=c(rep(0, u), lambdamultvecold))
    pred.start <- as.numeric(predict(fit.start, newx=x, type="response"))
    startparam <- est_param(xr, xu, kappa, m, n, p, pred.start, phi, 
                            rep(phi, r), lambda2, lambdamultvecold, 
                            lambdamultvecold, intercept, unpenalized, FALSE, 
                            FALSE, TRUE)
    dsigmaold <- as.numeric(startparam$dsigma)
    muold <- as.numeric(startparam$mu)
    ciold <- as.numeric(startparam$ci)
    chiold <- as.numeric(startparam$chi)
  } else {
    dsigmaold <- as.numeric(diag(init$sigma))
    muold <- as.numeric(init$mu)
    ciold <- as.numeric(init$ci)
    chiold <- as.numeric(init$chi)
  }
  
  # calculating the sums needed in optimisation routine
  sum1 <- sapply(uintsec, function(int) {
    ind <- which(intsec==int) + intercept;
    sum((dsigmaold[ind + u] + muold[ind + u]^2)*
          (1 + sqrt(phi/chiold[ind - intercept])))})
  
  # outer loop of algorithm:
  opt.iter <- vector(mode="list", length=0)
  vb.iter <- numeric(0)
  opt.conv <- logical(0)
  vb.conv <- logical(0)
  conv <- FALSE
  iter1 <- 0
  if(trace) {cat("\n", "Estimating penalty multipliers by empirical Bayes", "\n", sep="")}
  srt <- proc.time()[3]
  while(!conv & (iter1 < control$maxit)) {
    
    iter1 <- iter1 + 1
    
    # estimating new lambdag
    opt <- optim(par=c(s, log(unlist(lambdag, use.names=FALSE))), 
                 fn=fopt_groups, lambda2=lambda2, nparts=nparts, 
                 partsind=partsind, partsmat=partsmat, sizes=unlist(sizes), G=G,
                 sum1=sum1, method="BFGS", 
                 control=list(maxit=control$maxit.opt))
    opt.conv <- c(opt.conv, opt$convergence==0)
    opt.iter <- c(opt.iter, opt$counts)
    
    # assigning new penalty multipliers
    lambdagnew <- split(exp(opt$par[-1]), factor(rep(partnames, unlist(G)), 
                                                 levels=partnames))
    
    # estimating optional monotone penalties
    lambdagnew <- sapply(1:nparts, function(part) {
      if(monotone$monotone[part]) {
        return(as.numeric(pava(lambdagnew[[part]], sizes[[part]], 
                               monotone$decreasing[part])))
      } else {
        return(lambdagnew[[part]])
      }}, simplify=FALSE)
    if(any(monotone$monotone)) {
      lambdagnew <- sapply(1:nparts, function(part) {
        exp(log(lambdagnew[[part]]) - sum(sapply(1:nparts, function(part) {
          sum(sizes[[part]]*log(lambdagnew[[part]]))}))/(p*nparts))}, 
        simplify=FALSE)
    }
    s <- opt$par[1]
    
    # keeping track of penalty multipliers
    lambdagseq <- sapply(1:nparts, function(part) {
      cbind(lambdagseq[[part]], lambdagnew[[part]])}, simplify=FALSE)
    names(lambdagseq) <- partnames
    lambdamultvec <- apply(sapply(1:nparts, function(part) {
      lambdagnew[[part]][partitions[[part]]]}), 1, prod)
    
    # inner loop of algorithm:
    conv2 <- FALSE
    iter2 <- 0
    while(!conv2 & (iter2 < control$maxit.vb)) {
      iter2 <- iter2 + 1
      
      # estimating new model parameters
      newparam <- est_param(xr, xu, kappa, m, n, p, ciold, phi, chiold, 
                            lambda2, lambdamultvec, lambdamultvecold,
                            intercept, unpenalized, FALSE, FALSE, FALSE)
      
      dsigma <- as.numeric(newparam$dsigma)
      mu <- as.numeric(newparam$mu)
      ci <- as.numeric(newparam$ci)
      chi <- as.numeric(newparam$chi)
      
      # checking convergence of inner loop
      conv2 <- max(c(abs((mu - muold)/ifelse(muold==0, muold + 0.00001, muold)), 
                     abs((dsigma - dsigmaold)/
                           ifelse(dsigmaold==0, dsigmaold + 0.00001, 
                                  dsigmaold)))) < control$epsilon
      
      # updating vb parameters
      muold <- mu
      dsigmaold <- dsigma
      ciold <- ci
      chiold <- chi
      
    }
    
    # keeping track of convergence and number of iteration
    vb.conv <- c(vb.conv, conv2)
    vb.iter <- c(vb.iter, iter2)
    
    # sum is needed in optimisation routine
    sum1 <- sapply(uintsec, function(int) {
      ind <- which(intsec==int) + intercept;
      sum((dsigmaold[ind + u] + muold[ind + u]^2)*
            (1 + sqrt(phi/chiold[ind - intercept])))})
    
    # checking convergence of outer loop:
    conv <- max(abs((unlist(lambdagnew) - unlist(lambdag))/
                      unlist(lambdag))) < control$epsilon
    
    # updating lambdag for new iteration
    lambdagold <- lambdag
    lambdag <- lambdagnew
    lambdamultvecold <- lambdamultvec
    
    # printing progress
    if(trace) {
      cat("\r", "Estimated penalty multipliers for ", 
          paste(partnames, sapply(c(1:nparts), function(part) {
            paste(paste(round(lambdag[[part]], 2), " (", groupnames[[part]], ")", 
                        sep=""), collapse=", ")}), sep=": ", collapse=" and "), 
          "      ", sep="")
    }
  }
  
  # printing estimation time and final estimates
  eb.time <- proc.time()[3] - srt
  if(trace) {cat("\r", "Penalty multipliers estimated at ", 
                 paste(partnames, lapply(lambdag, function(part) {
                   paste(round(part, 2), collapse=", ")}), sep=": ", 
                   collapse=" and "), " in ", round(eb.time, 2), " seconds ", 
                 sep="")
  }
  
  # if the full vb posterior of beta is desired, estimated it here
  if(posterior) {
    newparam <- est_param(xr, xu, kappa, m, n, p, ci, phi, chi, lambda2, 
                          lambdamultvec, lambdamultvecold, intercept, 
                          !is.null(unpenalized), TRUE, FALSE, FALSE)
    dsigma <- newparam$sigma
    mu <- as.numeric(newparam$mu)
    ci <- as.numeric(newparam$ci)
    chi <- as.numeric(newparam$chi)
  }
  
  # objective function to minimise to find closest model size
  fsel <- function(lambda, maxselec, alpha, groupreg) {
    if(groupreg) {
      fit.sel <- glmnet(x, y, family="binomial", alpha=alpha, lambda=lambda,
                        standardize=FALSE, intercept=intercept,
                        penalty.factor=c(rep(0, u), lambdamultvec))
    } else {
      fit.sel <- glmnet(x, y, family="binomial", alpha=alpha, lambda=lambda,
                        standardize=FALSE, intercept=intercept,
                        penalty.factor=c(rep(0, u), rep(1, r)))
    }
    return(fit.sel$df - u - maxselec)
  }
  
  # variable selection and frequentist elastic net model estimation
  # if no specific model size is requested, let glmnet decide lambdas,
  # otherwise we try to find the closest model sizes
  if(all(psel==TRUE) | is.numeric(psel)) {
    fit.final <- glmnet(x, y=ymat, family="binomial", alpha=alpha, 
                        standardize=FALSE, intercept=intercept, 
                        penalty.factor=c(rep(0, u), lambdamultvec))
    if(is.numeric(psel)) {
      lambdamax <- max(fit.final$lambda)
      lambdamin <- 1e-10
    
      # calculate lambda sequence for desired model sizes
      sel.out <- sapply(psel, function(parsel) {
        lambda.sel <- uniroot(fsel, interval=c(lambdamin, lambdamax), 
                              maxiter=50, maxselec=parsel, alpha=alpha, 
                              groupreg=TRUE)$root
        return(lambda.sel)}, simplify=TRUE)
    
      # estimate final models using lambda sequence
      fit.final <- glmnet(x, y=ymat, family="binomial", alpha=alpha, 
                          lambda=sel.out, standardize=FALSE, intercept=intercept, 
                          penalty.factor=c(rep(0, u), lambdamultvec))
    }
  }
  
  # fit the regular elastic net model if compare is true
  if(compare) {
    fit.regul <- glmnet(x, y, family="binomial", alpha=alpha,
                        standardize=FALSE, intercept=intercept,
                        penalty.factor=c(rep(0, u), rep(1, r)))
    
    # if we want a specif number of features
    if(is.numeric(psel)) {
      lambdamax <- max(fit.regul$lambda)
      lambdamin <- 1e-10
      
      # calculate lambda sequence for desired model sizes
      sel.out <- sapply(psel, function(parsel) {
        lambda.sel <- uniroot(fsel, interval=c(lambdamin, lambdamax), 
                              maxiter=50, maxselec=parsel, alpha=alpha, 
                              groupreg=FALSE)$root
        return(lambda.sel)}, simplify=TRUE)
      fit.regul <- glmnet(x, y=ymat, family="binomial", alpha=alpha, 
                          lambda=sel.out, standardize=FALSE, 
                          intercept=intercept, 
                          penalty.factor=c(rep(0, u), rep(1, r)))
    }
  }
  
  # if no frequentist models estimated, set to NULL
  if(all(psel==FALSE)) {
    groupreg <- NULL
  } else {
    groupreg <- fit.final
  }
  # if no regular elastic net estimated, set to NULL
  if(compare) {
    regular <- fit.regul
  } else {
    regular <- NULL
  }
  freq <- list(groupreg=groupreg, regular=regular)
  iter <- list(lambda.iter=iter1, opt.iter=opt.iter, vb.iter=vb.iter)
  conv <- list(lambda.conv=conv, opt.conv=opt.conv, vb.conv=vb.conv)
  vb.post <- list(mu=mu, sigma=dsigma, ci=ci, chi=chi)
  lambdag.est <- sapply(lambdagseq, function(lg) {
    lg[, iter$lambda.iter + 1]}, simplify=FALSE)
  names(lambdag.est) <- partnames
  for(l in 1:nparts) {
    names(lambdag.est[[l]]) <- groupnames[[l]]
    rownames(lambdagseq[[l]]) <- groupnames[[l]]
  }
  
  # output list
  out <- list(call=fit.call, alpha=alpha, lambda=lambda, lambdag.seq=lambdagseq, 
              lambdag=lambdag.est, vb.post=vb.post, freq.model=freq, iter=iter, 
              conv=conv, args=argum)
  class(out) <- "gren"
  return(out)
  
}