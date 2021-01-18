setClass(
  Class = 'tune',
  slots = c(
    est = 'matrix',
    weight = 'vector'
  )
)


setClass(
  Class = 'cv',
  slots = c(
    input = 'matrix',
    weight = 'vector',
    lambda1 = 'vector',
    lambda2 = 'vector',
    lambda1_max = 'numeric',
    lambda1_1se = 'numeric',
    lambda2_max = 'numeric',
    lambda2_1se = 'numeric',
    cvm = 'matrix',
    cvse = 'matrix',
    cvn = 'matrix',
    cvmax = 'numeric',
    cv1se = 'numeric',
    n_max = 'numeric',
    n_1se = 'numeric',
    intercept_max = 'numeric',
    intercept_1se = 'numeric',
    coef_max = 'vector',
    coef_1se = 'vector', 
    fam = 'character',
    measure = 'character'
  )
)




##' @title Get optimal cut points for binary or right censored phenotype
##' @description Obtain optimal cut point based on Youden index for binary phenotype and log rank test for right censored phenotype.
##' @param pre predicted value
##' @param act actual values (class for binary phenotype and Surv object for right censored phenotype)
##' @param fam the family of the outcome, can be "Gaussian", "Logistic" or "Cox"
##' @return optimal cut point
##' @importFrom OptimalCutpoints optimal.cutpoints
##' @importFrom maxstat maxstat.test
##' @examples x <- rnorm(100)
##' y <- as.numeric(x + rlogis(100) > 0)
##' getcut(x, y)
##' @export
getcut <- function(pre, act, fam = 'Logistic') {
  if(fam == 'Logistic') {
    summary(optimal.cutpoints(X = 'pred', status = 'y', tag.healthy = 0, methods = "Youden", data = data.frame(pred = pre, y = act)))$p.table$Global$Youden[[1]][1]
  } else {
    maxstat.test(act ~ pre, data.frame(pre = pre, act = act), smethod = 'LogRank')$estimate
  }
}

##' @importFrom stats dnorm qnorm sd
biserial <- function(x, y) {
  tab <- as.vector(table(y))/length(y)
  if (length(tab) < 2) {
    r <- NA
  } else {
    zp <- dnorm(qnorm(tab[2]))
    hi <- mean(x[y == 1])
    lo <- mean(x[y == 0])
    r <- (hi - lo)*(prod(tab))/(zp * sd(x))
    if (!is.na(r)) {
      if (r > 1) {
        r <- 1
      } else if (r < -1) {
        r <- -1
      }
    }
  }
  return(r)
}


##' @title Evaluate prediction
##' @description Evaluate goodness of prediction.
##' @param y_pre predicted value
##' @param y actual values (class for binary phenotype and Surv object for right censored phenotype)
##' @param cutpoint cutpoints for binary phenotype, default to be 0.5
##' @param fam family of the phenotype, can be "continous", "binary", or "Cox"
##' @return goodness of prediction
##' @importFrom stats cor
##' @importFrom pROC auc
##' @importFrom survival survConcordance
##' @examples x <- rnorm(100)
##' y <- rnorm(100)
##' evaluate(x, y)
##' @export
evaluate <- function(y_pre, y, cutpoint = .5, fam = 'Gaussian') {
  if(fam == 'Gaussian') {
    res <- c('MAE' = mean(abs(y_pre - y)), 'MSE' = mean((y_pre - y) ^ 2), 'Pearson' = cor(y_pre, y), 'Spearman' = cor(y_pre, y, method = 'spearman'))
  } else if(fam == 'Logistic') {
    AUC <- auc(y, y_pre)
    y_cat <- as.numeric(y_pre >= cutpoint)
    tp <- sum(y_cat == 1 & y == 1)
    tn <- sum(y_cat == 0 & y == 0)
    fn <- sum(y_cat == 0 & y == 1)
    fp <- sum(y_cat == 1 & y == 0)
    pos <- tp + fn
    neg <- tn + fp
    tol <- tp + tn + fn + fp
    Sn <- tp / pos
    Sp <- tn / neg
    J <- Sn + Sp - 1
    ACC <- (tp + tn) / tol
    MCC <- (1 - (fn / pos +fp / neg)) / sqrt((1 + (fp - fn) / pos)*(1 + (fn - fp) / neg))
    bis <- biserial(y_pre, y)
    res <- c('AUC' = AUC, 'ACC' = ACC, 'MCC' = MCC, 'J' = J, 'Sn' = Sn, 'Sp' = Sp, 'Biserial' = bis)
  } else {
    res <- survConcordance(y ~ y_pre)$concordance
  }
  return(res)
}


##' @title Prediction visualization
##' @description Sample plots for prediction evaluation (scatter plot for Gaussian, ROC curve for logistic, and Kaplan Meier curve for Cox)
##' @param y_pre predicted value
##' @param y_test actual value
##' @param fam type of predicted outcome, can be "Gaussian" (default), "Logisitc", and "Cox"
##' @param mod fitted glmagarph model, must be available for Cox if cutpoint not provided
##' @param y_train the training outcome to obtain
##' @param cutp cutpoint for Cox model
##' @importFrom ggplot2 qplot annotate theme_bw xlab ylab aes_
##' @importFrom plotROC geom_roc style_roc
##' @importFrom survminer ggsurvplot surv_fit
##' @importFrom survival Surv
##' @examples x <- rnorm(100)
##' y <- x + rnorm(100)
##' evaluate_plot(x, y)
##' @return plots
##' @export
evaluate_plot <- function(y_pre, y_test, fam = 'Gaussian', mod, y_train, cutp) {
  if (fam == 'Gaussian') {
    qplot(y_test, y_pre) + xlab('actual value') + ylab('predicted value') + theme_bw()
  } else if (fam == 'Logistic') {
    aucp <- auc(y_test ~ y_pre)
    datroc <- data.frame(test = y_test, pre = y_pre)
    ggplot(datroc, aes_(d =~ test, m =~ pre)) + geom_roc(n.cuts = 0) + style_roc() +  annotate("text", x = .5, y = .5, label = paste("AUC", round(aucp, 3), sep = ' = '))
  } else {
    if (missing(cutp)) {
      if (missing(mod) | missing(y_train)) {
        stop('Please input the model and training label!')
      }
      datsurv <- data.frame(train = y_train, pre = predict(mod, type = 'response'))
      cutp <- maxstat.test(train ~ pre, smethod = 'LogRank', datsurv)$estimate
    }
    test1 <- test2 <- NULL
    datsurvtest <- data.frame(test1 = y_test[, 1], test2 = y_test[, 2], group = as.numeric(y_pre < cutp))
    ggsurvplot(surv_fit(Surv(test1, test2) ~ group, data = datsurvtest), data = datsurvtest, risk.table = T, conf.int = T, legend.labs = c('high risk', 'low risk'), pval = T, ggtheme = theme_bw())
  }
}


##' @title Standardized Laplacian matrix
##' @description Obtain standardized Laplacian matrix given adjacency matrix
##' @param A adjacency matrix
##' @return Laplacian matrix
##' @examples a <- matrix(0, 2, 2)
##' la <- laps(a)
##' @export
laps <- function(A){
  dd <- rowSums(A)
  ddD <- dd ^ -.5
  ddD[is.infinite(ddD)] <- 0
  D <- diag(ddD)
  L0 <- diag(nrow(A)) - D %*% A %*% D
  diag(L0) <- 1
  return(L0)
}


##' @title Estimate standardized Laplacian matrix
##' @description Estimate standardized Laplacian matrix given data using gene co-expression network method
##' @param x data
##' @param sparse estimate a sparse network or not, default to be T, but may be slow
##' @return standardized laplacian matrix
##' @importFrom huge huge huge.select
##' @importFrom stats cor
##' @references Ucar D, Neuhaus I, Ross-MacDonald P, Tilford C, Parthasarathy S, et al. (2007) Construction of a reference gene association network from multiple profiling data: application to data analysis. Bioinformatics 23: 2716-2724.
##' @references Meinshausen, N., & B{\"u}hlmann, P. (2006). High-dimensional graphs and variable selection with the lasso. The annals of statistics, 1436-1462.
##' @examples xx <- matrix(rnorm(12), 3, 4)
##' ss <- getS(xx, FALSE)
##' @export
getS <- function(x, sparse = T) {
  x <- scale(as.matrix(x))
  p <- ncol(x)
  IXI <- apply(-abs(cor(x, use = 'pairwise.complete.obs')), 2, rank)
  w <- 1 / (IXI * t(IXI))
  diag(w) <- 0
  a0 <- w / sqrt(tcrossprod(rowSums(w), colSums(w)))
  L <- diag(p) - a0
  if (sparse) {
    out <- as.matrix(huge.select(huge(x, method = 'mb', verbose = F), verbose = F)$refit)
    diag(out) <- 1
    L <- L*out
  }
  if(min(eigen(L)$va) < .001) {
    L <- diag(p) - .95*a0
  }
  return(L)
}



##' @title tune two network
##' @description Tune two network for better prediction.
##' @param y outcome
##' @param x predictors matrix
##' @param L1 Laplacian matrix for the first network
##' @param L2 Laplacian matrix for the second network
##' @param adaptl2 whether to adapt the sign for quadratic penalty, default to be TRUE
##' @param nfolds number of folds used in cross validation, default to be five
##' @param cvwhich fold assignment, start from zero, if missing do random cross validation
##' @param foldseed the random seed for cross validation design
##' @param stratify whether to do stratified cross validation for Logistic or Cox model, default to be TRUE
##' @param lam0 The tuning parameters for quadratic penalty. If not defined, tuned by default
##' @param bets The candidate weight for the first network, must be between 0 and 1, default to be 0, 0.1,..., 1
##' @param intercept whether to include intercept. Ignore for Cox regression
##' @param standardize whether to standardize predictors
##' @param fam family for the outcome, can be "Gaussian", "Logistic", and "Cox"
##' @param type1se whether to use one standard error or maximum rule, default to be one standard error rule
##' @param measdev Whether to use deviance to tune, default to be deviance. If not, use mean absolue error, area under ROC curve, or concordance index for Gaussian, Logistic, and Cox
##' @param maxiter maximum number of iterations, default to be 500
##' @param cri stoppint criterion, default to be 0.001
##' @param parallel whether to do parallel computing at each fold
##' @return \item{est}{estimated mixed Laplacian matrix}
##' @return \item{weight}{weights for the two Laplacian matrix}
##' @useDynLib glmaag, .registration = TRUE
##' @importFrom Rcpp evalCpp
##' @importFrom stats coef
##' @importFrom foreach foreach %dopar%
##' @examples data(sampledata)
##' data(L0)
##' data(L1)
##' y <- sampledata$Y_Gau
##' x <- sampledata[, -(1:3)]
##' Ltune <- tune_network(y, x, L0, L1, adaptl2 = FALSE)
##' weight <- Ltune@weight
##' Lest <- Ltune@est
##' @export
tune_network <- function(y, x, L1, L2, adaptl2 = T, nfolds = 5, cvwhich, foldseed, stratify = T, lam0, bets, intercept = T, standardize = T, fam = 'Gaussian', type1se = T, measdev = T, maxiter = 10000, cri = .001, parallel = F) {
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  if (fam == 'Gaussian') {
    y <- as.vector(y)
    if(intercept) {
      ycen <- scale(y, scale = F)
    } else {
      ycen <- y
    }
    if (missing(cvwhich)) {
      if (!missing(foldseed)) {
        set.seed(foldseed)
      }
      cvwhich <- sample(rep(0:(nfolds - 1), length.out = n))
    }
    
  } else if(fam == 'Logistic') {
    if (sum(y == 0) < nfolds | sum(y) < nfolds) {
      stop('Too unbalance')
    }
    y <- as.vector(y)
    
    if (missing(cvwhich)) {
      if (!missing(foldseed)) {
        set.seed(foldseed)
      }
      if (stratify) {
        pos <- which(y == 1)
        neg <- which(y == 0)
        poswhich <- sample(rep_len(0:(nfolds - 1), length(pos)))
        negwhich <- sample(rep_len(0:(nfolds - 1), length(neg)))
        cvwhich <- c()
        cvwhich[pos] <- poswhich
        cvwhich[neg] <- negwhich
      } else {
        cvwhich <- sample(rep(0:(nfolds - 1), length.out = n))
      }
    }
  } else {
    if(nnzero(y[, 2]) < nfolds) {
      stop('Too many censoring!')
    }
    if (intercept) {
      message('No intercept for Cox')
      intercept <- F
    }
    yord <- order(y)
    y <- y[yord]
    x <- x[yord, ]
    if (which(y[, 2] == 1)[1] > 1) {
      y <- y[-(1:(which(y[, 2] == 1)[1] - 1))]
      x <- x[-(1:(which(y[, 2] == 1)[1] - 1)), ]
    }
    tie <- any(duplicated(y[, 1]))
    if(tie) {
      tie1 <- which(y[, 1] %in% y[, 1][duplicated(y[, 1])])
      tie2 <- rank(y[, 1], ties.method = 'min')[tie1] - 1
      tie3 <- rank(y[, 1], ties.method = 'max')[tie1] - 1
      tie1 <- tie1 - 1
      ntie <- length(tie1)
    } else {
      tie1 <- tie2 <- tie3 <- 0
      ntie <- 0
    }
    if (missing(cvwhich)) {
      if (!missing(foldseed)) {
        set.seed(foldseed)
      }
      obs <- which(y[, 2] == 1)
      cen <- which(y[, 2] == 0)
      if (length(cen) < nfolds | !stratify) {
        cvwhich <- sample(rep(0:(nfolds - 1), length.out = n))
      } else {
        obswhich <- sample(rep_len(0:(nfolds - 1), length(obs)))
        cenwhich <- sample(rep_len(0:(nfolds - 1), length(cen)))
        cvwhich <- c()
        cvwhich[obs] <- obswhich
        cvwhich[cen] <- cenwhich
      }
    }
    for (i in 0:(nfolds - 1)) {
      if (which(cumsum(y[cvwhich == i, 2]) == 1)[1] > 1) {
        cvwhich[which(cvwhich == i)[1:(which(cumsum(y[cvwhich == i, 2]) == 1)[1] - 1)]] <- nfolds + 1
      }
    }
  }
  
  if (standardize & intercept) {
    meanx <- colMeans(x)
    sdx <- apply(x,2,sd)
    xstd <- scale(x)
  } else if(standardize & !intercept) {
    meanx <- rep(0, p)
    sdx <- apply(x, 2, sd)
    xstd <- sweep(x, 2, sdx, '/')
  } else {
    meanx <- rep(0, p)
    sdx <- rep(1, p)
    xstd <- sweep(x, 2, sdx, '/')
  }
  
  if (missing(bets)) {
    bets <- seq(0, 1, .1)
  }
  names(xstd) <- NULL
  if (missing(lam0)) {
    if (n < p) {
      lam0 <- .01*2^(7:0)
    } else {
      lam0 <- c(.01*2^(7:0), 0)
    }
  }
  
  if(any(eigen(L1)$valu < .001) | any(eigen(L2)$valu < .001)){
    dl <- rep(1, p)
    Lad <- matrix(0, p, p)
    L <- diag(p)
    stop('Networks is not positive definite, please change to elastic net!')
  }
  
  
  if (adaptl2) {
    II <- diag(p)
    if (parallel) {
      cvmod0 <- foreach(i = 0:(nfolds - 1), .combine = 'cbind') %dopar% {
        if (fam == 'Gaussian') {
          cvlrnet2_pal(xstd[cvwhich != i, ], xstd[cvwhich == i, ], ycen[cvwhich != i], ycen[cvwhich == i], II, lam0, measdev)
        } else if (fam == 'Logistic') {
          cvloginet2_pal(xstd[cvwhich != i, ], xstd[cvwhich == i, ], y[cvwhich != i], y[cvwhich == i], rep(0, p), II, lam0, intercept, measdev, maxiter, cri)
        } else {
          ytr <- y[cvwhich != i]
          tietr <- any(duplicated(ytr[, 1]))
          if(tietr) {
            tie1tr <- which(ytr[, 1] %in% ytr[, 1][duplicated(ytr[, 1])])
            tie2tr <- rank(ytr[, 1], ties.method = 'min')[tie1tr] - 1
            tie3tr <- rank(ytr[, 1], ties.method = 'max')[tie1tr] - 1
            tie1tr <- tie1tr - 1
            ntietr <- length(tie1tr)
          } else {
            tie1tr <- tie2tr <- tie3tr <- 0
            ntietr <- 0
          }
          cvCoxnet2_pal(rep(0, p), xstd, xstd[cvwhich != i, ], xstd[cvwhich == i, ], y[cvwhich == i, 1], tie, ntie, tie1, tie2, tie3, tietr, ntietr, tie1tr, tie2tr, tie3tr, y[, 2], y[cvwhich != i, 2], y[cvwhich == i, 2], II, lam0, measdev, maxiter, cri)
        }
      }
      
    } else {
      cvmod0 <- switch (fam,
                        Gaussian = cvlrnet2(nfolds, xstd, ycen, II, lam0, cvwhich, measdev),
                        Logistic = cvloginet2(nfolds, xstd, y, rep(0, p), II, lam0, intercept, cvwhich, measdev, maxiter, cri),
                        Cox = cvCoxnet2(nfolds, rep(0, p), xstd, y[, 1], tie, ntie, tie1, tie2, tie3, y[, 2], II, lam0, cvwhich, measdev, maxiter, cri)
      )
    }
    
    cvmean0 <- rowMeans(cvmod0, na.rm = T)
    
    if (type1se) {
      cvse0 <- apply(cvmod0, 1, sd, na.rm = T)/sqrt(nfolds)
      b <- as.vector(switch (fam,
                             Gaussian = lrnet(xstd, ycen, lam0[which(cvmean0 >= max(cvmean0, na.rm = T) - cvse0[which.max(cvmean0)])[1]]*II),
                             Logistic = loginet(0, rep(0, p), xstd, y, lam0[which(cvmean0 >= max(cvmean0, na.rm = T) - cvse0[which.max(cvmean0)])[1]]*II, intercept, maxiter, cri)[-1],
                             Cox = Coxnet(rep(0, p), xstd, tie, tie1, tie2, tie3, y[, 2], lam0[which(cvmean0 >= max(cvmean0, na.rm = T) - cvse0[which.max(cvmean0)])[1]]*II, maxiter, cri)
      ))
    } else {
      b <- as.vector(switch (fam,
                             Gaussian = lrnet(xstd, ycen, lam0[which.max(cvmean0)]*II),
                             Logistic = loginet(0, rep(0, p), xstd, y, lam0[which.max(cvmean0)]*II, intercept, maxiter, cri)[-1],
                             Cox = Coxnet(rep(0, p), xstd, tie, tie1, tie2, tie3, y[, 2], lam0[which.max(cvmean0)]*II, maxiter, cri)
      ))
    }
    
    ss <- diag(sign(b))
    L1 <- ss %*% L1 %*% ss
    L2 <- ss %*% L2 %*% ss
  }
  
  if (parallel) {
    outlist <- foreach(i = 0:(nfolds - 1)) %dopar% {
      if (fam == 'Gaussian') {
        cvlrnet1_pal(xstd[cvwhich != i, ], xstd[cvwhich == i, ], ycen[cvwhich != i], ycen[cvwhich == i], L1, L2, lam0, bets, measdev)
      } else if (fam == 'Logistic') {
        cvloginet1_pal(xstd[cvwhich != i, ], xstd[cvwhich == i, ], y[cvwhich != i], y[cvwhich == i], rep(0, p), L1, L2, lam0, bets, intercept, measdev, maxiter, cri)
      } else {
        ytr <- y[cvwhich != i]
        tietr <- any(duplicated(ytr[, 1]))
        if(tietr) {
          tie1tr <- which(ytr[, 1] %in% ytr[, 1][duplicated(ytr[, 1])])
          tie2tr <- rank(ytr[, 1], ties.method = 'min')[tie1tr] - 1
          tie3tr <- rank(ytr[, 1], ties.method = 'max')[tie1tr] - 1
          tie1tr <- tie1tr - 1
          ntietr <- length(tie1tr)
        } else {
          tie1tr <- tie2tr <- tie3tr <- 0
          ntietr <- 0
        }
        cvCoxnet1_pal(rep(0, p), xstd, xstd[cvwhich != i, ], xstd[cvwhich == i, ], y[cvwhich == i, 1], tie, tie1, tie2, tie3, tietr, tie1tr, tie2tr, tie3tr, y[, 2], y[cvwhich != i, 2], y[cvwhich == i, 2], L1, L2, lam0, bets, measdev, maxiter, cri)
      }
    }
    mod0 <- simplify2array(outlist)
  } else {
    mod0 <- switch (fam,
                    Gaussian = cvlrnet1(nfolds, xstd, ycen, L1, L2, lam0, bets, cvwhich, measdev),
                    Logistic = cvloginet1(nfolds, xstd, y, rep(0, p), L1, L2, lam0, bets, intercept, cvwhich, measdev, maxiter, cri),
                    Cox = cvCoxnet1(nfolds, rep(0, p), xstd, y[, 1], tie, tie1, tie2, tie3, y[, 2], L1, L2, lam0, bets, cvwhich, measdev, maxiter, cri)
    )
  }
  
  
  mod <- apply(mod0, 1:2, mean, na.rm = T)
  maxcv <- bets[which(mod == max(mod, na.rm = T), arr.ind = T)[1, 1]]
  return(new('tune', est = maxcv*L1 + (1 - maxcv)*L2, weight = c(maxcv, 1 - maxcv)))
}



##' @title Cross validation for glmaag
##' @description Do k-fold cross-validation for glmaag
##' @param y outcome
##' @param x predictors matrix
##' @param L Laplacian matrix for the first network
##' @param tune whether to tune the input network with estimated network or identity matrix, ignored if no input network
##' @param est when there is no input network whether to use estimated network or identiy matrix (elastic net) or mixed the network with estimated network or identity matrix, default to be estimated network
##' @param nfolds number of folds used in cross validation, default to be five
##' @param cvwhich fold assignment, start from zero, if missing do random cross validation
##' @param foldseed the random seed for cross validation design
##' @param stratify whether to do stratified cross validation for Logistic or Cox model, default to be TRUE
##' @param gam The power of weights of L1 penalty, default to be ones
##' @param lam1 The tuning parameters for L1 penalty. If not defined, searched by default
##' @param lam2 The tuning parameters for quadratic penalty. If not defined, searched by default
##' @param dfmax maximum number of parameters allowed in the model, default to be p/2
##' @param w0 Weights for L1 penalty. If not defined, estimated via quadratic penalyzed regression
##' @param adaptl1 whether to adapt the L1 penalty, default to be TRUE
##' @param adaptl2 whether to adapt the sign for quadratic penalty, default to be TRUE
##' @param pind indicator vector whether to put L1 penalty on the feature, 1 means penalyzed while 0 means not penalyzed, default to be all ones (all penalyzed)
##' @param intercept whether to include intercept. Ignore for Cox regression
##' @param standardize whether to standardize predictors
##' @param fam family for the outcome, can be "Gaussian", "Logistic", and "Cox"
##' @param type1se whether to use one standard error or maximum rule, default to be one standard error rule
##' @param measdev Whether to use deviance to tune, default to be deviance. If not, use mean absolue error, area under ROC curve, or concordance index for Gaussian, Logistic, and Cox
##' @param maxiter maximum number of iterations, default to be 500
##' @param cri stoppint criterion, default to be 0.001
##' @param parallel whether to do parallel computing at each fold, need to set up parallel first, default to be FALSE
##' @return \item{input}{input predictor matrix}
##' @return \item{inputweight}{estimated weights if mixing network}
##' @return \item{lambda1}{lambda1 path that has been searched}
##' @return \item{lambda1}{lambda1 path that has been searched}
##' @return \item{lambda1_max}{selected lambda1 based on maximum rule}
##' @return \item{lambda2_max}{selected lambda2 based on maximum rule}
##' @return \item{lambda1_1se}{selected lambda1 based on one standard error rule}
##' @return \item{lambda2_1se}{selected lambda2 based on one standard error rule}
##' @return \item{cvm}{the mean cross validation accuracy}
##' @return \item{cv1se}{the standard error of cross validation accuracy}
##' @return \item{cvn}{the mean number of parameter estimated among folds}
##' @return \item{n_max}{number of selected features based on maximum rule}
##' @return \item{n_1se}{number of selected features based on one standard error rule}
##' @return \item{intercept_max}{estimated intercept based on maximum rule}
##' @return \item{intercept_1se}{estimated intercept based on one standard error rule}
##' @return \item{coef_max}{estimated coefficients based on maximum rule}
##' @return \item{coef_1se}{estimated coefficients based on one standard error rule}
##' @return \item{fam}{family of outcome}
##' @return \item{measure}{measure in cross validation}
##' @useDynLib glmaag, .registration = TRUE
##' @importFrom Rcpp evalCpp
##' @importFrom foreach foreach %dopar%
##' @importFrom Matrix nnzero
##' @importFrom stats sd residuals predict glm lm resid coef binomial
##' @importFrom methods new
##' @importFrom survival coxph
##' @examples data(sampledata)
##' data(L0)
##' y <- sampledata$Y_Gau
##' x <- sampledata[, -(1:3)]
##' cvwhich <- sample(rep(0:4, length.out = length(y)))
##' mod <- cv_glmaag(y, x, L0, cvwhich = cvwhich)
##' @export
cv_glmaag <- function(y, x, L, nfolds = 5, cvwhich, foldseed, stratify = T, gam = 1, tune = F, est = T, lam1, lam2, dfmax, w0 , adaptl1 = T, adaptl2 = T, pind, intercept = T, standardize = T, maxiter = 10000, cri = .001, fam = 'Gaussian', measdev = T, type1se = T, parallel = F) {
  x <- as.matrix(x)
  p <- ncol(x)
  n <- nrow(x)
  input <- x
  if (fam == 'Gaussian') {
    if (intercept) {
      ycen <- scale(y, scale = F)
      my <- mean(y)
    } else {
      ycen <- y
      my <- 0
    }
    cvwhich <- sample(rep(0:(nfolds - 1), length.out = n))
  } else if (fam == 'Logistic') {
    if (sum(y == 1) < nfolds | sum(y == 0) < nfolds) {
      stop('Too unbalance!')
    }
    if (missing(cvwhich)) {
      if (!missing(foldseed)) {
        set.seed(foldseed)
      }
      if (stratify) {
        pos <- which(y == 1)
        neg <- which(y == 0)
        poswhich <- sample(rep_len(0:(nfolds - 1), length(pos)))
        negwhich <- sample(rep_len(0:(nfolds - 1), length(neg)))
        cvwhich <- c()
        cvwhich[pos] <- poswhich
        cvwhich[neg] <- negwhich
      } else {
        cvwhich <- sample(rep(0:(nfolds - 1), length.out = n))
      }
    }
    my <- mean(y)
    if(intercept) {
      b0 <- log(my / (1 - my))
    } else {
      b0 <- 0
    }
  } else {
    if(nnzero(y[, 2]) < nfolds) {
      stop('Too many censoring!')
    }
    if (intercept) {
      message('No intercept for Cox')
      intercept <- F
    }
    yord <- order(y)
    y <- y[yord]
    x <- x[yord, ]
    if (which(y[, 2] == 1)[1] > 1) {
      y <- y[-(1:(which(y[, 2] == 1)[1] - 1))]
      x <- x[-(1:(which(y[, 2] == 1)[1] - 1)), ]
    }
    tie <- any(duplicated(y[, 1]))
    if(tie) {
      tie1 <- which(y[, 1] %in% y[, 1][duplicated(y[, 1])])
      tie2 <- rank(y[, 1], ties.method = 'min')[tie1] - 1
      tie3 <- rank(y[, 1], ties.method = 'max')[tie1] - 1
      tie1 <- tie1 - 1
      ntie <- length(tie1)
    } else {
      tie1 <- tie2 <- tie3 <- 0
      ntie <- 0
    }
    if (missing(cvwhich)) {
      if (!missing(foldseed)) {
        set.seed(foldseed)
      }
      obs <- which(y[, 2] == 1)
      cen <- which(y[, 2] == 0)
      if (length(cen) < nfolds | !stratify) {
        cvwhich <- sample(rep(0:(nfolds - 1), length.out = n))
      } else {
        obswhich <- sample(rep_len(0:(nfolds - 1), length(obs)))
        cenwhich <- sample(rep_len(0:(nfolds - 1), length(cen)))
        cvwhich <- c()
        cvwhich[obs] <- obswhich
        cvwhich[cen] <- cenwhich
      }
    }
    for (i in 0:(nfolds - 1)) {
      if (which(cumsum(y[cvwhich == i, 2]) == 1)[1] > 1) {
        cvwhich[which(cvwhich == i)[1:(which(cumsum(y[cvwhich == i, 2]) == 1)[1] - 1)]] <- nfolds + 1
      }
    }
  }
  
  measure <- ifelse(measdev, 'deviance', switch (fam,
                                                 Gaussian = 'MAE',
                                                 Logistic = 'AUC',
                                                 Cox = 'C'
  ))
  
  if (missing(dfmax)) {
    dfmax <- round(.5*p)
    minlam10 <- F
  } else if (dfmax >= p) {
    minlam10 <- T
  } else {
    minlam10 <- F
  }
  
  if(standardize & intercept) {
    meanx <- colMeans(x)
    sdx <- apply(x,2,sd)
    xstd <- scale(x)
  } else if(standardize & !intercept) {
    meanx <- rep(0, p)
    sdx <- apply(x, 2, sd)
    xstd <- sweep(x, 2, sdx, '/')
  } else {
    meanx <- rep(0, p)
    sdx <- rep(1, p)
    xstd <- sweep(x, 2, sdx, '/')
  }
  
  
  
  names(xstd) <- NULL
  
  
  
  if(missing(pind)) {
    pind <- rep(1, p)
  }
  
  
  if (missing(lam1)) {
    if (fam == 'Gaussian'){
      if (all(pind == 1)) {
        lammax <- max(abs(crossprod(xstd, ycen))) / n
      } else {
        lammax <- max(abs(residuals(lm(ycen ~ xstd[, -which(pind == 0)])) %*% xstd[, -which(pind == 0)])) / n
      }
      
    } else if (fam == 'Logistic'){
      if (all(pind == 1)) {
        lammax <- max(abs(crossprod(xstd, y - my))) / n
      } else {
        lammax <- max(abs(residuals(glm(y ~ xstd[, -which(pind == 0)], family = binomial()), 'response') %*% xstd[, -which(pind == 0)])) / n
      }
    } else {
      if (all(pind == 1)) {
        lammax <- max(abs(resid(coxph(y ~ 1)) %*% xstd)) / n
      } else {
        lammax <- max(abs(resid(coxph(y ~ xstd[, -which(pind == 0)])) %*% xstd)) / n
      }
    }
    if (minlam10) {
      if (n < p) {
        lam1 <- c(exp(seq(log(lammax), log(.01*lammax), len = 100)), 0)
      } else {
        lam1 <- c(exp(seq(log(lammax), log(.0001*lammax), len = 100)), 0)
      }
    } else {
      if (n < p) {
        lam1 <- exp(seq(log(lammax), log(.01*lammax), len = 100))
      } else {
        lam1 <- exp(seq(log(lammax), log(.0001*lammax), len = 100))
      }
    }
    
  }
  if (missing(lam2)) {
    lam2 <- c(.01*2^(7:0), 0)
  }
  
  if (n < p) {
    lam20 <- lam2[lam2 != 0]
  } else {
    lam20 <- lam2
  }
  
  
  
  if(missing(L)) {
    if (tune) {
      message('Cannot tune without an input network!')
    }
    gototune <- F
    if (est) {
      message('Use estimated network')
      L <- getS(xstd)
      dl <- rep(1, p)
      Lad <- L
      diag(Lad) <- 0
    } else {
      message('Use elastic net')
      dl <- rep(1, p)
      Lad <- matrix(0, p, p)
      L <- diag(p)
    }
  } else if(any(eigen(L)$valu < .001)){
    dl <- rep(1, p)
    Lad <- matrix(0, p, p)
    L <- diag(p)
    warning('Networks is not positive definite, change to identity matrix!')
    if (tune) {
      if (est) {
        L2 <- getS(xstd)
        gototune <- T
      } else {
        gototune <- F
      }
    } else {
      gototune <- F
    }
  } else {
    dl <- diag(L)
    Lad <- L
    diag(Lad) <- 0
    gototune <- tune
    if (tune) {
      if (est) {
        L2 <- getS(xstd)
      } else {
        L2 <- diag(p)
      }
    }
  }
  
  tuneweight <- c(1, 0)
  
  if (gototune) {
    modL <- switch (fam,
                    Gaussian = tune_network(ycen, xstd, L, L2, adaptl2 = adaptl2, nfolds = nfolds, intercept = F, standardize = F, fam = fam, type1se = type1se, measdev = measdev, maxiter = maxiter, cri = cri, parallel = parallel),
                    Logistic = tune_network(y, xstd, L, L2, adaptl2 = adaptl2, nfolds = nfolds, intercept = intercept, standardize = F, fam = fam, type1se = type1se, measdev = measdev, maxiter = maxiter, cri = cri, parallel = parallel),
                    Cox = tune_network(y, xstd, L, L2, adaptl2 = adaptl2, nfolds = nfolds, intercept = F, standardize = F, fam = fam, type1se = type1se, measdev = measdev, maxiter = maxiter, cri = cri, parallel = parallel)
    )
    L <- modL@est
    tuneweight <- modL@weight
    dl <- diag(L)
    Lad <- L
    diag(Lad) <- 0
  } else if(adaptl2) {
    II <- diag(p)
    if (parallel) {
      cvmod0 <- foreach(i = 0:(nfolds - 1), .combine = 'cbind') %dopar% {
        if (fam == 'Gaussian') {
          cvlrnet2_pal(xstd[cvwhich != i, ], xstd[cvwhich == i, ], ycen[cvwhich != i], ycen[cvwhich == i], II, lam20, measdev)
        } else if (fam == 'Logistic') {
          cvloginet2_pal(xstd[cvwhich != i, ], xstd[cvwhich == i, ], y[cvwhich != i], y[cvwhich == i], rep(0, p), II, lam20, intercept, measdev, maxiter, cri)
        } else {
          ytr <- y[cvwhich != i]
          tietr <- any(duplicated(ytr[, 1]))
          if(tietr) {
            tie1tr <- which(ytr[, 1] %in% ytr[, 1][duplicated(ytr[, 1])])
            tie2tr <- rank(ytr[, 1], ties.method = 'min')[tie1tr] - 1
            tie3tr <- rank(ytr[, 1], ties.method = 'max')[tie1tr] - 1
            tie1tr <- tie1tr - 1
            ntietr <- length(tie1tr)
          } else {
            tie1tr <- tie2tr <- tie3tr <- 0
            ntietr <- 0
          }
          cvCoxnet2_pal(rep(0, p), xstd, xstd[cvwhich != i, ], xstd[cvwhich == i, ], y[cvwhich == i, 1], tie, ntie, tie1, tie2, tie3, tietr, ntietr, tie1tr, tie2tr, tie3tr, y[, 2], y[cvwhich != i, 2], y[cvwhich == i, 2], II, lam20, measdev, maxiter, cri)
        }
      }
      
    } else {
      cvmod0 <- switch (fam,
                        Gaussian = cvlrnet2(nfolds, xstd, ycen, II, lam20, cvwhich, measdev),
                        Logistic = cvloginet2(nfolds, xstd, y, rep(0, p), II, lam20, intercept, cvwhich, measdev, maxiter, cri),
                        Cox = cvCoxnet2(nfolds, rep(0, p), xstd, y[, 1], tie, ntie, tie1, tie2, tie3, y[, 2], II, lam20, cvwhich, measdev, maxiter, cri)
      )
    }
    
    cvmean0 <- rowMeans(cvmod0, na.rm = T)
    
    if (type1se) {
      cvse0 <- apply(cvmod0, 1, sd, na.rm = T)/sqrt(nfolds)
      b <- as.vector(switch (fam,
                             Gaussian = lrnet(xstd, ycen, lam20[which(cvmean0 >= max(cvmean0, na.rm = T) - cvse0[which.max(cvmean0)])[1]]*II),
                             Logistic = loginet(0, rep(0, p), xstd, y, lam20[which(cvmean0 >= max(cvmean0, na.rm = T) - cvse0[which.max(cvmean0)])[1]]*II, intercept, maxiter, cri)[-1],
                             Cox = Coxnet(rep(0, p), xstd, tie, tie1, tie2, tie3, y[, 2], lam20[which(cvmean0 >= max(cvmean0, na.rm = T) - cvse0[which.max(cvmean0)])[1]]*II, maxiter, cri)
      ))
    } else {
      b <- as.vector(switch (fam,
                             Gaussian = lrnet(xstd, ycen, lam20[which.max(cvmean0)]*II),
                             Logistic = loginet(0, rep(0, p), xstd, y, lam20[which.max(cvmean0)]*II, intercept, maxiter, cri)[-1],
                             Cox = Coxnet(rep(0, p), xstd, tie, tie1, tie2, tie3, y[, 2], lam20[which.max(cvmean0)]*II, maxiter, cri)
      ))
    }
    
    ss <- diag(sign(b))
    Lad <- ss %*% -abs(Lad) %*% ss
    L <- ss %*% L %*% ss
  }
  
  
  
  
  if(minlam10 | adaptl1) {
    if (parallel) {
      cvmod0 <- foreach(i = 0:(nfolds - 1), .combine = 'cbind') %dopar% {
        if (fam == 'Gaussian') {
          cvlrnet2_pal(xstd[cvwhich != i, ], xstd[cvwhich == i, ], ycen[cvwhich != i], ycen[cvwhich == i], L, lam20, measdev)
        } else if (fam == 'Logistic') {
          cvloginet2_pal(xstd[cvwhich != i, ], xstd[cvwhich == i, ], y[cvwhich != i], y[cvwhich == i], rep(0, p), L, lam20, intercept, measdev, maxiter, cri)
        } else {
          ytr <- y[cvwhich != i]
          tietr <- any(duplicated(ytr[, 1]))
          if(tietr) {
            tie1tr <- which(ytr[, 1] %in% ytr[, 1][duplicated(ytr[, 1])])
            tie2tr <- rank(ytr[, 1], ties.method = 'min')[tie1tr] - 1
            tie3tr <- rank(ytr[, 1], ties.method = 'max')[tie1tr] - 1
            tie1tr <- tie1tr - 1
            ntietr <- length(tie1tr)
          } else {
            tie1tr <- tie2tr <- tie3tr <- 0
            ntietr <- 0
          }
          cvCoxnet2_pal(rep(0, p), xstd, xstd[cvwhich != i, ], xstd[cvwhich == i, ], y[cvwhich == i, 1], tie, ntie, tie1, tie2, tie3, tietr, ntietr, tie1tr, tie2tr, tie3tr, y[, 2], y[cvwhich != i, 2], y[cvwhich == i, 2], L, lam20, measdev, maxiter, cri)
        }
      }
      
    } else {
      cvmod0 <- switch (fam,
                        Gaussian = cvlrnet2(nfolds, xstd, ycen, L, lam20, cvwhich, measdev),
                        Logistic = cvloginet2(nfolds, xstd, y, rep(0, p), L, lam20, intercept, cvwhich, measdev, maxiter, cri),
                        Cox = cvCoxnet2(nfolds, rep(0, p), xstd, y[, 1], tie, ntie, tie1, tie2, tie3, y[, 2], L, lam20, cvwhich, measdev, maxiter, cri)
      )
    }
    
    cvmean0 <- rowMeans(cvmod0, na.rm = T)
    cvse0 <- apply(cvmod0, 1, sd, na.rm = T)/sqrt(nfolds)
    
    
    if (type1se) {
      b <- as.vector(switch (fam,
                             Gaussian = lrnet(xstd, ycen, lam20[which(cvmean0 >= max(cvmean0, na.rm = T) - cvse0[which.max(cvmean0)])[1]]*L),
                             Logistic = loginet(0, rep(0, p), xstd, y, lam20[which(cvmean0 >= max(cvmean0, na.rm = T) - cvse0[which.max(cvmean0)])[1]]*L, intercept, maxiter, cri)[-1],
                             Cox = Coxnet(rep(0, p), xstd, tie, tie1, tie2, tie3, y[, 2], lam20[which(cvmean0 >= max(cvmean0, na.rm = T) - cvse0[which.max(cvmean0)])[1]]*L, maxiter, cri)
      ))
    } else {
      b <- as.vector(switch (fam,
                             Gaussian = lrnet(xstd, ycen, lam20[which.max(cvmean0)]*L),
                             Logistic = loginet(0, rep(0, p), xstd, y, lam20[which.max(cvmean0)]*L, intercept, maxiter, cri)[-1],
                             Cox = Coxnet(rep(0, p), xstd, tie, tie1, tie2, tie3, y[, 2], lam20[which.max(cvmean0)]*L, maxiter, cri)
      ))
    }
  }
  
  
  
  if (adaptl1) {
    w <- (abs(b)^(-gam))*pind
  } else {
    w <- rep(1, p)*pind
  }
  
  if (!missing(w0)) {
    w <- w0
  }
  lam1 <- lam1/min(w)
  lam10 <- lam1[lam1 != 0]
  if (parallel) {
    outlist <- foreach(i = 0:(nfolds - 1)) %dopar% {
      if (fam == 'Gaussian') {
        cvlraagg_pal(dfmax, rep(0, p), w, xstd[cvwhich != i, ], xstd[cvwhich == i, ], ycen[cvwhich != i], ycen[cvwhich == i], Lad, dl, lam10, lam2, measdev, maxiter, cri)
      } else if (fam == 'Logistic') {
        cvlogiaagg_pal(dfmax, rep(0, p), xstd[cvwhich != i, ], xstd[cvwhich == i, ], y[cvwhich != i], y[cvwhich == i], Lad, dl, w, lam10, lam2, intercept, measdev, maxiter, cri)
      } else {
        ytr <- y[cvwhich != i]
        tietr <- any(duplicated(ytr[, 1]))
        if(tietr) {
          tie1tr <- which(ytr[, 1] %in% ytr[, 1][duplicated(ytr[, 1])])
          tie2tr <- rank(ytr[, 1], ties.method = 'min')[tie1tr] - 1
          tie3tr <- rank(ytr[, 1], ties.method = 'max')[tie1tr] - 1
          tie1tr <- tie1tr - 1
          ntietr <- length(tie1tr)
        } else {
          tie1tr <- tie2tr <- tie3tr <- 0
          ntietr <- 0
        }
        cvCoxaagg_pal(dfmax, ntie, ntietr, rep(0, p), w, xstd, xstd[cvwhich != i, ], xstd[cvwhich == i, ], y[cvwhich == i, 1], tie, tie1, tie2, tie3, tietr, tie1tr, tie2tr, tie3tr, y[, 2], y[cvwhich != i, 2], y[cvwhich == i, 2], Lad, dl, lam10, lam2, measdev, maxiter, cri)
      }
    }
    cvmod <- list(simplify2array(lapply(outlist, '[[', 1)), simplify2array(lapply(outlist, '[[', 2)))
    
  } else {
    cvmod <- switch (fam,
                     Gaussian = cvlraagg(nfolds, dfmax, rep(0, p), w, xstd, ycen, Lad, dl, lam10, lam2, cvwhich, measdev, maxiter, cri),
                     Logistic = cvlogiaagg(nfolds, dfmax, rep(0, p), xstd, y, Lad, dl, w, lam10, lam2, cvwhich, intercept, measdev, maxiter, cri),
                     Cox = cvCoxaagg(nfolds, dfmax, ntie, rep(0, p), w, xstd, y[, 1], tie, tie1, tie2, tie3, y[, 2], Lad, dl, lam10, lam2, cvwhich, measdev, maxiter, cri)
    )
  }
  
  
  
  
  if (minlam10) {
    if (min(lam20) == min(lam2)) {
      cvmean <- rbind(apply(cvmod[[1]], 1:2, mean, na.rm = T), cvmean0)
      cvse <- rbind(apply(cvmod[[1]], 1:2, sd, na.rm = T) / sqrt(nfolds), cvse0)
      cvn <- rbind(apply(cvmod[[2]], 1:2, mean, na.rm = T), rep(p, length(lam20)))
    } else {
      cvmean <- rbind(apply(cvmod[[1]], 1:2, mean, na.rm = T), c(cvmean0, NA))
      cvse <- rbind(apply(cvmod[[1]], 1:2, sd, na.rm = T) / sqrt(nfolds), c(cvse0, NA))
      cvn <- rbind(apply(cvmod[[2]], 1:2, mean, na.rm = T), rep(p, length(lam20) + 1))
    }
  } else {
    cvmean <- apply(cvmod[[1]], 1:2, mean, na.rm = T)
    cvse <- apply(cvmod[[1]], 1:2, sd, na.rm = T) / sqrt(nfolds)
    cvn <- apply(cvmod[[2]], 1:2, mean, na.rm = T)
    cvmean[cvn > dfmax] <- cvse[cvn > dfmax]  <- NA
  }
  
  cvmax <- max(cvmean, na.rm = T)
  maxloc <- which(cvmean == cvmax, arr.ind = T)[1, ]
  lambda1_max <- lam1[maxloc[1]]
  lambda2_max <- lam2[maxloc[2]]
  
  if (fam == 'Gaussian') {
    if (lam1[maxloc[1]] != 0){
      modmax <- as.vector(lraagg(rep(0, p), lambda1_max, lambda2_max, w, xstd, ycen, Lad, dl, maxiter, cri))
    } else {
      modmax <- as.vector(lrnet(xstd, ycen, lambda2_max*L))
    }
    if (intercept) {
      my <- mean(y)
    } else {
      my <- 0
    }
    b0_max <- my-sum(meanx/sdx*modmax)
    b_max <- modmax/sdx
  } else if (fam == 'Logistic'){
    if (lam1[maxloc[1]] != 0) {
      modmax <- as.vector(logiaagg(b0, rep(0, p), lambda1_max, lambda2_max, w, xstd, y, Lad, dl, intercept, maxiter, cri))
    } else {
      modmax <- as.vector(loginet(b0, rep(0, p), xstd, y, lambda2_max*L, intercept, maxiter, cri))
    }
    bbb0 <- modmax[1]
    bbb <- modmax[-1]
    b0_max <- bbb0 - sum(meanx/sdx*bbb)
    b_max <- bbb/sdx
  } else {
    if (lam1[maxloc[1]] != 0) {
      modmax <- as.vector(Coxaagg(rep(0, p), lambda1_max, lambda2_max, ntie, w, xstd, tie, tie1, tie2, tie3, y[, 2], Lad, dl, maxiter, cri))
      
    } else {
      modmax <- as.vector(Coxnet(rep(0, p), xstd, tie, tie1, tie2, tie3, y[, 2], lambda2_max*L, maxiter, cri))
    }
    b_max <- modmax / sdx
    b0_max <- 0
  }
  
  n_max <- nnzero(b_max)
  
  cv1se <- cvmax - cvse[maxloc[1], maxloc[2]]
  
  id1se <- matrix(which(cvmean >= cv1se, arr.ind = T), ncol = 2)
  
  colnames(id1se) <- c('lambda1','lambda2')
  if (nrow(id1se) == 1) {
    mod1se <- modmax
    lambda1_1se <- lambda1_max
    lambda2_1se <- lambda2_max
    b0_1se <- b0_max
    b_1se <- b_max
  } else {
    id1se <- matrix(id1se[!duplicated(id1se[, 2]), ], ncol = 2)
    id1se <- matrix(id1se[!duplicated(id1se[, 1]), ], ncol = 2)
    n1se <- nrow(id1se)
    if (n1se == 1) {
      if (lam1[id1se[1, 1]] == 0) {
        if (fam == 'Gaussian') {
          mod1se <- as.vector(lrnet(xstd, ycen, lam2[id1se[1, 2]]*L))
          b0_1se <- my - sum(meanx/sdx*mod1se)
          b_1se <- mod1se/sdx
        } else if (fam == 'Logistic') {
          mod1se <- as.vector(loginet(b0, rep(0, p), xstd, y, lam2[id1se[1, 2]]*L, intercept, maxiter, cri))
          bbb0 <- mod1se[1]
          bbb <- mod1se[-1]
          b0_1se <- bbb0-sum(meanx/sdx*bbb)
          b_1se <- bbb/sdx
        } else {
          mod1se <- as.vector(Coxnet(rep(0, p), xstd, tie, tie1, tie2, tie3, y[, 2], lam2[id1se[1, 2]]*L, maxiter, cri))
          b0_1se <- 0
          b_1se <- mod1se/sdx
        }
        lambda1_1se <- 0
        lambda2_1se <- lam2[id1se[1, 2]]
      } else {
        if(fam == 'Gaussian') {
          mod1se <- as.vector(lraagg(rep(0, p), lam1[id1se[1, 1]], lam2[id1se[1, 2]], w, xstd, ycen, Lad, dl, maxiter, cri))
          b0_1se <- my-sum(meanx/sdx*mod1se)
          b_1se <- mod1se/sdx
        } else if(fam == 'Logistic') {
          mod1se <- as.vector(logiaagg(b0, rep(0, p), lam1[id1se[1, 1]], lam2[id1se[1, 2]], w, xstd, y, Lad, dl, intercept, maxiter, cri))
          bbb0 <- mod1se[1]
          bbb <- mod1se[-1]
          b0_1se <- bbb0-sum(meanx/sdx*bbb)
          b_1se <- bbb/sdx
        } else {
          mod1se <- as.vector(Coxaagg(rep(0, p), lam1[id1se[1, 1]], lam2[id1se[1, 2]], ntie, w, xstd, tie, tie1, tie2, tie3, y[, 2], Lad, dl, maxiter, cri))
          b0_1se <- 0
          b_1se <- mod1se/sdx
          
        }
        lambda1_1se <- lam1[id1se[1, 1]]
        lambda2_1se <- lam2[id1se[1, 2]]
      }
    } else {
      if(any(lam1[id1se[, 1]] == 0)) {
        id1se <- matrix(id1se[lam1[id1se[, 1]] != 0, ], ncol = 2)
      }
      n1se <- nrow(id1se)
      if (n1se == 1) {
        if(fam == 'Gaussian') {
          mod1se <- as.vector(lraagg(rep(0, p), lam1[id1se[1, 1]], lam2[id1se[1, 2]], w, xstd, ycen, Lad, dl, maxiter, cri))
          b0_1se <- my-sum(meanx/sdx*mod1se)
          b_1se <- mod1se/sdx
        } else if(fam == 'Logistic') {
          mod1se <- as.vector(logiaagg(b0, rep(0, p), lam1[id1se[1, 1]], lam2[id1se[1, 2]], w, xstd, y, Lad, dl, intercept, maxiter, cri))
          bbb0 <- mod1se[1]
          bbb <- mod1se[-1]
          b0_1se <- bbb0-sum(meanx/sdx*bbb)
          b_1se <- bbb/sdx
        } else {
          mod1se <- as.vector(Coxaagg(rep(0, p), lam1[id1se[1, 1]], lam2[id1se[1, 2]], ntie, w, xstd, tie, tie1, tie2, tie3, y[, 2], Lad, dl, maxiter, cri))
          b0_1se <- 0
          b_1se <- mod1se/sdx
          
        }
        lambda1_1se <- lam1[id1se[1, 1]]
        lambda2_1se <- lam2[id1se[1, 2]]
      } else {
        param <- cbind(lam1[id1se[, 1]], lam2[id1se[, 2]])
        if(fam == 'Gaussian'){
          mod1semod <- findlr1se(rep(0, p), w, xstd, ycen, Lad, dl, maxiter, cri, param)
          mod1se <- mod1semod$b
          b0_1se <- my-sum(meanx/sdx*mod1se)
          b_1se <- mod1se/sdx
        }else if (fam == 'Logistic'){
          mod1semod <- findlogi1se(b0, rep(0, p), w, xstd, y, Lad, dl, intercept, maxiter, cri, param)
          b0 <- mod1semod$b0
          b <- mod1semod$b
          b0_1se <- b0-sum(meanx/sdx*b)
          b_1se <- b/sdx
        } else {
          mod1semod <- findCox1se(n, p, ntie, rep(0, p), w, xstd, tie, tie1, tie2, tie3, y[, 2], Lad, dl,maxiter, cri, param)
          b0_1se <- 0
          b <- as.vector(mod1semod$b)
          b_1se <- b/sdx
        }
        lambda1_1se <- mod1semod$tune[1]
        lambda2_1se <- mod1semod$tune[2]
      }
    }
  }
  n_1se <- nnzero(b_1se)
  rownames(cvmean) <- rownames(cvse) <- rownames(cvn) <- lam1
  colnames(cvmean) <- colnames(cvse) <- colnames(cvn) <- lam2
  
  res <- new('cv', input = input, weight = tuneweight, lambda1 = lam1, lambda2 = lam2, lambda1_max = lambda1_max, lambda1_1se = lambda1_1se,
             lambda2_max = lambda2_max, lambda2_1se = lambda2_1se, cvm = cvmean, cvse = cvse, cvn = cvn, cvmax = cvmax, cv1se = cv1se, n_max = n_max,
             n_1se = n_1se, intercept_max = b0_max, intercept_1se = b0_1se, coef_max = b_max, coef_1se = b_1se, fam = fam, measure = measure)
  structure(res, class = 'cv_glmaag')
}



##' @title Coefficients
##' @description Get the coefficients estimated by the cv_glmaag model
##' @param object the estimated cv_glmaag model
##' @param type1se whether or not used 1 SE error (default to be TRUE)
##' @param ... \dots
##' @return estimated coefficient included intercept (Cox model does not return intercept)
##' @method coef cv_glmaag
##' @examples data(sampledata)
##' data(L0)
##' y <- sampledata$Y_Gau
##' x <- sampledata[, -(1:3)]
##' cvwhich <- sample(rep(0:4, length.out = length(y)))
##' mod <- cv_glmaag(y, x, L0, cvwhich = cvwhich)
##' cc <- coef(mod)
##' @export
coef.cv_glmaag <- function(object, type1se = T, ...) {
  if(object@fam != 'Cox') {
    if(type1se) {
      return(c(object@intercept_1se, object@coef_1se))
    } else {
      return(c(object@intercept_max, object@coef_max))
    }
  } else {
    if(type1se) {
      return(object@coef_1se)
    } else {
      return(object@coef_max)
    }
  }
}

##' @title Predict
##' @description Prediction for cv_glmaag model
##' @param object the estimated cv_glmaag model
##' @param x the new dataset for prediction, if omitted returns the training prediction
##' @param type1se whether or not using the coefficients by one standard error ruld, default to be TRUE
##' @param type can be either "link", or "response", link returns linear predicted score, For Gaussian model this option can be ingnored, for logistic model "response" returns predicted probability, for Cox model "reponse" returns relative risk
##' @param ... \dots
##' @return the predicted value
##' @method predict cv_glmaag
##' @examples data(sampledata)
##' data(L0)
##' y <- sampledata$Y_Gau
##' x <- sampledata[, -(1:3)]
##' cvwhich <- sample(rep(0:4, length.out = length(y)))
##' mod <- cv_glmaag(y, x, L0, cvwhich = cvwhich)
##' pp <- predict(mod)
##' @export
predict.cv_glmaag <- function(object, x, type1se = T, type = 'link', ...) {
  if (missing(x)) {
    x <- object@input
  }
  x <- as.matrix(x)
  if(type1se) {
    ypre <- as.vector(object@intercept_1se + x %*% object@coef_1se)
  } else {
    ypre <- as.vector(object@intercept_max + x %*% object@coef_max)
  }
  if(type == 'link' | object@fam == 'Gaussian') {
    return(ypre)
  } else {
    if (object@fam == 'Logistic') {
      return(1 / (1 + exp(-ypre)))
    } else {
      return(exp(ypre))
    }
  }
}

##' @importFrom ggplot2 ggplot aes_ geom_line geom_ribbon element_blank element_text theme theme_bw ylab ggtitle
plotpath <- function(x, meas, datoplot) {
  lambda2 <- NULL
  ggplot(datoplot[lambda2 == x], aes_(x =~ lambda1, y =~ measure)) + geom_line() + geom_ribbon(aes_(ymin =~ measure - SE, ymax =~ measure + SE), alpha = .1) + ylab(meas) + theme_bw() + ggtitle(substitute(lambda[2] == p, list(p = x))) + xlab(expression(lambda[1])) + theme(plot.title = element_text(hjust = .5), axis.title.y = element_blank())
}



##' @title Cross validation plot
##' @description plot cross validation performance paths
##' @param x the cv_glmaag object
##' @param col_count number of columns in the plots
##' @param SE whether or not plot the standared error curves (when SE = TRUE)
##' @param ... \dots
##' @return plot generated by the model
##' @importFrom data.table as.data.table melt setnames :=
##' @importFrom gridExtra grid.arrange
##' @importFrom ggplot2 ggplot geom_line labs theme theme_bw aes_ element_blank element_text
##' @method plot cv_glmaag
##' @examples data(sampledata)
##' data(L0)
##' y <- sampledata$Y_Gau
##' x <- sampledata[, -(1:3)]
##' cvwhich <- sample(rep(0:4, length.out = length(y)))
##' mod <- cv_glmaag(y, x, L0, cvwhich = cvwhich)
##' gg <- plot(mod, SE = FALSE)
##' @export
plot.cv_glmaag <- function(x, col_count = 3, SE = T, ...) {
  lambda1 <- lambda2 <- NULL
  meanmat <- as.data.table(x@cvm)
  setnames(meanmat, as.character(x@lambda2))
  meanmat[, lambda1 := x@lambda1]
  meandat <- melt(meanmat, 'lambda1')
  setnames(meandat, c('lambda1', 'lambda2', 'measure'))
  if (SE) {
    semat <- as.data.table(x@cvse)
    setnames(semat, as.character(x@lambda2))
    semat[, lambda1 := x@lambda1]
    sedat <- melt(semat, 'lambda1')
    setnames(sedat, c('lambda1', 'lambda2', 'SE'))
    datoplot <- merge(meandat, sedat, c('lambda1', 'lambda2'))
    datoplot[, lambda2 := as.numeric(as.character(lambda2))]
    plotlist <- lapply(x@lambda2, plotpath, meas = x@measure, datoplot = datoplot)
    grid.arrange(grobs = plotlist, ncol = col_count, left = x@measure)
  } else {
    ggplot(meandat, aes_(x =~ lambda1, y =~ measure, col =~ lambda2)) + geom_line() + labs(title = x@measure, color = expression(lambda[2]~''), x = expression(lambda[1])) + theme_bw() + theme(axis.title.y = element_blank(), plot.title = element_text(hjust = .5), legend.position = 'bottom')
  }
}

##' @title the results of the cross validation model
##' @description print fitted information
##' @param x the fitted cv_glmaag object
##' @param ... \dots
##' @method print cv_glmaag
##' @examples data(sampledata)
##' data(L0)
##' y <- sampledata$Y_Gau
##' x <- sampledata[, -(1:3)]
##' cvwhich <- sample(rep(0:4, length.out = length(y)))
##' mod <- cv_glmaag(y, x, L0, cvwhich = cvwhich)
##' print(mod)
##' @export
print.cv_glmaag <- function(x, ...) {
  cat('The maximum solution is: lambda1 = ', x@lambda1_max, ', lambda2 = ', x@lambda2_max, ', ', x@measure, ' = ', x@cvmax, '.\n', sep = '')
  cat('The one standard error solution is: lambda1 = ', x@lambda1_1se, ', lambda2 = ', x@lambda2_1se, ', ', x@measure, ' = ', x@cv1se, '.', sep = '')
}

