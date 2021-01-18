## Transform inpute data into log compositional data.
proc.comp <- function(x) {
  if(any(x <= 0)) stop("Data must be positive")
  if( is.null(dim(x)) ) x <- matrix(x, nrow = 1) else x <- as.matrix(x)

  if( any(abs(rowSums(x) - 1) > 1e-10) ) {
    message("Data is transformed into compositional data by deviding rowSums")
    x <- x / rowSums(x)
  }
  x <- log(x)
  return(x)
}


## Point linear interpolation
point.interp <- function(sseq_obs, sseq_full) {
  ## formula: sfrac*right+(1-sfrac)*left
  if (length(sseq_obs) == 1) {
    nums <- length(sseq_full)
    left <- rep(1, nums)
    right <- left
    sfrac <- rep(1, nums)
  } else {
    sseq_full[sseq_full > max(sseq_obs)] <- max(sseq_obs)
    sseq_full[sseq_full < min(sseq_obs)] <- min(sseq_obs)
    k <- length(sseq_obs)
    # sfrac <- (sseq_full - sseq_obs[1])/(sseq_obs[k] - sseq_obs[1])
    # sseq_obs <- (sseq_obs - sseq_obs[1])/(sseq_obs[k] - sseq_obs[1])
    # coord <- approx(sseq_obs, seq(sseq_obs), sfrac, method = "linear")$y
    coord <- approx(sseq_obs, seq(sseq_obs), sseq_full, method = "linear")$y
    left <- floor(coord)
    right <- ceiling(coord)
    sfrac <- (sseq_full - sseq_obs[left])/(sseq_obs[right] - sseq_obs[left])
    sfrac[left == right] <- 1
  }

  list(left = left, right = right, frac = sfrac)
}


## Map coefficients \code{beta} for specified value(s) of \code{lambda} from
## fitted coefficient matrix for \code{lam} sequence via linear interpolation.
lamtobeta <- function(lambda, beta, s) {
  lamlist <- point.interp(lambda, s)
  if(length(s) == 1) {
    beta = beta[, lamlist$left, drop=FALSE] * (1 - lamlist$frac) +
      beta[, lamlist$right, drop=FALSE] * lamlist$frac
  } else {
    beta = beta[, lamlist$left, drop=FALSE] %*% diag(1 - lamlist$frac) +
      beta[, lamlist$right, drop=FALSE] %*% diag(lamlist$frac)
  }
  colnames(beta) <- paste(seq(along = s))
  return(beta)
}


## Linear prediction
predict.linear <- function(object, newX, s = NULL) {
  beta <- object$beta

  if(is.null(dim(newX))) newX = matrix(newX, nrow = 1)

  if( (nrow(beta) - 1) != ncol(newX) ) {
    stop("numbers of variables in data is not consistant with estimated coefficients")
  }

  if(!is.null(s)) {
    beta <- lamtobeta(lambda = object$lam, beta = beta, s = s)
  }

  fitting <- as.matrix(cbind2(newX, 1)) %*% beta  #fitting <- cbind2(newX, 1) %*% beta
  return(fitting)
}


## Get cross-validation result
cv.test <- function(outlist, y, X, foldid, lam, trim = 0, keep = FALSE) {
  nlam <- length(lam)
  n <- length(y)
  predmat <- matrix(NA, n, nlam) # prediction matrix
  nfolds <- max(foldid)

  for (i in seq(nfolds)) {
    which <- foldid == i
    pred <- predict.linear(outlist[[i]], X[which, , drop = FALSE], s = NULL)
    nlami <- length(outlist[[i]]$lam)
    predmat[which, seq(nlami)] <- pred
  }

  cvraw <- (y - predmat)^2
  N <- length(y) - apply(is.na(predmat), 2, sum)
  cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
  cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(N - 1))

  if(trim > 0) {
    cv.trim <- apply(cvraw, 2, function(x) {
      x <- x[!is.na(x)]
      x.boundary <- quantile(x, probs = c(trim / 2, 1 - trim / 2))
      x <- x[x < x.boundary[2]]
      x <- x[x >= x.boundary[1]]
      x.mean <- mean(x)
      x.sd <- sd(x)
      return(c(MEAN = x.mean, SD = x.sd))
    })
  } else {cv.trim <- NULL}

  output <- list()
  output$cvm <- cvm
  output$cvsd <- cvsd
  output$cvmtrim <- cv.trim[1, ]
  output$cvsdtrim <- cv.trim[2, ]
  if(keep) output$fit.preval <- predmat
  return(output)
}


## Get optimal value of lambda from cross-validation tuned compCL model
getmin <- function(lam, cvm, cvsd, digits = 5) {
  idx <- !is.na(cvm)
  cvm1 <- cvm[idx]
  cvsd1 <- cvsd[idx]
  cvm1 <- round(cvm1, digits = digits)
  cvmin <- min(cvm1, na.rm = TRUE)
  idmin <- cvm1 <= cvmin
  lam.min <- max(lam[idmin])
  idmin <- match(lam.min, lam)
  semin <- round((cvm1 + cvsd1), digits = digits)[idmin]
  idmin <- cvm1 <= semin
  lam.1se <- max(lam[idmin])
  output <- list(lam.min = lam.min, lam.1se = lam.1se)
  return(output)
}

## Get optimal value of lam and k
## from cross-validation-tuned or GIC-tuned FuncompCGL model.
ggetmin <- function(lam, cvm, cvsd, digits = 10, k_list) {
  cvm1 <- round(cvm, digits = digits)
  cvmin <- min(cvm1, na.rm = TRUE)
  lam.min <- apply(cvm1 <= cvmin, 1, function(x, lam)
                   ifelse(length(lam[x]) == 0, 0, max(lam[x], na.rm = TRUE)),
                   lam=lam)
  lam.min_k <- which.max(lam.min)
  lam.min <- max(lam.min)
  output <- list(lam.min = c("lam" = lam.min, "df" = k_list[lam.min_k]))
  if(!missing(cvsd)) {
    idmin <- match(lam.min, lam)
    semin <- round((cvm + cvsd), digits = digits)[lam.min_k, idmin]
    lam.1se <- apply(cvm1 <= semin, 1, function(x, lam)
                     ifelse(length(lam[x]) == 0, 0, max(lam[x], na.rm = TRUE)),
                     lam=lam)
    lam.1se_k <- which.max(lam.1se)
    lam.1se <- max(lam.1se)
    output <- c(output,
                list(lam.1se = c("lam" = lam.1se, "df" = k_list[lam.1se_k])))
  }
  return(output)
}

## Upper and lower curvers for cross-validation plot
error.bars <- function(x, upper, lower, width = 0.02, ...)
{
  xlim <- range(x)
  barw <- diff(xlim) * width
  segments(x, upper, x, lower, ...)
  segments(x - barw, upper, x + barw, upper, ...)
  segments(x - barw, lower, x + barw, lower, ...)
  range(upper, lower)
}


### Initial value for lam in FuncompCGL method.
### Maximum lam value is specified.
### Minmun lam value is computed as lambda.factor times the maximum lam value.
lam.ini <- function(y, Z, ix, iy, pf) {
  lam_max <- 0
  n <- length(y)
  X <- t(Z) %*% y / n
  for(i in 1:length(ix)) {
    a <- sqrt(sum(X[ix[i]:iy[i]]^2))
    if(pf[i] > 0) lam_max <- max(lam_max / pf[i], a)
  }
  return(lam_max)
}


vet <- function(X, p, k) {
  b <- X[-(1: (p * k))]
  C <- matrix(X[1 : (p * k)], byrow = TRUE, nrow = p)
  rownames(C) <- paste0("X", 1:p)
  return(list(b = b, C = C ))
}

Nzero <- function(X, p, k) {
  Non.zero <- apply(vet(X, p, k)$C, 1, function(x) {
    test <- max(abs(x)) == 0 #all.equal(x, rep(0, k))
    #ifelse(test %in% "TRUE", TRUE, FALSE)
  })
  Non.zero <- (1:p)[ !Non.zero]
  return(Non.zero)
}

