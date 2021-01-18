#' Tool for vector to matrix
#'
#' @param ZOneDim a vector.
#' @param m the number of cluster.
#' @param n sample size.
#' @return adjacency matrix
#' @examples
#' m <- 20
#' n <- 500
#' ZOneDim <- sample(seq_len(m), n, replace = TRUE)
#' #'
#' \donttest{
#' getZmat(ZOneDim, m, n)
#' }
getZmat <- function(ZOneDim, m, n) {
  Zmat <- matrix(NA, m, n)
  for (i in 1:m) {
    Zmat[i, ] <- as.numeric(ZOneDim == i)
  }
  Zmat
}

#' Log scale ratio calculation
#'
#' @param deno denominator.
#' @param nume numerator.
#' @return result of ratio
#'
#' @examples
#' deno <- log(1)
#' nume <- log(2)
#' #'
#' \donttest{
#' calculateRatio(deno, nume)
#' }
calculateRatio <- function(deno, nume) {
  ## deno nume both in log scale
  maxNume <- max(nume)
  transDeno <- deno - maxNume
  transNume <- nume - maxNume
  res <- exp(transDeno) / (sum(exp(transNume)))
  res
}


#' Convert list of string to vector of string
#'
#' @param stringList list of string
#' @return vector of string
#'
#' @examples
#' stringList <- list("abc")
#' #'
#' \donttest{
#' listToStrVec(stringList)
#' }
listToStrVec <- function(stringList) {
  for (i in 1:length(stringList)) {
    stringList[[i]] <-
      paste0("(", paste0(stringList[[i]], collapse = ""), ")")
  }
  unlist(stringList)
}



#' likelihood
#'
#' @param thetaYList thetaYList
#' @param ZOneDim ZOneDim
#' @param qqVec qqVec
#' @param muBar muBar
#' @param X X
#'
#' @examples
#' set.seed(100)
#' n <- 10
#' p <- 2
#' q <- 1
#' K <- 2
#' m <- 1
#' muBar <- c(0, 0)
#' qVec <- c(1, 1)
#' constraint <- c(0, 0, 0)
#' X <- t(
#'   fabMix::simData(
#'     sameLambda = TRUE,
#'     sameSigma = TRUE,
#'     K.true = K,
#'     n = n,
#'     q = q,
#'     p = p,
#'     sINV_values = 1 / ((1:p))
#'   )$data
#' )
#' hparam <- new(
#'   "Hparam",
#'   alpha1 = 0.567755037123148,
#'   alpha2 = 1.1870201935945,
#'   delta = 2,
#'   ggamma = 2,
#'   bbeta = 3.39466184520673
#' )
#' ZOneDim <- sample(seq_len(m), n, replace = TRUE)
#' thetaYList <-
#'   new(
#'     "ThetaYList",
#'     tao = 0.366618687752634,
#'     psy = list(structure(
#'       c(
#'         4.18375613018654,
#'         0, 0, 5.46215996830771
#'       ),
#'       .Dim = c(2L, 2L)
#'     )),
#'     M = list(structure(
#'       c(
#'         3.27412045866392,
#'         -2.40544145363349
#'       ),
#'       .Dim = 1:2
#'     )),
#'     lambda = list(structure(
#'       c(
#'         2.51015961514781,
#'         -0.0741189919182549
#'       ),
#'       .Dim = 2:1
#'     )),
#'     Y = list(structure(
#'       c(
#'         -0.244239011725104,
#'         -0.26876172736886,
#'         0.193431511203083,
#'         0.41624466812811,
#'         -0.54581548068437,
#'         -0.0479517628308146,
#'         -0.633383997203325,
#'         0.856855296613208,
#'         0.792850576988512,
#'         0.268208848994559
#'       ),
#'       .Dim = c(1L, 10L)
#'     ))
#'   )
#' #'
#' \donttest{
#' likelihood(thetaYList, ZOneDim, qVec, muBar, X)
#' }
likelihood <- function(thetaYList, ZOneDim, qqVec, muBar, X) {
  m <- length(qqVec)
  n <- ncol(X)
  xvec <- c()
  ## 2.1: X
  xVal <- 0
  for (i in 1:n) {
    ## for subject i
    k <- ZOneDim[i]

    # if(is.vector(X[, ZOneDim == k])){
    #   meanx = X[, ZOneDim == k]
    # }else{
    #   meanx = apply(X[, ZOneDim == k], 1, mean)
    # }

    meanx <- thetaYList@lambda[[k]] %*% thetaYList@Y[[k]][, i] + c(thetaYList@M[[k]])
    varx <- thetaYList@psy[[k]]
    xvec[i] <- (mvtnorm::dmvnorm(x = X[, i], mean = meanx, sigma = varx, log = T))
    xVal <- xVal + mvtnorm::dmvnorm(x = X[, i], mean = meanx, sigma = varx, log = T)
  }

  xVal
  # xvec
}

#' sumerizeZ
#'
#' @param Zlist  Zlist
#' @param index  index
#' @examples
#' Zlist <- list(c(1, 2, 3), c(3, 2, 1), c(2, 2, 2))
#' #'
#' \donttest{
#' sumerizeZ(Zlist)
#' }
sumerizeZ <- function(Zlist, index = 1:length(Zlist)) {
  sampleSize <- length(Zlist[[1]])
  res <- c()
  for (i in 1:sampleSize) {
    temp <- c()
    for (j in index) {
      temp[j] <- Zlist[[j]][i]
    }
    res[i] <- getmode(temp)
  }
  return(res)
}



#' getmode
#'
#' @param v  v
#' @examples
#' #'
#' \donttest{
#' getmode(c(1, 1, 2, 3))
#' }
getmode <- function(v) {
  v <- v[!is.na(v)]
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


#' calculateVarList
#'
#' @param psyList psyList
#' @param lambdaList lambdaList
#' @examples
#' set.seed(100)
#' n <- 10
#' p <- 2
#' q <- 1
#' K <- 2
#' m <- 1
#' muBar <- c(0, 0)
#' qVec <- c(1, 1)
#' constraint <- c(0, 0, 0)
#' X <- t(
#'   fabMix::simData(
#'     sameLambda = TRUE,
#'     sameSigma = TRUE,
#'     K.true = K,
#'     n = n,
#'     q = q,
#'     p = p,
#'     sINV_values = 1 / ((1:p))
#'   )$data
#' )
#' hparam <- new(
#'   "Hparam",
#'   alpha1 = 0.567755037123148,
#'   alpha2 = 1.1870201935945,
#'   delta = 2,
#'   ggamma = 2,
#'   bbeta = 3.39466184520673
#' )
#' ZOneDim <- sample(seq_len(m), n, replace = TRUE)
#' thetaYList <-
#'   new(
#'     "ThetaYList",
#'     tao = 0.366618687752634,
#'     psy = list(structure(
#'       c(
#'         4.18375613018654,
#'         0, 0, 5.46215996830771
#'       ),
#'       .Dim = c(2L, 2L)
#'     )),
#'     M = list(structure(
#'       c(
#'         3.27412045866392,
#'         -2.40544145363349
#'       ),
#'       .Dim = 1:2
#'     )),
#'     lambda = list(structure(
#'       c(
#'         2.51015961514781,
#'         -0.0741189919182549
#'       ),
#'       .Dim = 2:1
#'     )),
#'     Y = list(structure(
#'       c(
#'         -0.244239011725104,
#'         -0.26876172736886,
#'         0.193431511203083,
#'         0.41624466812811,
#'         -0.54581548068437,
#'         -0.0479517628308146,
#'         -0.633383997203325,
#'         0.856855296613208,
#'         0.792850576988512,
#'         0.268208848994559
#'       ),
#'       .Dim = c(1L, 10L)
#'     ))
#'   )
#' #'
#' \donttest{
#' calculateVarList(thetaYList@psy, thetaYList@lambda)
#' }
calculateVarList <- function(psyList, lambdaList) {
  m <- length(psyList)
  varList <- list()
  for (i in 1:m) {
    tempCov <- psyList[[i]] + lambdaList[[i]] %*% t(lambdaList[[i]])
    varList[[i]] <- tempCov
  }
  return(varList)
}


#' clearCurrentThetaYlist
#'
#' @param thetaYList thetaYList
#' @param clusInd clusInd
#' @param mMax mMax
#' @examples
#' set.seed(100)
#' n <- 10
#' p <- 2
#' q <- 1
#' K <- 2
#' m <- 1
#' muBar <- c(0, 0)
#' qVec <- c(1, 1)
#' constraint <- c(0, 0, 0)
#' X <- t(
#'   fabMix::simData(
#'     sameLambda = TRUE,
#'     sameSigma = TRUE,
#'     K.true = K,
#'     n = n,
#'     q = q,
#'     p = p,
#'     sINV_values = 1 / ((1:p))
#'   )$data
#' )
#' hparam <- new(
#'   "Hparam",
#'   alpha1 = 0.567755037123148,
#'   alpha2 = 1.1870201935945,
#'   delta = 2,
#'   ggamma = 2,
#'   bbeta = 3.39466184520673
#' )
#' ZOneDim <- sample(seq_len(m), n, replace = TRUE)
#' thetaYList <-
#'   new(
#'     "ThetaYList",
#'     tao = 0.366618687752634,
#'     psy = list(structure(
#'       c(
#'         4.18375613018654,
#'         0, 0, 5.46215996830771
#'       ),
#'       .Dim = c(2L, 2L)
#'     )),
#'     M = list(structure(
#'       c(
#'         3.27412045866392,
#'         -2.40544145363349
#'       ),
#'       .Dim = 1:2
#'     )),
#'     lambda = list(structure(
#'       c(
#'         2.51015961514781,
#'         -0.0741189919182549
#'       ),
#'       .Dim = 2:1
#'     )),
#'     Y = list(structure(
#'       c(
#'         -0.244239011725104,
#'         -0.26876172736886,
#'         0.193431511203083,
#'         0.41624466812811,
#'         -0.54581548068437,
#'         -0.0479517628308146,
#'         -0.633383997203325,
#'         0.856855296613208,
#'         0.792850576988512,
#'         0.268208848994559
#'       ),
#'       .Dim = c(1L, 10L)
#'     ))
#'   )
#' clusInd <- rep(1, m)
#' mMax <- 1
#' #'
#' \donttest{
#' clearCurrentThetaYlist(thetaYList, clusInd, mMax)
#' }
clearCurrentThetaYlist <- function(thetaYList, clusInd, mMax) {
  resThetaYList <- thetaYList

  for (i in 1:mMax) {
    if (clusInd[i] == 1) {
      resThetaYList@tao[i] <- thetaYList@tao[i]
      resThetaYList@psy[[i]] <- thetaYList@psy[[i]]
      resThetaYList@M[[i]] <- thetaYList@M[[i]]
      resThetaYList@lambda[[i]] <- thetaYList@lambda[[i]]
      resThetaYList@Y[[i]] <- thetaYList@Y[[i]]
    } else {
      resThetaYList@tao[i] <- 0
      resThetaYList@psy[[i]] <- NA
      resThetaYList@M[[i]] <- NA
      resThetaYList@lambda[[i]] <- NA
      resThetaYList@Y[[i]] <- NA
    }
  }
  return(resThetaYList)
}


#' combineClusterPara
#'
#' @param oldList oldList
#' @param newList newList
#' @param ind ind
#'
#' @examples
#' set.seed(100)
#' n <- 10
#' p <- 2
#' q <- 1
#' K <- 2
#' m <- 1
#' muBar <- c(0, 0)
#' qVec <- c(1, 1)
#' constraint <- c(0, 0, 0)
#' X <- t(
#'   fabMix::simData(
#'     sameLambda = TRUE,
#'     sameSigma = TRUE,
#'     K.true = K,
#'     n = n,
#'     q = q,
#'     p = p,
#'     sINV_values = 1 / ((1:p))
#'   )$data
#' )
#' hparam <- new(
#'   "Hparam",
#'   alpha1 = 0.567755037123148,
#'   alpha2 = 1.1870201935945,
#'   delta = 2,
#'   ggamma = 2,
#'   bbeta = 3.39466184520673
#' )
#' ZOneDim <- sample(seq_len(m), n, replace = TRUE)
#' newList <- oldList <-
#'   new(
#'     "ThetaYList",
#'     tao = 0.366618687752634,
#'     psy = list(structure(
#'       c(
#'         4.18375613018654,
#'         0, 0, 5.46215996830771
#'       ),
#'       .Dim = c(2L, 2L)
#'     )),
#'     M = list(structure(
#'       c(
#'         3.27412045866392,
#'         -2.40544145363349
#'       ),
#'       .Dim = 1:2
#'     )),
#'     lambda = list(structure(
#'       c(
#'         2.51015961514781,
#'         -0.0741189919182549
#'       ),
#'       .Dim = 2:1
#'     )),
#'     Y = list(structure(
#'       c(
#'         -0.244239011725104,
#'         -0.26876172736886,
#'         0.193431511203083,
#'         0.41624466812811,
#'         -0.54581548068437,
#'         -0.0479517628308146,
#'         -0.633383997203325,
#'         0.856855296613208,
#'         0.792850576988512,
#'         0.268208848994559
#'       ),
#'       .Dim = c(1L, 10L)
#'     ))
#'   )
#' #'
#' \donttest{
#' combineClusterPara(oldList, newList, 1)
#' }
combineClusterPara <- function(oldList, newList, ind) {
  resList <- oldList
  resList@tao <- oldList@tao * (1 - newList@tao)
  resList@tao[ind] <- newList@tao

  resList@psy[ind] <- newList@psy
  resList@M[ind] <- newList@M
  resList@lambda[ind] <- newList@lambda
  resList@Y[ind] <- newList@Y
  return(resList)
}




#' getIndThetaY
#'
#' @param thetaYList  thetaYList
#' @param Ind Ind
#' @examples
#' set.seed(100)
#' n <- 10
#' p <- 2
#' q <- 1
#' K <- 2
#' m <- 2
#' muBar <- c(0, 0)
#' qVec <- c(1, 1)
#' constraint <- c(0, 0, 0)
#' X <- t(
#'   fabMix::simData(
#'     sameLambda = TRUE,
#'     sameSigma = TRUE,
#'     K.true = K,
#'     n = n,
#'     q = q,
#'     p = p,
#'     sINV_values = 1 / ((1:p))
#'   )$data
#' )
#' hparam <- new(
#'   "Hparam",
#'   alpha1 = 0.567755037123148,
#'   alpha2 = 1.1870201935945,
#'   delta = 2,
#'   ggamma = 2,
#'   bbeta = 3.39466184520673
#' )
#' ZOneDim <- sample(seq_len(m), n, replace = TRUE)
#' thetaYList <-
#'   new(
#'     "ThetaYList",
#'     tao = 0.366618687752634,
#'     psy = list(structure(
#'       c(
#'         4.18375613018654,
#'         0, 0, 5.46215996830771
#'       ),
#'       .Dim = c(2L, 2L)
#'     )),
#'     M = list(structure(
#'       c(
#'         3.27412045866392,
#'         -2.40544145363349
#'       ),
#'       .Dim = 1:2
#'     )),
#'     lambda = list(structure(
#'       c(
#'         2.51015961514781,
#'         -0.0741189919182549
#'       ),
#'       .Dim = 2:1
#'     )),
#'     Y = list(structure(
#'       c(
#'         -0.244239011725104,
#'         -0.26876172736886,
#'         0.193431511203083,
#'         0.41624466812811,
#'         -0.54581548068437,
#'         -0.0479517628308146,
#'         -0.633383997203325,
#'         0.856855296613208,
#'         0.792850576988512,
#'         0.268208848994559
#'       ),
#'       .Dim = c(1L, 10L)
#'     ))
#'   )
#' #'
#' \donttest{
#' getIndThetaY(thetaYList, 1)
#' }
getIndThetaY <- function(thetaYList, Ind) {
  new("ThetaYList",
    tao = thetaYList@tao[Ind],
    psy = thetaYList@psy[Ind],
    M = thetaYList@M[Ind],
    lambda = thetaYList@lambda[Ind],
    Y = thetaYList@Y[Ind]
  )
}

#' getRemovedIndThetaY
#'
#' @param thetaYList thetaYList
#' @param Ind Ind
#'
#' @examples
#' set.seed(100)
#' n <- 10
#' p <- 2
#' q <- 1
#' K <- 2
#' m <- 2
#' muBar <- c(0, 0)
#' qVec <- c(1, 1)
#' constraint <- c(0, 0, 0)
#' X <- t(
#'   fabMix::simData(
#'     sameLambda = TRUE,
#'     sameSigma = TRUE,
#'     K.true = K,
#'     n = n,
#'     q = q,
#'     p = p,
#'     sINV_values = 1 / ((1:p))
#'   )$data
#' )
#' hparam <- new(
#'   "Hparam",
#'   alpha1 = 0.567755037123148,
#'   alpha2 = 1.1870201935945,
#'   delta = 2,
#'   ggamma = 2,
#'   bbeta = 3.39466184520673
#' )
#' ZOneDim <- sample(seq_len(m), n, replace = TRUE)
#' thetaYList <-
#'   new("ThetaYList", tao = c(0.90162050961987, 0.0983794903801295),
#'   psy = list(structure(c(3.68472841602225, 0, 0, 8.34691978354054),
#'   .Dim = c(2L, 2L)), structure(c(0.785011896130842, 0, 0, 1.19022383323437),
#'   .Dim = c(2L, 2L))), M = list(structure(c(
#'     2.96424305287004,
#'     1.08454861414306
#'   ), .Dim = 1:2), structure(c(
#'     -0.232625450433964,
#'     0.984505960868685
#'   ), .Dim = 1:2)), lambda = list(structure(c(
#'     -0.964026624054337,
#'     0.89378616732449
#'   ), .Dim = 2:1), structure(c(
#'     0.533334148228635,
#'     -1.80033696090263
#'   ), .Dim = 2:1)), Y = list(structure(c(
#'     -0.15346475266988,
#'     1.6584112693271, 0.409294936277862, -1.46628591247549, -0.532753243163142,
#'     -0.332143130316749, 0.307558110800446, -0.525374243612587, 0.527667526535661,
#'     0.748193650431916
#'   ), .Dim = c(1L, 10L)), structure(c(
#'     0.571325118638535,
#'     0.542462985882966, 0.559971315637159, -1.73905343105432, -0.583549598471542,
#'     1.71264245945391, -0.327119395945831, 1.02464651767821, -1.11462280255215,
#'     0.81095592501554
#'   ), .Dim = c(1L, 10L))))
#' Ind <- 1
#' #'
#' \donttest{
#' getRemovedIndThetaY(thetaYList, Ind)
#' }
getRemovedIndThetaY <- function(thetaYList, Ind) {
  res <- thetaYList

  res@tao[Ind] <- NA
  res@tao <- res@tao / sum(res@tao, na.rm = T)
  res@psy[[Ind]] <- NA
  res@M[[Ind]] <- NA
  res@lambda[[Ind]] <- NA
  res@Y[[Ind]] <- NA

  res
}

#' changeConstraintFormat
#'
#' @param strNum strNum
#' @examples
#' #'
#' \donttest{
#' changeConstraintFormat(c(0, 0, 0))
#' }
changeConstraintFormat <- function(strNum) {
  res <- ""
  res <- gsub(pattern = "0", replacement = "U", x = strNum)
  res <- gsub(pattern = "1", replacement = "C", x = res)
  res <- gsub(pattern = "\\(", replacement = "", x = res)
  res <- gsub(pattern = "\\)", replacement = "", x = res)
  return(res)
}
