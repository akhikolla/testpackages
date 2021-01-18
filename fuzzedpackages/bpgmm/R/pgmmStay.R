#' stayMCMCupdate
#'
#' @param X X
#' @param thetaYList thetaYList
#' @param ZOneDim ZOneDim
#' @param hparam hparam
#' @param qVec qVec
#' @param qnew qnew
#' @param dVec dVec
#' @param sVec sVec
#' @param constraint constraint
#' @param clusInd clusInd
#' @export
#'
#' @examples
#' #'
#' set.seed(110)
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
#' new("ThetaYList", tao = c(0.90162050961987, 0.0983794903801295),
#' psy = list(structure(c(3.68472841602225, 0, 0, 8.34691978354054),
#' .Dim = c(2L, 2L)), structure(c(0.785011896130842, 0, 0, 1.19022383323437),
#' .Dim = c(2L, 2L))), M = list(structure(c(
#'   2.96424305287004,
#'   1.08454861414306
#' ), .Dim = 1:2), structure(c(
#'   -0.232625450433964,
#'   0.984505960868685
#' ), .Dim = 1:2)), lambda = list(structure(c(
#'   -0.964026624054337,
#'   0.89378616732449
#' ), .Dim = 2:1), structure(c(
#'   0.533334148228635,
#'   -1.80033696090263
#' ), .Dim = 2:1)), Y = list(structure(c(
#'   -0.15346475266988,
#'   1.6584112693271, 0.409294936277862, -1.46628591247549, -0.532753243163142,
#'   -0.332143130316749, 0.307558110800446, -0.525374243612587, 0.527667526535661,
#'   0.748193650431916
#' ), .Dim = c(1L, 10L)), structure(c(
#'   0.571325118638535,
#'   0.542462985882966, 0.559971315637159, -1.73905343105432, -0.583549598471542,
#'   1.71264245945391, -0.327119395945831, 1.02464651767821, -1.11462280255215,
#'   0.81095592501554
#' ), .Dim = c(1L, 10L))))
#' qnew <- 1
#' dVec <- c(1, 1, 1)
#' sVec <- c(1, 1, 1)
#' constraint <- c(0, 0, 0)
#' clusInd <- rep(1, m)
#' \donttest{
#' stayMCMCupdate(
#'   X,
#'   thetaYList,
#'   ZOneDim,
#'   hparam,
#'   qVec,
#'   qnew,
#'   dVec,
#'   sVec,
#'   constraint,
#'   clusInd
#' )
#' }
stayMCMCupdate <- function(X,
                           thetaYList,
                           ZOneDim,
                           hparam,
                           qVec,
                           qnew,
                           dVec,
                           sVec,
                           constraint,
                           clusInd) {
  ggamma <- hparam@ggamma
  m <- length(which(clusInd == 1))
  p <- nrow(X)
  n <- ncol(X)
  ## transfor to non empty obj
  sampleVec <- which(clusInd == 1)
  NEthetaYList <- getIndThetaY(thetaYList, sampleVec)

  NEZOneDim <- ZOneDim
  for (i in 1:length(sampleVec)) {
    NEZOneDim[NEZOneDim == sampleVec[i]] <- i
  }
  ## update

  NEthetaYList <- updatePostThetaY(m, n, p, hparam, NEthetaYList, NEZOneDim, qVec[qVec != 0], constraint, X, ggamma)
  NEZOneDim <- updatePostZ(X, m, n, NEthetaYList)
  hparam <- update_Hyperparameter(m, p, qnew, hparam, NEthetaYList, dVec, sVec)


  ZOneDim <- sampleVec[NEZOneDim]

  thetaYList <- getThetaYWithEmpty(NEthetaYList, clusInd)
  return(list(thetaYList = thetaYList, ZOneDim = ZOneDim, hparam = hparam))
}



#' getThetaYWithEmpty
#'
#' @param NEthetaYList NEthetaYList
#' @param clusInd clusInd
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
#' \donttest{
#' getThetaYWithEmpty(thetaYList, clusInd)
#' }
#'
getThetaYWithEmpty <- function(NEthetaYList, clusInd) {
  thetaYList <- new("ThetaYList")
  j <- 1
  for (i in 1:length(clusInd)) {
    if (clusInd[i] == 1) {
      thetaYList@tao[i] <- NEthetaYList@tao[j]
      thetaYList@psy[[i]] <- NEthetaYList@psy[[j]]
      thetaYList@M[[i]] <- NEthetaYList@M[[j]]
      thetaYList@lambda[[i]] <- NEthetaYList@lambda[[j]]
      thetaYList@Y[[i]] <- NEthetaYList@Y[[j]]
      j <- j + 1
    } else {
      thetaYList@tao[i] <- 0
      thetaYList@psy[[i]] <- NA
      thetaYList@M[[i]] <- NA
      thetaYList@lambda[[i]] <- NA
      thetaYList@Y[[i]] <- NA
    }
  }
  thetaYList
}

#' toNEthetaYlist
#'
#' @param thetaYList thetaYList
#' @param ZOneDim ZOneDim
#' @param qVec qVec
#' @param clusInd clusInd
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
#' \donttest{
#' toNEthetaYlist(thetaYList, ZOneDim, qVec, clusInd)
#' }
#'
toNEthetaYlist <- function(thetaYList, ZOneDim, qVec, clusInd) {
  ## transfor to non empty obj
  sampleVec <- which(clusInd == 1)
  NEthetaYList <- getIndThetaY(thetaYList, sampleVec)

  NEZOneDim <- ZOneDim
  for (i in 1:length(sampleVec)) {
    NEZOneDim[NEZOneDim == sampleVec[i]] <- i
  }
  qnew <- max(qVec)
  qVec <- qVec[qVec == qnew]
  return(list(thetaYList = NEthetaYList, ZOneDim = NEZOneDim, qVec = qVec))
}

#' Title
#'
#' @param NEthetaYList NEthetaYList
#' @param NEZOneDim NEZOneDim
#' @param qnew qnew
#' @param clusInd clusInd
#'
#' @export
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
#' qnew <- 1
#' toEthetaYlist(thetaYList, ZOneDim, qnew, clusInd)
toEthetaYlist <- function(NEthetaYList, NEZOneDim, qnew, clusInd) {
  ## transfor back
  sampleVec <- which(clusInd == 1)
  ZOneDim <- NEZOneDim
  for (i in length(sampleVec):1) {
    ZOneDim[ZOneDim == i] <- sampleVec[i]
  }
  qVec <- clusInd
  qVec[clusInd == 1] <- qnew
  thetaYList <- getThetaYWithEmpty(NEthetaYList, clusInd)
  return(list(thetaYList = thetaYList, ZOneDim = ZOneDim, qVec = qVec))
}
