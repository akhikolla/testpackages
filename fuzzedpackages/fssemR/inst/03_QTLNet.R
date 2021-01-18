##' test QTLnet and Qgraph v.s FSSEM
library(fssemR)
library(qtlnet)
source("./inst/00_SparsemaximumLiklihood.R")

##' @title implQDG
##' @description implementor function of FSSEM solver
##' @param data Data archive of experiment measurements, includeing
##' eQTL matrices, Gene expression matrices of different conditions,
##' marker of eQTLs and data generation SEM model
##' @return List of TPR and FDR
##' @export
##' @importFrom qtlnet qdg
##' @importFrom stringr str_remove
implQDG = function(data = NULL) {
  require(qtlnet)
  B = vector("list", 2)
  for (i in 1:2) {
    out = qdg(
      cross = data$QTL$Cross[[i]],
      phenotype.names = data$QTL$Pheno,
      marker.names = data$QTL$marker,
      QTL = data$QTL$eQTLs,
      alpha = 0.10,
      n.qdg.random.starts = 1,
      skel.method = "pcskel"
    )
    DG = apply(out$DG, 1, function(x) {
      if (x[2] == "---->") {
        c(x[1], x[3])
      } else {
        c(x[3], x[1])
      }
    })
    DG = t(DG)
    DG = apply(DG, 2, function(x) {
      as.numeric(str_remove(x, "g"))
    })
    B[[i]] = matrix(0, nrow = data$Vars$p, ncol = data$Vars$p)
    apply(DG, 1, function(x) {
      B[[i]][x[2], x[1]] <<- TRUE
    })
  }
  TPR4GRN = (TPR(B[[1]], data$Vars$B[[1]], PREC = 1e-3) + TPR(B[[2]], data$Vars$B[[2]], PREC = 1e-3)) / 2
  FDR4GRN = (FDR(B[[1]], data$Vars$B[[1]], PREC = 1e-3) + FDR(B[[2]], data$Vars$B[[2]], PREC = 1e-3)) / 2
  TPR4DiffGRN = TPR(B[[1]] - B[[2]], data$Vars$B[[1]] - data$Vars$B[[2]], PREC = 1e-3)
  FDR4DiffGRN = FDR(B[[1]] - B[[2]], data$Vars$B[[1]] - data$Vars$B[[2]], PREC = 1e-3)

  data.frame(
    TPR = TPR4GRN,
    FDR = FDR4GRN,
    TPRofDiffGRN = TPR4DiffGRN,
    FDRofDiffGRN = FDR4DiffGRN
  )
}


##' @title implSML
##' @description implementor function of Sparse-Maximum-Likelihood algorithm
##' @param data Data archive of experiment measurements, includeing
##' eQTL matrices, Gene expression matrices of different conditions,
##' marker of eQTLs and data generation SEM model
##' @return List of TPR and FDR
implSML = function(data = NULL) {
  gamma = cv.multiRegression(
    data$Data$X,
    data$Data$Y,
    data$Data$Sk,
    ngamma = 50,
    nfold = 5,
    data$Vars$n,
    data$Vars$p,
    data$Vars$k
  )
  fit   = multiRegression(
    data$Data$X,
    data$Data$Y,
    data$Data$Sk,
    gamma,
    data$Vars$n,
    data$Vars$p,
    data$Vars$k,
    trans = FALSE
  )
  Xs = transx(data)
  Ys = data$Data$Y
  Sk = data$Data$Sk
  cv.fit = cv_SMLasso(
    Bs = fit$Bs,
    fs = fit$fs,
    Ys = Ys,
    Xs = Xs,
    sigma2 = fit$sigma2,
    Ng = data$Vars$p,
    nlambda = 20
  )
  cv.lambda = optLasso_cv(cv.fit)
  fitc = vector("list", 2)
  for (i in 1:2) {
    fitc[[i]] = sparse_maximum_likehood_iPALM(
      B = fit$Bs[[i]],
      f = fit$fs[[i]],
      Y = Ys[[i]],
      X = Xs,
      sigma2 = fit$sigma2,
      N = data$Vars$n,
      Ng = data$Vars$p,
      lambda = cv.lambda,
      maxit = 50
    )
  }
  TPR4GRN = (TPR(fitc[[1]]$B, data$Vars$B[[1]], PREC = 1e-3) + TPR(fitc[[1]]$B, data$Vars$B[[2]], PREC = 1e-3)) / 2
  FDR4GRN = (FDR(fitc[[1]]$B, data$Vars$B[[1]], PREC = 1e-3) + FDR(fitc[[1]]$B, data$Vars$B[[2]], PREC = 1e-3)) / 2
  TPR4DiffGRN = TPR(fitc[[1]]$B - fitc[[2]]$B, data$Vars$B[[1]] - data$Vars$B[[2]], PREC = 1e-3)
  FDR4DiffGRN = FDR(fitc[[1]]$B - fitc[[2]]$B, data$Vars$B[[1]] - data$Vars$B[[2]], PREC = 1e-3)

  data.frame(
    TPR = TPR4GRN,
    FDR = FDR4GRN,
    TPRofDiffGRN = TPR4DiffGRN,
    FDRofDiffGRN = FDR4DiffGRN
  )
}
