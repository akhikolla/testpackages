Cpredict <- function (focal, Krigobject) {
    CKrigidx <- Krigobject$CKrigidx
    if (is.null(CKrigidx)) {
      #localCKrigptr->fillaxialFocal(focal_C):
      if (is.null(Krigobject$unique_x)) Krigobject$unique_x <- as.matrix(unique(Krigobject$x))
      dX <- sweep(Krigobject$unique_x,2L,focal,FUN=`-`)
      dX2 <- dX^2
      #localCKrigptr->filleuclFocal():
      if (is.null(Krigobject$covTheta2)){
        CovFnParam=blackbox.getOption("CovFnParam")
        Krigobject$smoothness=CovFnParam["smoothness"]
        Krigobject$covTheta2 <- CovFnParam[colnames(Krigobject$x)]^2
      }
      eucl <- sweep(dX2,2L,Krigobject$covTheta2,FUN=`/`)
      eucl <- sqrt(rowSums(eucl))
      #localCKrigptr->fillcovFocal<covTypedef>():
      covfocal <- MaternCorr(eucl,smoothness=Krigobject$smoothness)
    } else {
      if (CKrigidx < 0) {
        stop.redef(paste("(!) Invalid 'CKrigidx' value", CKrigidx,
                         "in Cpredict(...)"))
      }
      if (CKrigidx > 100) {
        message.redef(paste("(!) Suspicious 'CKrigidx' value",
                            CKrigidx, "in Cpredict(...)"))
        message.redef(paste("(!) Might be correct, but then still suggests a poor organization of code."))
      }
      covfocal <- CcovFocal(focal,CKrigidx)
    }
    return(Krigobject$d + sum(covfocal * Krigobject$c))
}
