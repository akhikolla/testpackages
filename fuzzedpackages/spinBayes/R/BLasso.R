
BLasso <- function(xx, y, CLC, EX, ZX, s, max.steps, hat.m, hat.r, hat.clc, hat.zeta, hyper, debugging){

  invTAUsq.star = rep(0.1, s)
  invTAUsq.zeta = rep(0.1, length(hat.zeta))
  invSig.clc0 = 10^-3 #1 # 10^-3
  invSigM0 = 10^-3
  hat.sigma.sq = 1

  lambda.star = 1
  lambda.zeta = 1

  aE = ifelse(is.null(hyper$a.e), 1, hyper$a.e);
  bE = ifelse(is.null(hyper$b.e), 1.5, hyper$b.e);
  a.star = ifelse(is.null(hyper$a.v), 1, hyper$a.v);
  b.star = ifelse(is.null(hyper$b.v), 1.5, hyper$b.v);
  alpha = ifelse(is.null(hyper$s), 0.2, hyper$s);
  gamma = ifelse(is.null(hyper$h), 0.1, hyper$h);       # for sigma.sq

  EXZX = cbind(EX, ZX)
  progress = ifelse(debugging, 10^(floor(log10(max.steps))-1), 0)
  noE = FALSE

  # if(debugging) message("No. of CLC: ", ncol(CLC), "\n")
  nclc = ncol(CLC)
  invSigCLC0 = diag(invSig.clc0, nrow=nclc, ncol=nclc)

  out = BayesLasso(xx, y, CLC, EXZX, s, max.steps, hat.m, hat.r, hat.clc, hat.zeta, invSigM0, invTAUsq.star, invSigCLC0, invTAUsq.zeta, lambda.star, lambda.zeta, hat.sigma.sq, aE, bE, a.star, b.star, alpha, gamma, progress)

  if(is.null(EX)){
    GS.zeta = NULL
	GS.rs = out$GS.zeta
  }else{
	GS.zeta = out$GS.zeta[,1:s]
	GS.rs = out$GS.zeta[,-c(1:s)]
  }

  if(nclc==1){
	GS.Z = out$GS.alpha
	GS.clc = NULL
  }else{
	GS.clc = out$GS.alpha[,-nclc,drop=FALSE]
	GS.Z = out$GS.alpha[,nclc,drop=FALSE]
  }

  BVC = list(posterior = list( GS.m = out$GS.m,
                               GS.clc = GS.clc,
							   GS.Z = GS.Z,
                               GS.zeta = GS.zeta,
							   GS.r0 = out$GS.r0,
                               GS.rs = GS.rs,
                               GS.invTAUsq.zeta = out$GS.invTAUsq.zeta,
                               GS.invTAUsq.star = out$GS.invTAUsq.star,
                               GS.lambdaSq.zeta = out$GS.lambda.sq.zeta,
                               GS.lambdaSq.star = out$GS.lambda.sq.star,
                               GS.sigma.sq = out$GS.sigma.sq))

  class(BVC) = c("LinOnly", "BVCNonSparse")
  BVC
}
