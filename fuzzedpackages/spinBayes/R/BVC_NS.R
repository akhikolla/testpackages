
BVC_NS <- function(xx, y, CLC0, EX, s, q, max.steps, hat.m, hat.r, hat.clc, hat.zeta, sparse, hyper, debugging){

  invTAUsq.star = rep(0.1, s)
  invTAUsq.zeta = rep(0.1, s)
  invSig.clc0 = 10^-3 #1 # 10^-3
  invSigM0 = rep(10^-3, q)
  hat.sigma.sq = 1

  lambda.star = 1 #length(hat.r.star) * sqrt(hat.sigma.sq) / sum(abs(hat.r.star))
  lambda.zeta = 1

   hat.pi.s = 0.9
  hat.pi.ze = 0.9

  aE = ifelse(is.null(hyper$a.e), 1, hyper$a.e);
  bE = ifelse(is.null(hyper$b.e), 1.5, hyper$b.e);
  a.star = ifelse(is.null(hyper$a.v), 1, hyper$a.v);
  b.star = ifelse(is.null(hyper$b.v), 1.5, hyper$b.v);
  alpha = ifelse(is.null(hyper$s), 0.2, hyper$s);
  gamma = ifelse(is.null(hyper$h), 0.1, hyper$h);

  mu.star = ifelse(is.null(hyper$r.v), 1, hyper$r.v);
  nu.star = ifelse(is.null(hyper$w.v), 1, hyper$w.v);
  muE = ifelse(is.null(hyper$r.e), 1, hyper$r.e);
  nuE = ifelse(is.null(hyper$w.e), 1, hyper$w.e);

  progress = ifelse(debugging, 10^(floor(log10(max.steps))-1), 0)
  noE = TRUE

  if(!is.null(EX)){
    noE = FALSE
    # if(debugging) message("No. of CLC0: ", ncol(CLC0), "\n")
    CLC = CLC0[,-ncol(CLC0), drop=FALSE]
    hat.clin = utils::head(hat.clc, -1)
    invSigCLC0 = diag(invSig.clc0, nrow=ncol(CLC), ncol=ncol(CLC))

	if(sparse){
		out = BVCPointMassNonStr(xx, y, CLC, EX, s, q, max.steps, hat.m, hat.r, hat.clin, hat.zeta, invSigM0, invTAUsq.star, invSigCLC0, invTAUsq.zeta, hat.pi.s, hat.pi.ze, lambda.star, lambda.zeta, hat.sigma.sq, aE, bE, a.star, b.star, alpha, gamma, muE, mu.star, nuE, nu.star, progress)
	}else{
		out = BVCNoStr(xx, y, CLC, EX, s, q, max.steps, hat.m, hat.r, hat.clin, hat.zeta, invSigM0, invTAUsq.star, invSigCLC0, invTAUsq.zeta, lambda.star, lambda.zeta, hat.sigma.sq, aE, bE, a.star, b.star, alpha, gamma, progress)
	}

  }else{
    # if(debugging) message("No. of CLC0: ", ncol(CLC0), "\n")
    if(ncol(CLC0)==1){
      CLC = matrix(0, nrow = nrow(CLC0), ncol = 1)
      CLIN = FALSE
      hat.clin = 0
    }else{
      CLC = CLC0[,-ncol(CLC0),drop=FALSE]
      hat.clin = utils::head(hat.clc, -1)
      CLIN = TRUE
    }
    invSigCLC0 = diag(invSig.clc0, nrow=ncol(CLC), ncol=ncol(CLC))

	if(sparse){
		out = BVCPointMassNonStr_NoE(xx, y, CLC, CLIN, s, q, max.steps, hat.m, hat.r, hat.clin, invSigM0, invTAUsq.star, invSigCLC0, hat.pi.s, lambda.star, hat.sigma.sq, a.star, b.star, alpha, gamma, mu.star, nu.star, progress)
	}else{
		out = BVCNoStr_NoE(xx, y, CLC, CLIN, s, q, max.steps, hat.m, hat.r, hat.clin, invSigM0, invTAUsq.star, invSigCLC0, lambda.star, hat.sigma.sq, a.star, b.star, alpha, gamma, progress)
	}

  }

  BVC = list(posterior = list( GS.m = out$GS.m,
                               GS.clc = out$GS.alpha,
                               GS.zeta = out$GS.zeta,
                               GS.rs = out$GS.rs,
                               # GS.phi = out$GS.tRsRs,
                               GS.invTAUsq.zeta = out$GS.invTAUsq.zeta,
                               GS.invTAUsq.star = out$GS.invTAUsq.star,
                               GS.lambdaSq.zeta = out$GS.lambda.sq.zeta,
                               GS.lambdaSq.star = out$GS.lambda.sq.star,
                               GS.sigma.sq = out$GS.sigma.sq))

  if(sparse){
    BVC$posterior$GS.phi = out$GS.tRsRs
    BVC$posterior$GS.phi[BVC$posterior$GS.phi!=0] = 1
    class(BVC) = "BVCSparse"
    if(debugging){
      debugList = list(	GS.tRsRs = out$GS.tRsRs,
                        GS.l0 = out$GS.l0,
                        GS.lz = out$GS.lZ,
                        GS.ls = out$GS.lS,
                        GS.pi0 = out$GS.pi0,
                        GS.pi.zeta = out$GS.pi.zeta,
                        GS.pi.star = out$GS.pi.star,
						hyper = c(r.v=mu.star, r.e=muE, w.v=nu.star, w.e=nuE))
      BVC$debugList = debugList
    }
  }else{
    class(BVC) = "BVCNonSparse"
  }

  if(noE){
    BVC$posterior$GS.zeta = BVC$posterior$GS.invTAUsq.zeta = BVC$posterior$GS.lambdaSq.zeta = NULL
    BVC$debugList$GS.lz = BVC$debugList$GS.pi.zeta = NULL
    if(!CLIN) BVC$posterior$GS.clc = NULL
    class(BVC) = c("VarOnly", class(BVC))
  }else{
    class(BVC) = c("VarLin", class(BVC))
  }
  BVC
}
