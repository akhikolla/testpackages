
Robust <- function(xx, y, CLC, s, L, max.steps, hatAlpha, hatBeta, sparse, structure, hyper, debugging){
  nclc = ncol(CLC)
  n = nrow(xx)

  quant = 0.5
  xi1 = (1-2*quant)/(quant*(1-quant))
  xi2 = sqrt(2/(quant*(1-quant)))
  hatTau = 1
  hatV = rep(1,n) # rgamma(n, shape=1, rate=hatTau)
  hatEtaSq = 1
  hatSg = rep(1, s)
  hatPi = 0.5
  invSigAlpha0 = diag(10^-3, nclc)

  sh0_1 = ifelse(is.null(hyper$a0), 1, hyper$a0)
  sh0_0 = ifelse(is.null(hyper$b0), 1, hyper$b0)
  sh1_1 = ifelse(is.null(hyper$a1), 1, hyper$a1)
  sh1_0 = ifelse(is.null(hyper$b1), 1, hyper$b1)

  a = ifelse(is.null(hyper$c1), 1, hyper$c1)
  b = ifelse(is.null(hyper$c2), 1, hyper$c2)

  s1 = s2 = ifelse(is.null(hyper$d1), 1, hyper$d1)
  r = ifelse(is.null(hyper$d2), 1, hyper$d2)

  hatPi0 = hatPi1 = hatPi
  hatEta1Sq = hatEta2Sq = hatEtaSq
  r1 = r2 = r


  hatGamma = matrix(1, nrow = L, ncol = s)

  progress = ifelse(debugging, 10^(floor(log10(max.steps))-1), 0)

  if(sparse){
    fit=switch (structure,
                "sparsegroup" = BRSGL_SS(xx, y, CLC, s=s, L=L, max.steps, hatAlpha, hatBg=hatBeta, hatTau, hatV, hatGamma, invSigAlpha0,
                                         hatSsq=1, hatPi0=hatPi, hatPi1=hatPi, xi1, xi2, hatT=1, a, b, sh0_1, sh0_0, sh1_1, sh1_0, cutoff=0.05, progress),
                "group" = BRGL_SS(xx, y, CLC, s, L, max.steps, hatAlpha, hatBeta, hatTau, hatV, hatSg, invSigAlpha0, hatPi, hatEtaSq, xi1, xi2, r, a, b, sh0_1, sh0_0, progress),
                "individual" = BRL_SS(xx, y, CLC, max.steps, hatAlpha, c(hatBeta), hatTau, hatV, hatSg=rep(1, s*L), invSigAlpha0, hatPi, hatEtaSq, xi1, xi2, r, a, b, sh1_1, sh1_0, progress)
    )
  }else{
    hatBeta = hatBeta + 10^-5
    fit=switch (structure,
                "sparsegroup" = BRSGL(xx, y, CLC, s, L, max.steps, hatAlpha, hatBeta, hatTau, hatV, hatSg, hatGamma, invSigAlpha0, hatEta1Sq, hatEta2Sq, xi1, xi2, s1, s2, r1, r2, a, b, progress),
                "group" = BRGL(xx, y, CLC, s, L, max.steps, hatAlpha, hatBeta, hatTau, hatV, hatSg, invSigAlpha0, hatEtaSq, xi1, xi2, r, a, b, progress),
                "individual" = BRL(xx, y, CLC, max.steps, hatAlpha, c(hatBeta), hatTau, hatV, hatSg=rep(1, s*L), invSigAlpha0, hatEtaSq, xi1, xi2, r, a, b, progress)
    )
  }

  # out = list(posterior = list( GS.alpha = fit$GS.alpha,
  #                              GS.beta = fit$GS.beta))
  out = list( GS.alpha = fit$GS.alpha,
              GS.beta = fit$GS.beta)

  if(sparse){
    class(out)=c("Sparse", "RBVS")
  }else{
    class(out)=c("NonSparse", "RBVS")
  }

  out

}
