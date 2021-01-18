#############################################################################
#############################################################################
#### Functions in Univariate Autoregressive Cokriging Models


#############################################################################
#############################################################################

margin.posterior= function(param, input, output, level, H, dist, cov.model="matern_5_2",
  prior="JR", hyperparam=NULL){

Dim = dim(dist)[3]
p.x = Dim


S = length(output)

if(is.null(hyperparam)){
  al = 0.5-p.x
  bl = 1 * dim(dist)[1]^(-1/p.x) * (al + p.x)
  nugget.UB = 1
}else{
  al = hyperparam$a
  bl = hyperparam$b * dim(dist)[1]^(-1/p.x) * (al + p.x)
  nugget.UB = hyperparam$nugget.UB
}

  # param contains phi and nugget (maybe)
  if(length(param)==Dim){ #no nugget
    phi = exp(-param)
    is.nugget=FALSE
  }else{
    phi = exp(-param[1:Dim]) 
    nugget = nugget.UB*exp(param[Dim+1]) / (1 + exp(param[Dim+1]))
    is.nugget = TRUE
    phi = c(phi, nugget)
  }

q = rep(NA, S)
for(t in 1:S){
  q[t] = dim(H[[t]])[2]
}

n = sapply(output, length)

# input.max = apply(input[[level]], 2, max)
# input.min = apply(input[[level]], 2, min)
# Cl =  abs(input.max-input.min)
Cl = hyperparam$Cl



t = level

if(level==1){
  R = buildcov(phi, dist, covmodel=cov.model, nugget=is.nugget)
  U = chol(R)
  RInv = chol2inv(U)
  X.aug = H[[t]]
  y.aug = output[[t]]
    HRH = t(X.aug)%*%RInv%*%X.aug #+ diag(1e-6,dim(X.aug)[2])
    HRHchol = chol(HRH)
    HRHInv = chol2inv(HRHchol)
    betahat = HRHInv%*%(t(X.aug)%*%(RInv%*%y.aug))
    res = y.aug - X.aug%*%betahat
    SSE = c(t(res)%*%RInv%*%res)

    logdetG = -sum(log(diag(U)))

    g = logdetG - sum(log(diag(HRHchol))) + (-0.5)*(n[t]-q[t])*log(SSE)
    if(prior=="Ind_Jeffreys"){
      g = g + 0.5*(q[t])*log(SSE)
    }

}else{

    R = buildcov(phi, dist, covmodel=cov.model, nugget=is.nugget)
    U = chol(R)
    RInv = chol2inv(U)

    IB = match.input(input[[t]], input[[t-1]])$IB
    y_t1 = as.matrix(output[[t-1]][IB])

    X.aug = cbind(H[[t]], y_t1)
    y.aug = output[[t]]

    HRH = t(X.aug)%*%RInv%*%X.aug #+ diag(1e-6,dim(X.aug)[2])
    HRHchol = chol(HRH)
    HRHInv = chol2inv(HRHchol)
    betahat = HRHInv%*%(t(X.aug)%*%(RInv%*%y.aug))
    res = y.aug - X.aug%*%betahat
    SSE = c(t(res)%*%RInv%*%res) 

    logdetG = -sum(log(diag(U)))

    g = logdetG - sum(log(diag(HRHchol))) + (-0.5)*(n[t]-q[t]-1)*log(SSE) 
    if(prior=="Ind_Jeffreys"){
      g = g + 0.5*(q[t]+1)*log(SSE)
    }

}

X = X.aug 

## compute priors

  if(is.nugget){
    lnJacobian = sum(param[1:Dim]) + log(nugget) + log(nugget.UB - nugget) - 
             log(abs(nugget.UB+(nugget.UB-1)*nugget))  # - sum(log(phi)) + log(nugget) 

    if(prior=="JR"){
      temp = sum(1/phi[1:Dim] * Cl)  + nugget
      ln_prior = lnJacobian + al*log(temp) - bl*temp 
    }else{
      lnIJ_prior = log_objective_prior(c(1/phi, nugget), dist, RInv, X, cov.model, is.nugget, prior)
      ln_prior = lnIJ_prior + lnJacobian
    }

  }else{
    lnJacobian = sum(param)
    if(prior=="JR"){
        temp = sum(1/phi * Cl) 
        ln_prior = lnJacobian + al*log(temp) - bl*temp      
    }else{
      lnIJ_prior = log_objective_prior(c(1/phi), dist, RInv, X, cov.model, is.nugget, prior)
      ln_prior = lnIJ_prior + lnJacobian 
    }

  }

  g = ln_prior + g 


  return(g)
}

#############################################################################
#############################################################################






#############################################################################
#############################################################################
condsim.ND.univariate = function(formula, output, input, input.new, phi, cov.model="matern_5_2",
              nsample){


Dim = dim(input[[1]])[2]
p.x = Dim
if(dim(phi)[1]==Dim){
  is.nugget=FALSE
}else{
  is.nugget=TRUE
}

S = length(output)
n = sapply(output, length)
np = dim(input.new)[1] 

y = output

## add new inputs to missing inputs
# input.miss = list()
# input.union = list()
# for(t in 1:S){
#   input.miss[[t]] = input.new
#   input.union[[t]] = rbind(input.new, input[[t]])
# }


input.miss = list()
input.union = list()
pred.ID = list()
ID.orig = list()
index.full = 1:np

indB = list()
for(t in 1:(S)){
  ind.list = match.input(input[[t]], input.new)
  # indlist = ismember(input.new, input[[t]])
  if(!is.null(ind.list$IA)){
    indA = ind.list$IA
    indB[[t]] = ind.list$IB
    pred.ID.exist = indA
    input.exist = input[[t]][indA, ,drop=FALSE] # common inputs 
    input.added = input.new[-indB[[t]], ,drop=FALSE] # inputs in input.new but not in input[[t]]
    n.added = dim(input.added)[1]
    pred.ID.added = seq(1, n.added, by=1)
    ID.orig[[t]] = c(index.full[-indB[[t]]], indB[[t]])
    pred.ID[[t]] = c(pred.ID.added, pred.ID.exist+n.added)
    input.miss[[t]] = input.added
    input.union[[t]] = rbind(input.added, input[[t]])
  }else{
    ID.orig[[t]] = 1:np
    pred.ID[[t]] = 1:np
    input.miss[[t]] = input.new
    input.union[[t]] = rbind(input.new, input[[t]])
  }
}

n.m = rep(NA, S)
for(t in 1:S){
  n.m[t] = dim(input.miss[[t]])[1]
}

###################################################################
#### create covariates
###################################################################
H = list()
Hm = list()
for(t in 1:S){
  colnames(input[[t]]) = paste0("x", 1:p.x)
  df = data.frame(input[[t]])
  H[[t]] = model.matrix(formula[[t]], df)

  
  colnames(input.miss[[t]]) = paste0("x", 1:p.x)
  df = data.frame(input.miss[[t]])
  Hm[[t]] = model.matrix(formula[[t]], df)
  
}


###################################################################
#### compute intermediate quantities
###################################################################
dist.o = list()
dist.m = list()
dist.mo = list()
distlist = list()

for(t in 1:S){
  dist.o[[t]] = compute_distance(input[[t]], input[[t]])

  #if(t<S){
    dist.m[[t]] = compute_distance(input.miss[[t]], input.miss[[t]])
    dist.mo[[t]] = compute_distance(input.miss[[t]], input[[t]])
  #}

  distlist[[t]] = compute_distance(input.union[[t]], input.union[[t]])
}






ym.hat = list()
q = rep(NA, S)
for(t in 1:S){
  q[t] = dim(H[[t]])[2]
  ym.hat[[t]] = matrix(NA, nsample, n.m[[t]])
}


  #################################################################################
  #### Sampling from predictive distribution
  #################################################################################
  t = 1
  R = buildcov(phi[ ,t], dist.o[[t]], covmodel=cov.model, nugget=is.nugget)
  U = chol(R)
  RInv = chol2inv(U)
  #UH = backsolve(U, H[[t]], transpose=TRUE)
  #HRHInv = solve(crossprod(UH))
  HRHInv = solve(t(H[[t]])%*%RInv%*%H[[t]])
  betahat = HRHInv%*%t(H[[t]])%*%(RInv%*%y[[t]])
  res = y[[t]] - H[[t]]%*%betahat
  SSE = c(t(res)%*%RInv%*%res)

  Rm = buildcov(phi[ ,t], dist.m[[t]], covmodel=cov.model, nugget=FALSE)
  Rmo = buildcov(phi[ ,t], dist.mo[[t]], covmodel=cov.model, nugget=is.nugget)
  
  RmoU = t(backsolve(U, t(Rmo), transpose=TRUE))

  RmoRInv = Rmo%*%RInv

  # XXRR = t(Hm[[t]]) - t(H[[t]])%*%RInv%*%t(Rmo)
  XXRR = Hm[[t]] - RmoRInv%*%H[[t]]
  
  # Sig_ymymy = Rm - Rmo%*%RInv%*%t(Rmo) + (XXRR)%*%HRHInv%*%t(XXRR)
  Sig_ymymy = Rm - tcrossprod(RmoU) + (XXRR)%*%HRHInv%*%t(XXRR)

  Sig = Sig_ymymy * SSE/(n[t] - q[t])
  #Sig = Sig_ymymy * SSE/(n[t])

  mu_ymy = Hm[[t]]%*%betahat + RmoRInv%*%res
  
  #ym.hat[[t]] = c(mu_ymy)
  #krige[[t]] = ym.hat[[t]][pred.ID[[t]]]

  ym.hat[[t]] = mvtnorm::rmvt(nsample, sigma=Sig, df=n[t]-q[t], delta=mu_ymy, type="shifted")

  #Sig[Sig<0] = 0
  #krigeSE[[t]] = sqrt(diag(Sig)[pred.ID[[t]]])


  # for(k in 1:nsample){
    for(t in 2:(S)){
    
      ############################################################################
      #### estimate missing data ym
      ############################################################################

      R = buildcov(phi[ ,t], dist.o[[t]], covmodel=cov.model, nugget=is.nugget)
      U = chol(R)
      RInv = chol2inv(U)
      IB = match.input(input[[t]], input[[t-1]])$IB
      y_t1 = y[[t-1]][IB]
      
      X = cbind(H[[t]], y_t1)
      XRXInv = solve(t(X)%*%RInv%*%X)
      betahat = XRXInv%*%t(X)%*%(RInv%*%y[[t]])
      res = y[[t]] - X%*%betahat
      SSE = c(t(res)%*%RInv%*%res)



      Rm = buildcov(phi[ ,t], dist.m[[t]], covmodel=cov.model, nugget=FALSE)
      Rmo = buildcov(phi[ ,t], dist.mo[[t]], covmodel=cov.model, nugget=is.nugget)
      
      # IB = match.input(input.miss[[t]], input.miss[[t-1]])$IB
      # ym_t1 = ym.hat[[t-1]][ ,IB]

      ym_t1 = matrix(NA, nsample, n.m[t])
      for(k in 1:nsample){
        ym_t1[k, ] = create.w.new(t=t, input=input[[t-1]], input.miss=input.miss,
              y=y[[t-1]], ym=as.matrix(ym.hat[[t-1]][k,]))
        Xm = cbind(Hm[[t]], as.matrix(ym_t1[k, ]))

        XXRR = t(Xm) - t(X)%*%RInv%*%t(Rmo)
        Sig_ymymy = Rm - Rmo%*%RInv%*%t(Rmo) + t(XXRR)%*%XRXInv%*%XXRR
        Sig = Sig_ymymy * SSE/(n[t]-q[t]-1)
        mu_ymy = Xm%*%betahat + Rmo%*%(RInv%*%res)
        ym.hat[[t]][k, ] = c(mvtnorm::rmvt(1, sigma=Sig, df=n[t]-q[t]-1, delta=mu_ymy, type="shifted"))
      }

    }
  # }



## get summary statistics
krige = list()
krigeSE = list()
krige.lower95 = list()
krige.upper95 = list()
for(t in 1:S){
  # yhat[[t]] = ym.hat[[t]][ ,pred.ID[[t]], ]

  krige[[t]] = apply(ym.hat[[t]], 2, mean)
  krigeSE[[t]] = apply(ym.hat[[t]], 2, sd)
  krige.lower95[[t]] = apply(ym.hat[[t]], 2, quantile, 0.025)
  krige.upper95[[t]] = apply(ym.hat[[t]], 2, quantile, 0.975)
}


pred.mu = list()
pred.SE = list()
pred.lower95 = list()
pred.upper95 = list()
for(t in 1:S){
  pred.mu[[t]] = rep(NA, np)
  pred.SE[[t]] = rep(0, np)
  pred.lower95[[t]] = rep(NA, np)
  pred.upper95[[t]] = rep(NA, np)

  ind.list = ismember(input.new, input.miss[[t]])
  pred.mu[[t]][ind.list$IIA] = krige[[t]][ind.list$IA]
  pred.SE[[t]][ind.list$IIA] = krigeSE[[t]][ind.list$IA]
  pred.lower95[[t]][ind.list$IIA] = krige.lower95[[t]][ind.list$IA]
  pred.upper95[[t]][ind.list$IIA] = krige.upper95[[t]][ind.list$IA]
  if(length(ind.list$IIA)<np){
    ind.input = ismember(input.new, input[[t]])
    pred.mu[[t]][ind.input$IIA] = y[[t]][ind.input$IA]
    pred.lower95[[t]][ind.input$IIA] = y[[t]][ind.input$IA]
    pred.upper95[[t]][ind.input$IIA] = y[[t]][ind.input$IA]
  }
}
names(pred.mu) = paste0("Level", seq(1:S), "")
names(pred.SE) = paste0("Level", seq(1:S), "")
names(pred.lower95) = paste0("Level", seq(1:S), "")
names(pred.upper95) = paste0("Level", seq(1:S), "")

pred = list(mu=pred.mu, SE=pred.SE, 
       lower95=pred.lower95,upper95=pred.upper95)

return(pred)

}


#############################################################################
#############################################################################


#############################################################################
#############################################################################
compute.g.univ = function(param, input.list, level, y, H, ym, Hm, dist, hyper, cov.model){
  ## compute g function at fidelity level t
  ## 

  # param is a vector 
  
  Dim = dim(dist)[3]
  p.x = Dim
  N = ncol(y[[1]])
  nsample = dim(ym[[1]])[1] 

  if(is.null(hyper)){
      al = 0.5-p.x
      bl = dim(dist)[1]^(-1/p.x) * (al + p.x)
      nugget.UB = 1
  }else{
    al = hyper$a
    bl = hyper$b * dim(dist)[1]^(-1/p.x) * (al + p.x)
    nugget.UB = hyper$nugget.UB
  }

  # param contains phi and nugget (maybe)
  if(length(param)==Dim){ #no nugget
    phi = exp(-param)
    # phi = 1/param
    is.nugget=FALSE
  }else{
    phi = exp(-param[1:Dim]) 
    nugget = nugget.UB*exp(param[Dim+1]) / (1 + exp(param[Dim+1]))
    is.nugget = TRUE
    phi = c(phi, nugget)
  }

  # inputlist = augment.input(input)
  # input.miss = inputlist$miss
  # input.union = inputlist$union  
  input = input.list$input
  input.miss = input.list$input.miss
  
  S = length(y)
  n = rep(NA, S)
  n.aug = rep(NA, S)
  for(t in 1:S){
    n[t] = dim(y[[t]])[1]
    if(t<S){
      n.aug[t] = n[t] + dim(Hm[[t]])[1]
    }else{
      n.aug[t] = n[t]
    }
  }

  q = rep(NA, S)
  for(t in 1:S){
    q[t] = dim(H[[t]])[2]
  }

  # input.max = apply(input[[level]], 2, max)
  # input.min = apply(input[[level]], 2, min)
  # Cl =  abs(input.max-input.min)
  Cl = hyper$Cl



  if(is.nugget){
    lnJacobian = sum(param[1:Dim]) + log(nugget) + log(nugget.UB - nugget) - 
                 log(abs(nugget.UB+(nugget.UB-1)*nugget))
    temp = sum(1/phi[1:Dim] * Cl)  + nugget 

  }else{
    lnJacobian = sum(param)
    temp = sum(1/phi * Cl) 
  }

  # if(is.nugget){
  #   lnJacobian = sum(log(param[1:Dim])) + log(nugget) + log(nugget.UB - nugget) - 
  #                log(abs(nugget.UB+(nugget.UB-1)*nugget))
  #   temp = sum(1/phi[1:Dim] * Cl) + nugget 

  # }else{
  #   lnJacobian = sum(log(param))
  #   temp = sum(1/phi * Cl) 
  # }
  
  lnJacobian = lnJacobian + al*log(temp) - bl*temp 

  #print(phi)

  t = level 


  if(level==1){
    R = buildcov(phi, dist, covmodel=cov.model, nugget=is.nugget)
    U = chol(R)
    RInv = chol2inv(U)
    X = H[[t]]
    Xm = Hm[[t]]

    X.aug = rbind(X, Xm)
    RInvX = RInv%*%X.aug
    HRH = t(X.aug)%*%RInvX 
    HRHchol = chol(HRH)
    HRHInv = chol2inv(HRHchol)

    # compute Q
    Q = RInv - RInvX%*%HRHInv%*%t(RInvX)

    S_2_log = 0
    for(k in 1:nsample){
      y.aug = rbind(y[[t]], as.matrix(ym[[t]][k, , ]))
      S_2_log = S_2_log + compute_S(y.aug, Q)
    }
    S_2_log = S_2_log / nsample
    
    # S_2_log = compute_S3D(y[[t]], ym[[t]], Q)  # too slow


    g = -N*sum(log(diag(U))) - N*sum(log(diag(HRHchol))) - 0.5*(n.aug[t]-q[t])*S_2_log

  }else if(level>1 & level<S){
    R = buildcov(phi, dist, covmodel=cov.model, nugget=is.nugget)
    U = chol(R)
    RInv = chol2inv(U)

    ## compute QH
    H.aug = rbind(H[[t]], Hm[[t]])
    RInvH = RInv%*%H.aug
    HRH = t(H.aug)%*%RInvH
    HRHchol = chol(HRH)
    HRHInv = chol2inv(HRHchol)

    K = RInvH%*%HRHInv%*%t(RInvH)
    #QH = RInv - K

    y_t1 = array(NA, dim=c(nsample, n[t], N))
    ym_t1 = array(NA, dim=c(nsample, nrow(Hm[[t]]), N))
    
    S_2_log_sum = 0
    for(k in 1:nsample){
      # IB = match.input(input[[t]], input.miss[[t-1]])$IB
      # y_t1 = as.matrix(ym[[t-1]][IB])
      y_t1[k, , ] = create.w(t=t, input=input, input.miss=input.miss[[t-1]], 
                  y=y[[t-1]], ym=as.matrix(ym[[t-1]][k, , ]))

      # IB = match.input(input.miss[[t]], input.miss[[t-1]])$IB
      # ym_t1 = as.matrix(ym[[t-1]][IB])
      ym_t1[k, , ] = create.w(t=t, input=input.miss, input.miss=input.miss[[t-1]], 
                  y=as.matrix(ym[[t-1]][k, , ]), ym=as.matrix(ym[[t-1]][k, , ]))


      y.aug = rbind(y[[t]], as.matrix(ym[[t]][k, , ]))
      y_t1.aug = rbind(as.matrix(y_t1[k, , ]), as.matrix(ym_t1[k, , ]))

      S_2_log_sum = S_2_log_sum + compute_S_sum(y.aug, H.aug, y_t1.aug, RInv, K)

    }
    S_2_log_sum = S_2_log_sum / nsample

    g = - N*sum(log(diag(U)))  - 0.5*S_2_log_sum


  }else{

    R = buildcov(phi, dist, covmodel=cov.model, nugget=is.nugget) 
    U = chol(R)
    RInv = chol2inv(U)
    RInvH = RInv%*%H[[t]]
    HRH = t(H[[t]])%*%RInvH
    HRHchol = chol(HRH)
    HRHInv = chol2inv(HRHchol)
    K = RInvH%*%HRHInv%*%t(RInvH)
    #QH = RInv - K 

    y_t1 = array(NA, dim=c(nsample, n[t], N))

    logdetX = matrix(0, nsample, N)
    S_2_log = matrix(0, nsample, N)

    S_2_log_sum = 0
    for(k in 1:nsample){
      y_t1[k, , ] = create.w(t=t, input=input, input.miss=input.miss[[t-1]], 
                  y=y[[t-1]], ym=as.matrix(ym[[t-1]][k, , ]))
      
      S_2_log_sum = S_2_log_sum + compute_S_sum(y[[t]], H[[t]], as.matrix(y_t1[k, ,]), RInv, K)
      # out = compute_S_sum2(y[[t]], H[[t]], y_t1[k, ,], RInv)

    }

    S_2_log_sum = S_2_log_sum / nsample

    g = -N*sum(log(diag(U))) - 0.5*S_2_log_sum






  }


  g = lnJacobian + g 

  return(g)

}
#############################################################################
#############################################################################


#############################################################################
#############################################################################
condsim.NN.univariate <- function(formula,output,input,input.new,phi,cov.model,
                    nsample){



Dim = dim(input[[1]])[2]
p.x = Dim
if(dim(phi)[1]==Dim){
  is.nugget=FALSE
}else{
  is.nugget=TRUE
}

###################################################################
#### augment input
###################################################################
S = length(output)   # number of code
out = augment.input(input)
input.union = out$union
input.miss = out$miss

np = dim(input.new)[1] 

n = sapply(output, length)
n.aug = rep(NA, S)

y = output

## add new inputs to missing inputs
# for(t in 1:(S-1)){
#   input.miss[[t]] = rbind(input.new, input.miss[[t]])
# }

pred.ID = list()
ID.org = list()
index.full = 1:np 
for(t in 1:(S-1)){
  ind.list = match.input(input.union[[t]], input.new)
  if(!is.null(ind.list$IA)){
      indA = ind.list$IA 
      indB = ind.list$IB 
      pred.ID.exist = indA 
      input.exist = input.union[[t]][indA, ,drop=FALSE]
      input.added = input.new[-indB, ,drop=FALSE]
      n.added = dim(input.added)[1]
      pred.ID.added = seq(1,n.added,by=1)
      ID.org[[t]] = c(index.full[-indB], indB)
      pred.ID[[t]] = c(pred.ID.added, pred.ID.exist+n.added)
      input.miss[[t]] = rbind(input.added, input.miss[[t]])
      input.union[[t]] = rbind(input.miss[[t]], input[[t]])
  }else{
    ID.org[[t]] = 1:np 
    pred.ID[[t]] = 1:np 
    input.miss[[t]] = rbind(input.new, input.miss[[t]])
    input.union[[t]] = rbind(input.new, input.union[[t]])
  }

}

t = S
ind.list = match.input(input[[t]], input.new)
if(!is.null(ind.list$IA)){
  indA = ind.list$IA
  indB = ind.list$IB
  pred.ID.exist = indA
  input.exist = input.union[[t]][indA, , drop=FALSE]
  input.added = input.new[-indB, , drop=FALSE]
  n.added = dim(input.added)[1]
  pred.ID.added = seq(1, n.added, by=1)
  ID.org[[t]] = c(index.full[-indB], indB)
  pred.ID[[t]] = c(pred.ID.added, pred.ID.exist+n.added)
  input.miss[[t]] = input.added
}else{
  input.miss[[t]] = input.new
  ID.org[[t]] = 1:np
  pred.ID[[t]] = 1:np
}
input.union[[S]] = rbind(input.new, input.union[[S]])


## add new inputs to union of inputs
# for(t in 1:S){
#   input.union[[t]] = rbind(input.new, input.union[[t]])
# }


###################################################################
#### create covariates
###################################################################
H = list()
Hm = list()
for(t in 1:S){
  colnames(input[[t]]) = paste0("x", 1:p.x)
  df = data.frame(input[[t]])
  H[[t]] = model.matrix(formula[[t]], df)

  colnames(input.miss[[t]]) = paste0("x", 1:p.x)
  df = data.frame(input.miss[[t]])
  Hm[[t]] = model.matrix(formula[[t]], df)
  
}



n.m = rep(NA, S)
for(t in 1:S){
  n.m[t] = dim(input.miss[[t]])[1]
}





###################################################################
#### compute intermediate quantities
###################################################################
dist.o = list()
dist.m = list()
dist.mo = list()
distlist = list()

for(t in 1:S){
  dist.o[[t]] = compute_distance(input[[t]], input[[t]])

  #if(t<S){
    dist.m[[t]] = compute_distance(input.miss[[t]], input.miss[[t]])
    dist.mo[[t]] = compute_distance(input.miss[[t]], input[[t]])
  #}

  distlist[[t]] = compute_distance(input.union[[t]], input.union[[t]])
}



##
ym.hat = list()
q = rep(NA, S)
for(t in 1:S){
  q[t] = dim(H[[t]])[2]
  ym.hat[[t]] = matrix(NA, nsample, dim(Hm[[t]])[1])
}

 

  #################################################################################
  #### Sampling from predictive distribution
  #################################################################################
  t = 1
  R = buildcov(phi[ ,t], dist.o[[t]], covmodel=cov.model, nugget=is.nugget)
  U = chol(R)
  RInv = chol2inv(U)
  #UH = backsolve(U, H[[t]], transpose=TRUE)
  #HRHInv = solve(crossprod(UH))
  HRHInv = solve(t(H[[t]])%*%RInv%*%H[[t]])
  betahat = HRHInv%*%t(H[[t]])%*%(RInv%*%y[[t]])
  res = y[[t]] - H[[t]]%*%betahat
  SSE = c(t(res)%*%RInv%*%res)

  Rm = buildcov(phi[ ,t], dist.m[[t]], covmodel=cov.model, nugget=FALSE)
  Rmo = buildcov(phi[ ,t], dist.mo[[t]], covmodel=cov.model, nugget=is.nugget)
  
  #RmoU = t(backsolve(U, t(Rmo), transpose=TRUE))

  XXRR = t(Hm[[t]]) - t(H[[t]])%*%RInv%*%t(Rmo)
  #XXRR = t(Hm[[t]]) - t(RmoU%*%UH)
  Sig_ymymy = Rm - Rmo%*%RInv%*%t(Rmo) + t(XXRR)%*%HRHInv%*%XXRR
  #Sig_ymymy = Rm - tcrossprod(RmoU) + t(XXRR)%*%HRHInv%*%XXRR
  Sig = Sig_ymymy * SSE/(n[t] - q[t])
  #Sig = Sig_ymymy * SSE/(n[t])

  mu_ymy = Hm[[t]]%*%betahat + Rmo%*%(RInv%*%res)
  
  ym.hat[[t]] = mvtnorm::rmvt(nsample, sigma=Sig, df=n[t]-q[t], delta=mu_ymy, type="shifted")


  
  for(t in 2:(S)){
    
    ############################################################################
    #### estimate missing data ym
    ############################################################################

    R = buildcov(phi[ ,t], dist.o[[t]], covmodel=cov.model, nugget=is.nugget)
    U = chol(R)
    RInv = chol2inv(U)

    Rm = buildcov(phi[ ,t], dist.m[[t]], covmodel=cov.model, nugget=FALSE)
    Rmo = buildcov(phi[ ,t], dist.mo[[t]], covmodel=cov.model, nugget=is.nugget)
    
    R_sk = Rm - Rmo%*%RInv%*%t(Rmo)



    for(k in 1:nsample){

      y_t1 = create.w.univ(t=t, input=input, input.miss=input.miss[[t-1]],
             y=y[[t-1]], ym=(ym.hat[[t-1]][k, ]))

      ym_t1 = create.w.new(t=t, input=input[[t-1]], input.miss=input.miss,
              y=y[[t-1]], ym=(ym.hat[[t-1]][k, ]))

      # ym.hat[[t]][k, ] = conditional_simulation(y[[t]], H[[t]], as.matrix(y_t1),  
      #               RInv, Hm[[t]], as.matrix(ym_t1), Rmo, R_sk)


      X = cbind(H[[t]], y_t1)
      XRXInv = solve(t(X)%*%RInv%*%X)
      betahat = XRXInv%*%t(X)%*%(RInv%*%y[[t]])
      res = y[[t]] - X%*%betahat
      SSE = c(t(res)%*%RInv%*%res)

      Xm = cbind(Hm[[t]], ym_t1)

      XXRR = t(Xm) - t(X)%*%RInv%*%t(Rmo)
      Sig_ymymy = R_sk + t(XXRR)%*%XRXInv%*%XXRR
      Sig = Sig_ymymy * SSE/(n[t]-q[t]-1)

      mu_ymy = Xm%*%betahat + Rmo%*%(RInv%*%res)

      ym.hat[[t]][k, ] = c(mvtnorm::rmvt(1, sigma=Sig, df=n[t]-q[t]-1, 
                          delta=mu_ymy, type="shifted"))

    }
  }

################################################################################

## get summary statistics
krige = list()
krigeSE = list()
krige.lower95 = list()
krige.upper95 = list()
for(t in 1:S){
  # yhat[[t]] = ym.hat[[t]][ ,pred.ID[[t]], ]

  krige[[t]] = apply(ym.hat[[t]], 2, mean)
  krigeSE[[t]] = apply(ym.hat[[t]], 2, sd)
  krige.lower95[[t]] = apply(ym.hat[[t]], 2, quantile, 0.025)
  krige.upper95[[t]] = apply(ym.hat[[t]], 2, quantile, 0.975)
}


pred.mu = list()
pred.SE = list()
pred.lower95 = list()
pred.upper95 = list()
for(t in 1:S){
  pred.mu[[t]] = rep(NA, np)
  pred.SE[[t]] = rep(0, np)
  pred.lower95[[t]] = rep(NA, np)
  pred.upper95[[t]] = rep(NA, np)

  ind.list = ismember(input.new, input.miss[[t]])
  pred.mu[[t]][ind.list$IIA] = krige[[t]][ind.list$IA]
  pred.SE[[t]][ind.list$IIA] = krigeSE[[t]][ind.list$IA]
  pred.lower95[[t]][ind.list$IIA] = krige.lower95[[t]][ind.list$IA]
  pred.upper95[[t]][ind.list$IIA] = krige.upper95[[t]][ind.list$IA]
  if(length(ind.list$IIA)<np){
    ind.input = ismember(input.new, input[[t]])
    pred.mu[[t]][ind.input$IIA] = y[[t]][ind.input$IA]
    pred.lower95[[t]][ind.input$IIA] = y[[t]][ind.input$IA]
    pred.upper95[[t]][ind.input$IIA] = y[[t]][ind.input$IA]
  }
}

names(pred.mu) = paste0("Level", seq(1:S), "")
names(pred.SE) = paste0("Level", seq(1:S), "")
names(pred.lower95) = paste0("Level", seq(1:S), "")
names(pred.upper95) = paste0("Level", seq(1:S), "")

pred = list(mu=pred.mu, SE=pred.SE, 
       lower95=pred.lower95,upper95=pred.upper95)

return(pred)

}


#############################################################################
#############################################################################


#############################################################################
#############################################################################
#### Functions in Multivariate Autoregressive Cokriging Models


#############################################################################
#############################################################################
margin.posterior.mv = function(param, input, output, level, H, dist, cov.model="matern_5_2",
  hyper=NULL){

Dim = dim(dist)[3]
p.x = Dim

t = level 
n = dim(output[[t]])[1] # number of model runs
N = dim(output[[t]])[2] # number of spatial lomessageions


if(is.null(hyper)){
  al = 0.5-p.x
  bl = 1 * dim(dist)[1]^(-1/p.x) * (al + p.x)
  nugget.UB = 1
}else{
  al = hyper$a
  bl = hyper$b * dim(dist)[1]^(-1/p.x) * (al + p.x)
  nugget.UB = hyper$nugget.UB
}

  # param contains phi and nugget (maybe)
  if(length(param)==Dim){ #no nugget
    phi = exp(-param)
    is.nugget=FALSE
  }else{
    phi = exp(-param[1:Dim]) 
    nugget = nugget.UB*exp(param[Dim+1]) / (1 + exp(param[Dim+1]))
    is.nugget = TRUE
    phi = c(phi, nugget)
  }

q = dim(H)[2]

#Cl = dim(dist)[1]^(-1/p.x)

# input.max = apply(input[[level]], 2, max)
# input.min = apply(input[[level]], 2, min)
# Cl =  abs(input.max-input.min)
Cl = hyper$Cl

  if(is.nugget){
    lnJacobian = sum(param[1:Dim]) + log(nugget) + log(nugget.UB - nugget) - 
                 log(abs(nugget.UB+(nugget.UB-1)*nugget))
    temp = sum(1/phi[1:Dim] * Cl)  + nugget 

  }else{
    lnJacobian = sum(param)
    temp = sum(1/phi * Cl) 
  }
  
  lnJacobian = lnJacobian + al*log(temp) - bl*temp 



if(level==1){
    R = buildcov(phi, dist, covmodel=cov.model, nugget=is.nugget)
    U = chol(R)
    RInv = chol2inv(U)
    
    X.aug = H
    HRH = t(X.aug)%*%RInv%*%X.aug #+ diag(1e-6,dim(X.aug)[2])
    HRHchol = chol(HRH)
    HRHInv = chol2inv(HRHchol)
    # compute Q 
    Q = RInv%*%(diag(n)-X.aug%*%HRHInv%*%t(X.aug)%*%RInv)
    S_2_log = compute_S(output[[t]], Q)

    g = -N*sum(log(diag(U))) - N*sum(log(diag(HRHchol))) - 0.5*(n-q)*S_2_log

}else{
    R = buildcov(phi, dist, covmodel=cov.model, nugget=is.nugget)
    U = chol(R)
    RInv = chol2inv(U)
    RInvH = RInv%*%H
    HRH = t(H)%*%RInvH
    HRHchol = chol(HRH)
    HRHInv = chol2inv(HRHchol)

    
    K = RInvH%*%HRHInv%*%t(RInvH)
    #QH = RInv - K 

    IB = match.input(input[[t]], input[[t-1]])$IB
    y_t1 = output[[t-1]][IB, ]

    S_2_log_sum = compute_S_sum(output[[t]], H, y_t1, RInv, K)

    g = -N*sum(log(diag(U))) - 0.5*S_2_log_sum 

}

  g = lnJacobian + g 


  return(g)
}

#############################################################################
#############################################################################


#############################################################################
#############################################################################
compute.g = function(param, input.list, level, y, H, ym, Hm, dist, hyper, cov.model="matern_5_2"){
  ## compute g function at fidelity level t
  ## 

  # param is a vector 
  
  Dim = dim(dist)[3]
  p.x = Dim
  N = ncol(y[[1]])
  nsample = dim(ym[[1]])[1] 

  if(is.null(hyper)){
      al = 0.5-p.x
      bl = dim(dist)[1]^(-1/p.x) * (al + p.x)
      nugget.UB = 1
  }else{
    al = hyper$a
    bl = hyper$b * dim(dist)[1]^(-1/p.x) * (al + p.x)
    nugget.UB = hyper$nugget.UB
  }

  # param contains phi and nugget (maybe)
  if(length(param)==Dim){ #no nugget
    phi = exp(-param)
    # phi = 1/param
    is.nugget=FALSE
  }else{
    phi = exp(-param[1:Dim]) 
    nugget = nugget.UB*exp(param[Dim+1]) / (1 + exp(param[Dim+1]))
    is.nugget = TRUE
    phi = c(phi, nugget)
  }

  # inputlist = augment.input(input)
  # input.miss = inputlist$miss
  # input.union = inputlist$union  
  input = input.list$input
  input.miss = input.list$input.miss
  
  S = length(y)
  n = rep(NA, S)
  n.aug = rep(NA, S)
  for(t in 1:S){
    n[t] = dim(y[[t]])[1]
    if(t<S){
      n.aug[t] = n[t] + dim(Hm[[t]])[1]
    }else{
      n.aug[t] = n[t]
    }
  }

  q = rep(NA, S)
  for(t in 1:S){
    q[t] = dim(H[[t]])[2]
  }

  # input.max = apply(input[[level]], 2, max)
  # input.min = apply(input[[level]], 2, min)
  # Cl =  abs(input.max-input.min)
  Cl = hyper$Cl



  if(is.nugget){
    lnJacobian = sum(param[1:Dim]) + log(nugget) + log(nugget.UB - nugget) - 
                 log(abs(nugget.UB+(nugget.UB-1)*nugget))
    temp = sum(1/phi[1:Dim] * Cl)  + nugget 

  }else{
    lnJacobian = sum(param)
    temp = sum(1/phi * Cl) 
  }

  # if(is.nugget){
  #   lnJacobian = sum(log(param[1:Dim])) + log(nugget) + log(nugget.UB - nugget) - 
  #                log(abs(nugget.UB+(nugget.UB-1)*nugget))
  #   temp = sum(1/phi[1:Dim] * Cl) + nugget 

  # }else{
  #   lnJacobian = sum(log(param))
  #   temp = sum(1/phi * Cl) 
  # }
  
  lnJacobian = lnJacobian + al*log(temp) - bl*temp 

  #print(phi)

  t = level 


  if(level==1){
    R = buildcov(phi, dist, covmodel=cov.model, nugget=is.nugget)
    U = chol(R)
    RInv = chol2inv(U)
    X = H[[t]]
    Xm = Hm[[t]]

    X.aug = rbind(X, Xm)
    RInvX = RInv%*%X.aug
    HRH = t(X.aug)%*%RInvX 
    HRHchol = chol(HRH)
    HRHInv = chol2inv(HRHchol)

    # compute Q
    Q = RInv - RInvX%*%HRHInv%*%t(RInvX)

    S_2_log = 0
    for(k in 1:nsample){
      y.aug = rbind(y[[t]], as.matrix(ym[[t]][k, , ]))
      S_2_log = S_2_log + compute_S(y.aug, Q)
    }
    S_2_log = S_2_log / nsample
    
    # S_2_log = compute_S3D(y[[t]], ym[[t]], Q)  # too slow


    g = -N*sum(log(diag(U))) - N*sum(log(diag(HRHchol))) - 0.5*(n.aug[t]-q[t])*S_2_log

  }else if(level>1 & level<S){
    R = buildcov(phi, dist, covmodel=cov.model, nugget=is.nugget)
    U = chol(R)
    RInv = chol2inv(U)

    ## compute QH
    H.aug = rbind(H[[t]], Hm[[t]])
    RInvH = RInv%*%H.aug
    HRH = t(H.aug)%*%RInvH
    HRHchol = chol(HRH)
    HRHInv = chol2inv(HRHchol)

    K = RInvH%*%HRHInv%*%t(RInvH)
    #QH = RInv - K

    y_t1 = array(NA, dim=c(nsample, n[t], N))
    ym_t1 = array(NA, dim=c(nsample, nrow(Hm[[t]]), N))
    
    S_2_log_sum = 0
    for(k in 1:nsample){
      # IB = match.input(input[[t]], input.miss[[t-1]])$IB
      # y_t1 = as.matrix(ym[[t-1]][IB])
      y_t1[k, , ] = create.w(t=t, input=input, input.miss=input.miss[[t-1]], 
                  y=y[[t-1]], ym=as.matrix(ym[[t-1]][k, , ]))

      # IB = match.input(input.miss[[t]], input.miss[[t-1]])$IB
      # ym_t1 = as.matrix(ym[[t-1]][IB])
      ym_t1[k, , ] = create.w(t=t, input=input.miss, input.miss=input.miss[[t-1]], 
                  y=as.matrix(ym[[t-1]][k, , ]), ym=as.matrix(ym[[t-1]][k, , ]))


      y.aug = rbind(y[[t]], as.matrix(ym[[t]][k, , ]))
      y_t1.aug = rbind(as.matrix(y_t1[k, , ]), as.matrix(ym_t1[k, , ]))

      S_2_log_sum = S_2_log_sum + compute_S_sum(y.aug, H.aug, y_t1.aug, RInv, K)

    }
    S_2_log_sum = S_2_log_sum / nsample

    g = - N*sum(log(diag(U)))  - 0.5*S_2_log_sum


  }else{

    R = buildcov(phi, dist, covmodel=cov.model, nugget=is.nugget) 
    U = chol(R)
    RInv = chol2inv(U)
    RInvH = RInv%*%H[[t]]
    HRH = t(H[[t]])%*%RInvH
    HRHchol = chol(HRH)
    HRHInv = chol2inv(HRHchol)
    K = RInvH%*%HRHInv%*%t(RInvH)
    #QH = RInv - K 

    y_t1 = array(NA, dim=c(nsample, n[t], N))

    logdetX = matrix(0, nsample, N)
    S_2_log = matrix(0, nsample, N)

    S_2_log_sum = 0
    for(k in 1:nsample){
      y_t1[k, , ] = create.w(t=t, input=input, input.miss=input.miss[[t-1]], 
                  y=y[[t-1]], ym=as.matrix(ym[[t-1]][k, , ]))
      
      S_2_log_sum = S_2_log_sum + compute_S_sum(y[[t]], H[[t]], y_t1[k, ,], RInv, K)
      # out = compute_S_sum2(y[[t]], H[[t]], y_t1[k, ,], RInv)

    }

    S_2_log_sum = S_2_log_sum / nsample

    g = -N*sum(log(diag(U))) - 0.5*S_2_log_sum


  }


  g = lnJacobian + g 

  return(g)

}

#############################################################################
#############################################################################



#############################################################################
#############################################################################
fit.ND=function(formula, output, input, phi, cov.model,
  prior, opt){

hyperparam = prior$hyper

Dim = dim(input[[1]])[2]
p.x = Dim
if(dim(phi)[1]==Dim){
  is.nugget=FALSE
}else{
  is.nugget=TRUE
}


S = length(output)

###################################################################
#### create covariates
###################################################################
H = list()
Hm = list()
for(t in 1:S){
  colnames(input[[t]]) = paste0("x", 1:p.x)
  df = data.frame(input[[t]])
  H[[t]] = model.matrix(formula[[t]], df)

}



###################################################################
#### compute intermediate quantities
###################################################################

distlist = list()
for(t in 1:S){
  distlist[[t]] = compute_distance(input[[t]], input[[t]])
}

Cl = list()
for(t in 1:S){
  input.max = apply(input[[t]], 2, max)
  input.min = apply(input[[t]], 2, min)
  Cl[[t]] =  abs(input.max-input.min)
}

for(t in 1:S){
  hyperparam[[t]]$Cl = Cl[[t]]
}
###################################################################
#### begin optimization algorithm
###################################################################

phi.new = phi 

for(t in 1:S){

  if(is.nugget){
    nu = log(phi[p.x+1, t]) - log(hyperparam[[t]]$nugget.UB-phi[p.x+1, t])   # logit of nugget
    init.val = c(-log(phi[1:p.x, t]), nu) 
  }else{
    init.val = -log(phi[ ,t])
  }

  fit = try(optim(init.val, margin.posterior.mv, input=input, output=output, level=t,
              H=H[[t]], dist=distlist[[t]], cov.model=cov.model,
              hyper=hyperparam[[t]], 
              control=list(fnscale=-1, maxit=opt$maxit),
              method=opt$method, lower=opt$lower, upper=opt$upper),
      silent=T)

  if(inherits(fit, "try-error")){
    phi.new[ ,t] = phi[ ,t]
    message("\n optimization error, skip t=", t, "\n")
    print(fit)
  }else{
    if(is.nugget){
      nugget = hyperparam[[t]]$nugget.UB*exp(fit$par[p.x+1]) / (1+exp(fit$par[p.x+1]))
      phi.new[ ,t] = c(exp(-fit$par[1:p.x]), nugget)
    }else{
      phi.new[ ,t] = exp(-fit$par)
    } 
  }

}

colnames(phi.new) = paste0("Level", seq(1:S), "")

return(list(par=phi.new))




}


#############################################################################
#############################################################################




#############################################################################
#############################################################################
fit.NN <- function(formula,output,input,phi,cov.model,prior,
               opt, MCEM){

hyperparam = prior$hyperparam

maxit = MCEM$maxit 
tol = MCEM$tol
n.sample = MCEM$n.sample
verbose = MCEM$verbose

Dim = dim(input[[1]])[2]
p.x = Dim
if(dim(phi)[1]==Dim){
  is.nugget=FALSE
}else{
  is.nugget=TRUE
}

###################################################################
#### augment input
###################################################################
S = length(output)   # number of code
out = augment.input(input)
input.union = out$union
input.miss = out$miss
input.list = list(input=input, input.miss=input.miss)

Cl = list()
for(t in 1:S){
  input.max = apply(input.list$input[[t]], 2, max)
  input.min = apply(input.list$input[[t]], 2, min)
  Cl[[t]] =  abs(input.max-input.min)
}

for(t in 1:S){
  hyperparam[[t]]$Cl = Cl[[t]]
}

###################################################################
#### create covariates
###################################################################
H = list()
Hm = list()
for(t in 1:S){
  colnames(input[[t]]) = paste0("x", 1:p.x)
  df = data.frame(input[[t]])
  H[[t]] = model.matrix(formula[[t]], df)

  if(t<S){
    colnames(input.miss[[t]]) = paste0("x", 1:p.x)
    df = data.frame(input.miss[[t]])
    Hm[[t]] = model.matrix(formula[[t]], df)
  }

}




###################################################################
#### compute intermediate quantities
###################################################################
dist.o = list()
dist.m = list()
dist.mo = list()
distlist = list()

for(t in 1:S){
  dist.o[[t]] = compute_distance(input[[t]], input[[t]])

  if(t<S){
    dist.m[[t]] = compute_distance(input.miss[[t]], input.miss[[t]])
    dist.mo[[t]] = compute_distance(input.miss[[t]], input[[t]])
  }

  distlist[[t]] = compute_distance(input.union[[t]], input.union[[t]])
}



###################################################################
#### begin MCEM algorithm
###################################################################
phi.new = phi
conv = FALSE
iter = 1
while(!conv){

  ###############################################################
  #### generate M Monte Carlo samples for missing data
  ###############################################################
  # y.m = list()
  #system.time(
  # for(k in 1:n.sample){
  #   y.m[[k]] = sample.ym(y=output,input=input,param=phi,Ho=H,Hm=Hm,dist.o=dist.o,
  #           dist.m=dist.m,dist.mo=dist.mo,cov.model=cov.model)
  # }
  y.m = sample.ym(y=output,input=input,param=phi,Ho=H,Hm=Hm,dist.o=dist.o,
            dist.m=dist.m,dist.mo=dist.mo,cov.model=cov.model, 
            nsample=n.sample)
  #)
  ###############################################################
  #### compute and maximize the Q function at each fidelity
  ###############################################################
  
  for(t in 1:S){

    if(is.nugget){
      nu = log(phi[p.x+1, t]) - log(hyperparam[[t]]$nugget.UB-phi[p.x+1, t])   # logit of nugget
      init.val = c(-log(phi[1:p.x, t]), nu) 
    }else{
      init.val = -log(phi[ ,t])
    }

    fit = try(optim(init.val, compute.g, input.list=input.list, level=t, y=output, H=H, ym=y.m, Hm=Hm,
                dist=distlist[[t]], hyper=hyperparam[[t]], cov.model=cov.model,
                  control=list(fnscale=-1, maxit=opt$maxit),
                  method=opt$method, lower=opt$lower, upper=opt$upper),
        silent=T)

    # fit = try(optimr(init.val, compute.Q.default, input=input, level=t, y=output, H=H, y.m=y.m, Hm=Hm,
    #           distlist=distlist, cov.model=cov.model,
    #             control=list(fnscale=-1, maxit=opt$maxit),
    #             method=opt$method, lower=opt$lower, upper=opt$upper),
    #   silent=T)

    if(inherits(fit, "try-error")){
      phi.new[ ,t] = phi[ ,t]
      message("\n optimization error, skip t=", t, "\n")
      print(fit)
    }else{
      if(is.nugget){
        nugget = hyperparam[[t]]$nugget.UB*exp(fit$par[p.x+1]) / (1+exp(fit$par[p.x+1]))
        phi.new[ ,t] = c(exp(-fit$par[1:p.x]), nugget)
      }else{
        phi.new[ ,t] = exp(-fit$par)
      } 
    }

  }

  ###############################################################
  #### check convergence
    if(inherits(fit, "try-error")){
      diff = tol + 1
    }else{
      diff = mean((phi.new - phi)^2)
    }

  if(verbose){
    message("iter=", iter, "\n")
  }
  


  if(iter>maxit){
    conv = TRUE
  }else{
    if(diff<tol){
      conv = TRUE
    }
  }

  phi = phi.new

  iter = iter + 1

}



colnames(phi) = paste0("Level", seq(1:S), "")

return(list(par=phi, eps=diff, iter=iter)) 

}

#############################################################################
#############################################################################




#############################################################################
#############################################################################
predict.ND = function(formula, output, input, input.new, phi, cov.model){



Dim = dim(input[[1]])[2]
p.x = Dim
if(dim(phi)[1]==Dim){
  is.nugget=FALSE
}else{
  is.nugget=TRUE
}

S = length(output)
n = rep(NA, S)
for(t in 1:S){
  n[t] = dim(output[[t]])[1]
}

np = dim(input.new)[1] 

y = output

## add new inputs to missing inputs
input.miss = list()
input.union = list()
for(t in 1:S){
  input.miss[[t]] = input.new
  input.union[[t]] = rbind(input.new, input[[t]])
}


###################################################################
#### create covariates
###################################################################
H = list()
Hm = list()
for(t in 1:S){
  colnames(input[[t]]) = paste0("x", 1:p.x)
  df = data.frame(input[[t]])
  H[[t]] = model.matrix(formula[[t]], df)

  
  colnames(input.miss[[t]]) = paste0("x", 1:p.x)
  df = data.frame(input.miss[[t]])
  Hm[[t]] = model.matrix(formula[[t]], df)
  
}


###################################################################
#### compute intermediate quantities
###################################################################
dist.o = list()
dist.m = list()
dist.mo = list()
distlist = list()

for(t in 1:S){
  dist.o[[t]] = compute_distance(input[[t]], input[[t]])

  #if(t<S){
    dist.m[[t]] = compute_distance(input.miss[[t]], input.miss[[t]])
    dist.mo[[t]] = compute_distance(input.miss[[t]], input[[t]])
  #}

  distlist[[t]] = compute_distance(input.union[[t]], input.union[[t]])
}


ym.hat = list()
q = rep(NA, S)
for(t in 1:S){
  q[t] = dim(H[[t]])[2]
}

krige = list()
krige.var = list()

  #################################################################################
  #### Get predictive mean and predictive variance
  #################################################################################
  t = 1
  R = buildcov(phi[ ,t], dist.o[[t]], covmodel=cov.model, nugget=is.nugget)
  U = chol(R)
  RInv = chol2inv(U)
  #UH = backsolve(U, H[[t]], transpose=TRUE)
  #HRHInv = solve(crossprod(UH))
  HRHInv = solve(t(H[[t]])%*%RInv%*%H[[t]])

  Rm = buildcov(phi[ ,t], dist.m[[t]], covmodel=cov.model, nugget=FALSE)
  Rmo = buildcov(phi[ ,t], dist.mo[[t]], covmodel=cov.model, nugget=is.nugget)
  
  # compute predictive mean
  RmoRInv = Rmo%*%RInv 
  KW = (Hm[[t]] - RmoRInv%*%H[[t]])%*%HRHInv%*%t(H[[t]])%*%RInv + RmoRInv
  krige[[t]] = KW%*%y[[t]]

  # compute predictive variance
  XXRR = t(Hm[[t]]) - t(H[[t]])%*%RInv%*%t(Rmo)
  c_star = diag(Rm - Rmo%*%RInv%*%t(Rmo) + t(XXRR)%*%HRHInv%*%XXRR)
  c_star[c_star<0] = 0  # avoid numerical problems

  Q = RInv - (RInv%*%H[[t]])%*%HRHInv%*%t(RInv%*%H[[t]])
  sigma2.hat = compute_Svec(y[[t]], Q) / (n[t]-q[t])

  constant = (n[t]-q[t])/(n[t]-q[t]-2)
  krige.var[[t]] =  constant*c_star %*% t(sigma2.hat)  # m-by-N matrix


  for(t in 2:S){
  
    ############################################################################
    #### estimate missing data ym
    ############################################################################

    R = buildcov(phi[ ,t], dist.o[[t]], covmodel=cov.model, nugget=is.nugget)
    U = chol(R)
    RInv = chol2inv(U)
    # RInvH = RInv%*%H[[t]]
    # HRH = t(H[[t]])%*%RInvH
    # HRHchol = chol(HRH)
    # HRHInv = chol2inv(HRHchol)

    # K = RInvH%*%HRHInv%*%t(RInvH)
    # QH = RInv - K 


    IB = match.input(input[[t]], input[[t-1]])$IB
    y_t1 = y[[t-1]][IB, ]
    IB = match.input(input.miss[[t]], input.miss[[t-1]])$IB
    ym_t1 = krige[[t-1]][IB, ]
    # Xm = cbind(Hm[[t]], ym_t1)

    
    Rm = buildcov(phi[ ,t], dist.m[[t]], covmodel=cov.model, nugget=FALSE)
    Rmo = buildcov(phi[ ,t], dist.mo[[t]], covmodel=cov.model, nugget=is.nugget)
    
    R_sk_diag = diag(Rm - Rmo%*%RInv%*%t(Rmo)) 
    R_sk_diag[R_sk_diag<0] = 0 

    pred.list = compute_prediction(y[[t]], H[[t]], y_t1, krige[[t-1]],krige.var[[t-1]],
              RInv, Hm[[t]], ym_t1, Rmo, R_sk_diag)

    krige[[t]] = pred.list$krige
    krige.var[[t]] = pred.list$krige.var

  }

################################################################################

# for(t in 1:(S-1)){
#   ind = sort(ID.org[[t]], index.return=TRUE)$ix 
#   krige[[t]][ind] = krige[[t]]
#   krigeSE[[t]][ind] = krigeSE[[t]]
# }

krigeSE = list()
lower95 = list()
upper95 = list()
for(t in 1:S){
  krige.var[[t]][krige.var[[t]]<0] = 0
  krigeSE[[t]] = sqrt(krige.var[[t]])
  # degree = ifelse(t==1, n[t]-q[t], n[t]-q[t]-1)
  lower95[[t]] = krige[[t]] - 2*krigeSE[[t]]
  upper95[[t]] = krige[[t]] + 2*krigeSE[[t]]
}

names(krige) = paste0("Level", seq(1:S), "")
names(krigeSE) = paste0("Level", seq(1:S), "")
names(lower95) = paste0("Level", seq(1:S), "")
names(upper95) = paste0("Level", seq(1:S), "")


out = list(mu=krige, SE=krigeSE, lower95=lower95, upper95=upper95)




return(out)



}



#############################################################################
#############################################################################



#############################################################################
#############################################################################
predict.NN <- function(formula,output,input,input.new,phi,cov.model="matern_5_2",
                    nsample=30){


  
  
  Dim = dim(input[[1]])[2]
  N = dim(output[[1]])[2]
  
  p.x = Dim
  if(dim(phi)[1]==Dim){
    is.nugget=FALSE
  }else{
    is.nugget=TRUE
  }
  
  ###################################################################
  #### augment input
  ###################################################################
  S = length(output)   # number of code
  out = augment.input(input)
  input.union = out$union
  input.miss = out$miss
  
  np = dim(input.new)[1] 
  
  n = rep(NA, S)
  for(t in 1:S){
    n[t] = dim(input[[t]])[1]
  }
  
  y = output
  
  ## add new inputs to missing inputs
  # for(t in 1:(S-1)){
  #   input.miss[[t]] = rbind(input.new, input.miss[[t]])
  # }
  
  pred.ID = list()
  ID.org = list()
  index.full = 1:np 
  for(t in 1:(S-1)){
    ind.list = match.input(input.union[[t]], input.new)
    if(!is.null(ind.list$IA)){
      indA = ind.list$IA 
      indB = ind.list$IB 
      pred.ID.exist = indA 
      input.exist = input.union[[t]][indA, ,drop=FALSE]
      input.added = input.new[-indB, ,drop=FALSE]
      n.added = dim(input.added)[1]
      pred.ID.added = seq(1,n.added,by=1)
      ID.org[[t]] = c(index.full[-indB], indB)
      pred.ID[[t]] = c(pred.ID.added, pred.ID.exist+n.added)
      input.miss[[t]] = rbind(input.added, input.miss[[t]])
      input.union[[t]] = rbind(input.miss[[t]], input[[t]])
    }else{
      ID.org[[t]] = 1:np 
      pred.ID[[t]] = 1:np 
      input.miss[[t]] = rbind(input.new, input.miss[[t]])
      input.union[[t]] = rbind(input.new, input.union[[t]])
    }
    
  }
  
  t = S
  ind.list = match.input(input[[t]], input.new)
  if(!is.null(ind.list$IA)){
    indA = ind.list$IA
    indB = ind.list$IB
    pred.ID.exist = indA
    input.exist = input.union[[t]][indA, , drop=FALSE]
    input.added = input.new[-indB, , drop=FALSE]
    n.added = dim(input.added)[1]
    pred.ID.added = seq(1, n.added, by=1)
    ID.org[[t]] = c(index.full[-indB], indB)
    pred.ID[[t]] = c(pred.ID.added, pred.ID.exist+n.added)
    input.miss[[t]] = input.added
  }else{
    input.miss[[t]] = input.new
    ID.org[[t]] = 1:np
    pred.ID[[t]] = 1:np
  }
  input.union[[S]] = rbind(input.new, input.union[[S]])
  
  
  
  ###################################################################
  #### create covariates
  ###################################################################
  H = list()
  Hm = list()
  for(t in 1:S){
    colnames(input[[t]]) = paste0("x", 1:p.x)
    df = data.frame(input[[t]])
    H[[t]] = model.matrix(formula[[t]], df)
    
    colnames(input.miss[[t]]) = paste0("x", 1:p.x)
    df = data.frame(input.miss[[t]])
    Hm[[t]] = model.matrix(formula[[t]], df)
    
  }
  
  
  n.m = rep(NA, S)
  for(t in 1:S){
    n.m[t] = dim(input.miss[[t]])[1]
  }
  
  
  
  
  ###################################################################
  #### compute intermediate quantities
  ###################################################################
  dist.o = list()
  dist.m = list()
  dist.mo = list()
  distlist = list()
  
  for(t in 1:S){
    dist.o[[t]] = compute_distance(input[[t]], input[[t]])
    
    #if(t<S){
    dist.m[[t]] = compute_distance(input.miss[[t]], input.miss[[t]])
    dist.mo[[t]] = compute_distance(input.miss[[t]], input[[t]])
    #}
    
    distlist[[t]] = compute_distance(input.union[[t]], input.union[[t]])
  }
  
  
  
  
  ##
  krige = list()
  krigeSE = list()
  ym.hat = list()
  q = rep(NA, S)
  for(t in 1:S){
    q[t] = dim(H[[t]])[2]
    ym.hat[[t]] = array(NA, dim=c(nsample, dim(Hm[[t]])[1], N))
  }
  
  
  
  
  #################################################################################
  #### get predictive mean and predictive variance
  #################################################################################
  t = 1
  R = buildcov(phi[ ,t], dist.o[[t]], covmodel=cov.model, nugget=is.nugget)
  U = chol(R)
  RInv = chol2inv(U)
  HRHInv = solve(t(H[[t]])%*%RInv%*%H[[t]])
  
  Rm = buildcov(phi[ ,t], dist.m[[t]], covmodel=cov.model, nugget=FALSE)
  Rmo = buildcov(phi[ ,t], dist.mo[[t]], covmodel=cov.model, nugget=is.nugget)
  
  
  # compute conditional mean
  RmoRInv = Rmo%*%RInv
  KW = (Hm[[t]]-RmoRInv%*%H[[t]]) %*% HRHInv %*% t(H[[t]]) %*% RInv + RmoRInv
  mu_ymy = KW %*% y[[t]]  # n.m-by-N matrix
  
  # compute predictive variance
  XXRR = t(Hm[[t]]) - t(H[[t]])%*%RInv%*%t(Rmo)
  c_star = Rm - Rmo%*%RInv%*%t(Rmo) + t(XXRR)%*%HRHInv%*%XXRR
  ## c_star is not positive definite!!! So, there is no cholesky decomposition available
  
  # compute S2
  Q = RInv - RInv%*%H[[t]]%*%HRHInv%*%t(H[[t]])%*%t(RInv)
  sigma2.hat = compute_Svec(y[[t]], Q) / (n[t]-q[t])    
  
  
  L = t(chol(c_star))
  ym.hat[[t]] = sample_mvt(mu_ymy, L=L, sigma=sigma2.hat, df=n[t]-q[t], nsample)
  
  
  #for(k in 1:nsample){
  for(t in 2:(S)){
    
    ############################################################################
    #### simulating missing data ym
    ############################################################################
    
    R = buildcov(phi[ ,t], dist.o[[t]], covmodel=cov.model, nugget=is.nugget)
    U = chol(R)
    RInv = chol2inv(U)
    
    
    Rm = buildcov(phi[ ,t], dist.m[[t]], covmodel=cov.model, nugget=FALSE)
    Rmo = buildcov(phi[ ,t], dist.mo[[t]], covmodel=cov.model, nugget=is.nugget)
    R_sk = Rm - Rmo%*%RInv%*%t(Rmo)
    # RmoRInv = Rmo%*%RInv 
    
    
    y_t1 = array(NA, dim=c(nsample, n[t], N))
    ym_t1 = array(NA, dim=c(nsample, n.m[t], N))
    # IB = match.input(input[[t]], input.miss[[t-1]])$IB
    # y_t1 = ym.hat[[t-1]][IB]
    
    
    # b = matrix(NA, q[t]+1, N)
    # mu_y = matrix(NA, n.m[t], N)
    for(k in 1:nsample){
      
      y_t1[k, , ] = create.w(t=t, input=input, input.miss=input.miss[[t-1]],
                             y=y[[t-1]], ym=as.matrix(ym.hat[[t-1]][k, , ]))
      
      ym_t1[k, , ] = create.w.pred(t=t, input=input[[t-1]], input.miss=input.miss,
                                   y=y[[t-1]], ym=as.matrix(ym.hat[[t-1]][k, ,]))
      
      ym.hat[[t]][k, , ] = conditional_simulation(y[[t]], H[[t]], y_t1[k, , ],  
                                                  RInv, Hm[[t]], ym_t1[k, , ], Rmo, R_sk)
      
    }
    
    
    
  }
  #}
  
  ################################################################################
  
  krige = list()
  krigeSE = list()
  krige.lower95 = list()
  krige.upper95 = list()
  for(t in 1:S){
    # yhat[[t]] = ym.hat[[t]][ ,pred.ID[[t]], ]
    
    krige[[t]] = apply(ym.hat[[t]], c(2,3), mean)
    krigeSE[[t]] = apply(ym.hat[[t]], c(2,3), sd)
    krige.lower95[[t]] = apply(ym.hat[[t]], c(2,3), quantile, 0.025)
    krige.upper95[[t]] = apply(ym.hat[[t]], c(2,3), quantile, 0.975)
  }
  
  pred.mu = list()
  pred.SE = list()
  pred.lower95 = list()
  pred.upper95 = list()
  for(t in 1:S){
    pred.mu[[t]] = matrix(NA, np, N)
    pred.SE[[t]] = matrix(0, np, N)
    pred.lower95[[t]] = matrix(NA, np, N)
    pred.upper95[[t]] = matrix(NA, np, N)
    
    ind.list = ismember(input.new, input.miss[[t]])
    pred.mu[[t]][ind.list$IIA, ] = krige[[t]][ind.list$IA, ]
    pred.SE[[t]][ind.list$IIA, ] = krigeSE[[t]][ind.list$IA, ]
    pred.lower95[[t]][ind.list$IIA, ] = krige.lower95[[t]][ind.list$IA, ]
    pred.upper95[[t]][ind.list$IIA, ] = krige.upper95[[t]][ind.list$IA, ]
    if(length(ind.list$IIA)<np){
      ind.input = ismember(input.new, input[[t]])
      pred.mu[[t]][ind.input$IIA, ] = y[[t]][ind.input$IA, ]
      pred.lower95[[t]][ind.input$IIA, ] = y[[t]][ind.input$IA, ]
      pred.upper95[[t]][ind.input$IIA, ] = y[[t]][ind.input$IA, ]
    }
  }
  
  names(pred.mu) = paste0("Level", seq(1:S), "")
  names(pred.SE) = paste0("Level", seq(1:S), "")
  names(pred.lower95) = paste0("Level", seq(1:S), "")
  names(pred.upper95) = paste0("Level", seq(1:S), "")
  
  pred = list(mu=pred.mu, SE=pred.SE, 
              lower95=pred.lower95,upper95=pred.upper95)
  
  return(pred)
  
}

#############################################################################
#############################################################################



#############################################################################
#############################################################################
condsim.ND = function(formula, output, input, input.new, phi, cov.model="matern_5_2", 
              nsample=30){


Dim = dim(input[[1]])[2]
N = dim(output[[1]])[2]

p.x = Dim
if(dim(phi)[1]==Dim){
  is.nugget=FALSE
}else{
  is.nugget=TRUE
}

S = length(output)
np = dim(input.new)[1] 

n = rep(NA, S)
for(t in 1:S){
  n[t] = dim(input[[t]])[1]
}

y = output

## add new inputs to missing inputs
# input.miss = list()
# input.union = list()
# for(t in 1:S){
#   input.miss[[t]] = input.new
#   input.union[[t]] = rbind(input.new, input[[t]])
# }
# input.miss = list()
# for(t in 1:S){
#   input.miss[[t]] = input.new
# }

input.miss = list()
input.union = list()
pred.ID = list()
ID.orig = list()
index.full = 1:np

indB = list()
for(t in 1:(S)){
  ind.list = match.input(input[[t]], input.new)
  # indlist = ismember(input.new, input[[t]])
  if(!is.null(ind.list$IA)){
    indA = ind.list$IA
    indB[[t]] = ind.list$IB
    pred.ID.exist = indA
    input.exist = input[[t]][indA, ,drop=FALSE] # common inputs 
    input.added = input.new[-indB[[t]], ,drop=FALSE] # inputs in input.new but not in input[[t]]
    n.added = dim(input.added)[1]
    pred.ID.added = seq(1, n.added, by=1)
    ID.orig[[t]] = c(index.full[-indB[[t]]], indB[[t]])
    pred.ID[[t]] = c(pred.ID.added, pred.ID.exist+n.added)
    input.miss[[t]] = input.added
    input.union[[t]] = rbind(input.added, input[[t]])
  }else{
    ID.orig[[t]] = 1:np
    pred.ID[[t]] = 1:np
    input.miss[[t]] = input.new
    input.union[[t]] = rbind(input.new, input[[t]])
  }
}


# 
n.m = rep(NA, S)
for(t in 1:S){
  n.m[t] = dim(input.miss[[t]])[1]
}


###################################################################
#### create covariates
###################################################################
H = list()
Hm = list()
for(t in 1:S){
  colnames(input[[t]]) = paste0("x", 1:p.x)
  df = data.frame(input[[t]])
  H[[t]] = model.matrix(formula[[t]], df)

  
  colnames(input.miss[[t]]) = paste0("x", 1:p.x)
  df = data.frame(input.miss[[t]])
  Hm[[t]] = model.matrix(formula[[t]], df)
  
}


###################################################################
#### compute intermediate quantities
###################################################################
dist.o = list()
dist.m = list()
dist.mo = list()
distlist = list()

for(t in 1:S){
  dist.o[[t]] = compute_distance(input[[t]], input[[t]])

  #if(t<S){
    dist.m[[t]] = compute_distance(input.miss[[t]], input.miss[[t]])
    dist.mo[[t]] = compute_distance(input.miss[[t]], input[[t]])
  #}

  distlist[[t]] = compute_distance(input.union[[t]], input.union[[t]])
}




ym.hat = list()
q = rep(NA, S)
for(t in 1:S){
  q[t] = dim(H[[t]])[2]
  ym.hat[[t]] = array(NA, dim=c(nsample, dim(Hm[[t]])[1], N))
}


  #################################################################################
  #### Sampling from predictive distribution
  #################################################################################
  t = 1
  R = buildcov(phi[ ,t], dist.o[[t]], covmodel=cov.model, nugget=is.nugget)
  U = chol(R)
  RInv = chol2inv(U)
  HRHInv = solve(t(H[[t]])%*%RInv%*%H[[t]])

  Rm = buildcov(phi[ ,t], dist.m[[t]], covmodel=cov.model, nugget=FALSE)
  Rmo = buildcov(phi[ ,t], dist.mo[[t]], covmodel=cov.model, nugget=is.nugget)
  
  # compute conditional mean
  RmoRInv = Rmo%*%RInv
  KW = (Hm[[t]]-RmoRInv%*%H[[t]]) %*% HRHInv %*% t(H[[t]]) %*% RInv + RmoRInv
  mu_ymy = KW %*% y[[t]]

  RmoU = t(backsolve(U, t(Rmo), transpose=TRUE))

  # compute predictive variance
  XXRR = t(Hm[[t]]) - t(H[[t]])%*%RInv%*%t(Rmo)
  c_star = Rm - Rmo%*%RInv%*%t(Rmo) + t(XXRR)%*%HRHInv%*%XXRR
  ## c_star is not positive definite!!! So, there is no cholesky decomposition available

  # compute S2
  Q = RInv - RInv%*%H[[t]]%*%HRHInv%*%t(H[[t]])%*%t(RInv)
  sigma2.hat = compute_Svec(y[[t]], Q) / (n[t]-q[t]) 

  # for(j in 1:N){
  #   Sig = c_star * sigma2.hat[j] 
  #   ym.hat[[t]][ , ,j] = mvtnorm::rmvt(nsample, sigma=Sig, df=n[t]-q[t], delta=mu_ymy[ ,j], type="shifted")
  # }

  L = t(chol(c_star))
  ym.hat[[t]] = sample_mvt(mu_ymy, L=L, sigma=sigma2.hat, df=n[t]-q[t], nsample)


  # ym.hatC = array(NA, dim=c(nsample, dim(Hm[[t]])[1], N))


  for(t in 2:(S)){
  
    ############################################################################
    #### estimate missing data ym
    ############################################################################

    R = buildcov(phi[ ,t], dist.o[[t]], covmodel=cov.model, nugget=is.nugget)
    U = chol(R)
    RInv = chol2inv(U)
    IB = match.input(input[[t]], input[[t-1]])$IB
    y_t1 = y[[t-1]][IB, ]


    # for(k in 1:nsample){
    #   ym_t1[k, ,] = create.w.pred(t=t, input=input[[t-1]], input.miss=input.miss,
    #                 y=y[[t-1]], ym=ym.hat[k, , ])
    # }
    ym_t1 = array(NA, dim=c(nsample, n.m[t], N))

    Rm = buildcov(phi[ ,t], dist.m[[t]], covmodel=cov.model, nugget=FALSE)
    Rmo = buildcov(phi[ ,t], dist.mo[[t]], covmodel=cov.model, nugget=is.nugget)
    
    R_sk = Rm - Rmo%*%RInv%*%t(Rmo)
    RmoRInv = Rmo%*%RInv

    # b = matrix(NA, q[t]+1, N)
    # mu_y = matrix(NA, nrow(Hm[[t]]), N)
    # ts = proc.time()
    for(k in 1:nsample){

      ym_t1[k, ,] = create.w.pred(t=t, input=input[[t-1]], input.miss=input.miss,
                    y=y[[t-1]], ym=as.matrix(ym.hat[[t-1]][k, , ]))
      # ym.hatC = array(NA, dim=c(nsample, dim(Hm[[t]])[1], N))

      ym.hat[[t]][k, , ] = conditional_simulation(y[[t]], H[[t]], y_t1,  
                    RInv, Hm[[t]], ym_t1[k, , ], Rmo, R_sk)
      # ym.hatC[[t]][k, , ] = pred.list$ym

      # sigma = pred.list$sigma

      # for(j in 1:N){
      #   X = cbind(H[[t]], y_t1[ ,j])
      #   XRXInv = solve(t(X)%*%RInv%*%X)
      #   b[ ,j] = XRXInv%*%t(X)%*%RInv%*%y[[t]][ ,j]
      #   Xp = cbind(Hm[[t]], ym_t1[k, ,j])
      #   mu_y[ ,j] = Xp%*%b[ ,j] + RmoRInv%*%(y[[t]][ ,j]-X%*%b[ ,j])
      #   temp = Xp - RmoRInv%*%X
      #   c_star = R_sk + temp%*%XRXInv%*%t(temp)

      #   ym.hat[[t]][k, ,j] = mvtnorm::rmvt(1, sigma=sigma[j]*c_star, df=n[t]-q[t]-1, 
      #                       delta=mu_y[ ,j], type="shifted")   
      # }

    }
    # te = proc.time() - ts


  }



## get summary statistics
krige = list()
krigeSE = list()
krige.lower95 = list()
krige.upper95 = list()
for(t in 1:S){
  # yhat[[t]] = ym.hat[[t]][ ,pred.ID[[t]], ]

  krige[[t]] = apply(ym.hat[[t]], c(2,3), mean)
  krigeSE[[t]] = apply(ym.hat[[t]], c(2,3), sd)
  krige.lower95[[t]] = apply(ym.hat[[t]], c(2,3), quantile, 0.025)
  krige.upper95[[t]] = apply(ym.hat[[t]], c(2,3), quantile, 0.975)
}


pred.mu = list()
pred.SE = list()
pred.lower95 = list()
pred.upper95 = list()
for(t in 1:S){
  pred.mu[[t]] = matrix(NA, np, N)
  pred.SE[[t]] = matrix(0, np, N)
  pred.lower95[[t]] = matrix(NA, np, N)
  pred.upper95[[t]] = matrix(NA, np, N)

  ind.list = ismember(input.new, input.miss[[t]])
  pred.mu[[t]][ind.list$IIA, ] = krige[[t]][ind.list$IA, ]
  pred.SE[[t]][ind.list$IIA, ] = krigeSE[[t]][ind.list$IA, ]
  pred.lower95[[t]][ind.list$IIA, ] = krige.lower95[[t]][ind.list$IA, ]
  pred.upper95[[t]][ind.list$IIA, ] = krige.upper95[[t]][ind.list$IA, ]
  if(length(ind.list$IIA)<np){
    ind.input = ismember(input.new, input[[t]])
    pred.mu[[t]][ind.input$IIA, ] = y[[t]][ind.input$IA, ]
    pred.lower95[[t]][ind.input$IIA, ] = y[[t]][ind.input$IA, ]
    pred.upper95[[t]][ind.input$IIA, ] = y[[t]][ind.input$IA, ]
  }
}

names(pred.mu) = paste0("Level", seq(1:S), "")
names(pred.SE) = paste0("Level", seq(1:S), "")
names(pred.lower95) = paste0("Level", seq(1:S), "")
names(pred.upper95) = paste0("Level", seq(1:S), "")

pred = list(mu=pred.mu, SE=pred.SE, 
       lower95=pred.lower95,upper95=pred.upper95)

return(pred)

}

#############################################################################
#############################################################################



#############################################################################
#############################################################################
condsim.NN <- function(formula,output,input,input.new,phi,cov.model="matern_5_2",
              nsample=30){



Dim = dim(input[[1]])[2]
N = dim(output[[1]])[2]

p.x = Dim
if(dim(phi)[1]==Dim){
  is.nugget=FALSE
}else{
  is.nugget=TRUE
}

###################################################################
#### augment input
###################################################################
S = length(output)   # number of code
out = augment.input(input)
input.union = out$union
input.miss = out$miss

np = dim(input.new)[1] 

n = rep(NA, S)
for(t in 1:S){
  n[t] = dim(input[[t]])[1]
}

y = output

## add new inputs to missing inputs
# for(t in 1:(S-1)){
#   input.miss[[t]] = rbind(input.new, input.miss[[t]])
# }

pred.ID = list()
ID.org = list()
index.full = 1:np 
for(t in 1:(S-1)){
  ind.list = match.input(input.union[[t]], input.new)
  if(!is.null(ind.list$IA)){
      indA = ind.list$IA 
      indB = ind.list$IB 
      pred.ID.exist = indA 
      input.exist = input.union[[t]][indA, ,drop=FALSE]
      input.added = input.new[-indB, ,drop=FALSE]
      n.added = dim(input.added)[1]
      pred.ID.added = seq(1,n.added,by=1)
      ID.org[[t]] = c(index.full[-indB], indB)
      pred.ID[[t]] = c(pred.ID.added, pred.ID.exist+n.added)
      input.miss[[t]] = rbind(input.added, input.miss[[t]])
      input.union[[t]] = rbind(input.miss[[t]], input[[t]])
  }else{
    ID.org[[t]] = 1:np 
    pred.ID[[t]] = 1:np 
    input.miss[[t]] = rbind(input.new, input.miss[[t]])
    input.union[[t]] = rbind(input.new, input.union[[t]])
  }

}

t = S
ind.list = match.input(input[[t]], input.new)
if(!is.null(ind.list$IA)){
  indA = ind.list$IA
  indB = ind.list$IB
  pred.ID.exist = indA
  input.exist = input.union[[t]][indA, , drop=FALSE]
  input.added = input.new[-indB, , drop=FALSE]
  n.added = dim(input.added)[1]
  pred.ID.added = seq(1, n.added, by=1)
  ID.org[[t]] = c(index.full[-indB], indB)
  pred.ID[[t]] = c(pred.ID.added, pred.ID.exist+n.added)
  input.miss[[t]] = input.added
}else{
  input.miss[[t]] = input.new
  ID.org[[t]] = 1:np
  pred.ID[[t]] = 1:np
}
input.union[[S]] = rbind(input.new, input.union[[S]])



###################################################################
#### create covariates
###################################################################
H = list()
Hm = list()
for(t in 1:S){
  colnames(input[[t]]) = paste0("x", 1:p.x)
  df = data.frame(input[[t]])
  H[[t]] = model.matrix(formula[[t]], df)

  colnames(input.miss[[t]]) = paste0("x", 1:p.x)
  df = data.frame(input.miss[[t]])
  Hm[[t]] = model.matrix(formula[[t]], df)
  
}


n.m = rep(NA, S)
for(t in 1:S){
  n.m[t] = dim(input.miss[[t]])[1]
}




###################################################################
#### compute intermediate quantities
###################################################################
dist.o = list()
dist.m = list()
dist.mo = list()
distlist = list()

for(t in 1:S){
  dist.o[[t]] = compute_distance(input[[t]], input[[t]])

  #if(t<S){
    dist.m[[t]] = compute_distance(input.miss[[t]], input.miss[[t]])
    dist.mo[[t]] = compute_distance(input.miss[[t]], input[[t]])
  #}

  distlist[[t]] = compute_distance(input.union[[t]], input.union[[t]])
}




##
krige = list()
krigeSE = list()
ym.hat = list()
q = rep(NA, S)
for(t in 1:S){
  q[t] = dim(H[[t]])[2]
  ym.hat[[t]] = array(NA, dim=c(nsample, dim(Hm[[t]])[1], N))
}


 

  #################################################################################
  #### get predictive mean and predictive variance
  #################################################################################
  t = 1
  R = buildcov(phi[ ,t], dist.o[[t]], covmodel=cov.model, nugget=is.nugget)
  U = chol(R)
  RInv = chol2inv(U)
  HRHInv = solve(t(H[[t]])%*%RInv%*%H[[t]])

  Rm = buildcov(phi[ ,t], dist.m[[t]], covmodel=cov.model, nugget=FALSE)
  Rmo = buildcov(phi[ ,t], dist.mo[[t]], covmodel=cov.model, nugget=is.nugget)
  

  # compute conditional mean
  RmoRInv = Rmo%*%RInv
  KW = (Hm[[t]]-RmoRInv%*%H[[t]]) %*% HRHInv %*% t(H[[t]]) %*% RInv + RmoRInv
  mu_ymy = KW %*% y[[t]]  # n.m-by-N matrix

  # compute predictive variance
  XXRR = t(Hm[[t]]) - t(H[[t]])%*%RInv%*%t(Rmo)
  c_star = Rm - Rmo%*%RInv%*%t(Rmo) + t(XXRR)%*%HRHInv%*%XXRR
  ## c_star is not positive definite!!! So, there is no cholesky decomposition available

  # compute S2
  Q = RInv - RInv%*%H[[t]]%*%HRHInv%*%t(H[[t]])%*%t(RInv)
  sigma2.hat = compute_Svec(y[[t]], Q) / (n[t]-q[t])    


  L = t(chol(c_star))
  ym.hat[[t]] = sample_mvt(mu_ymy, L=L, sigma=sigma2.hat, df=n[t]-q[t], nsample)


  #for(k in 1:nsample){
    for(t in 2:(S)){
      
      ############################################################################
      #### simulating missing data ym
      ############################################################################

      R = buildcov(phi[ ,t], dist.o[[t]], covmodel=cov.model, nugget=is.nugget)
      U = chol(R)
      RInv = chol2inv(U)


      Rm = buildcov(phi[ ,t], dist.m[[t]], covmodel=cov.model, nugget=FALSE)
      Rmo = buildcov(phi[ ,t], dist.mo[[t]], covmodel=cov.model, nugget=is.nugget)
      R_sk = Rm - Rmo%*%RInv%*%t(Rmo)
      # RmoRInv = Rmo%*%RInv 


      y_t1 = array(NA, dim=c(nsample, n[t], N))
      ym_t1 = array(NA, dim=c(nsample, n.m[t], N))
      # IB = match.input(input[[t]], input.miss[[t-1]])$IB
      # y_t1 = ym.hat[[t-1]][IB]

    
      # b = matrix(NA, q[t]+1, N)
      # mu_y = matrix(NA, n.m[t], N)
      for(k in 1:nsample){

        y_t1[k, , ] = create.w(t=t, input=input, input.miss=input.miss[[t-1]],
                y=y[[t-1]], ym=as.matrix(ym.hat[[t-1]][k, , ]))
    
        ym_t1[k, , ] = create.w.pred(t=t, input=input[[t-1]], input.miss=input.miss,
              y=y[[t-1]], ym=as.matrix(ym.hat[[t-1]][k, ,]))

        ym.hat[[t]][k, , ] = conditional_simulation(y[[t]], H[[t]], y_t1[k, , ],  
                    RInv, Hm[[t]], ym_t1[k, , ], Rmo, R_sk)

      }



    }
  #}

################################################################################

krige = list()
krigeSE = list()
krige.lower95 = list()
krige.upper95 = list()
for(t in 1:S){
  # yhat[[t]] = ym.hat[[t]][ ,pred.ID[[t]], ]

  krige[[t]] = apply(ym.hat[[t]], c(2,3), mean)
  krigeSE[[t]] = apply(ym.hat[[t]], c(2,3), sd)
  krige.lower95[[t]] = apply(ym.hat[[t]], c(2,3), quantile, 0.025)
  krige.upper95[[t]] = apply(ym.hat[[t]], c(2,3), quantile, 0.975)
}

pred.mu = list()
pred.SE = list()
pred.lower95 = list()
pred.upper95 = list()
for(t in 1:S){
  pred.mu[[t]] = matrix(NA, np, N)
  pred.SE[[t]] = matrix(0, np, N)
  pred.lower95[[t]] = matrix(NA, np, N)
  pred.upper95[[t]] = matrix(NA, np, N)

  ind.list = ismember(input.new, input.miss[[t]])
  pred.mu[[t]][ind.list$IIA, ] = krige[[t]][ind.list$IA, ]
  pred.SE[[t]][ind.list$IIA, ] = krigeSE[[t]][ind.list$IA, ]
  pred.lower95[[t]][ind.list$IIA, ] = krige.lower95[[t]][ind.list$IA, ]
  pred.upper95[[t]][ind.list$IIA, ] = krige.upper95[[t]][ind.list$IA, ]
  if(length(ind.list$IIA)<np){
    ind.input = ismember(input.new, input[[t]])
    pred.mu[[t]][ind.input$IIA, ] = y[[t]][ind.input$IA, ]
    pred.lower95[[t]][ind.input$IIA, ] = y[[t]][ind.input$IA, ]
    pred.upper95[[t]][ind.input$IIA, ] = y[[t]][ind.input$IA, ]
  }
}

names(pred.mu) = paste0("Level", seq(1:S), "")
names(pred.SE) = paste0("Level", seq(1:S), "")
names(pred.lower95) = paste0("Level", seq(1:S), "")
names(pred.upper95) = paste0("Level", seq(1:S), "")

pred = list(mu=pred.mu, SE=pred.SE, 
       lower95=pred.lower95,upper95=pred.upper95)

  return(pred)

}

#############################################################################
#############################################################################






#############################################################################
#############################################################################

sample.ym <- function(y, input, param, Ho, Hm, dist.o, dist.m, dist.mo, cov.model="matern_5_2", nsample=30){
  
  S = length(y)


  Dim = dim(dist.o[[1]])[3]

  N = dim(y[[1]])[2]
  # param contains phi and nugget (maybe)
  if(length(param[ ,1])==Dim){ #no nugget
    #phi = exp(-param)
    is.nugget=FALSE
  }else{
    #phi = exp(-param[1:Dim, ,drop=FALSE])
    #nugget = exp(param[Dim+1, ,drop=FALSE]) / (1 + exp(param[Dim+1, ,drop=FALSE]))
    is.nugget = TRUE
    #phi = rbind(phi, nugget)
  }

  phi = param 

  inputlist = augment.input(input)
  input.miss = inputlist$miss

  

  nm = rep(NA, S-1)
  for(t in 1:(S-1)){
    nm[t] = dim(dist.m[[t]])[1]
  }
  n = rep(NA, S)
  q = rep(NA, S)
  for(t in 1:S){
    n[t] = dim(y[[t]])[1]
    q[t] = dim(Ho[[t]])[2]
  }
  ym = list()
  for(t in 1:(S-1)){
    ym[[t]] = array(NA, dim=c(nsample, nm[t], N))
  }
  
  names(ym) = paste0("Level", seq(1:(S-1)), "")


  t=1

  R = buildcov(phi[ ,t], dist.o[[t]], covmodel=cov.model, nugget=is.nugget)
  U = chol(R)
  RInv = chol2inv(U)
  HRHInv = solve(t(Ho[[t]])%*%RInv%*%Ho[[t]])

  Rm = buildcov(phi[ ,t], dist.m[[t]], covmodel=cov.model, nugget=is.nugget)
  Rmo = buildcov(phi[ ,t], dist.mo[[t]], covmodel=cov.model, nugget=FALSE)
  
  # compute conditional mean
  RmoRInv = Rmo%*%RInv
  KW = (Hm[[t]]-RmoRInv%*%Ho[[t]]) %*% HRHInv %*% t(Ho[[t]]) %*% RInv + RmoRInv
  mu_ymy = KW %*% y[[t]]  # n.m-by-N matrix

  RmoU = t(backsolve(U, t(Rmo), transpose=TRUE))

  # compute predictive variance 
  XXRR = t(Hm[[t]]) - t(Ho[[t]])%*%RInv%*%t(Rmo)
  c_star = Rm - Rmo%*%RInv%*%t(Rmo) + t(XXRR)%*%HRHInv%*%XXRR

  # compute S2
  Q = RInv - RInv%*%Ho[[t]]%*%HRHInv%*%t(Ho[[t]])%*%t(RInv)
  sigma2.hat = compute_Svec(y[[t]], Q) / (n[t]-q[t])

  # for(j in 1:N){
  #   Sig = c_star * sigma2.hat[j] 
  #   ym[[t]][ ,j] = c(mvtnorm::rmvt(1, sigma=Sig, df=n[t]-q[t], delta=mu_ymy[ ,j], type="shifted"))
  # }
  L = t(chol(c_star))
  ym[[t]] = sample_mvt(mu=mu_ymy, L=L, sigma=sigma2.hat, df=n[t]-q[t], nsample)


  if(S>2){
    for(t in 2:(S-1)){
      R = buildcov(phi[ ,t], dist.o[[t]], covmodel=cov.model, nugget=is.nugget)
      U = chol(R)
      RInv = chol2inv(U)
      # IB = match.input(input[[t]], input.miss[[t-1]])$IB
      # y_t1 = ym[[t-1]][IB]
      y_t1 = create.w(t=t, input=input, input.miss=input.miss[[t-1]], 
                y=y[[t-1]], ym=ym[[t-1]])

      # IB = match.input(input.miss[[t]], input.miss[[t-1]])$IB
      # ym_t1 = ym[[t-1]][IB]
      ym_t1 = create.w(t=t, input=input.miss, input.miss=input.miss[[t-1]], 
                y=ym[[t-1]], ym=ym[[t-1]])


      Rm = buildcov(phi[ ,t], dist.m[[t]], covmodel=cov.model, nugget=is.nugget)
      Rmo = buildcov(phi[ ,t], dist.mo[[t]], covmodel=cov.model, nugget=FALSE)
      
      R_sk = Rm - Rmo%*%RInv%*%t(Rmo)

      pred.list = conditional_simulation(y[[t]], Ho[[t]], y_t1,  
                  RInv, Hm[[t]], ym_t1, Rmo, R_sk)


      # sample ym at t
      #ym[[t]] = c(mvtnorm::rmvt(1, sigma=Sig, df=n[t]-q[t]-1, delta=mu_ymy, type="shifted"))

    }
  }



  
return(ym)
  
}

#############################################################################
#############################################################################



#############################################################################
#############################################################################




#############################################################################
#############################################################################





#############################################################################
#############################################################################




#############################################################################
#############################################################################


