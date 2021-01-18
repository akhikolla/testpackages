#' @title Prediction at new inputs in the autoregressive cokriging model
#' @description This function makes prediction in
#'  autogressive cokriging models. If a nested design is used, the predictive mean and predictive variance are
#'  computed exactly; otherwise, Monte Carlo simulation from the predictive distribution is used to approximate
#' the predictive mean and predictive variance. 
#'  
#' @param obj a \code{\link{cokm}} object construted via the function \code{\link{cokm}} in 
#' this package
#' @param input.new a matrix including new inputs for making prediction
#' @author Pulong Ma <mpulong@gmail.com>
#' 
#' @export
#' 
#' @seealso \code{\link{cokm}}, \code{\link{cokm.fit}}, \code{\link{cokm.condsim}}, \code{\link{ARCokrig}}
#' 
#' @examples 
#' Funcc = function(x){
#'   return(0.5*(6*x-2)^2*sin(12*x-4)+10*(x-0.5)-5)
#' }
#' 
#' Funcf = function(x){
#'   z1 = Funcc(x)
#'   z2 = 2*z1-20*x+20 + sin(10*cos(5*x))
#'   return(z2)
#' }
#' 
#' #####################################################################
#' ###### Nested design 
#' #####################################################################
#' Dc <- seq(-1,1,0.1)
#' indDf <- c(1, 3, 6, 8, 10, 13, 17, 21)
#' zc <- Funcc(Dc)
#' Df <- Dc[indDf]
#' zf <- Funcf(Df)
#' 
#' input.new = as.matrix(seq(-1,1,length.out=200))
#' 
#' 
#' ## create the cokm object
#' prior = list(name="Reference")
#' obj = cokm(formula=list(~1,~1+x1), output=list(c(zc), c(zf)),
#'               input=list(as.matrix(Dc), as.matrix(Df)),
#'               prior=prior, cov.model="matern_5_2")
#'
#' ## update model parameters in the cokm object
#' 
#' obj = cokm.fit(obj)
#' 
#' cokrige = cokm.predict(obj, input.new)
#' 
#' 

cokm.predict <- function(obj, input.new){

formula = obj@formula
output = obj@output 
input = obj@input 
param = obj@param
cov.model = obj@cov.model 
nugget.est = obj@nugget.est
NestDesign = obj@NestDesign

phi = do.call(cbind, param)

if(!is.matrix(input.new)){
  stop("input.new should be a matrix.")
}

if(NestDesign){

  Dim = dim(input[[1]])[2]
  p.x = Dim
  is.nugget = nugget.est

  S = length(output)
  n = sapply(output, length)
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
    betahat = HRHInv%*%t(H[[t]])%*%(RInv%*%y[[t]])
    res = y[[t]] - H[[t]]%*%betahat
    SSE = c(t(res)%*%RInv%*%res) / (n[t]-q[t])    # \hat{sigma}^2

    Rm = buildcov(phi[ ,t], dist.m[[t]], covmodel=cov.model, nugget=FALSE)
    Rmo = buildcov(phi[ ,t], dist.mo[[t]], covmodel=cov.model, nugget=is.nugget)
    

    #RmoU = t(backsolve(U, t(Rmo), transpose=TRUE))

    XXRR = t(Hm[[t]]) - t(H[[t]])%*%RInv%*%t(Rmo)
    #XXRR = t(Hm[[t]]) - t(RmoU%*%UH)
    Sig_ymymy = Rm - Rmo%*%RInv%*%t(Rmo) + t(XXRR)%*%HRHInv%*%XXRR
    #Sig_ymymy = Rm - tcrossprod(RmoU) + t(XXRR)%*%HRHInv%*%XXRR
    Sig = Sig_ymymy * SSE
    #Sig = Sig_ymymy * SSE/(n[t]-q[t])

    # mu_ymy = Hm[[t]]%*%betahat + Rmo%*%(RInv%*%res)
    krige[[t]] = c(Hm[[t]]%*%betahat + Rmo%*%(RInv%*%res))

    const = (n[t]-q[t])/(n[t]-q[t]-2)
    #krigeSE[[t]] = sqrt(diag(Sig)[pred.ID[[t]]])
    krige.var[[t]] = const*diag(Sig)
    krige.var[[t]][krige.var[[t]]<0] = 0 

    #for(k in 1:nsample){
      for(t in 2:(S)){
      
        ############################################################################
        #### estimate missing data ym
        ############################################################################

        R = buildcov(phi[ ,t], dist.o[[t]], covmodel=cov.model, nugget=is.nugget)
        U = chol(R)
        RInv = chol2inv(U)
        IB = match.input(input[[t]], input[[t-1]])$IB
        y_t1 = y[[t-1]][IB]
        IB = match.input(input.miss[[t]], input.miss[[t-1]])$IB
        ym_t1 = as.matrix(krige[[t-1]][IB])

        X = cbind(H[[t]], y_t1)
        XRXInv = solve(t(X)%*%RInv%*%X)
        betahat = XRXInv%*%t(X)%*%(RInv%*%y[[t]])
        gamma = betahat[q[t]+1]

        res = y[[t]] - X%*%betahat
        SSE = c(t(res)%*%RInv%*%res) / (n[t]-q[t]-1)

        Xm = cbind(Hm[[t]], ym_t1)
        Rm = buildcov(phi[ ,t], dist.m[[t]], covmodel=cov.model, nugget=FALSE)
        Rmo = buildcov(phi[ ,t], dist.mo[[t]], covmodel=cov.model, nugget=is.nugget)
        
        RInvH = RInv %*% H[[t]]
        HRHInv = solve(t(H[[t]])%*%RInvH)
        Q = RInv - RInvH %*% HRHInv %*% t(RInvH)
        XXRR = t(Xm) - t(X)%*%RInv%*%t(Rmo)
        #RmoU = t(backsolve(U, t(Rmo), transpose=TRUE))
        Sig_ymymy = diag(Rm - Rmo%*%RInv%*%t(Rmo) + t(XXRR)%*%XRXInv%*%XXRR)  + 
                    krige.var[[t-1]] / c(crossprod(y_t1, Q %*% y_t1))

        #Sig_ymymy = Rm - tcrossprod(RmoU) + t(XXRR)%*%XRXInv%*%XXRR
        const = (n[t]-q[t]-1) / (n[t]-q[t]-3)
        Sig = gamma^2*krige.var[[t-1]] + const*Sig_ymymy * SSE
        #Sig = Sig_ymymy * SSE/(n[t])

        krige[[t]] = c(Xm%*%betahat + Rmo%*%(RInv%*%res))

        krige.var[[t]] = Sig
        krige.var[[t]][krige.var[[t]]<0] = 0 

      }
    #}

  ################################################################################

  # for(t in 1:(S-1)){
  #   ind = sort(ID.org[[t]], index.return=TRUE)$ix 
  #   krige[[t]][ind] = krige[[t]]
  #   krigeSE[[t]][ind] = krigeSE[[t]]
  # }


  confint = list()
  krigeSE = list()
  lower95 = list()
  upper95 = list()
  for(t in 1:S){
    krigeSE[[t]] = sqrt(krige.var[[t]])
    # degree = ifelse(t==1, n[t]-q[t], n[t]-q[t]-1)
    # confint[[t]] = cbind(krige[[t]] + qt(0.025, df=degree)*krigeSE[[t]], 
    #                      krige[[t]] + qt(0.975, df=degree)*krigeSE[[t]])
    lower95[[t]]=krige[[t]] - 2*krigeSE[[t]] 
    upper95[[t]]=krige[[t]] + 2*krigeSE[[t]]
  }

  names(krige) = paste0("Level", seq(1:S), "")
  names(krigeSE) = paste0("Level", seq(1:S), "")
  names(lower95) = paste0("Level", seq(1:S), "")
  names(upper95) = paste0("Level", seq(1:S), "")


  out = list(mu=krige, SE=krigeSE, lower95=lower95, upper95=upper95)


}else{
  
  nsample = obj@tuning$n.sample
  out = condsim.NN.univariate(formula=formula, output=output, input=input,
                input.new=input.new, phi=phi, cov.model=cov.model,
                nsample=nsample)
}



return(out)




}
