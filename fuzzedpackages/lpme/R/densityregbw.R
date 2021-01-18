"densityregbw" <- function(Y, W, h1=NULL, h2=NULL, sig = NULL,
                           xinterval = quantile(W, probs=c(0.025, 0.975), names = FALSE),
                           K1 = "Gauss", K2 = "Gauss", mean.estimate = NULL, spline.df = 5){
  #########################################################################################
  # data structure
  #########################################################################################
  ## initial normal reference bandwidth;
  #Ker = function(x) .Call("Kern_2nd_order", x, PACKAGE = "lpme"); #Ker = function(x) dnorm(x);
  #intK2 = integrate(function(x) Ker(x)^2, -50, 50)$value;
  #intx2K = integrate(function(x) x^2*Ker(x), -50, 50)$value;
  #const = (8*sqrt(pi)*intK2/(3*intx2K^2))^(1/5); const;
  const = 0.427398;
  n = length(Y);
  sdW = sd(W); sdY = sd(Y);
  ## weights
  fWhat=function(w) {
    resp = ifelse((w<=max(xinterval)&w>=min(xinterval)), 1, 0);
    resp[resp<0]=0; resp
  }
  pW = fWhat(W);
  if(is.null(mean.estimate)){
    if(is.null(sig)){
      if((K1=="Gauss")&(K2=="Gauss")){
        hxyhat = c(sdW*1.06, sdY*1.06)*n^(-1/5);
        if(is.null(h1)) h1 = hxyhat[1]*seq(0.2, 1.5, length.out = 20 )
        if(is.null(h2)) h2 = hxyhat[2]*seq(0.2, 1.5, length.out = 15 )
        hh = as.matrix( expand.grid(h1, h2) );
        resCV = .Call("CVdens_LCfit", W, Y, pW, h1, h2, PACKAGE = "lpme")$CV;
        hhxy = as.vector(hh[which.min(as.vector(resCV)),]); 
      }else if((K1=="SecOrder")&(K2=="Gauss")){
        K1=="SecOrder"; K2=="Gauss"
        hxyhat = c(sdW*const, sdY*1.06)*n^(-1/5);
        if(is.null(h1)) h1 = hxyhat[1]*seq(0.2, 1.5, length.out = 20 )
        if(is.null(h2)) h2 = hxyhat[2]*seq(0.2, 1.5, length.out = 15 )
        hh = as.matrix( expand.grid(h1, h2) );
        resCV = .Call("CVdens_LCfit2", W, Y, pW, h1, h2, PACKAGE = "lpme")$CV;
        hhxy = as.vector(hh[which.min(as.vector(resCV)),]);
      }else{
        stop("For mean.estimate = NULL without measurement error: it currently only supports 
             (K1='SecOrder' and K2='Gauss') or (K1='Gauss' and K2='Gauss')")
      }
    }else{
      if((K1=="SecOrder")&(K2=="Gauss")){
        ## reliability estimate:
        lamhat = ((var(W)-sig^2))/var(W);
        corYW = abs(cor(Y,W))
        h_scale = (1+sqrt(1-lamhat)*corYW);
        hxyhat = c(sdW*const, sdY*1.06)*n^(-1/5);
        if(is.null(h1)) h1 = hxyhat[1]*seq(0.2, 1.5, length.out = 20 )*h_scale
        if(is.null(h2)) h2 = hxyhat[2]*seq(0.2, 1.5, length.out = 15 )
        hh = as.matrix( expand.grid(h1, h2) );
        resCV = .Call("CVdens_LCfit2", W, Y, pW, h1/h_scale, h2, PACKAGE = "lpme")$CV;
        hhxy = as.vector(hh[which.min(as.vector(resCV)),]); 
      }else{
        stop("For mean.estimate = NULL with measurement error: it currently only supports:
             (K1='SecOrder' and K2='Gauss')")
      }
    }
  }else{
    if(is.null(sig)){
      ## naive mean curve
      if(mean.estimate=="spline"){
        lm.sp = lm(Y~splines::ns(W, df=spline.df), data=data.frame(Y, W));
        newY0 = lm.sp$residuals;
      }else if(mean.estimate=="kernel"){
        bw0 = try(locpol::pluginBw(W, Y,deg=1,kernel=locpol::gaussK), TRUE);
        if('try-error' %in% class(bw0)){
          message("there is an error for in the function 'pluginBw'; the function 'thumbBw' was used instead.")
          bw0 = locpol::thumbBw(W, Y,deg=1,kernel=locpol::gaussK)
        }
        fit0 = lpme::meanreg(Y, W, bw0, method="naive", xgrid=seq(min(W), max(W), length.out=500));
        meanhat0 = splinefun(fit0$xgrid, fit0$yhat, method="natural");
        newY0 = Y-meanhat0(W);
      }else{
        stop("It currently only supports mean.estimate='spline' or 'kernel'.")
      }
      if((K1=="SecOrder")&(K2=="SecOrder")){
        hxyhat = c(sdW*const, sd(newY0)*const)*n^(-1/5);
        if(is.null(h1)) h1 = hxyhat[1]*seq(0.2, 1.5, length.out = 20 )
        if(is.null(h2)) h2 = hxyhat[2]*seq(0.2, 1.5, length.out = 15 )
        hh = as.matrix( expand.grid(h1, h2) );
        resCV = .Call("CVdens_LCfit2", W, newY0, pW, h1, h2*1.06/const, PACKAGE = "lpme")$CV;
        hhxy = as.vector(hh[which.min(as.vector(resCV)),]);
      }else if((K1=="SecOrder")&(K2=="Gauss")){
        hxyhat = c(sdW*const, sd(newY0)*1.06)*n^(-1/5);
        if(is.null(h1)) h1 = hxyhat[1]*seq(0.2, 1.5, length.out = 20 )
        if(is.null(h2)) h2 = hxyhat[2]*seq(0.2, 1.5, length.out = 15 )
        hh = as.matrix( expand.grid(h1, h2) );
        resCV = .Call("CVdens_LCfit2", W, newY0, pW, h1, h2, PACKAGE = "lpme")$CV;
        hhxy = as.vector(hh[which.min(as.vector(resCV)),]);
      }else if((K1=="Gauss")&(K2=="Gauss")){
        hxyhat = c(sdW*1.06, sd(newY0)*1.06)*n^(-1/5);
        if(is.null(h1)) h1 = hxyhat[1]*seq(0.2, 1.5, length.out = 20 )
        if(is.null(h2)) h2 = hxyhat[2]*seq(0.2, 1.5, length.out = 15 )
        hh = as.matrix( expand.grid(h1, h2) );
        resCV = .Call("CVdens_LCfit", W, newY0, pW, h1, h2, PACKAGE = "lpme")$CV;
        hhxy = as.vector(hh[which.min(as.vector(resCV)),]);
      }else{
        stop("Currently it only supports: (K1='SecOrder' and K2='SecOrder') or 
             (K1='SecOrder' and K2='Gauss') or (K1='Gauss' and K2='Gauss')")
      }
    }else{
      ## naive mean curve
      if(mean.estimate=="spline"){
        lm.sp = lm(Y~splines::ns(W, df=spline.df), data=data.frame(Y, W));
        newY0 = lm.sp$residuals;
      }else if(mean.estimate=="kernel"){
        bw0 = try(locpol::pluginBw(W, Y,deg=1,kernel=locpol::gaussK), TRUE);
        if('try-error' %in% class(bw0)){
          message("there is an error for in the function 'pluginBw'; the function 'thumbBw' was used instead.")
          bw0 = locpol::thumbBw(W, Y,deg=1,kernel=locpol::gaussK)
        }
        fit0 = lpme::meanreg(Y, W, bw0, method="naive", xgrid=seq(min(W), max(W), length.out=500));
        meanhat0 = splinefun(fit0$xgrid, fit0$yhat, method="natural");
        newY0 = Y-meanhat0(W);
      }else{
        stop("It currently only supports mean.estimate='spline' or 'kernel'.")
      }
      if((K1=="SecOrder")&(K2=="SecOrder")){
        ## reliability estimate:
        lamhat = ((var(W)-sig^2))/var(W);
        corYW = abs(cor(newY0,W))
        h_scale = (1+sqrt(1-lamhat)*corYW);
        hxyhat = c(sdW*const, sd(newY0)*const)*n^(-1/5);
        if(is.null(h1)) h1 = hxyhat[1]*seq(0.2, 1.5, length.out = 20 )*h_scale
        if(is.null(h2)) h2 = hxyhat[2]*seq(0.2, 1.5, length.out = 15 )
        hh = as.matrix( expand.grid(h1, h2) );
        resCV = .Call("CVdens_LCfit2", W, newY0, pW, h1/h_scale, h2*1.06/const, PACKAGE = "lpme")$CV;
        hhxy = as.vector(hh[which.min(as.vector(resCV)),]);
      }else{
        stop("Currently it only supports: (K1='SecOrder' and K2='SecOrder')")
      }
    }
  }
  output <- list(bw = hhxy,
                 h1 = h1,
                 h2 = h2);
  class(output) <- c("densityregbw")
  output
}