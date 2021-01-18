"densityreg" <- function(Y, W, bw, xgrid=NULL, ygrid = NULL, sig=NULL, K1 = "Gauss", 
                         K2 = "Gauss", mean.estimate = NULL, spline.df = 5){
  #########################################################################################
  # data structure
  #########################################################################################
  n = length(Y);
  sdW = sd(W); sdY = sd(Y);
  #if(is.null(bw)) stop("please specify the bandwidth vector bw");
  if(is.null(xgrid)){
    xgrid = seq(quantile(W, probs=0.025), quantile(W, probs=0.975), length.out = 200);
  }
  if(is.null(ygrid)){
    ygrid = seq(quantile(Y, probs=0), quantile(Y, probs=1), length.out=200);
  }
  nx	= length(xgrid)
  ny	= length(ygrid)
  dt = 0.001
  tt = seq(-1,1,dt)
  if(is.null(mean.estimate)){
    if(is.null(sig)){
      if((K1=="Gauss")&(K2=="Gauss")){
        fitxy = .Call("LCfitDensityReg", W, Y, xgrid, ygrid, bw[1], bw[2], rep(0,nx), PACKAGE = "lpme")$fitxy;
      }else if((K1=="SecOrder")&(K2=="Gauss")){
        fitxy = .Call("LCfitDensityRegK", W, Y, xgrid, ygrid, bw[1], bw[2], rep(0,nx), PACKAGE = "lpme")$fitxy;
      }else{
        stop("For mean.estimate = NULL without measurement error: it currently only supports 
             (K1='SecOrder' and K2='Gauss') or (K1='Gauss' and K2='Gauss')")
      }
    }else{
      if((K1=="SecOrder")&(K2=="Gauss")){
        fitxy = .Call("LCfitDensityRegLap", W, Y, xgrid, ygrid, bw[1], bw[2], rep(0,nx), 
                      dt, tt, sig, PACKAGE = "lpme")$fitxy;
      }else{
        stop("For mean.estimate = NULL with measurement error: it currently only supports: 
             (K1='SecOrder' and K2='Gauss')")
      }
    }
  }else{
    if(is.null(sig)){
      ## naive mean curve
      if(mean.estimate=="spline"){
        lm.sp = lm(Y~ns(W, df=spline.df), data=data.frame(Y, W));
        mean0 = predict(lm.sp, data.frame(W=xgrid));
        newY0 = lm.sp$residuals;
      }else if(mean.estimate=="kernel"){
        bw0 = try(locpol::pluginBw(W, Y,deg=1,kernel=locpol::gaussK), TRUE);
        if('try-error' %in% class(bw0)){
          message("there is an error for in the function 'pluginBw'; the function 'thumbBw' was used instead.")
          bw0 = locpol::thumbBw(W, Y,deg=1,kernel=locpol::gaussK)
        }
        fit0 = lpme::meanreg(Y, W, bw0, method="naive", xgrid=seq(min(W), max(W), length.out=500));
        meanhat0 = splinefun(fit0$xgrid, fit0$yhat, method="natural");
        mean0 = meanhat0(xgrid);
        newY0 = Y-meanhat0(W);
      }else{
        stop("It currently only supports mean.estimate='spline' or 'kernel'.")
      }
      if((K1=="SecOrder")&(K2=="SecOrder")){
        fitxy = .Call("LCfitDensityRegKK", W, newY0, xgrid, ygrid, bw[1], bw[2], mean0, PACKAGE = "lpme")$fitxy;
      }else if((K1=="SecOrder")&(K2=="Gauss")){
        fitxy = .Call("LCfitDensityRegK", W, newY0, xgrid, ygrid, bw[1], bw[2], mean0, PACKAGE = "lpme")$fitxy;
      }else if((K1=="Gauss")&(K2=="Gauss")){
        fitxy = .Call("LCfitDensityReg", W, newY0, xgrid, ygrid, bw[1], bw[2], mean0, PACKAGE = "lpme")$fitxy;
      }else{
        stop("Currently it only supports: (K1='SecOrder' and K2='SecOrder') or 
           (K1='SecOrder' and K2='Gauss') or (K1='Gauss' and K2='Gauss')")
      }
    }else{
      ## deconvolution setup
      m  = 2^16; m_mid = m/2+1;
      beta  = sqrt((2*pi)/m); 
      beta2  = (2*pi)/(m*beta); 
      input  = seq(-m*beta/2, m*beta/2-beta, by=beta)
      output= seq(-pi/beta, pi/beta-beta2, by=beta2)
      mconst= (-1)^(0:(m-1));
      ## Kernel for which CF is (1-t^2)^3 with Normal errors
      FKsup  = function(t) { ifelse( (t<=1 & t>=-1), (1-t^2)^3, 0) }
      FKoutput = rep(0, m);
      FKoutput[m_mid] = FKsup(output[m_mid]);
      i=1; indicator1 =1; indicator2 =1;
      while ( ( abs(indicator1)>1e-30 | abs(indicator2)>1e-30 ) & (i<(m/2)) ) {
        indicator1 = FKsup(output[m_mid-i]);
        indicator2 = FKsup(output[m_mid+i]);
        FKoutput[m_mid-i] = indicator1;
        FKoutput[m_mid+i] = indicator2;
        i=i+1;
      }
      ## inverse FFT to get K
      Kinput  = Re( (-1)^(0:(m-1))/beta*fft( (-1)^(0:(m-1))*FKoutput)/m );
      ## naive mean curve
      if(mean.estimate=="spline"){
        lm.sp = lm(Y~ns(W, df=spline.df), data=data.frame(Y, W));
        mean0decon = predict(lm.sp, data.frame(W=input));
        newY0 = lm.sp$residuals;
      }else if(mean.estimate=="kernel"){
        bw0 = try(locpol::pluginBw(W, Y,deg=1,kernel=locpol::gaussK), TRUE);
        if('try-error' %in% class(bw0)){
          message("there is an error for in the function 'pluginBw'; the function 'thumbBw' was used instead.")
          bw0 = locpol::thumbBw(W, Y,deg=1,kernel=locpol::gaussK)
        }
        fit0 = lpme::meanreg(Y, W, bw0, method="naive", xgrid=input);
        meanhat0 = splinefun(fit0$xgrid, fit0$yhat, method="natural");
        newY0 = Y-meanhat0(W);
        mean0decon = meanhat0(input);
      }else{
        stop("It currently only supports mean.estimate='spline' or 'kernel'.")
      }
      if((K1=="SecOrder")&(K2=="SecOrder")){
        fitxy = .Call("LCfitDensityRegLap2", W, newY0, xgrid, ygrid, bw[1], bw[2], mean0decon, 
                      input, output, beta, beta2, mconst, Kinput, sig, PACKAGE = "lpme")$fitxy;
      }else{
        stop("Currently it only supports: (K1='SecOrder' and K2='SecOrder')")
      }
    }
  }
  
  #### Save to a list
  output <- list(xgrid = xgrid, 
                 ygrid = ygrid,
                 fitxy = fitxy);
  class(output) <- c("densityreg")
  output
}
