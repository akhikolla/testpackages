"moderegbwSIMEX" <- function(Y, W, method="CV-density", p.order=0, sig, B=5, h1=NULL, h2=NULL,
                              length.h=10, CIregion = 0.95, nstart = 4){
  #########################################################################################
  # data structure
  #########################################################################################
  n = length(Y);
  dt=0.001;
  tt = seq(-1, 1, dt);
  Ws = matrix(0, n, B); Wss=Ws;
  for (i in 1:B){
    Ws[,i] = W+.rlaplace(n,0,sig/sqrt(2));
    Wss[,i] = Ws[,i]+.rlaplace(n,0,sig/sqrt(2));
  }
  maxWs = rep(0, B); maxWss = rep(0, B);
  for(b in 1:B){
    maxWs[b] = max(abs(max(Ws[,b])-min(W)), abs(max(Ws[,b])-max(W)), 
                   abs(min(Ws[,b])-min(W)), abs(min(Ws[,b])-max(W)))
    maxWss[b] = max(abs(max(Wss[,b])-min(Ws[,b])), abs(max(Wss[,b])-max(Ws[,b])), 
                    abs(min(Wss[,b])-min(Ws[,b])), abs(min(Wss[,b])-max(Ws[,b])))
  }
  maxW = max(maxWs,maxWss); 
  ## approximated kernels;
  delta=0.005;
  xgrid1 = seq(0, 5.1, delta);
  K2y = dnorm(xgrid1, 0, 1);
  xgrid2 = seq(0, 7.1, delta);
  intK2y = exp(-xgrid2^2/4)/(2*sqrt(pi));
  ## default bandwidth range
  if(is.null(h2)){
    h2ref = c(1.06*sd(Y)*n^(-0.2));
  }else{
    h2ref = h2;
  }
  if(is.null(h1)) {
    h1sB = rep(0, B); h1ssB = rep(0, B);
    for(b in 1:B) {
      h1sB[b] = decon::bw.dmise(Ws[,b], sig=sig, error="laplacian")
      h1ssB[b] = decon::bw.dmise(Wss[,b], sig=sig, error="laplacian")
    }
    bw1s = mean(h1sB); bw1ss = mean(h1ssB);
    if(p.order==0){
      h1s = seq(bw1s*0.2, bw1s*2, length.out=length.h)
      h1ss = h1s
    }else{
      h1s = seq(bw1s*0.2, bw1s*2, length.out=length.h)
      h1ss = h1s
    }
  }else{
    h1s = h1; h1ss = h1;
  }
  ## density of W
  dd=density(W)
  dens=splinefun(dd$x, dd$y);
  winterval = quantile(W, probs=c(0.5-CIregion/2, 0.5+CIregion/2), names = FALSE);
  fWhat=function(w) {
    resp = ifelse((w<=winterval[2]&w>=winterval[1]), 1, 0);
    resp[resp<0]=0; 
    resp
  }
  pW = fWhat(W);
  ## density of Ws
  pWs= matrix(1, n, B);
  for (b in 1:B){
    dd = density(Ws[,b]);
    dens = splinefun(dd$x, dd$y);
    wsinterval = quantile(Ws[,b], probs=c(0.5-CIregion/2, 0.5+CIregion/2), names = FALSE);
    fWhats = function(w){
      resp = ifelse((w<=wsinterval[2]&w>=wsinterval[1]), 1, 0);
      resp[resp<0]=0; 
      resp
    }
    pWs[,b] = fWhats(Ws[,b]);
  }
  ## estimate bandwidth
  if(method=="CV-density"){
    if(p.order==0){
      xgrid3 = seq(0, round(maxW/min(h1s)/delta)*delta, delta);
      Kux = .Call("Ku0_sec_order", length(xgrid3), delta, h1s, sig, tt, dt, PACKAGE = "lpme");
      Ku0x = Kux$Ku0;
      CV1 = matrix(0, length(h1s), B);
      for(b in 1:B){
        CV = .Call("CVdens_LCfitLap", Ws[,b], W, Y, pW, h1s, h2ref,
                   Ku0x, K2y, intK2y, delta, PACKAGE = "lpme")$CV;
        CV1[,b] = as.vector(CV);
      }
      CV1mean = rowMeans(CV1); hhx1 = as.vector(h1s[which.min(CV1mean)]);
      #h1ss = seq(hhx1*1, hhx1*3, length.out=10);
      if( any(h1s!=h1ss) ){
        xgrid3 = seq(0, round(maxW/min(h1ss)/delta)*delta, delta);
        Kux = .Call("Ku0_sec_order", length(xgrid3), delta, h1ss, sig, tt, dt, PACKAGE = "lpme");
        Ku0x = Kux$Ku0;
      }
      CV2 = matrix(0, length(h1ss), B);
      for(b in 1:B){
        CV = .Call("CVdens_LCfitLap", Wss[,b], Ws[,b], Y, pWs[,b], h1ss, h2ref,
                   Ku0x, K2y, intK2y, delta, PACKAGE = "lpme")$CV;
        CV2[,b] = as.vector(CV);
      }
      CV2mean = rowMeans(CV2); hhx2 = as.vector(h1ss[which.min(CV2mean)]);
      hhxy = c(hhx1^2/hhx2,h2ref);
    }else{
      xgrid3 = seq(0, round(maxW/min(h1s)/delta)*delta, delta);
      Kux = .Call("Ku_sec_order", length(xgrid3), delta, h1s, sig, tt, dt, PACKAGE = "lpme");
      Ku0x = Kux$Ku0; Ku1x = Kux$Ku1; Ku2x = Kux$Ku2;
      CV1 = matrix(0, length(h1s), B);
      for(b in 1:B){
        CV = .Call("CVdens_LLfitLap", Ws[,b], W, Y, pW, h1s, h2ref,
                   Ku0x, Ku1x, Ku2x, K2y, intK2y, delta, PACKAGE = "lpme")$CV;
        CV1[,b] = as.vector(CV);
      }
      CV1mean = rowMeans(CV1); hhx1 = as.vector(h1s[which.min(CV1mean)]);
      #h1ss = seq(hhx1*1, hhx1*3, length.out=10);
      if( any(h1s!=h1ss) ){
        xgrid3 = seq(0, round(maxW/min(h1ss)/delta)*delta, delta);
        Kux = .Call("Ku_sec_order", length(xgrid3), delta, h1ss, sig, tt, dt, PACKAGE = "lpme");
        Ku0x = Kux$Ku0; Ku1x = Kux$Ku1; Ku2x = Kux$Ku2;
      }
      CV2 = matrix(0, length(h1ss), B);
      for(b in 1:B){
        CV = .Call("CVdens_LLfitLap", Wss[,b], Ws[,b], Y, pWs[,b], h1ss, h2ref,
                   Ku0x, Ku1x, Ku2x, K2y, intK2y, delta, PACKAGE = "lpme")$CV;
        CV2[,b] = as.vector(CV);
      }
      CV2mean = rowMeans(CV2); hhx2 = as.vector(h1ss[which.min(CV2mean)]);
      hhxy = c(hhx1^2/hhx2,h2ref);
    }
  }else{
    dd = cbind(W, Y, pW); Worder = order(dd[,1]); 
    dd = dd[order(dd[,1]),]; ddW = dd;
    x = dd[,1];
    nx = length(x);
    nyx = nstart;
    x.num = rep(nyx, nx);
    yindx = c(0, cumsum(x.num));
    meshyW = rep(0, yindx[nx+1]);
    quanprobs = seq(0.05, 0.95, length=nyx);
    for(ii in 1:nx){
      indxl=yindx[ii]+1;
      indxr=yindx[ii+1];
      xdiff= dd[-1,1]-dd[-n,1];
      xdiff = xdiff[xdiff!=0];
      win0=quantile(abs(xdiff), 0.1);
      indx = which(dd[,1]>(x[ii]-win0*nyx)&dd[,1]<(x[ii]+win0*nyx));
      win = win0;
      while(length(indx)<10*nyx){
        win=win+win0;
        indx = which(dd[,1]>(x[ii]-win)&dd[,1]<(x[ii]+win));
      }
      yii = dd[indx, 2];
      meshyW[indxl:indxr] = sort(yii)[round(length(indx)*quanprobs)];
    }
    if(p.order==0){
      xgrid3 = seq(0, round(maxW/min(h1s)/delta)*delta, delta);
      Kux = .Call("Ku0_sec_order", length(xgrid3), delta, h1s, sig, tt, dt, PACKAGE = "lpme");
      Ku0x = Kux$Ku0;
      CV1 = matrix(0, length(h1s), B);
      for(b in 1:B){
        CV = .Call("CVmode_LCfitLap", Ws[Worder,b], ddW[,1], ddW[,2], ddW[,3], meshyW, yindx, h1s, h2ref,
                   200, 1e-6, Ku0x, delta, PACKAGE = "lpme")$CV;
        CV1[,b] = as.vector(CV);
      }
      CV1mean = rowMeans(CV1); hhx1 = as.vector(h1s[which.min(CV1mean)]);
      #h1ss = seq(hhx1*1, hhx1*3, length.out=10);
      if( any(h1s!=h1ss) ){
        xgrid3 = seq(0, round(maxW/min(h1ss)/delta)*delta, delta);
        Kux = .Call("Ku0_sec_order", length(xgrid3), delta, h1ss, sig, tt, dt, PACKAGE = "lpme");
        Ku0x = Kux$Ku0;
      }
      CV2 = matrix(0, length(h1ss), B);
      for(b in 1:B){
        dd = cbind(Ws[,b], Y, pWs[,b]); Wsorder = order(dd[,1]); 
        dd = dd[order(dd[,1]),]; ddWs = dd;
        x = dd[,1];
        nx = length(x);
        nyx = nstart;
        x.num = rep(nyx, nx);
        yindx = c(0, cumsum(x.num));
        meshyWs = rep(0, yindx[nx+1]);
        quanprobs = seq(0.05, 0.95, length=nyx);
        for(ii in 1:nx){
          indxl=yindx[ii]+1;
          indxr=yindx[ii+1];
          xdiff= dd[-1,1]-dd[-n,1];
          xdiff = xdiff[xdiff!=0];
          win0=quantile(abs(xdiff), 0.1);
          indx = which(dd[,1]>(x[ii]-win0*nyx)&dd[,1]<(x[ii]+win0*nyx));
          win = win0;
          while(length(indx)<10*nyx){
            win=win+win0;
            indx = which(dd[,1]>(x[ii]-win)&dd[,1]<(x[ii]+win));
          }
          yii = dd[indx, 2];
          meshyWs[indxl:indxr] = sort(yii)[round(length(indx)*quanprobs)];
        }
        CV = .Call("CVmode_LCfitLap", Wss[Wsorder,b], ddWs[,1], ddWs[,2], ddWs[,3], meshyWs, yindx, h1ss, h2ref,
                   200, 1e-6, Ku0x, delta, PACKAGE = "lpme")$CV;
        CV2[,b] = as.vector(CV);
      }
      CV2mean = rowMeans(CV2); hhx2 = as.vector(h1ss[which.min(CV2mean)]);
      hhxy = c(hhx1^2/hhx2,h2ref);
    }else{
      xgrid3 = seq(0, round(maxW/min(h1s)/delta)*delta, delta);
      Kux = .Call("Ku_sec_order", length(xgrid3), delta, h1s, sig, tt, dt, PACKAGE = "lpme");
      Ku0x = Kux$Ku0; Ku1x = Kux$Ku1; Ku2x = Kux$Ku2;
      CV1 = matrix(0, length(h1s), B);
      for(b in 1:B){
        CV = .Call("CVmode_LLfitLap", Ws[Worder,b], ddW[,1], ddW[,2], ddW[,3], meshyW, yindx, h1s, h2ref,
                   200, 1e-6, Ku0x, Ku1x, Ku2x, delta, PACKAGE = "lpme")$CV;
        CV1[,b] = as.vector(CV);
      }
      CV1mean = rowMeans(CV1); hhx1 = as.vector(h1s[which.min(CV1mean)]);
      #h1ss = seq(hhx1*1, hhx1*3, length.out=10);
      if( any(h1s!=h1ss) ){
        xgrid3 = seq(0, round(maxW/min(h1ss)/delta)*delta, delta);
        Kux = .Call("Ku_sec_order", length(xgrid3), delta, h1ss, sig, tt, dt, PACKAGE = "lpme");
        Ku0x = Kux$Ku0; Ku1x = Kux$Ku1; Ku2x = Kux$Ku2;
      }
      CV2 = matrix(0, length(h1ss), B);
      for(b in 1:B){
        dd = cbind(Ws[,b], Y, pWs[,b]); Wsorder = order(dd[,1]); 
        dd = dd[order(dd[,1]),]; ddWs = dd;
        x = dd[,1];
        nx = length(x);
        nyx = nstart;
        x.num = rep(nyx, nx);
        yindx = c(0, cumsum(x.num));
        meshyWs = rep(0, yindx[nx+1]);
        quanprobs = seq(0.05, 0.95, length=nyx);
        for(ii in 1:nx){
          indxl=yindx[ii]+1;
          indxr=yindx[ii+1];
          xdiff= dd[-1,1]-dd[-n,1];
          xdiff = xdiff[xdiff!=0];
          win0=quantile(abs(xdiff), 0.1);
          indx = which(dd[,1]>(x[ii]-win0*nyx)&dd[,1]<(x[ii]+win0*nyx));
          win = win0;
          while(length(indx)<10*nyx){
            win=win+win0;
            indx = which(dd[,1]>(x[ii]-win)&dd[,1]<(x[ii]+win));
          }
          yii = dd[indx, 2];
          meshyWs[indxl:indxr] = sort(yii)[round(length(indx)*quanprobs)];
        }
        CV = .Call("CVmode_LLfitLap", Wss[Wsorder,b], ddWs[,1], ddWs[,2], ddWs[,3], meshyWs, yindx, h1ss, h2ref,
                   200, 1e-6, Ku0x, Ku1x, Ku2x, delta, PACKAGE = "lpme")$CV;
        CV2[,b] = as.vector(CV);
      }
      CV2mean = rowMeans(CV2); hhx2 = as.vector(h1ss[which.min(CV2mean)]);
      hhxy = c(hhx1^2/hhx2,h2ref);
    }
  }
  output <- list(bw = hhxy,
                 hx1 = hhx1,
                 hx2 = hhx2,
                 h1s = h1s,
                 h1ss = h1ss,
                 CV1 = CV1, 
                 CV2 = CV2, 
                 CV1mean=CV1mean,
                 CV2mean=CV2mean);
  class(output) <- c("moderegbwSIMEX")
  output
}