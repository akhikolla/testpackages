"meanreg" <- function(Y, W, bw, xgrid=NULL, method="HZ", sig=NULL, error="laplace", FT_fu){
  #########################################################################################
  # data structure
  #########################################################################################
  if(is.null(xgrid)){
    xgrid = seq(quantile(W, probs=0.025), quantile(W,probs=0.975), length.out = 100);
  }
  xmin  = min(xgrid);
  xmax  = max(xgrid);
  m  = 2^16; m_mid = m/2+1;
  beta  = sqrt((2*pi)/m); 
  beta2  = (2*pi)/(m*beta); 
  input  = seq(-m*beta/2, m*beta/2-beta, by=beta)
  output= seq(-pi/beta, pi/beta-beta2, by=beta2)
  nlower= round(xmin/beta+m/2+1); nupper = round(xmax/beta+m/2+1); xindex = nlower:nupper; 
  mconst= (-1)^(0:(m-1)); xindexC = xindex-1;
  x  = input[xindex]; 
  wide	= beta;
  nx	= length(x);
  ## Kernel for which CF is (1-t^2)^8 with Normal errors
  FKsup  = function(t) { ifelse( (t<=1 & t>=-1), (1-t^2)^8, 0) }
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
  dt = 0.0001;
  tt = seq(-1,1,dt);
  
  #########################################################################################
  # calling the c++ code and # output
  #########################################################################################
  if(method=="HZ"){
    if(error=="laplace"){
      if(is.null(sig)) stop("please specify sig if non naive method is used")
      foo <- .Call("fitnewLap", 
                   x_ = x, 
                   input_ = input, 
                   output_ = output, 
                   beta_ = beta, 
                   beta2_ = beta2, 
                   mconst_ = mconst, 
                   Kinput_ = Kinput, 
                   W_ = W, 
                   Y_ = Y, 
                   sigU_ = sig, 
                   h_ = bw,
                   PACKAGE = "lpme");
    }else if(error=="normal") {
      if(is.null(sig)) stop("please specify sig if non naive method is used")
      foo <- .Call("fitnewGau", 
                   x_ = x, 
                   input_ = input, 
                   output_ = output, 
                   beta_ = beta, 
                   beta2_ = beta2, 
                   mconst_ = mconst, 
                   Kinput_ = Kinput, 
                   W_ = W, 
                   Y_ = Y, 
                   sigU_ = sig, 
                   h_ = bw,
                   PACKAGE = "lpme");
    }else{
      if(!is.function(FT_fu)) stop("pleaes specify the characteristic function of U as a function");
      FfU = FT_fu(output);
      foo <- .Call("fitNEW", 
                   x_ = x, 
                   input_ = input, 
                   output_ = output, 
                   beta_ = beta, 
                   beta2_ = beta2, 
                   mconst_ = mconst, 
                   Kinput_ = Kinput, 
                   W_ = W, 
                   Y_ = Y, 
                   FfU_ = FfU, 
                   h_ = bw,
                   PACKAGE = "lpme");
    }
  } else if (method=="DFC") {
    if(is.null(sig)) stop("please specify sig if non naive method is used")
    if(error=="laplace"){
      foo <- .Call("fitjasaLap", 
                   x_ = x, 
                   h_ = bw,
                   W_ = W, 
                   Y_ = Y, 
                   sigU_ = sig, 
                   dt_ = dt, 
                   t_ = tt,                     
                   PACKAGE = "lpme");
    } else {
      foo <- .Call("fitjasaGau", 
                   x_ = x, 
                   h_ = bw,
                   W_ = W, 
                   Y_ = Y, 
                   sigU_ = sig, 
                   dt_ = dt, 
                   t_ = tt,
                   PACKAGE = "lpme");
    }
  } else {
    foo <- .Call("fitlocpoly", 
                 x_ = x,
                 beta_ = beta,
                 Kinput_ = Kinput,
                 W_ = W, 
                 Y_ = Y, 
                 h_ = bw,
                 PACKAGE = "lpme");
  }
  
  #########################################################################################
  # save state
  #########################################################################################
  model.name <- "local polynomial estimation";
  ngrid = length(xgrid);
  yhat = rep(0, ngrid);
  for(i in 1:ngrid){
    yhat[i] = (foo$ghat)[which.min(abs(x-xgrid[i]))]
  }
  
  #### Save to a list
  output <- list(modelname=model.name,
                 xgrid = xgrid, 
                 yhat = yhat);
  class(output) <- c("meanreg")
  output
}